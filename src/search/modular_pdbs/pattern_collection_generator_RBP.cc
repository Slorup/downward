//#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_generator_RBP.h"

//#include "../causal_graph.h"
//#include "../globals.h"
//#include "../task_proxy.h"

//#include "../utils/markup.h"
#include "../utils/math.h"
#include "../utils/rng.h"
//#include "../utils/timer.h"
//#include "../heuristic.h"

//#include <algorithm>
//#include <cassert>
//#include <iostream>
//#include <unordered_set>
//#include <vector>
//#include <math.h>
//Hack to use SS get_type, it needs heuristic object in constructor
//#include "../heuristics/blind_search_heuristic.h"
//#include "../heuristics/lm_cut_heuristic.h"
//#include "../successor_generator.h"
//#include "../utils/countdown_timer.h"
//#include "pdb_factory.h"
//#include "pattern_database_interface.h"
#include "../utils/debug_macros.h"
//#include <random>
//#include "../sampling.h"
#include "../task_tools.h"
#include "../options/option_parser.h"
#include "../options/plugin.h"
#include "../causal_graph.h"

    
using namespace std;
//PatternCollectionGeneratorComplementary is to be the driver for 8 PDB-based options
//RandomCollectionGeneration: CBP, RBP, CGamer, the one used by *Pommerening et al. 
//Local Search: iPDB, gaPDB, CGamer, VPN, *none, *changing order of patterns in gaPDB mutation. (this one matters to 0-1 greedy cost partitioning).
//GenPDB: Symbolic, Explicit, Online, expressed on pdb_factory class
//PDBEval: AvgH, Random sampling, Stratified sampling, *original iPDB method
//CombPDBs->Canonical, hPO, Max
//CostPartition->*None, Saturated, 0-1 greedy 
//Learning: UCB1 to choose bin packing, pdb size. 
//Re-evaluate: None, RemovedPDBsDominated, Run GHS (note: I think if we run GHS here, it may help to generate PDBs that are complementary to LM-cut)
// And the pseudo-code logic for it:
//While (time < 900 seconds)
//     PC<-RandomCollectionGeneration(MaxSize,CostPartition)
//     setInterestingPCs<-LocalSearch(P,PDBEval,MaxSize,CostPartition)
//     selectedPCs<-SubsetSelection(setInterestingPCs,PDBEval,ComPDBs)
//     generatedPDBs <-GenPDB(selectedPCs)
//     H<-Re-evaluate (H¿generatedPCs)
//     Learning(BinPackingRewards)
//ComPDBs (Hl)
namespace pdbs3 {

/*std::ostream& operator<< (std::ostream &out, const set<int> x)
{
    // Since operator<< is a friend of the Point class, we can access Point's members directly.
    for(auto element : x)
      out<<element<<",";
 
    return out;
}*/
  PatternCollectionGeneratorRBP::PatternCollectionGeneratorRBP(const options::Options & opts) :
    time_limit (opts.get<int>("time_limit")),
    pdb_factory (opts.get<shared_ptr<PDBFactory>>("pdb_factory")){
      cout<<"hello modular_RBP"<<endl;
    //num_vars=task->get_num_variables();
    }
    
  void PatternCollectionGeneratorRBP::initialize(std::shared_ptr<AbstractTask> task) {
    cout<<"Calling initialize GeneratorRBP-style,InSituCausalCheck:"<<InSituCausalCheck<<flush<<endl;
    num_vars= task->get_num_variables();
    cout<<"PatternCollectionGeneratorRBP::num_vars:"<<num_vars<<endl;
    //TaskProxy task_proxy_temp(*task);
    //task_proxy=make_shared<TaskProxy>(task_proxy_temp);
    task_proxy=make_shared<TaskProxy>(*task);
    VariablesProxy variables = task_proxy->get_variables();
    Pattern pattern;
    for(size_t i=0;i<num_vars;i++){
      //cout<<"i:"<<i<<flush<<endl;
      pattern.push_back(i);
    }
    overall_problem_size=get_pattern_size(pattern);
    cout<<"Overall problem size:"<<overall_problem_size<<endl;
    //Store goal variables for multiple use
    for (FactProxy goal : task_proxy->get_goals()) {
      initial_remaining_goal_vars.insert(goal.get_variable().get_id());
    }
    cout<<"initial_remaining_goal_vars:";for( auto element : initial_remaining_goal_vars) cout<<element<<",";cout<<endl;
    for (size_t i = 0; i < variables.size(); ++i) {
          initial_remaining_vars.insert(i);
    }
    //Some problems have such large sizes that even 
    //symbolic can not deal with them, limiting it
    //to 20 orders of magnitude for now
    //PDB SIZE BEING CHOSEN BY UBC
    //max_single_PDB_size = min(20.0,log10 (overall_problem_size));
    if(pdb_factory->name().find("symbolic")!=string::npos){
      //initial_max_target_size=max_target_size;
      if(max_single_PDB_size>10){//start with a cautious max_target_size for symbolic
        max_single_PDB_size=max_single_PDB_size/2.0;
      }
    }
    else{
      max_single_PDB_size=8;
      //initial_max_target_size=max_target_size;
      cout<<"initial time_limit="<<time_limit<<endl;
    }
      cout<<"starting with a max_single_PDB_size of:"<<max_single_PDB_size<<" because PDB_factory is"<<pdb_factory->name()<<endl;
      min_single_PDB_size=min(max_single_PDB_size-2.0,4.0);
      cout<<"starting with a min_single_PDB_size of:"<<min_single_PDB_size<<endl;
	  
      cout<<"Finished Calling initialize GeneratorRBP-style,InSituCausalCheck:"<<InSituCausalCheck<<flush<<endl;
    
      //VariablesProxy variables = task_proxy->get_variables();
      //const CausalGraph &causal_graph = task_proxy->get_causal_graph();
  }
  PatternCollectionContainer PatternCollectionGeneratorRBP::generate(){
    //cout<<"Calling GeneratorRBP-style,InSituCausalCheck:"<<InSituCausalCheck<<flush<<endl;exit(1);
    double limit=max_single_PDB_size;
    //double limit=10*pow(10,11);
    //cout<<"fixed_pdb_size limit to:"<<limit<<endl;
    //cout<<"single pattern only true"<<endl;
	  
    if(InSituCausalCheck)
      CBP_count++;
    else
      RBP_count++;
     
    PatternCollectionContainer PC;
    //First choose PDB size between min and max target sizes
		//int temp_max_single_PDB_size=rand()%(max_single_PDB_size-min_single_PDB_size+1);temp_max_single_PDB_size+=min_single_PDB_size;
    //cout<<"\t random_selection of max_single_PDB_size:,"<<temp_max_single_PDB_size<<",cout:,"<<RBP_count<<endl;
	
    //SHOULD MOVE THIS TO CLASS VARS and populate in initialize
    DEBUG_COMP(cout<<"SHOULD MOVE  variables & causal_graph to class-wide vars and populate in initialize, not here"<<endl;);
    VariablesProxy variables = task_proxy->get_variables();
    const CausalGraph &causal_graph = task_proxy->get_causal_graph();
    DEBUG_COMP(cout<<"RBP,max_pdb_size:"<<max_single_PDB_size<<",goals_to_add:"<<goals_to_add<<flush<<endl;);
   
    //Now generate initial PCs

    set<int> remaining_vars=initial_remaining_vars;
    set<int> remaining_goal_vars;
    remaining_goal_vars=initial_remaining_goal_vars;
    DEBUG_COMP(cout<<"\t\tgoals_to_add:"<<goals_to_add<<",out of:"<<remaining_goal_vars.size()<<endl;);
    //cout<<"total_vars:"<<remaining_vars.size()<<",goals:";for (auto var : remaining_goal_vars) cout<<var<<",";cout<<endl;exit(1);

    for (int i = 0; i < num_patterns; ++i) {

      DEBUG_COMP(cout<<"Working on pattern:"<<i<<endl;);
      //cout<<"Working on pattern:"<<i<<endl;
      //cout<<",remaining_vars:";for (auto var : remaining_vars) cout<<var<<",";cout<<endl;

        vector<int> vars_to_check;
        
        //cout<<"\tremaining_goal_vars:";for (auto var : remaining_goal_vars) cout<<var<<",";cout<<endl;

        vector<vector<bool>> pattern_collection;
        vector<bool> pattern(variables.size(), false);
        int var_id;
        double current_size = 1;
        
        vector<int> pattern_int;
        Pattern candidate_pattern;

        while(!remaining_vars.empty()){
          if(pattern_int.size()>0){
          candidate_pattern=pattern_int;
          sort(candidate_pattern.begin(), candidate_pattern.end());
          set<int> rel_vars_set;
          vector<int> relevant_vars;
          vector<int> relevant_vars_in_remaining;
	  //The only diference in variable selection
	  //between RBP and CBP is that CBP limits the number of variables
	  //To those remaining which also are causally connected to already
	  //selected variabels.  For RBP, a posteriori check is performed
	  //to remove any variables non casually related
	  if(InSituCausalCheck){
	    for (auto var : pattern_int){
	      const vector<int> &rel_vars = causal_graph.get_eff_to_pre(var);
	      for(auto var2 : rel_vars){
		rel_vars_set.insert(var2);
	      }
	    }
	    //cout<<"\tpattern_int:";print_vect_int(pattern_int);cout<<",rel_vars_set:";for(auto element : rel_vars_set) cout<<element<<",";cout<<endl;
	  }
	  else{
	    rel_vars_set=remaining_vars;
	  }
          set_difference(rel_vars_set.begin(), rel_vars_set.end(),
            candidate_pattern.begin(), candidate_pattern.end(),
            back_inserter(relevant_vars));
            //cout<<"relevant vars to current_pattern:";for (auto item : relevant_vars) cout<<item<<",";cout<<endl;
          set_intersection(relevant_vars.begin(), relevant_vars.end(),
            remaining_vars.begin(), remaining_vars.end(),
            back_inserter(relevant_vars_in_remaining));
            //cout<<"\t\tbefore shuffle,relevant vars in remaining:";for (auto item : relevant_vars_in_remaining) cout<<item<<",";cout<<flush<<endl;
            g_rng()->shuffle(relevant_vars_in_remaining);
            //cout<<"\t\tafter shuffle,relevant vars in remaining:";for (auto item : relevant_vars_in_remaining) cout<<item<<",";cout<<flush<<endl;
          while(relevant_vars_in_remaining.size()>0){
            var_id=relevant_vars_in_remaining.back();
            relevant_vars_in_remaining.pop_back();
            double next_var_size = variables[var_id].get_domain_size();
            if(utils::is_product_within_limit(current_size, next_var_size, limit)) {
              candidate_pattern.push_back(var_id);
              current_size *= next_var_size;
              pattern[var_id] = true;
              remaining_vars.erase(var_id);
              remaining_goal_vars.erase(var_id);
              //cout<<"\t\tadded to pattern var_id:"<<var_id<<",current_size:"<<current_size<<",max_single_PDB_size:"<<max_single_PDB_size<<",new_pattern:";print_vect_int(candidate_pattern);cout<<"relevant_vars left:,";for (auto element : relevant_vars_in_remaining) cout<<element<<",";cout<<endl;//",remaining vars:"<<remaining_vars.size()<<endl;
              break;
            }
          }

          if(candidate_pattern!=pattern_int){

            pattern_int=candidate_pattern;
            if(remaining_vars.size()==0){//one single pattern with all vars, add now!
	      //Pattern creaton finished, removed irrelevant vars and return to list of available vars
	      //if doing disjunctive patterns and a goal is remaining, otherwise not needed because for 
	      //every new pattern all vars are available anyhow
	      Pattern removed_vars;
	      remove_irrelevant_variables(pattern_int,removed_vars);

	      if(disjunctive_patterns&&remaining_goal_vars.size()>0){
		remaining_vars.insert(removed_vars.begin(), removed_vars.end());
		cout<<"\t\tinserted removed_vars:";for(auto element : removed_vars) cout<<element<<",";cout<<endl;
	      }
	      //cout<<"\t\t no remaining_vars,adding pc:";for (auto element : pattern_int) cout<<element<<",";cout<<endl;
              PC.add_pc(pattern_int);
              break;
            }
          }
          else{//no var is small enough to be added, or none left
            if(pattern_int.size()>0){
	      //Pattern creaton finished, removed irrelevant vars and return to list of available vars
	      //if doing disjunctive patterns and a goal is remaining, otherwise not needed because for 
	      //every new pattern all vars are available anyhow
	      Pattern removed_vars;
	      remove_irrelevant_variables(pattern_int,removed_vars);

	      if(disjunctive_patterns&&remaining_goal_vars.size()>0){
		remaining_vars.insert(removed_vars.begin(), removed_vars.end());
		//cout<<"\t\tinserted removed_vars:";for(auto element : removed_vars) cout<<element<<",";cout<<endl;
	      }

	      //cout<<"\t\t no var is small enough,adding pc:";for (auto element : pattern_int) cout<<element<<",";cout<<endl;
              PC.add_pc(pattern_int);
            }
            else{
              cout<<"candidate_pattern size is empty!"<<endl;exit(1);
            }

            if(pattern_int.size()>0){
              pattern_collection.push_back(pattern);
              vector<int> trans_pattern=transform_to_pattern_normal_form(pattern_collection.back());
              //cout<<"\tno more relevant vars can be added, final_pattern["<<i<<"]:";print_vect_int(candidate_pattern);
              //cout<<",size:"<<get_pattern_size(trans_pattern)<<flush<<endl;
              //cout<<"added pattern["<<pattern_collection.size()-1<<"]:"<<trans_pattern<<",size:"<<get_pattern_size(trans_pattern)<<flush<<endl;
              pattern_int.clear();
              pattern.clear();
              pattern.resize(variables.size(), false);
              current_size = 1;
            }
          }
          }
          else{//choose a remaining var at random, nothing selected yet for this pattern
            auto temp_it=remaining_goal_vars.begin();
            if(remaining_goal_vars.empty()){
              //no more goal vars, so no more patterns, as we can not start it with a goal variable
              break;
            }
            double next_var_size=0;
            for(unsigned goal_vars=0;goal_vars<goals_to_add;goal_vars++){
              if(remaining_goal_vars.size()==0){
                break;
              }
              temp_it=remaining_goal_vars.begin();
              advance(temp_it,rand()%remaining_goal_vars.size());
              var_id=*temp_it;
              remaining_goal_vars.erase(temp_it);
              remaining_vars.erase(var_id);
              pattern[var_id]=true;
              pattern_int.push_back(var_id);
              next_var_size = variables[var_id].get_domain_size();
              current_size *= next_var_size;
            }
            //cout<<"starting vars for pattern:";
            //for(auto id : pattern_int) cout<<","<<id;
            //cout<<",remaining_goal_vars:";
            //for(auto id : remaining_goal_vars) cout<<","<<id;
            //cout<<endl;
          //cout<<",remaining_vars:";
          //for(auto id : remaining_vars) cout<<","<<id;
          //cout<<endl;
        
        }
      }

      if(pattern_int.size()==initial_remaining_vars.size()){//So pattern had all vars, we are finished!){
        break;
      }
      if(!disjunctive_patterns){
        remaining_vars=initial_remaining_vars;
        remaining_goal_vars=initial_remaining_goal_vars;
      }
      else if(remaining_goal_vars.size()==0){
        break;
      }
    }
    //  cout<<"Generate PC:";PC.print();

    return PC;
  }
  void PatternCollectionGeneratorRBP::remove_irrelevant_variables(Pattern &pattern,Pattern &removed_vars){
	set<int> in_original_pattern(pattern.begin(), pattern.end());
	set<int> in_pruned_pattern;
	removed_vars.clear();
	
	vector<int> vars_to_check;
	for ( auto var_id : initial_remaining_goal_vars) {
	  if (in_original_pattern.count(var_id)) {
	    vars_to_check.push_back(var_id);
	    in_pruned_pattern.insert(var_id);
	  }
	}
	//cout<<"vars_to_check:";print_vect_int(vars_to_check);
	    
	const CausalGraph &causal_graph = task_proxy->get_causal_graph();
	while (!vars_to_check.empty()) {
	    int var = vars_to_check.back();
	    vars_to_check.pop_back();
	    const vector<int> &rel = causal_graph.get_eff_to_pre(var);
	    for (size_t i = 0; i < rel.size(); ++i) {
		int var_no = rel[i];
		if (in_original_pattern.count(var_no) &&
		    !in_pruned_pattern.count(var_no)) {
		    // Parents of relevant variables are causally relevant.
		    vars_to_check.push_back(var_no);
		    in_pruned_pattern.insert(var_no);
		}
	    }
	}
	pattern.assign(in_pruned_pattern.begin(), in_pruned_pattern.end());
	sort(pattern.begin(), pattern.end());
        //Now get the vars which were removed  
	set_difference(in_original_pattern.begin(), in_original_pattern.end(),
            in_pruned_pattern.begin(), in_pruned_pattern.end(),
            back_inserter(removed_vars));
	//cout<<"in_original_pattern:";for(auto element : in_original_pattern) cout<<element<<",";
	//cout<<"in_prunned_pattern:";for(auto element: in_pruned_pattern) cout <<element<<","; 
	//cout<<"removed_vars:";print_vect_int(removed_vars);cout<<endl;
	if(InSituCausalCheck)
	  if(in_original_pattern.size()>in_pruned_pattern.size()){
	    cerr<<"Bug, if we are doing CBP, we should never select causally disconnected variables!";
	    cerr<<"in_original_pattern:";for(auto element : in_original_pattern) cout<<element<<",";
	    cerr<<"in_prunned_pattern:";for(auto element: in_pruned_pattern) cout <<element<<","; 
	    cerr<<"removed_vars:";print_vect_int(removed_vars);cout<<endl;
	    exit(1);
	  }
  }

  static shared_ptr<PatternCollectionGeneratorComplementary>_parse(options::OptionParser &parser) {
    parser.add_option<int> ("time_limit", "If populated,stop construction on first node past boundary and time limit", "100");
    
    parser.add_option<shared_ptr<PDBFactory>>(
        "pdb_factory",
        "See detailed documentation for pdb factories. ",
	      "modular_symbolic");
    
    //parser.add_option<unsigned> ("goals_to_add", "How many goals it uses to seed each pattern","0");

    options::Options options = parser.parse();
    parser.document_synopsis(
        "Pattern Generator for both RBP and CBP pattern generation styles",
        "selection of variables to generate Pattern Collection");
    options::Options opts = parser.parse();
    if (parser.dry_run())
      return 0;
    return make_shared<PatternCollectionGeneratorRBP>(opts);
  }

  static options::PluginShared<PatternCollectionGeneratorComplementary> _plugin("modular_rbp", _parse);
}
