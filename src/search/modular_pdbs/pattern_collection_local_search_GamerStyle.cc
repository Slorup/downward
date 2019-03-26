#include "pattern_collection_local_search_GamerStyle.h"
#include "../sampling.h"
#include "../utils/timer.h"
#include "../task_tools.h"
#include "../utils/countdown_timer.h"
#include "../options/option_parser.h"
#include "../options/plugin.h"

#include "../causal_graph.h"
//#include "../globals.h"
//#include "../task_proxy.h"

//#include "../utils/markup.h"
//#include "../utils/math.h"
//#include "../utils/rng.h"
//#include "../utils/timer.h"
//#include "../heuristic.h"
//#include "../utils/collections.h"

#include "pdb_factory.h"
#include "../utils/debug_macros.h"
#include "pattern_collection_generator_complementary.h"
#include "pattern_database_interface.h"
#include<algorithm>

    
using namespace std;
using namespace utils;
namespace pdbs3 {
template < class T >
inline std::ostream& operator << (std::ostream& os, const std::vector<T>& v) 
{
    os << "[";
    //for (std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    //for (std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    for(auto i : v)
    {
        //os << *ii << ",";
        os << i << ",";
    }
    os << "]";
    return os;
}

  PatternCollectionLocalSearchGamerStyle::PatternCollectionLocalSearchGamerStyle(const options::Options & opts)
  :time_limit(opts.get<int>("time_limit")),
      verbose(opts.get<bool>("verbose")){
  //PatternCollectionLocalSearchGamerStyle::PatternCollectionLocalSearchGamerStyle() 
      cout<<"hello LocalSearchGamerStyle"<<endl;
    local_search_timer = make_unique<utils::CountdownTimer>(time_limit);
    //num_vars=task->get_num_variables();
   //Setting up costs for top_patterns
      TaskProxy task_proxy(*(g_root_task()));
      OperatorsProxy operators = task_proxy.get_operators();
      operator_costs.reserve(operators.size());
      for (OperatorProxy op : operators)
	operator_costs.push_back(op.get_cost());
  }
    
  
    
  bool PatternCollectionLocalSearchGamerStyle::do_local_search(shared_ptr<PatternCollectionInformation> current_result, shared_ptr<PatternCollectionEvaluator> evaluation_method,
     shared_ptr<PDBFactory> pdb_factory){
    int start_local_search_time=utils::g_timer();
    //bool improvement_found=false;
    TaskProxy task_proxy(*(g_root_task()));
    size_t num_vars = task_proxy.get_variables().size();
    PatternCollectionContainer new_candidate_local_search;
    const State &initial_state = task_proxy.get_initial_state();


    if(verbose){
      cout<<"eval_method threshold:"<<evaluation_method->get_threshold()<<endl;
      cout<<"Starting do_local_search:"<<get_name()<<",num_vars:"<<num_vars<<",local_episodes:"<<get_episodes()<<",time_limit:"<<local_search_timer->get_remaining_time()<<endl;
    }
    //FIRST GET TOP PDBs IN EACH MAX_ADDITIVE_SUBSETS, SO THEY HAVE FULL COSTS
    //std::shared_ptr<PatternCollection> current_patterns=make_shared<PatternCollection>(*current_result->get_patterns());
    std::shared_ptr<MaxAdditivePDBSubsets> current_subsets=current_result->get_max_additive_subsets();

    vector<Pattern> top_patterns;
    vector<PatternCollectionContainer> max_subsets_container;
    PatternCollectionContainer PC_temp;
    for(size_t subset=0;subset<current_subsets->size();subset++){
      auto pdb=current_subsets->at(subset).at(0);
      auto pattern = pdb->get_pattern();
      for (const shared_ptr<PatternDatabaseInterface> &pdb : current_subsets->at(subset)) {
	PC_temp.clear();PC_temp.add_pc(pdb->get_pattern());
      }
      max_subsets_container.push_back(PC_temp);
      top_patterns.push_back(pattern);
      cout<<"pattern["<<subset<<"]:";
      for(auto var : pattern){
	cout<<var<<",";
      }
      cout<<endl;
    }
    cout<<"top_patterns.size:"<<top_patterns.size()<<endl;

    //THEN GET EACH CANDIDATE PC BY ADDING A SINGLE VAR TO TOP PDB 
    //AND RECALCULATING COSTS FOR REMAINING PDBS IF ANY
    bool improvement_found=false;

    for(auto pattern : top_patterns){
      if(local_search_timer->is_expired()){
	return improvement_found;
      }
      if(impossible_to_update_pattern.find(pattern)!=impossible_to_update_pattern.end()){
	cout<<"do_local_search, impossible to update"<<pattern<<endl;
	continue;
      }
      if(do_local_search_pattern(pattern,current_result,evaluation_method,pdb_factory))
	  improvement_found=true;
    }
    

    //FINALLY ADD BEST VARS TO TOP PDB, AND RECALCULATE COSTS FOR REMAINING PDBs IN PC
    //NEED TO CODE IT!
    cout<<"Current samples:"<<evaluation_method->get_num_samples()<<endl;
    cout<<"Finished do_local_search:,"<<get_name()<<",time_spent:,"<<utils::g_timer()-start_local_search_time<<endl;
    
    return improvement_found;
  }
 bool PatternCollectionLocalSearchGamerStyle::do_local_search_pattern(Pattern old_pattern,shared_ptr<PatternCollectionInformation> current_result,
     shared_ptr<PatternCollectionEvaluator> evaluation_method, shared_ptr<PDBFactory> pdb_factory){
   if(pdb_factory->is_solved()){
     cout<<"calling do_local_search_pattern with problem already solved!"<<endl;
     return false;
   }
   //Pattern temp2={0,3,6,7,28,32,33,37,38,42,43,47,58,60,63,64,65,66,67,68,69,70,71,72,73,74,76,77,78,79};
   //old_pattern=temp2;
   vector<int> candidates;
   set<int> input_pattern(old_pattern.begin(), old_pattern.end());
   TaskProxy task_proxy(*(g_root_task()));
   const CausalGraph &cg = task_proxy.get_causal_graph();
    
   //shared_ptr<PatternDatabaseInterface> old_pdb = pdb_factory->compute_pdb(task_proxy, old_pattern,operator_costs);

   bool all_pdbs_finished=true;
   //double starting_best_score=best_score;
   //We get baseline to compare every time we do local_search_pattern
   //For some methods sample changes constantly, for others it is a function of 
   //whether there has been improvements, we should add a marker to determine whether 
   //a new sample is necessary
   evaluation_method->sample_states(current_result);
   double starting_best_score=evaluation_method->get_sample_score();
   double best_score=starting_best_score;
   //cout<<"starting_best_score:"<<best_score<<",eval_method threshold:"<<evaluation_method->get_threshold()<<endl;
   //First get connected variables 
   for (size_t var = 0; var < g_variable_domain.size(); ++var) {
        if (input_pattern.count(var)){ 
          DEBUG_MSG(cout<<"\t\tGamer_local_search,skipping exisiting var:"<<var<<endl;);
          continue;
        }
      for (int succ : cg.get_pre_to_eff(var)) {
        if (input_pattern.count(succ)) {
          DEBUG_MSG(cout<<"\t\tGamer,connected variables:"<<succ<<"to var:"<<var<<" added."<<endl;);
          candidates.push_back(var); 
          break;
        }
      }
    }
    if(candidates.size()==0){
      impossible_to_update_pattern.insert(old_pattern);
      if(verbose)
	cout<<"input_pattern:"<<old_pattern<<"cannot be updated, no connected variables"<<endl;
      return false;
    }
    if(verbose){
      cout<<"input_pattern:";for (auto i : old_pattern) cout<<i<<",";
      cout<<"candidates:";for (auto i : candidates) cout<<i<<",";cout<<endl;
    }
    assert(candidates.size()>0);
    std::random_shuffle(candidates.begin(), candidates.end());//to aboid biases on repeat runs with partial PDBs
    vector<pair<int,double> > improving_vars;//var,score
   
    vector<shared_ptr<PatternDatabaseInterface> > improving_pdbs;
    
    //SET PRECONSTRUCTION TIME AS A FUNCTION OF NUMER OF VARIABLES AND AIALABLE TIME
    //We set the minimum time to be 100 ms
    double pdb_preconst_time=1000.0*(double(local_search_timer->get_remaining_time())/double(candidates.size()))*0.6666;
    pdb_preconst_time=max(100,int(pdb_preconst_time));
    cout<<"local_search_timer:"<<local_search_timer->get_remaining_time()<<",pdb_preconst_time:"<<pdb_preconst_time<<",num_candidates:"<<candidates.size()<<endl;
    pdb_factory->set_new_max_time(int(pdb_preconst_time));

    //Now create candidate patterns to evaluate which vars to add
    while(candidates.size()){
      if(local_search_timer->is_expired()){
	all_pdbs_finished=false;
	if(verbose)
	  cout<<"Breaking out of candidate loop, local_search_timer is expired"<<endl;
	break;
      }
      last_var=candidates.back();
      candidates.pop_back();
      Pattern candidate_pattern=old_pattern;
      candidate_pattern.push_back(last_var);
      sort(candidate_pattern.begin(), candidate_pattern.end());
      int start_pdb_time=utils::g_timer();
      shared_ptr<PatternDatabaseInterface> candidate_pdb = pdb_factory->compute_pdb(task_proxy, candidate_pattern, operator_costs);
      if(verbose)
	cout<<"candidate:"<<candidate_pattern<<",construction time:"<<utils::g_timer()-start_pdb_time<<endl;
      PDBCollection temp;temp.push_back(candidate_pdb);
      
      if(pdb_factory->is_solved()){
	current_result->include_additive_pdbs(pdb_factory->terminate_creation(temp,0, 0, 0));
	current_result->set_dead_ends(pdb_factory->get_dead_ends());
	return false;
      }

      if(local_search_timer->is_expired()){
	all_pdbs_finished=false;
	if(verbose)
	  cout<<"Breaking out of candidate loop, local_search_timer is expired"<<endl;
	break;
      }

      //No calling terminate creation when using online_pdbs
      /*if(pdb_factory->name().find("online")==string::npos){
	if(verbose)
	  cout<<"calling terminate_creation for each candidate in do_local_search, adding 50 secs"<<endl;
      }*/
      temp.clear();
	
      if(candidate_pdb->is_finished()){
	if(verbose)
	  cout<<"time,"<<utils::g_timer()<<",candidate_pattern:"<<candidate_pattern<<", is finished"<<endl;
      }
      else{
	if(verbose)
	  cout<<"time:,"<<utils::g_timer()<<",candidate_pattern:"<<candidate_pattern<<", is unfinished"<<endl;
	all_pdbs_finished=false;
      }
      //Now evaluate PDB
      shared_ptr<ModularZeroOnePDBs> candidate_ptr=make_shared<ModularZeroOnePDBs>(candidate_pdb);
      evaluation_method->evaluate(candidate_ptr);
      if(evaluation_method->get_eval_score()>starting_best_score){
        //if(candidate_pdb->compute_mean_finite_h()>starting_best_score)//If we find 2 equal improvements to best known we add both
	if(verbose)
	  cout<<"Selected,sample_score:"<<starting_best_score<<",previous_best_score:,"<<best_score<<",eval_score:,"<<evaluation_method->get_eval_score()<<",last_var:"<<last_var<<endl;
	improvement_found=true;
	//improving_vars.push_back(make_pair<int,double>(last_var,evaluation_method->get_reward()));
	improving_vars.push_back(make_pair<int,double>(int(last_var),evaluation_method->get_eval_score()));
	improving_pdbs.push_back(candidate_pdb);
	//best_score=max(best_score,evaluation_method->get_reward());
	best_score=max(best_score,evaluation_method->get_eval_score());
      }
      else{
	if(verbose)
	  cout<<"Unselected,best_score:"<<best_score<<",eval_score:"<<evaluation_method->get_eval_score()<<",last_var:"<<last_var<<endl;
      }
    }

    if(improving_vars.size()==0){
	if(pdb_factory->is_solved()){
	  cerr<<"no_improving_vars,problem solved with candidate PDB!"<<endl;exit(1);
	}
      if(all_pdbs_finished){
	if(verbose)
	  cout<<"adding to impossible_to_update,pattern:"<<old_pattern<<",all_pdbs were finished"<<endl;
	impossible_to_update_pattern.insert(old_pattern);
      }
      return false;
    }
    
    //Now find those var(s) with the best improvement
    //We remove from improving vars all variables 
    //whose score is 0.999 or lower of best_score
    //Same as CGAMER
    cout<<"Improving vars before prune:"<<endl;
    for (auto var : improving_vars) cout<<"before_prune_var:,"<<var.first<<",score:,"<<var.second<<endl;
    improving_vars.erase(std::remove_if(
	  improving_vars.begin(), 
	  improving_vars.end(),
	  [best_score](pair<int,double> x){
	      return x.second < 0.999*best_score;
	  }), improving_vars.end());
    cout<<"Improving vars after prune:"<<endl;
    Pattern new_pattern=old_pattern;
    for (auto var : improving_vars) {
      cout<<"after_prune_var:,"<<var.first<<",score:,"<<var.second<<endl;
      new_pattern.push_back(var.first);
    }
    
    //Note: if we have gone past time limit, we simply add
    //the best improving PDB, we have run out of time
    //to do any merges
    if(local_search_timer->is_expired()){
      assert(improving_pdbs.size()>0);//we should at least have keept one improving_var or returned false previously
      cout<<"time:,"<<utils::g_timer()<<",local_search_time is expired,last best candidate PDB is Final PDB:"<<*improving_pdbs.back()<<endl;
      PDBCollection new_Coll;
      new_Coll.push_back(improving_pdbs.back());//last improving PDB will be have the best found avg_h_value(ignoring ties for now)
      //current_result->include_additive_pdbs(pdb_factory->terminate_creation(new_Coll,0,0,0));
      current_result->include_additive_pdbs(pdb_factory->terminate_creation(new_Coll,250000, 50000, 10000000));
      current_result->set_dead_ends(pdb_factory->get_dead_ends());
      //Recompute as usual
      const State &initial_state = task_proxy.get_initial_state();
      cout<<"time:"<<utils::g_timer()<<",initial_h before recompute:,"<<current_result->get_value(initial_state)<<endl;
      current_result->recompute_max_additive_subsets();
      cout<<"time:"<<utils::g_timer()<<",initial_h after recompute:,"<<current_result->get_value(initial_state)<<endl;
      return true;
    }






    sort(new_pattern.begin(),new_pattern.end());
    cout<<"time:,"<<utils::g_timer()<<",Final PDB:,"<<new_pattern<<endl;
    shared_ptr<PatternDatabaseInterface> new_pdb = pdb_factory->compute_pdb(task_proxy, new_pattern, operator_costs);


    //Repeat timer_check before going to terminate final PDB
    if(local_search_timer->is_expired()){
      PDBCollection new_Coll;
      shared_ptr<ModularZeroOnePDBs> temp_ptr=make_shared<ModularZeroOnePDBs>(new_pdb);
      if(evaluation_method->evaluate(temp_ptr)){
	new_Coll.push_back(new_pdb);//
	cout<<"time:,"<<utils::g_timer()<<",local_search_time is expired,last best candidate PDB is Final PDB:"<<*new_pdb<<endl;
      }
      else{
	new_Coll.push_back(improving_pdbs.back());
	cout<<"time:,"<<utils::g_timer()<<",local_search_time is expired before terminate, Final PDB is last improving PDB:"<<*improving_pdbs.back()<<endl;
      }
      //We will go past time limit but this is the last PDB
      current_result->include_additive_pdbs(pdb_factory->terminate_creation(new_Coll,250000, 50000, 10000000));
      //current_result->include_additive_pdbs(pdb_factory->terminate_creation(new_Coll,0,0,0));
      current_result->set_dead_ends(pdb_factory->get_dead_ends());
      //Recompute as usual
      const State &initial_state = task_proxy.get_initial_state();
      cout<<"time:"<<utils::g_timer()<<",initial_h before recompute:,"<<current_result->get_value(initial_state)<<endl;
      current_result->recompute_max_additive_subsets();
      cout<<"time:"<<utils::g_timer()<<",initial_h after recompute:,"<<current_result->get_value(initial_state)<<endl;
      return true;
    }


    PDBCollection temp;temp.push_back(new_pdb);
    pdb_factory->terminate_creation(temp, 250000, 50000, 10000000);
    //WE DO NOT ALLOW TO TRY THE SAME INPUT PATTERN TWICE AS LONG AS WE GOT TO THE STAGE
    //OF DOING A 250 SECS TERMINATE ONCE
    impossible_to_update_pattern.insert(old_pattern);
   
    temp.clear();
    if(!new_pdb->is_finished()){
	  cout<<"Final PDB was unfinished"<<endl;

    }
    cout<<"time:,"<<utils::g_timer()<<",Final PDB constructed"<<endl;
    //Now update the current result
    PDBCollection new_Coll;
    new_Coll.push_back(new_pdb);
    current_result->include_additive_pdbs(pdb_factory->terminate_creation(new_Coll));
    current_result->set_dead_ends(pdb_factory->get_dead_ends());
    //Now get rid of old PDBs
		
    const State &initial_state = task_proxy.get_initial_state();
    cout<<"time:"<<utils::g_timer()<<",initial_h before recompute:,"<<current_result->get_value(initial_state)<<endl;
    current_result->recompute_max_additive_subsets();
    cout<<"time:"<<utils::g_timer()<<",initial_h after recompute:,"<<current_result->get_value(initial_state)<<endl;
    
    return true;
 }


  PatternCollectionContainer PatternCollectionLocalSearchGamerStyle::generate_next_candidate(PatternCollectionContainer candidate_collection){
    DEBUG_MSG(cout<<"starting generate_next_candidate"<<endl;);
    //If more than one pattern, we only do local search on the first pattern
    Pattern candidate_pattern=candidate_collection.get_top_pattern();
    std::set<int> pattern;
    for (auto var : candidate_pattern) pattern.insert(var);


    TaskProxy task_proxy(*(g_root_task()));
    const CausalGraph &cg = task_proxy.get_causal_graph();
    vector<int> candidates;
    for (size_t var = 0; var < g_variable_domain.size(); ++var) {
        if (pattern.count(var)){ 
          DEBUG_MSG(cout<<"\t\tGamer_local_search,skipping exisiting var:"<<var<<endl;);
          continue;
        }
	else if(forbidden_vars[candidate_pattern].count(var)){//This pattern was proven to not improve search
	  continue;
	}

      for (int succ : cg.get_pre_to_eff(var)) {
        if (pattern.count(succ)) {
          DEBUG_MSG(cout<<"\t\tGamer,connected variables:"<<succ<<"to var:"<<var<<" added."<<endl;);
          candidates.push_back(var); 
          break;
        }
      }
    }

    if(candidates.size()==0){
      new_patterns.restart_pc(candidate_pattern);
      return new_patterns;
    }
    

    //cout<<"input_pattern:";for (auto i : candidate_pattern) cout<<i<<",";
    //cout<<"candidates:";for (auto i : candidates) cout<<i<<",";cout<<endl;
    assert(candidates.size()>0);
    //last_var=candidates.at(rand()%candidates.size());
    last_var=candidates[0];
    forbidden_vars[candidate_pattern].insert(last_var);
    //cout<<"adding random_var:"<<last_var<<endl;
    candidate_pattern.push_back(last_var);
    //All patterns have to be sorted! for duplication checks
    //and any other comparisons of patterns which are vectors,
    //so order matters!
    sort(candidate_pattern.begin(), candidate_pattern.end());
    new_patterns.restart_pc(candidate_pattern);
    //cout<<"new_patterns:";new_patterns.print();
    return new_patterns;//Not adding collection
  }
  bool PatternCollectionLocalSearchGamerStyle::impossible_to_improve(shared_ptr<PatternCollectionInformation> current_result){//with regards to adding one variable
    bool impossible_to_improve=true;
    std::shared_ptr<MaxAdditivePDBSubsets> current_subsets=current_result->get_max_additive_subsets();
    if(current_subsets->size()==0){//no PDBs in result, so impossible to improve!
      if(verbose) cout<<"no patterns in results, so impossible to improve!!!"<<endl;
      return true;
    }
    for(size_t subset=0;subset<current_subsets->size();subset++){
      auto pdb=current_subsets->at(subset).at(0);
      auto pattern = pdb->get_pattern();
      if(impossible_to_update_pattern.find(pattern)==impossible_to_update_pattern.end()){
	if(verbose) cout<<"pattern:"<<pattern<<"is possible to improve as far as we know"<<endl;
	impossible_to_improve=false;
	break;
      }
    }
    
    if(verbose)
      cout<<"do_local_search, posible_improvement:"<<impossible_to_improve<<endl;
    return impossible_to_improve;
  }

  static shared_ptr<PatternCollectionLocalSearch>_parse(options::OptionParser &parser) {
    parser.add_option<int> ("time_limit", "If populated,stop construction on first node past boundary and time limit", "100");
    parser.add_option<bool> ("verbose", "debug_mode from command line", "false");
    options::Options options = parser.parse();
    parser.document_synopsis(
        "Pattern Generator Local Search module",
        "Selection of variables to generate Pattern Collection");
    options::Options opts = parser.parse();
    if (parser.dry_run())
        return 0;

    return make_shared<PatternCollectionLocalSearchGamerStyle>(opts);
  }
  

  static options::PluginShared<PatternCollectionLocalSearch> _plugin("local_search_gamer", _parse);
}
