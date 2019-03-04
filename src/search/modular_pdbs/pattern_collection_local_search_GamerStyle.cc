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
  :time_limit (opts.get<int>("time_limit")){
  //PatternCollectionLocalSearchGamerStyle::PatternCollectionLocalSearchGamerStyle() 
      cout<<"hello LocalSearchGamerStyle"<<endl;
    //num_vars=task->get_num_variables();
  }
    
  
    
  bool PatternCollectionLocalSearchGamerStyle::do_local_search(shared_ptr<PatternCollectionInformation> current_result, shared_ptr<PatternCollectionEvaluator> evaluation_method,
     shared_ptr<PDBFactory> pdb_factory){
    int start_local_search_time=utils::g_timer();
    bool improvement_found=false;
    TaskProxy task_proxy(*(g_root_task()));
    size_t num_vars = task_proxy.get_variables().size();
    shared_ptr<ModularZeroOnePDBs> candidate_ptr;
    PatternCollectionContainer new_candidate_local_search;
    const State &initial_state = task_proxy.get_initial_state();
    
    cout<<"Starting do_local_search:"<<get_name()<<",num_vars:"<<num_vars<<",local_episodes:"<<get_episodes()<<",time_limit:"<<get_time_limit()<<endl;
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
    
    for(size_t pattern=0;pattern<top_patterns.size();pattern++){
      vector<size_t> improving_vars;
      PatternCollectionContainer PC_old=max_subsets_container.at(pattern);
      PatternCollectionContainer PC_new;
      Pattern old_pattern=top_patterns.at(pattern);
      if(impossible_to_update_pattern.find(old_pattern)!=impossible_to_update_pattern.end()){
	cout<<"do_local_search, impossible to update"<<old_pattern<<endl;
	continue;
      }

      //forbidden_vars.clear();
      for (auto var : old_pattern) forbidden_vars[old_pattern].insert(var);
      bool all_pdbs_finished=true;
      for(size_t var=0;var<num_vars-forbidden_vars[old_pattern].size();var++){
	PC_new=generate_next_candidate(PC_old);
	candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, PC_new.get_PC(), *pdb_factory);
	//DEBUG
	if(!candidate_ptr->is_finished()){
	  all_pdbs_finished=false;
	  cout<<",candidate PC not finished,";cout<<"PC_new:,";PC_new.print();
	}
	if(pdb_factory->is_solved()){
	      cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
	      current_result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
	      current_result->set_dead_ends(pdb_factory->get_dead_ends());
	      return true;
	}
	if(evaluation_method->evaluate(candidate_ptr)){
	  //UPDATE PC_old to improved PC&add New_PC to result&resample
	  improving_vars.push_back(last_var);
	  PC_old=PC_new;
	  cout<<"time:,"<<utils::g_timer()<<",improvement found while doing local_search,var:,"<<var<<",type:"<<get_name();cout<<"PC_new:,";PC_new.print();
	  improvement_found=true;
	  current_result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
	  current_result->set_dead_ends(pdb_factory->get_dead_ends());
	  //GAMER NEEDS TO PRUNE PREVIOUS PATTERN OR IT WOULD END UP WITH LOTS OF DOMINATED PDBs
	  cout<<"time:"<<utils::g_timer()<<",initial_h before recompute additive sets:,"<<current_result->get_value(initial_state)<<endl;
	  current_result->recompute_max_additive_subsets();
	  cout<<"time:"<<utils::g_timer()<<",initial_h after recompute:,"<<current_result->get_value(initial_state)<<endl;
	  evaluation_method->sample_states(current_result);
	}
	/*else{
	  cout<<"time:,"<<utils::g_timer()<<",no improvement found while doing local_search,var:,"<<var<<",type:"<<get_name();cout<<"PC_new:,";
	  PC_new.print();
	}*/

	if(start_local_search_time-utils::g_timer()>time_limit){
	  cout<<"Interrupting do_local_search:,"<<get_name()<<",time_spent:,"<<utils::g_timer()-start_local_search_time<<",time_limit:"<<time_limit<<endl;
	  return improvement_found;
	}
      }
      if(improving_vars.size()==0&&all_pdbs_finished==true){
	cout<<"Impossible for do_local_search to improve pattern:"<<old_pattern<<"by adding a single var"<<endl;
	impossible_to_update_pattern.insert(old_pattern);
      }
      cout<<"improving_vars:";for (auto var : improving_vars) cout<<var<<",";cout<<endl;
    }
    

    //FINALLY ADD BEST VARS TO TOP PDB, AND RECALCULATE COSTS FOR REMAINING PDBs IN PC
    //NEED TO CODE IT!
    cout<<"Current samples:"<<evaluation_method->get_num_samples()<<endl;
    cout<<"Finished do_local_search:,"<<get_name()<<",time_spent:,"<<utils::g_timer()-start_local_search_time<<endl;
    return improvement_found;
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
    //cout<<"\t forbidden_vars_size:"<<forbidden_vars.size()<<",candidate_vars_size:"<<candidates.size()<<endl;
    if(candidates.size()==0){
      new_patterns.restart_pc(candidate_pattern);
      return new_patterns;
    }
    //cout<<"input_pattern:";for (auto i : candidate_pattern) cout<<i<<",";
    //cout<<"candidates:";for (auto i : candidates) cout<<i<<",";cout<<endl;
    assert(candidates.size()>0);
    last_var=candidates.at(rand()%candidates.size());
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

  static shared_ptr<PatternCollectionLocalSearch>_parse(options::OptionParser &parser) {
    parser.add_option<int> ("time_limit", "If populated,stop construction on first node past boundary and time limit", "100");
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
