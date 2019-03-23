//#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_evaluator_RandWalk_Avg_h.h"
#include "../sampling.h"
#include "../utils/timer.h"
#include "../task_tools.h"
#include "../utils/countdown_timer.h"
#include "../options/option_parser.h"
#include "../options/plugin.h"

//#include "../causal_graph.h"
//#include "../globals.h"
//#include "../task_proxy.h"

//#include "../utils/markup.h"
//#include "../utils/math.h"
//#include "../utils/rng.h"
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
#include "pdb_factory.h"
//#include "pattern_database_interface.h"
#include "../utils/debug_macros.h"
//#include <random>
//#include "../sampling.h"
#include "../state_registry.h"

    
using namespace std;
//PatternCollectionEvaluator is to be the driver for 8 PDB-based options
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
PatterCollectionEvaluatorRandWalk_Avg_H::PatterCollectionEvaluatorRandWalk_Avg_H(const options::Options & opts) :
	time_limit (opts.get<int>("time_limit")){
    cout<<"hello EvaluatorRandWalk_Avg_H"<<flush<<endl;
  //num_vars=task->get_num_variables();
}
  void PatterCollectionEvaluatorRandWalk_Avg_H::initialize(std::shared_ptr<AbstractTask> task) {
    int num_vars= task->get_num_variables();
    cout<<"num_vars:"<<num_vars<<flush<<endl;
    TaskProxy task_proxy_temp(*task);
    task_proxy=make_shared<TaskProxy>(task_proxy_temp);
    successor_generator=utils::make_unique_ptr<SuccessorGenerator>(task);
    //result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
//    cout << "Manual pattern collection: " << *patterns << endl;
//    return PatternCollectionInformation(task, patterns);
  }
  bool PatterCollectionEvaluatorRandWalk_Avg_H::evaluate(std::shared_ptr<ModularZeroOnePDBs> candidate_PC){
    //cout<<"candidate_pc.size:"<<candidate_PC->get_size()<<endl;
    increased_states=0;
    double sum_h=0;
    size_t additions=0;
    for(auto state_pair : samples){
      if(candidate_PC->get_value(state_pair.first)>state_pair.second){
        DEBUG_MSG(cout<<"\th improved from "<<state_pair.second<<" to "<<candidate_PC->get_value(state_pair.first)<<endl;);
        increased_states++;
        if (candidate_PC->get_value(state_pair.first)!=numeric_limits<int>::max()){
	  sum_h+=candidate_PC->get_value(state_pair.first);
	  additions++;
	}
	else{//Found new dead_end, adding 1 for reflecting beneficial_impact
	  sum_h+=state_pair.second+1.0;
	  additions++;
	}
      }
      else{
	sum_h+=state_pair.second;
	additions++;
      }
    }
    set_eval_score(sum_h/double(additions));
    if(get_eval_score()>get_sample_score()){
      cout<<"time:"<<utils::g_timer()<<",Improving PC,increased_states:"<<increased_states<<",eval_score:"<<get_eval_score()<<",sample_score"<<sample_score<<endl;
      /*for(std::map<size_t,std::pair<State,int> >::iterator it=unique_samples.begin(); it!=unique_samples.end(); ++it){
        it->second.second=max(it->second.second,candidate_PC->get_value(it->second.first));
      }*/

      //Updating both unique_samples, used for pruning dominated heuristics and samples, and samples, used to evaluate heuristics
      //We should use samples to get a current view of the search, and unique_samples only to determine quasy-domination
      //using unique_samples to evaluate a heuristic can use too many nodes and distorts the value of the threshold
      //as the number of unique_samples keep growing
      return true;//Add collection
    }
    //else if(increased_states>0)
      //cout<<"time:"<<utils::g_timer()<<",Not_Improving PC,increased_states:"<<increased_states<<",sample_score:"<<get_sample_score()<<",eval_score:"<<get_eval_score()<<",additions:"<<additions<<"samples:,"<<samples.size()<<endl;
      if(additions==0){
          cerr<<"No additions for calculating avg_h_val, debug me!!!"<<endl;
	  exit(1);
      }
    //
    return false;//Not adding collection
  }

  void PatterCollectionEvaluatorRandWalk_Avg_H::sample_states(std::shared_ptr<PatternCollectionInformation> current_result){
    double sum_h=0;
    size_t additions=0;
    samples.clear();//We only use as samples the states, ordered in the open they would be fetched, in the current open list
    //Need to keep pointer to result or sample_states... function will complain current_result is not captured
    result=current_result;
    //evaluator_timer = new utils::CountdownTimer(time_limit);
    float start_time=utils::g_timer();
    //DEBUG_MSG(cout<<"adding to samples, unique_size prior:"<<unique_samples.size()<<",sampled_states:"<<samples.size()<<endl;);


// IF STATES_LOADED_FROM_OPEN_LIST IS EMPTY IT MEANS THIS IS FIRST CALL, REAL SEARCH HAS NOT STARTED
// FOR INITIALIZATION PURPOSES WE DO SAME RANDOM WALK AS IN RandWalk class

      evaluator_timer = new utils::CountdownTimer(time_limit);
      const State &initial_state = task_proxy->get_initial_state();
      int init_h=current_result->get_value(initial_state);
      double average_operator_cost=get_average_operator_cost(*task_proxy);
      cout<<"calling sample_states_with_random_walks"<<flush<<endl;
      vector<State> samples_just_states;
      try {
	  samples_just_states = sample_states_with_random_walks(
	      *task_proxy, *successor_generator, get_num_samples(), init_h,
	      average_operator_cost,
	      [this](const State &state) {
		  return result->is_dead_end(state);
	      },
	      evaluator_timer);
      } catch (SamplingTimeout &) {
	cout<<"We are finished,sampling_timeout in random_walk,random_walk_time:"<<utils::g_timer()-start_time<<endl;
      }

      for (auto state : samples_just_states){
	int h_val=current_result->get_value(state);
        if (h_val==numeric_limits<int>::max())//Skipping dead_ends when doing avg_h
	  continue;
	additions++;
	sum_h+=h_val;
	auto ret=unique_samples.insert(make_pair(state.hash(),make_pair(state,h_val)));
	if(ret.second==false){//Node already exist so update h_value
	  ret.first->second.second=max(h_val,ret.first->second.second);
	}
	samples.push_back(make_pair(state,h_val));
      }
      set_num_samples(samples.size());
      set_sample_score(sum_h/double(additions));
      //cout<<"Sampled_score:"<<get_sample_score()<<endl;
      cout<<"adding to samples using RAND_WALK, unique_size prior:"<<unique_samples.size()<<flush<<",num_samples:,"<<get_num_samples()<<",sample_score:,"<<get_sample_score()<<endl;
      return;
  }



  static shared_ptr<PatternCollectionEvaluator>_parse(options::OptionParser &parser) {
    parser.add_option<int> ("time_limit", "If populated,stop construction on first node past boundary and time limit", "100");
    options::Options options = parser.parse();
    parser.document_synopsis(
        "Pattern Generator RBP",
        "RBP-stype selection of variables to generate Pattern Collection");
    options::Options opts = parser.parse();
    if (parser.dry_run())
        return 0;

    return make_shared<PatterCollectionEvaluatorRandWalk_Avg_H>(opts);
  }

  static options::PluginShared<PatternCollectionEvaluator> _plugin("rand_walk_evaluator_avg_h", _parse);
}
