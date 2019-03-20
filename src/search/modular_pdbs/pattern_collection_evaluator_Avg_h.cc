//#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_evaluator_Avg_h.h"
#include "../sampling.h"
#include "../utils/timer.h"
#include "../task_tools.h"
#include "../utils/countdown_timer.h"
#include "../options/option_parser.h"
#include "../options/plugin.h"
#include "pattern_database.h"

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
PatterCollectionEvaluator_Avg_H::PatterCollectionEvaluator_Avg_H(const options::Options & opts) :
	time_limit (opts.get<int>("time_limit")){
    cout<<"hello Evaluator_Avg_H"<<flush<<endl;
}
  void PatterCollectionEvaluator_Avg_H::initialize(std::shared_ptr<AbstractTask> task) {
    int num_vars= task->get_num_variables();
    cout<<"num_vars:"<<num_vars<<flush<<endl;
    TaskProxy task_proxy_temp(*task);
    task_proxy=make_shared<TaskProxy>(task_proxy_temp);
    successor_generator=utils::make_unique_ptr<SuccessorGenerator>(task);
  }
  bool PatterCollectionEvaluator_Avg_H::evaluate(std::shared_ptr<ModularZeroOnePDBs> candidate_PC){
    set_eval_score(candidate_PC->compute_approx_mean_finite_h());
    if(get_eval_score()>get_sample_score()){
      cout<<"time:"<<utils::g_timer()<<",Improving PC,"<<",eval_score:,"<<get_eval_score()<<",sample_score:,"<<sample_score<<endl;
      return true;//Add collection
    }
    return false;//Not adding collection
  }

  void PatterCollectionEvaluator_Avg_H::sample_states(std::shared_ptr<PatternCollectionInformation> current_result){
    double max_avg_h=0;
    
    std::shared_ptr<MaxAdditivePDBSubsets> current_subsets=current_result->get_max_additive_subsets();
    for(size_t subset=0;subset<current_subsets->size();subset++){
      double temp_avg_h=0;
      for(auto pdb : current_subsets->at(subset)){
	  temp_avg_h+=pdb->compute_mean_finite_h();
      }
      max_avg_h=max(max_avg_h,temp_avg_h);
    }

    set_sample_score(max_avg_h);
    set_num_samples(0);
    set_threshold(0);
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

    return make_shared<PatterCollectionEvaluator_Avg_H>(opts);
  }

  static options::PluginShared<PatternCollectionEvaluator> _plugin("evaluator_avg_h", _parse);
}
