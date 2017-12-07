//#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_generator_RBP.h"

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
//#include "../successor_generator.h"
//#include "../utils/countdown_timer.h"
//#include "pdb_factory.h"
//#include "pattern_database_interface.h"
//#include "../utils/debug_macros.h"
//#include <random>
//#include "../sampling.h"
#include "../task_tools.h"
#include "../options/option_parser.h"
#include "../options/plugin.h"

    
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
PatternCollectionGeneratorRBP::PatternCollectionGeneratorRBP(const options::Options & opts) :
  num_vars(0),
	time_limit (opts.get<int>("time_limit")){
    cout<<"hello modular_RBP"<<endl;
  //num_vars=task->get_num_variables();
}
  void PatternCollectionGeneratorRBP::initialize(std::shared_ptr<AbstractTask> task) {
    num_vars= task->get_num_variables();
//    cout << "Manual pattern collection: " << *patterns << endl;
//    return PatternCollectionInformation(task, patterns);
  }
  PatternCollectionContainer PatternCollectionGeneratorRBP::generate(){
      PatternCollectionContainer PC;
      Pattern temp_pattern;
      unsigned vars_to_choose=rand()%num_vars;
      for (unsigned i=0;i< vars_to_choose;i++){
        temp_pattern.push_back(i);
      }
      PC.add_pc(temp_pattern);
      return PC;
    }

  static shared_ptr<PatternCollectionGeneratorComplementary>_parse(options::OptionParser &parser) {
    parser.add_option<int> ("time_limit", "If populated,stop construction on first node past boundary and time limit", "100");
    options::Options options = parser.parse();
    parser.document_synopsis(
        "Pattern Generator RBP",
        "RBP-stype selection of variables to generate Pattern Collection");
    options::Options opts = parser.parse();
    if (parser.dry_run())
        return 0;

    return make_shared<PatternCollectionGeneratorRBP>(opts);
  }

  static options::PluginShared<PatternCollectionGeneratorComplementary> _plugin("modular_rbp", _parse);
}
