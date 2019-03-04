//#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_generator_GamerStyle.h"

#include "../causal_graph.h"
#include "../globals.h"
//#include "../task_proxy.h"

//#include "../utils/markup.h"
//#include "../utils/math.h"
//#include "../utils/rng.h"
#include "../utils/timer.h"
//#include "../heuristic.h"

//#include <algorithm>
//#include <cassert>
//#include <iostream>
//#include <unordered_set>
#include <vector>
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
//RandomCollectionGeneration: CBP, Gamer, CGamer, the one used by *Pommerening et al. 
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
PatternCollectionGeneratorGamer::PatternCollectionGeneratorGamer(const options::Options & opts) :
	time_limit (opts.get<int>("time_limit")){
    cout<<"hello GamerStyle pattern selection!"<<endl;
    //num_vars=task->get_num_variables();
}
std::ostream & operator<<(std::ostream &os, set<int> pattern){
    for (int v : pattern) os << "," << v;  
    return os;
}
std::ostream & operator<<(std::ostream &os, vector<int> pattern){
    for (int v : pattern) os << "," << v;  
    return os;
}
  PatternCollectionContainer PatternCollectionGeneratorGamer::generate(){
    cout<<"Calling PatternCollectionGeneratorGamer, it just creates a PC with one pattern with all goals, from there it is up to local search algorithms how to update pattern."<<endl;
    pattern.clear();
    for (auto goal : g_goal) pattern.push_back(goal.first);
    cout<<"Gamer, initial pattern made of all goals:"<<pattern<<endl;
    PatternCollectionContainer PC;
    PC.restart_pc(pattern);
    return PC;
  }

  static shared_ptr<PatternCollectionGeneratorComplementary>_parse(options::OptionParser &parser) {
    parser.add_option<int> ("time_limit", "If populated,stop construction on first node past boundary and time limit", "100");
    options::Options options = parser.parse();
    parser.document_synopsis(
        "C-Gamer style pattern generation","");
    options::Options opts = parser.parse();
    if (parser.dry_run())
        return 0;

    return make_shared<PatternCollectionGeneratorGamer>(opts);
  }

  static options::PluginShared<PatternCollectionGeneratorComplementary> _plugin("gamer", _parse);
}
