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
  //void PatternCollectionGeneratorGamer::_initialize(std::shared_ptr<AbstractTask> task) 
  //void PatternCollectionGeneratorGamer::initialize(std::shared_ptr<AbstractTask> task) 
  void PatternCollectionGeneratorGamer::initialize() {
    //1) Get initial abstraction
    for (auto goal : g_goal) pattern.insert(goal.first);
    cout<<"Gamer, initial pattern made of all goals:"<<pattern<<endl;
    
    candidate_pattern=candidate_vars();
    cout<<"initial candidate_pattern for Gamer-Style selection, candidate_vars:"<<candidate_pattern.size()<<endl;
  }

  PatternCollectionContainer PatternCollectionGeneratorGamer::generate(){
    cout<<"Calling PatternCollectionGeneratorGamer"<<endl;
    PatternCollectionContainer PC;
    if(candidate_pattern.empty()){
      cout<<"No more possible candidates for Gamer-Style Selection"<<endl;
      return PC;
    }
    else{
      cout<<"time:"<<utils::g_timer()<<",pattern:"<<pattern<<",candidate list:"<<candidate_pattern<<endl;
      while(!candidate_pattern.empty()){
        child_pattern=pattern;
        child_pattern.insert(candidate_pattern.back());
        candidate_pattern.pop_back();
        
        vector<int> child_pattern_vect;
        std::copy(child_pattern.begin(), child_pattern.end(), std::back_inserter(child_pattern_vect));
        PC.add_pc(child_pattern_vect);
      }
      return PC;
    }
  }
  PatternCollectionContainer PatternCollectionGeneratorGamer::get_PC(){//for calculating avg_h_val
    PatternCollectionContainer PC;
    vector<int> pattern_vect;
    std::copy(pattern.begin(), pattern.end(), std::back_inserter(pattern_vect));
    PC.add_pc(pattern_vect);
      return PC;
    }




std::vector<int> PatternCollectionGeneratorGamer::candidate_vars() const {
    TaskProxy task_proxy(*(g_root_task()));
    const CausalGraph &cg = task_proxy.get_causal_graph();
    vector<int> candidates;
    for (size_t var = 0; var < g_variable_domain.size(); ++var) {
        if (pattern.count(var)){ 
          //cout<<"\t\tGamer,skipping exisiting var:"<<var<<endl;
          continue;
        }

      for (int succ : cg.get_pre_to_eff(var)) {
        //cout<<"\t\tGamer,checking connected variables to var:"<<var<<"trying succ:"<<succ<<endl;
        if (pattern.count(succ)) {
          //cout<<"\t\tGamer,connected variables:"<<succ<<"to var:"<<var<<" added."<<endl;
          candidates.push_back(var); 
          break;
        }
      }
    }
    return candidates;
}
void PatternCollectionGeneratorGamer::confirm(set<int> input_pattern){
  //Update pattern with new improved pattern and also generate a new list of candidates
  pattern=input_pattern;
  candidate_pattern=candidate_vars();
  
  cout<<"updated Gamer_pattern:"<<pattern<<",new Gamer-style candidate vars:"<<child_pattern<<endl;
}
void PatternCollectionGeneratorGamer::check_improv(double input_avg_h=0){
  if(avg_h_val==0){
    cout<<"First call to check_improv,initial avg_h:"<<input_avg_h<<endl;
    return;
  }
  else if(input_avg_h>avg_h_val){
    cout<<"Gamer-style,avg_h_val raised to:,"<<input_avg_h<<",from:,"<<avg_h_val<<",confirming new pattern"<<endl;
    avg_h_val=input_avg_h;
  }
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

  static options::PluginShared<PatternCollectionGeneratorComplementary> _plugin("2RandomPatterns", _parse);
}
