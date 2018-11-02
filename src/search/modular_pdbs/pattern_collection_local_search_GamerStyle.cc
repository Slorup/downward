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

#include "pdb_factory.h"
#include "../utils/debug_macros.h"

    
using namespace std;
namespace pdbs3 {
  PatternCollectionLocalSearchGamerStyle::PatternCollectionLocalSearchGamerStyle(const options::Options & opts)
  :time_limit (opts.get<int>("time_limit")){
  //PatternCollectionLocalSearchGamerStyle::PatternCollectionLocalSearchGamerStyle() 
      cout<<"hello LocalSearchGamerStyle"<<endl;
    //num_vars=task->get_num_variables();
  }

  PatternCollectionContainer PatternCollectionLocalSearchGamerStyle::generate_next_candidate(PatternCollectionContainer candidate_collection){
    //If more than one pattern, we only do local search on the first pattern
    Pattern candidate_pattern=candidate_collection.get_top_pattern();
    std::set<int> pattern;
    for (auto var : candidate_pattern) pattern.insert(var);


    TaskProxy task_proxy(*(g_root_task()));
    const CausalGraph &cg = task_proxy.get_causal_graph();
    vector<int> candidates;
    for (size_t var = 0; var < g_variable_domain.size(); ++var) {
        if (pattern.count(var)){ 
          //cout<<"\t\tGamer_local_search,skipping exisiting var:"<<var<<endl;
          continue;
        }
	else if(forbidden_vars.count(var)){//This pattern was proven to not improve search
	  continue;
	}

      for (int succ : cg.get_pre_to_eff(var)) {
        if (pattern.count(succ)) {
          //cout<<"\t\tGamer,connected variables:"<<succ<<"to var:"<<var<<" added."<<endl;
          candidates.push_back(var); 
          break;
        }
      }
    }
    cout<<"\t forbidden_vars_size:"<<forbidden_vars.size()<<",candidate_vars_size:"<<candidates.size()<<endl;
    if(candidates.size()==0){
      new_patterns.restart_pc(candidate_pattern);
      return new_patterns;
    }
    cout<<"input_pattern:";for (auto i : candidate_pattern) cout<<i<<",";
    cout<<"candidates:";for (auto i : candidates) cout<<i<<",";cout<<endl;
    assert(candidates.size()>0);
    last_var=candidates.at(rand()%candidates.size());
    //cout<<"adding random_var:"<<last_var<<endl;
    candidate_pattern.push_back(last_var);
    //All patterns have to be sorted! for duplication checks
    //and any other comparisons of patterns which are vectors,
    //so order matters!
    sort(candidate_pattern.begin(), candidate_pattern.end());
    new_patterns.restart_pc(candidate_pattern);
    cout<<"new_patterns:";new_patterns.print();
    return new_patterns;//Not adding collection
  }

  void PatternCollectionLocalSearchGamerStyle::forbid_last_var(){
    forbidden_vars.insert(last_var);
  }
  void PatternCollectionLocalSearchGamerStyle::reset_forbidden_vars(){
    forbidden_vars.clear();
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
