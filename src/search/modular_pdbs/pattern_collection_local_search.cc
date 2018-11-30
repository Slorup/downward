#include "pattern_collection_local_search.h"

#include "../option_parser.h"
#include "../plugin.h"
//#include "../task_proxy.h"

#include "../task_tools.h"
#include "../causal_graph.h"
#include "../utils/collections.h"


    
using namespace std;
namespace pdbs3 {
//    virtual void initialize(std::shared_ptr<AbstractTask> task) {
//    cout << "Manual pattern collection: " << *patterns << endl;
//    return PatternCollectionInformation(task, patterns);
//
//   }

static options::PluginTypePlugin<PatternCollectionLocalSearch> _type_plugin(
"PatterCollectionLocalSearch",
"The various local search algorithms for pattern(s) improvement usable by the complementary_modular_pdbs heurisitic.");
}
  
namespace pdbs3 {
  
  void PatternCollectionLocalSearch::initialize(std::shared_ptr<AbstractTask> task) {
    num_vars= task->get_num_variables();
    cout<<"ini,num_vars:"<<num_vars<<endl;
    task_proxy=make_shared<TaskProxy>(*task);
    for (FactProxy goal : task_proxy->get_goals()) {
      initial_remaining_goal_vars.insert(goal.get_variable().get_id());
    }
  }
    
  void PatternCollectionLocalSearch::remove_irrelevant_variables(Pattern &pattern,Pattern &removed_vars){
	set<int> in_original_pattern(pattern.begin(), pattern.end());
	set<int> in_pruned_pattern;
	removed_vars.clear();
	
	//cout<<"in_original_pattern:,";for (auto i : in_original_pattern) cout<<i<<",";cout<<endl;
	//cout<<"initial_remaining_goal_vars:,";for (auto i : initial_remaining_goal_vars) cout<<i<<",";cout<<endl;
	
	vector<int> vars_to_check;
	for ( auto var_id : initial_remaining_goal_vars) {
	  if (in_original_pattern.count(var_id)) {
	    vars_to_check.push_back(var_id);
	    in_pruned_pattern.insert(var_id);
	  }
	}
	//cout<<"vars_to_check:";for (auto i : vars_to_check) cout<<i<<",";cout<<flush<<endl;
	    
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
	//cout<<"in_original_pattern:";for(auto element : in_original_pattern) cout<<element<<",";cout<<endl;
	//cout<<"in_prunned_pattern:";for(auto element: in_pruned_pattern) cout <<element<<",";cout<<endl; 
	//cout<<"removed_vars:";for(auto element: removed_vars) cout <<element<<",";cout<<endl; 
  }
}
