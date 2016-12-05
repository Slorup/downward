#include "pattern_database_online_plus.h"
#include "pattern_database_symbolic.h"

using namespace std;

#include "../priority_queue.h"

#include "../utils/debug_macros.h"

namespace pdbs {




    PatternDatabaseOnlinePlus::PatternDatabaseOnlinePlus(const TaskProxy &task, 
							 const Pattern &pattern,
							 const std::vector<int> &operator_costs,
							 std::shared_ptr<AbstractTask> pdb_task_,
							 std::vector<Heuristic *> heuristics_, 
							 symbolic::SymController * engine, 
							 std::shared_ptr<symbolic::SymVariables> vars, 
							 std::shared_ptr<symbolic::SymStateSpaceManager> manager, 
							 const symbolic::SymParamsSearch & params, 
							 double generationTime, double generationMemoryGB)
	: PatternDatabaseInterface(task, pattern, operator_costs), pdb_task(pdb_task_),
	  heuristics(heuristics_), task_proxy(task), 
	successor_generator(pdb_task),
	symbolic_pdb(make_unique<PatternDatabaseSymbolic>(task, pattern, operator_costs, 
							  engine, vars, manager, params, 
							  generationTime,
							  generationMemoryGB)){
	
    }

    int PatternDatabaseOnlinePlus::get_value(const State & initial_state) const {

	int goal_cost = get_goal_cost(initial_state);
	if(goal_cost >= 0) {
	    return goal_cost;
	}

	int initial_h = compute_heuristic(initial_state); 
	if (initial_h == std::numeric_limits<int>::max()) {
	    return initial_h;
	}
	
        // (first implicit entry: priority,) second entry: index for an abstract state
	AdaptiveQueue<LocalStateID> open_list; 

	SearchInfo search_info; 

	LocalStateID initial_id = search_info.get_id(initial_state);
	SearchStateInfo & initial_node = search_info.get_state_info(initial_id);
	initial_node.h = initial_h;
	initial_node.g = 0;
	DEBUG_MSG(cout << "Initial: " << initial_id << " with h=" << initial_h << endl;);
	cout << "Initial: " << initial_id << " with h=" << initial_h << endl;
	open_list.push(initial_h, initial_id);
	int upper_bound = std::numeric_limits<int>::max();
	while (!open_list.empty()) {	    
	    pair<int, size_t> node = open_list.pop();
	    if (node.first > upper_bound) {
		return upper_bound;
	    }
	    SearchStateInfo & node_info = search_info.get_state_info(node.second);

	    if(node_info.closed){
		continue;
	    }
	    node_info.closed = true;

	    int parent_g = node_info.g;
	    const State & state = search_info.get_state(node.second);
	    int goal_cost = get_goal_cost(state);
	    if (goal_cost >= 0) { //Goal cost could be determined
		if (goal_cost < upper_bound - parent_g) {
		    upper_bound = min(upper_bound, parent_g + goal_cost) ;
		}
		continue;
	    }

	    vector<OperatorProxy> applicable_ops;
	    successor_generator.generate_applicable_ops(state, applicable_ops);

	    for (const auto & op : applicable_ops) {
		int succ_g = parent_g + op.get_cost(); // get_adjusted_cost(op)
		if (succ_g >= upper_bound) {
		    continue;
		}

		State succ_state = state.get_successor(op);
		int goal_cost = get_goal_cost(succ_state);
		if (goal_cost >= 0) { //Goal cost could be determined
		    if (goal_cost < upper_bound - succ_g) {
			upper_bound = min(upper_bound, succ_g + goal_cost);
		    }
		    continue;
		}

		LocalStateID succ_id = search_info.get_id(succ_state); 
	        SearchStateInfo & succ_node = search_info.get_state_info(succ_id);
		if(succ_node.h ==  std::numeric_limits<int>::max()) {
		    continue;
		}

		if (succ_node.g == -1) { //new node
		    succ_node.g = succ_g;
		    succ_node.h = compute_heuristic(succ_state);
		    if (succ_node.h == std::numeric_limits<int>::max()) {
			continue;
		    }
		    open_list.push(succ_node.f(), succ_id);
		} else if (succ_node.g > succ_g) {
		    // We found a new cheapest path to an open or closed state.
		    succ_node.closed = false;
		    open_list.push(succ_node.f(), succ_id);
		}
	    }
	}


	cout << "Upper bound" << upper_bound << endl;
	return upper_bound;
    }


    int PatternDatabaseOnlinePlus::compute_heuristic(const State & /*state*/) const {
	return symbolic_pdb->get_hvalue_unseen_states();
    }

    int PatternDatabaseOnlinePlus::get_goal_cost(const State & state) const {
	return symbolic_pdb->get_goal_cost(state);
    }

    double PatternDatabaseOnlinePlus::compute_mean_finite_h() const {
	return symbolic_pdb->compute_mean_finite_h();
    }



}
