#include "pattern_database_symbolic.h"

#include "../symbolic/uniform_cost_search.h"
#include "../symbolic/sym_controller.h"

#include "../utils/timer.h"
#include "../utils/debug_macros.h"

#include "../symbolic/sym_solution.h"



using namespace std;

using namespace symbolic;

using utils::Timer;

namespace pdbs {
    
    PatternDatabaseSymbolic::PatternDatabaseSymbolic(const TaskProxy &task_proxy, 
						     const Pattern &pattern, 
						     const std::vector<int> &operator_costs,
						     SymController * engine,
						     std::shared_ptr<SymVariables> vars_, 
						     std::shared_ptr<SymStateSpaceManager> manager_, 
						     const SymParamsSearch & params, 
						     int max_time_ms, int max_nodes, int global_limit_memory_MB) : 
	PatternDatabaseInterface(task_proxy, pattern, operator_costs), 
	vars (vars_), manager (manager_), heuristic(vars->getADD(0)), dead_ends(vars->zeroBDD()), 
	finished(false), hvalue_unseen_states(0), average(0) {

	create_pdb(engine,params, max_time_ms, max_nodes, global_limit_memory_MB);
    }


    void PatternDatabaseSymbolic::create_pdb(SymController * engine, const SymParamsSearch & params, 
					     int max_time_ms, int max_nodes, int global_limit_memory_MB) {
	//float start_time=utils::g_timer();
	//cout<<"start_time_create_pdb:"<<utils::g_timer()<<",";
	search= make_unique<symbolic::UniformCostSearch> (engine, params);
	search->set_limits(max_time_ms, max_nodes);
	//cout<<"UniformCostSearch.time:"<<utils::g_timer()-start_time<<",";
	search->init(manager, false);
	//cout<<"serach.init.time:"<<utils::g_timer()-start_time<<",";

	Timer time; 
	while (!search->finished() && 
	       time()*1000.0 < (double)max_time_ms &&
	       vars->totalMemoryGB()*1024 < global_limit_memory_MB &&
	       search->isSearchable()  && 
	       !engine->solved()) {
	    search->step();
	} 
	
	finished = search->finished();
	hvalue_unseen_states = search->getHNotClosed();
	average = search->getClosed()->average_hvalue();
	DEBUG_MSG(for (int v : pattern) cout << v << " ";);
	
	DEBUG_MSG(cout << "Solved: " << engine->solved() << " Finished: " << search->finished() <<  ", Average: " << average << endl;);

	if(engine->solved()) {
	    heuristic = engine->get_solution()->getADD();
	    search.reset();
	} else {
	    //cout<<"time before serch.getHeuristic(false):"<<time()<<endl;
	    heuristic = search->getHeuristic(false);
	    //cout<<"time after serch.getHeuristic(false):"<<time()<<endl;
	    if(finished) { 
		dead_ends += search->notClosed(); 
		search.reset();
	    }
	}
	//cout<<"Overall generationTime:,"<<utils::g_timer()-start_time<<endl;
    }

    int PatternDatabaseSymbolic::get_value(const State & state) const {
	return get_value (vars->getBinaryDescription(state.get_values()));
    }

    int PatternDatabaseSymbolic::get_value(const vector<int> &state) const {
	return get_value (vars->getBinaryDescription(state));
    }

    int PatternDatabaseSymbolic::get_value(int * inputs) const {
	if(!dead_ends.Eval(inputs).IsZero()){
	    return numeric_limits<int>::max();
	}

	ADD evalNode = heuristic.Eval(inputs);
	int abs_cost = Cudd_V(evalNode.getRegularNode());

	return (abs_cost == -1 ? numeric_limits<int>::max() : abs_cost);    
    }

    int PatternDatabaseSymbolic::get_goal_cost(const vector<int> & state_pattern, const State & state) const {
	assert(std::includes(pattern.begin(), pattern.end(), state_pattern.begin(), state_pattern.end()));
	auto bin = vars->getBinaryDescription(state_pattern, state.get_values());
	int value = get_value (bin);
	if(is_finished() || value < hvalue_unseen_states) {
	    return value;
	}else {
	    return -1;
	}
    }

    double PatternDatabaseSymbolic::compute_mean_finite_h() const {
	return average;
    }


    void PatternDatabaseSymbolic::terminate_creation (int max_time_ms, int max_nodes, 
						      int global_limit_memory_MB) {

	if(!search) {
	    return;
	}
	search->set_limits(max_time_ms, max_nodes);

	Timer time; 
	while (!search->finished() && 
	       time()*1000 < (double)max_time_ms &&
	       vars->totalMemoryGB()*1024 < global_limit_memory_MB &&
	       search->isSearchable()  && 
	       !search->getEngine()->solved()) {
	    search->step();
	} 
	
	finished = search->finished();
	hvalue_unseen_states = search->getHNotClosed();
	average = search->getClosed()->average_hvalue();

	if(search->getEngine()->solved()) {
	    heuristic = search->getEngine()->get_solution()->getADD();	    
	} else {
	    //cout<<"time before serch.getHeuristic(false):"<<time()<<endl;
	    heuristic = search->getHeuristic(false);
	    //cout<<"time after serch.getHeuristic(false):"<<time()<<endl;
	    if(finished) dead_ends += search->notClosed(); 
	}
    }
}

