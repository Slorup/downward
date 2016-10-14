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
						     int generationTime, 
						     double generationMemoryGB) : 
	PatternDatabaseInterface(task_proxy, pattern, operator_costs), 
	vars (vars_), manager (manager_), heuristic(vars->getADD(0)), dead_ends(vars->zeroBDD()), 
	finished(false), hvalue_unseen_states(0), average(0) {
	
	create_pdb(engine,params, generationTime, generationMemoryGB);
    }


    void PatternDatabaseSymbolic::create_pdb(SymController * engine, const SymParamsSearch & params, 
					     int generationTime, double generationMemoryGB) {
  	manager->init();
	symbolic::UniformCostSearch search (engine, params);
	search.init(manager, false);

	Timer time; 
	while (!search.finished() && 
	       time() < generationTime &&
	       vars->totalMemoryGB() < generationMemoryGB &&
	       search.isSearchable()  && 
	       !engine->solved()) {
	    search.step();
	} 
	
	finished = search.finished();
	hvalue_unseen_states = search.getHNotClosed();
	average = search.getClosed()->average_hvalue();
	DEBUG_MSG(for (int v : pattern) cout << v << " ";);
	
	//cout << " Finished: " << search.finished() <<  ", Average: " << average << endl;
	if(time()>generationTime){
	  cout<<"generationTimeLimit:"<<generationTime<<">GenTime:"<<time()<<",symbolic pdb interrupted"<<endl;
	}
	if(engine->solved()) {
	    heuristic = engine->get_solution()->getADD();	    
	} else {
	    heuristic = search.getHeuristic(false);
	    if(finished) dead_ends += search.notClosed(); 
	}
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

    double PatternDatabaseSymbolic::compute_mean_finite_h() const {
	return average;
    }

}
