#include "pattern_database_symbolic.h"

#include "../symbolic/uniform_cost_search.h"
#include "../utils/timer.h"

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
	vars (vars_), 
	manager (manager_), average(0) {
	
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
	       search.isSearchable() 
	       /*&& (!solution.solved() || originalsearch.getF() < solution.getCost())*/) {
	    search.step();
	} 
	
	average = search.getClosed()->average_hvalue();
	//cout << "Average: " << average << endl;
	heuristic = make_unique<ADD>(search.getHeuristic());
    }

int PatternDatabaseSymbolic::get_value(const State & state) const {
    return get_value (vars->getBinaryDescription(state.get_values()));
}

int PatternDatabaseSymbolic::get_value(const vector<int> &state) const {
    return get_value (vars->getBinaryDescription(state));
}

int PatternDatabaseSymbolic::get_value(int * inputs) const {
    // for(const BDD & bdd : notMutexBDDs){
    // 	if(bdd.Eval(inputs).IsZero()){
    // 	    return DEAD_END;
    // 	}
    // }

    if(!heuristic) return 0;

    ADD evalNode = heuristic->Eval(inputs);
    int abs_cost = Cudd_V(evalNode.getRegularNode());
    
    return (abs_cost == -1 ? -1 : abs_cost);    
}

double PatternDatabaseSymbolic::compute_mean_finite_h() const {
    
    return average;
}

}
