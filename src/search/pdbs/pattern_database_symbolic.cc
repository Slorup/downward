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
						     double generationTime, 
						     double generationMemoryGB) : 
	PatternDatabaseInterface(task_proxy, pattern, operator_costs), 
	vars (vars_), manager (manager_), heuristic(vars->getADD(0)), dead_ends(vars->zeroBDD()), 
	finished(false), hvalue_unseen_states(0), average(0) {
	  DEBUG_MSG(cout<<"start_time_pdb_constructor:"<<utils::g_timer()<<",";);
	  generationMemoryGB=(memory_limit/1024.0)-(utils::get_peak_memory_in_kb()/(1024.0*1024.0));
	  if(generationMemoryGB<=0){
	    cout<<"No more pdb gen, reached Memory limit!!!"<<endl;
	    return;
	  }
	  DEBUG_MSG(cout<<"generationMemoryGB for symbolic:"<<generationMemoryGB<<endl;);
	
	create_pdb(engine,params, generationTime, generationMemoryGB);
    }


    void PatternDatabaseSymbolic::create_pdb(SymController * engine, const SymParamsSearch & params, 
					     double generationTime, double generationMemoryGB) {
	//float start_time=utils::g_timer();
	//cout<<"start_time_create_pdb:"<<utils::g_timer()<<",";
	symbolic::UniformCostSearch search (engine, params);
	//cout<<"UniformCostSearch.time:"<<utils::g_timer()-start_time<<",";
	search.init(manager, false);
	//cout<<"serach.init.time:"<<utils::g_timer()-start_time<<",";


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
	
	DEBUG_MSG(cout << "Solved: " << engine->solved() << " Finished: " << search.finished() <<  ", Average: " << average << endl;);
	if(time()>=generationTime){
	  cout<<"generationTimeLimit:"<<generationTime<<">GenTime:"<<time()<<",symbolic pdb interrupted"<<endl;
	}
	if(vars->totalMemoryGB() >= generationMemoryGB){
	  cout<<"vars->totalMemoryGB():"<<vars->totalMemoryGB()<<">GenMemoryLimit(GB):"<<generationMemoryGB<<endl;
	}

	if(engine->solved()) {
	    heuristic = engine->get_solution()->getADD();	    
	} else {
	  //cout<<"time before serch.getHeuristic(false):"<<time()<<endl;
	    heuristic = search.getHeuristic(false);
	  //cout<<"time after serch.getHeuristic(false):"<<time()<<endl;
	    if(finished) dead_ends += search.notClosed(); 
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

}

