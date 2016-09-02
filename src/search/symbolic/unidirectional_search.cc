#include "unidirectional_search.h"

using namespace std;

namespace symbolic {
    UnidirectionalSearch::UnidirectionalSearch(const SymParamsSearch &params) : 
	SymSearch(params), fw(true) {}



void UnidirectionalSearch::statistics() const {
    cout << "Exp " << (fw ? "fw" : "bw") << " time: " << stats.step_time << "s (img:" <<
        stats.image_time << "s, heur: " << stats.time_heuristic_evaluation <<
        "s) in " << stats.num_steps_succeeded << " steps ";
}
}
