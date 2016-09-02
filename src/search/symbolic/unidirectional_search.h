#ifndef SYMBOLIC_UNIDIRECTIONAL_SEARCH_H
#define SYMBOLIC_UNIDIRECTIONAL_SEARCH_H

#include "sym_search.h"

#include "sym_estimate.h"
#include "sym_util.h"
#include <vector>
#include <map>
#include <memory>

namespace symbolic {

class SymExpStatistics {
public:
    double image_time, image_time_failed;
    double time_heuristic_evaluation;
    int num_steps_succeeded;
    double step_time;

    SymExpStatistics() :
        image_time(0),
        image_time_failed(0), time_heuristic_evaluation(0),
        num_steps_succeeded(0), step_time(0) {  }


    void add_image_time(double t) {
        image_time += t;
        num_steps_succeeded += 1;
    }

    void add_image_time_failed(double t) {
        image_time += t;
        image_time_failed += t;
        num_steps_succeeded += 1;
    }
};

class UnidirectionalSearch : public SymSearch {
protected:
    bool fw; //Direction of the search. true=forward, false=backward

    SymExpStatistics stats;

public:
    UnidirectionalSearch (const SymParamsSearch &params);

    inline bool isFW() const {
        return fw;
    }

    virtual void getHeuristic(std::vector<ADD> &heuristics,
			      std::vector <int> &maxHeuristicValues) const = 0;

    void statistics() const;
};
}
#endif // SYMBOLIC_EXPLORATION
