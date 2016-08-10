#include "pattern_database_interface.h"

#include "../task_tools.h" 

#include <algorithm>


using namespace std;

namespace pdbs {

PatternDatabaseInterface::PatternDatabaseInterface(
    const TaskProxy &task_proxy,
    const Pattern &pattern,
    const vector<int> &operator_costs_)
    : task_proxy(task_proxy),
      pattern(pattern), 
      operator_costs(operator_costs_) {
    verify_no_axioms(task_proxy);
    assert(operator_costs.empty() ||
           operator_costs.size() == task_proxy.get_operators().size());
    assert(utils::is_sorted_unique(pattern));
}


bool PatternDatabaseInterface::is_operator_relevant(const OperatorProxy &op) const {
    for (EffectProxy effect : op.get_effects()) {
        int var_id = effect.get_fact().get_variable().get_id();
        if (binary_search(pattern.begin(), pattern.end(), var_id)) {
            return true;
        }
    }
    return false;
}

}
