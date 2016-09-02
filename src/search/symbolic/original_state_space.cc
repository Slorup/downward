#include "original_state_space.h"

#include "../utils/debug_macros.h"
using namespace std;

namespace symbolic {
OriginalStateSpace::OriginalStateSpace(SymVariables *v,
                                       const SymParamsMgr &params,
                                       OperatorCost cost_type_) :
    SymStateSpaceManager(v, params, cost_type_) {
    for (size_t i = 0; i < g_variable_domain.size(); ++i) {
        fullVars.insert(i);
    }
}

void OriginalStateSpace::init_initial_state() {
    initialState = vars->getStateBDD(g_initial_state_data);
}

void OriginalStateSpace::init_goal() {
    goal = vars->getPartialStateBDD(g_goal);
}

void OriginalStateSpace::init_individual_trs() {
    if (!indTRs.empty())
        return;

    DEBUG_MSG(cout << "Initialize individual TRs of original state space" << endl;
              );
    for (size_t i = 0; i < g_operators.size(); i++) {
        const GlobalOperator *op = &(g_operators[i]);
        // Skip irrelevant operators
        /*if (op->is_dead()){
          continue;
          }*/
        int cost = get_adjusted_action_cost(*op, cost_type);
        DEBUG_MSG(cout << "Creating TR of op " << i << " of cost " << cost << endl;
                  );
        indTRs[cost].push_back(move(SymTransition(vars, op, cost)));
        if (p.mutex_type == MutexType::MUTEX_EDELETION) {
            indTRs[cost].back().edeletion(*this);
        }
    }
}
}
