#include "sym_solution.h"

#include <vector>       // std::vector
#include "../state_registry.h"

#include "sym_search.h"


using namespace std;

namespace symbolic {
void SymSolution::getPlan(vector <const GlobalOperator *> &path) const {
    if (path.empty()) { //This code should be modified to allow appending things to paths
        exp->getPlan(cut, g, h, path);
    }
}

ADD SymSolution::getADD() const {
    vector <const GlobalOperator *> path;
    exp->getPlan(cut, g, h, path);

    SymVariables *vars = exp->getStateSpace()->getVars();
    ADD hADD = vars->getADD(-1);
    int h_val = g + h;

    vector<int> s = g_initial_state_data;
    BDD sBDD = vars->getStateBDD(s);
    hADD += sBDD.Add() * (vars->getADD(h_val + 1));
    for (auto op : path) {
        h_val -= op->get_cost();
        for (const GlobalEffect &eff : op->get_effects()) {
            if (eff.does_fire(s)) {
                s[eff.var] = eff.val;
            }
        }
        sBDD = vars->getStateBDD(s);
        hADD += sBDD.Add() * (vars->getADD(h_val + 1));
    }
    return hADD;
}
}
