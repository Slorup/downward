#ifndef SYMBOLIC_SYM_PDB_H
#define SYMBOLIC_SYM_PDB_H

#include "sym_state_space_manager.h"
#include "sym_variables.h"
#include <set>

namespace symbolic {
class SymPDB : public SymStateSpaceManager {
    BDD nonRelVarsCube; //Cube BDD representing relevantVars
    BDD nonRelVarsCubeWithPrimes; //Cube BDD representing relevantVars
    std::string abstractionName;

    virtual void init_initial_state() override;
    virtual void init_goal() override;
    virtual void init_individual_trs() override;

public:
    SymPDB(SymVariables *bdd_vars, const SymParamsMgr &params, OperatorCost cost_type_); //Creates a BDD with all variables relevant
    SymPDB(SymVariables *bdd_vars, const SymParamsMgr &params, OperatorCost cost_type_,
           AbsTRsStrategy absTRsStrategy, const std::set<int> &relVars);

    SymPDB(std::shared_ptr<SymStateSpaceManager> &parent,
           AbsTRsStrategy absTRsStrategy, const std::set<int> &relVars);

    virtual ~SymPDB() {}
    virtual BDD shrinkExists(const BDD &bdd, int maxNodes) const override;
    virtual BDD shrinkForall(const BDD &bdd, int maxNodes) const override;
    virtual BDD shrinkTBDD(const BDD &tBDD, int maxNodes) const override;

    virtual ADD getExplicitHeuristicADD(bool fw) override;
    virtual void getExplicitHeuristicBDD(bool fw, std::map<int, BDD> &res) override;

    virtual std::string tag() const override;

    inline SymVariables *getVars() const {
        return vars;
    }

    virtual void print(std::ostream &os, bool fullInfo) const override;

    //virtual void init_mutex(const std::vector<MutexGroup> & mutex_groups) override;
};
}



#endif
