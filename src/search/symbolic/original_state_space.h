#ifndef SYMBOLIC_ORIGINAL_STATE_SPACE_H
#define SYMBOLIC_ORIGINAL_STATE_SPACE_H

#include "sym_state_space_manager.h"

namespace symbolic {
class OriginalStateSpace : public SymStateSpaceManager {
protected:
    virtual void init_initial_state() override;
    virtual void init_goal() override;
    virtual void init_individual_trs() override;

public:

    OriginalStateSpace(SymVariables *v, const SymParamsMgr &params, OperatorCost cost_type_);


    virtual std::string tag() const override {
        return "original";
    }


    virtual BDD shrinkExists(const BDD &bdd, int) const override {
        return bdd;
    }
    virtual BDD shrinkForall(const BDD &bdd, int) const override {
        return bdd;
    }
    virtual BDD shrinkTBDD(const BDD &bdd, int) const override  {
        return bdd;
    }

    virtual ADD getExplicitHeuristicADD(bool) override {
        return ADD();
    }

    virtual void getExplicitHeuristicBDD(bool, std::map<int, BDD> &) override {
    }
};
}
#endif
