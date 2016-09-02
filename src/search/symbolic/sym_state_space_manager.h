#ifndef SYMBOLIC_SYM_STATE_SPACE_MANAGER_H
#define SYMBOLIC_SYM_STATE_SPACE_MANAGER_H

#include "../operator_cost.h"
#include "sym_enums.h"
#include "sym_bucket.h"
#include "sym_variables.h"
#include "sym_util.h"

#include <vector>
#include <set>
#include <map>
#include <memory>
#include <cassert>

namespace options {
class OptionParser;
class Options;
}

namespace symbolic {
class SymVariables;
class SymTransition;

/*
 * All the methods may throw exceptions in case the time or nodes are exceeded.
 *
 */
class SymParamsMgr {
public:
    //Parameters to generate the TRs
    int max_tr_size, max_tr_time;

    //Parameters to generate the mutex BDDs
    MutexType mutex_type;
    int max_mutex_size, max_mutex_time;

    //Time and memory bounds for auxiliary operations
    int max_aux_nodes, max_aux_time;

    SymParamsMgr();
    SymParamsMgr(const options::Options &opts);
    static void add_options_to_parser(options::OptionParser &parser);
    void print_options() const;
};

class SymStateSpaceManager {
protected:
    SymVariables *vars;
    const SymParamsMgr p;
    const OperatorCost cost_type;

    // Hold a reference to the parent manager that can be used during
    // initialization. Uses weak_ptr in order to allow releasing the
    // resources used by the parent manager if necessary.
    std::weak_ptr<SymStateSpaceManager> parent_mgr;

    const AbsTRsStrategy abs_trs_strategy;

    //If the variable is fully/partially/not considered in the abstraction
    std::set <int> fullVars, absVars, nonRelVars;

    BDD initialState; // initial state
    BDD goal; // bdd representing the true (i.e. not simplified) goal-state

    std::map<int, std::vector <SymTransition>> transitions; //TRs
    int min_transition_cost; //minimum cost of non-zero cost transitions
    bool hasTR0; //If there is transitions with cost 0

    //Individual TRs: Useful for shrink and plan construction
    std::map<int, std::vector <SymTransition>> indTRs;

    bool mutexInitialized, mutexByFluentInitialized;

    //BDD representation of valid states (wrt mutex) for fw and bw search
    std::vector<BDD> notMutexBDDsFw, notMutexBDDsBw;

    //Dead ends for fw and bw searches. They are always removed in
    //filter_mutex (it does not matter which mutex_type we are using).
    std::vector<BDD> notDeadEndFw, notDeadEndBw;

    //notMutex relative for each fluent
    std::vector<std::vector<BDD>> notMutexBDDsByFluentFw, notMutexBDDsByFluentBw;
    std::vector<std::vector<BDD>> exactlyOneBDDsByFluent;

    void zero_preimage(const BDD &bdd, std::vector <BDD> &res, int maxNodes) const;
    void cost_preimage(const BDD &bdd, std::map <int, std::vector<BDD>> &res, int maxNodes) const;
    void zero_image(const BDD &bdd, std::vector <BDD> &res, int maxNodes) const;
    void cost_image(const BDD &bdd, std::map <int, std::vector<BDD>> &res, int maxNodes) const;

    virtual void init_initial_state() = 0;
    virtual void init_goal() = 0;

    virtual ADD getExplicitHeuristicADD(bool fw) = 0;
    virtual void getExplicitHeuristicBDD(bool fw, std::map<int, BDD> &res) = 0;

    virtual void getTransitions(const std::map<int, std::vector <SymTransition>> & /*individualTRs*/,
                                std::map<int, std::vector <SymTransition>> & /*res*/) const {
        std::cerr << "REBUILD TRs not supported by " << *this << std::endl;
        utils::exit_with(utils::ExitCode::UNSUPPORTED);
    }

    void shrinkTransitions(const std::map<int, std::vector <SymTransition>> &trs,
                           const std::map<int, std::vector <SymTransition>> &indTRs,
                           std::map<int, std::vector <SymTransition>> &res,
                           int maxTime, int maxNodes) const;

    BDD getRelVarsCubePre() const {
        return vars->getCubePre(fullVars) + vars->getCubePre(absVars);
    }

    BDD getRelVarsCubeEff() const {
        return vars->getCubeEff(fullVars) + vars->getCubeEff(absVars);
    }

    friend std::ostream &operator<<(std::ostream &os, const SymStateSpaceManager &state_space);


    virtual std::string tag() const = 0;

    void init_transitions_from_individual_trs();

    virtual void init_individual_trs() = 0;


    //Be careful of calling init_mutex and init_transitions before actually calling filter_mutex or image
    void init_mutex(const std::vector<MutexGroup> &mutex_groups,
                    bool genMutexBDDs, bool genMutexBDDsByFluent);
    void init_mutex(const std::vector<MutexGroup> &mutex_groups,
                    bool genMutexBDD, bool genMutexBDDByFluent, bool fw);


public:
    SymStateSpaceManager(SymVariables *v,
                         const SymParamsMgr &params,
                         OperatorCost cost_type_); //Original state space: All vars are relevant

    SymStateSpaceManager(std::shared_ptr<SymStateSpaceManager> &parent,
                         AbsTRsStrategy abs_trs_strategy_,
                         const std::set<int> &relevantVars); //Abstract state space (PDBs)

    virtual void init_mutex(const std::vector<MutexGroup> &mutex_groups);

    void init() {
        init_mutex(g_mutex_groups);
        init_transitions();
    }

    void init_transitions();

    inline bool isAbstracted() const {
        return !(absVars.empty() && nonRelVars.empty());
    }

    inline bool isOriginal() const {
        return !isAbstracted();
    }

    virtual BDD shrinkExists(const BDD &bdd, int maxNodes) const = 0;
    virtual BDD shrinkForall(const BDD &bdd, int maxNodes) const = 0;
    virtual BDD shrinkTBDD(const BDD &bdd, int maxNodes) const = 0;

    void filterMutex(Bucket &bucket, bool fw, bool initialization);
    void mergeBucket(Bucket &bucket) const;
    void mergeBucketAnd(Bucket &bucket) const;

    void shrinkBucket(Bucket &bucket, int maxNodes);

    void addDeadEndStates(bool fw, BDD bdd);

    void addDeadEndStates(const std::vector<BDD> &fw_dead_ends,
                          const std::vector<BDD> &bw_dead_ends);


    inline SymVariables *getVars() const {
        return vars;
    }

    inline const std::set <int> &getFullVars() const {
        return fullVars;
    }

    inline const std::set <int> &getAbsVars() const {
        return absVars;
    }

    inline const std::set <int> &getNonRelVars() const {
        return nonRelVars;
    }

    inline bool isRelevantVar(int var) const {
        return fullVars.count(var) > 0 || absVars.count(var);
    }

    int numVariablesToAbstract() const {
        return fullVars.size();
    }

    int numVariablesAbstracted() const {
        return absVars.size() + nonRelVars.size();
    }



    double stateCount(const Bucket &bucket) const {
        return vars->numStates(bucket);
    }

    inline BDD shrinkForall(const BDD &bdd) {
        setTimeLimit(p.max_aux_time);
        try{
            BDD res = shrinkForall(bdd, p.max_aux_nodes);
            unsetTimeLimit();
            return res;
        }catch (BDDError e) {
            unsetTimeLimit();
        }
        return zeroBDD();
    }


    inline long totalNodes() const {
        return vars->totalNodes();
    }

    inline unsigned long totalMemory() const {
        return vars->totalMemory();
    }

    inline const BDD &getGoal() {
        if (goal.IsZero()) {
            init_goal();
            assert(!goal.IsZero());
        }
        return goal;
    }

    const std::map<int, std::vector <SymTransition>> &getIndividualTRs() {
        if (indTRs.empty())
            init_individual_trs();
        return indTRs;
    }

    inline const BDD &getInitialState() {
        if (initialState.IsZero()) {
            init_initial_state();
            assert(!initialState.IsZero());
        }
        return initialState;
    }

    //Update binState
    inline int *getBinaryDescription(const GlobalState &state) const {
        return vars->getBinaryDescription(state);
    }

    inline BDD getBDD(int variable, int value) const {
        return vars->preBDD(variable, value);
    }

    inline Cudd *mgr() const {
        return vars->mgr();
    }

    inline BDD zeroBDD() const {
        return vars->zeroBDD();
    }

    inline BDD oneBDD() const {
        return vars->oneBDD();
    }

    inline const std::vector<BDD> &getNotMutexBDDs(bool fw) {
        init_mutex(g_mutex_groups);
        return fw ? notMutexBDDsFw : notMutexBDDsBw;
    }

    inline const std::vector<int> &vars_index_pre(int variable) const {
        return vars->vars_index_pre(variable);
    }

    inline const std::vector<int> &vars_index_eff(int variable) const {
        return vars->vars_index_eff(variable);
    }

    inline const std::vector<int> &vars_index_abs(int variable) const {
        return vars->vars_index_abs(variable);
    }

    bool mergeBucket(Bucket &bucket, int maxTime, int maxNodes) const {
        auto mergeBDDs = [] (BDD bdd, BDD bdd2, int maxNodes) {
                             return bdd.Or(bdd2, maxNodes);
                         };
        merge(vars, bucket, mergeBDDs, maxTime, maxNodes);
        removeZero(bucket); //Be sure that we do not contain only the zero BDD

        return bucket.size() <= 1;
    }

    bool mergeBucketAnd(Bucket &bucket, int maxTime, int maxNodes) const {
        auto mergeBDDs = [] (BDD bdd, BDD bdd2, int maxNodes) {
                             return bdd.And(bdd2, maxNodes);
                         };
        merge(vars, bucket, mergeBDDs, maxTime, maxNodes);
        removeZero(bucket); //Be sure that we do not contain only the zero BDD

        return bucket.size() <= 1;
    }

    void dumpMutexBDDs(bool fw) const;

    //Methods that require of TRs initialized
    inline int getMinTransitionCost() {
        assert(!transitions.empty());
        return min_transition_cost;
    }

    inline bool hasTransitions0() {
        assert(!transitions.empty());
        return hasTR0;
    }

    inline void zero_image(bool fw,
                           const BDD &bdd, std::vector<BDD> &res,
                           int maxNodes) {
        init_transitions();
        if (fw)
            zero_image(bdd, res, maxNodes);
        else
            zero_preimage(bdd, res, maxNodes);
    }

    inline void cost_image(bool fw,
                           const BDD &bdd, std::map <int, std::vector<BDD>> &res,
                           int maxNodes) {
        init_transitions();
        if (fw) {
            cost_image(bdd, res, maxNodes);
        } else {
            cost_preimage(bdd, res, maxNodes);
        }
    }

    //Methods that require of mutex initialized
    inline const BDD &getNotMutexBDDFw(int var, int val) {
        init_mutex(g_mutex_groups, false, true);
        return notMutexBDDsByFluentFw[var][val];
    }

    //Methods that require of mutex initialized
    inline const BDD &getNotMutexBDDBw(int var, int val) {
        init_mutex(g_mutex_groups, false, true);
        return notMutexBDDsByFluentBw[var][val];
    }

    //Methods that require of mutex initialized
    inline const BDD &getExactlyOneBDD(int var, int val) {
        init_mutex(g_mutex_groups, false, true);
        return exactlyOneBDDsByFluent[var][val];
    }

    BDD filter_mutex(const BDD &bdd,
                     bool fw, int maxNodes,
                     bool initialization);

    int filterMutexBucket(std::vector<BDD> &bucket, bool fw,
                          bool initialization, int maxTime, int maxNodes);


    inline void setTimeLimit(int maxTime) {
        vars->setTimeLimit(maxTime);
    }

    inline void unsetTimeLimit() {
        vars->unsetTimeLimit();
    }

    virtual void print(std::ostream &os, bool /*fullInfo*/) const {
        os << tag() << " (" << fullVars.size() << ")";
    }
};



/* class AbstractStateSpace : public SymStateSpaceManager { */

/*     AbsTRsStrategy absTRsStrategy; */

/*  public: */
/*     virtual std::string tag() const = 0; */
/* }; */
}
#endif
