#include "sym_state_space_manager.h"

#include "sym_enums.h"
#include "../utils/debug_macros.h"
#include <queue>
#include <limits>
#include <algorithm>

#include "sym_util.h"
#include "../globals.h"
#include "../global_operator.h"
#include "../mutex_group.h"
#include "../utils/timer.h"
#include "../options/options.h"
#include "../options/option_parser.h"
#include "../abstract_task.h"
#include "../global_operator.h"

using namespace std;

namespace symbolic {
    SymStateSpaceManager::SymStateSpaceManager(shared_ptr<SymStateSpaceManager> parent,
					       AbsTRsStrategy abs_trs_strategy_,
					       const std::set<int> &relevantVars) : SymStateSpaceManager(parent,
													 abs_trs_strategy_, 
													 relevantVars, 
													 parent->cost_type) {
}

    bool is_relevant_op(const GlobalOperator & op, const std::set<int> &relevantVars) {
	for(auto eff: op.get_effects ()) {
	    if(relevantVars.count(eff.var)) {
		return true;
	    }
	} 
	return false;
    }

    SymStateSpaceManager::SymStateSpaceManager(shared_ptr<SymStateSpaceManager> parent,
					       AbsTRsStrategy abs_trs_strategy_,
					       const std::set<int> &relevantVars, 
					       shared_ptr<OperatorCostFunction> cost_type_):
      vars(parent->vars), p(parent->p), cost_type(cost_type_),
      parent_mgr(parent), abs_trs_strategy(abs_trs_strategy_),
      fullVars(relevantVars),
      initialState(vars->zeroBDD()), goal(vars->zeroBDD()),
      min_transition_cost(parent->min_transition_cost),
      hasTR0(parent->hasTR0), mutexInitialized(false),
	mutexByFluentInitialized(false) {


	for (size_t i = 0; i < g_operators.size(); ++i) {
	    if (!is_relevant_op(g_operators[i], relevantVars)) continue;     
	    if (min_transition_cost == 0 || min_transition_cost > cost_type->get_adjusted_cost(i)) {
		min_transition_cost = cost_type->get_adjusted_cost(i);
	    }
	    if (cost_type->get_adjusted_cost(i) == 0) {
		hasTR0 = true;
	    }
	}
    

        

	for (size_t i = 0; i < g_variable_name.size(); i++) {
	    if (!fullVars.count(i)) {
		nonRelVars.insert(i);
	    }
	}
    }



SymStateSpaceManager::SymStateSpaceManager(SymVariables *v,
					   const SymParamsMgr &params,
					   shared_ptr<OperatorCostFunction> cost_type_) :
	vars(v), p(params), cost_type(cost_type_),
	abs_trs_strategy(AbsTRsStrategy::REBUILD_TRS),
	initialState(v->zeroBDD()), goal(v->zeroBDD()),
	min_transition_cost(0), hasTR0(false),
	mutexInitialized(false),
	mutexByFluentInitialized(false) {
    for (size_t i = 0; i < g_operators.size(); ++i) {
	if (min_transition_cost == 0 || min_transition_cost > cost_type->get_adjusted_cost(i)) {
	    min_transition_cost = cost_type->get_adjusted_cost(i);
	}
	if (cost_type->get_adjusted_cost(i) == 0) {
	    hasTR0 = true;
	}
    }
}



void SymStateSpaceManager::init_transitions_from_individual_trs() {
    if (!transitions.empty())
        return;                          //Already initialized!
    DEBUG_MSG(cout << "Init transitions" << endl;);

    DEBUG_MSG(cout << "Generate individual TRs" << endl;);
    transitions = map<int, vector <TransitionRelation>> (indTRs);     //Copy
    DEBUG_MSG(cout << "Individual TRs generated" << endl;);
    min_transition_cost = 0;
    hasTR0 = transitions.count(0) > 0;

    for (map<int, vector<TransitionRelation>>::iterator it = transitions.begin();
         it != transitions.end(); ++it) {
        merge(vars, it->second, mergeTR, p.max_tr_time, p.max_tr_size);

        if (min_transition_cost == 0 || min_transition_cost > it->first) {
            min_transition_cost = it->first;
        }

        cout << "TRs cost=" << it->first << " (" << it->second.size() << ")" << endl;
    }
}


void SymStateSpaceManager::init_mutex(const std::vector<MutexGroup> &mutex_groups,
                                      bool genMutexBDD, bool genMutexBDDByFluent) {
    //Check if I should initialize something and return
    if (mutexInitialized)
        genMutexBDD = false;
    if (mutexByFluentInitialized)
        genMutexBDDByFluent = false;
    if (!genMutexBDD && !genMutexBDDByFluent)
        return;
    if (genMutexBDD)
        mutexInitialized = true;
    if (genMutexBDDByFluent)
        mutexByFluentInitialized = true;

    if (genMutexBDDByFluent) {
        //Initialize structure for exactlyOneBDDsByFluent (common to both init_mutex calls)
        exactlyOneBDDsByFluent.resize(g_variable_domain.size());
        for (size_t i = 0; i < g_variable_domain.size(); ++i) {
            exactlyOneBDDsByFluent[i].resize(g_variable_domain[i]);
            for (int j = 0; j < g_variable_domain[i]; ++j) {
                exactlyOneBDDsByFluent[i][j] = oneBDD();
            }
        }
    }

    init_mutex(mutex_groups, genMutexBDD, genMutexBDDByFluent, false);
    init_mutex(mutex_groups, genMutexBDD, genMutexBDDByFluent, true);
}

void SymStateSpaceManager::init_mutex(const std::vector<MutexGroup> &mutex_groups,
                                      bool genMutexBDD, bool genMutexBDDByFluent, bool fw) {
    DEBUG_MSG(cout << "Init mutex BDDs " << (fw ? "fw" : "bw") << ": "
                   << genMutexBDD << " " << genMutexBDDByFluent << endl;);

    vector<vector<BDD>> &notMutexBDDsByFluent =
        (fw ? notMutexBDDsByFluentFw : notMutexBDDsByFluentBw);

    vector<BDD> &notMutexBDDs =
        (fw ? notMutexBDDsFw : notMutexBDDsBw);

    //BDD validStates = vars->oneBDD();
    int num_mutex = 0;
    int num_invariants = 0;

    if (genMutexBDDByFluent) {
        //Initialize structure for notMutexBDDsByFluent
        notMutexBDDsByFluent.resize(g_variable_domain.size());
        for (size_t i = 0; i < g_variable_domain.size(); ++i) {
            notMutexBDDsByFluent[i].resize(g_variable_domain[i]);
            for (int j = 0; j < g_variable_domain[i]; ++j) {
                notMutexBDDsByFluent[i][j] = oneBDD();
            }
        }
    }

    //Initialize mBDDByVar and invariant_bdds_by_fluent
    vector<BDD>  mBDDByVar;
    mBDDByVar.reserve(g_variable_domain.size());
    vector<vector<BDD>> invariant_bdds_by_fluent(g_variable_domain.size());
    for (size_t i = 0; i < invariant_bdds_by_fluent.size(); i++) {
        mBDDByVar.push_back(oneBDD());
        invariant_bdds_by_fluent[i].resize(g_variable_domain[i]);
        for (size_t j = 0; j < invariant_bdds_by_fluent[i].size(); j++) {
            invariant_bdds_by_fluent[i][j] = oneBDD();
        }
    }

    for (auto &mg : mutex_groups) {
        if (mg.pruneFW() != fw)
            continue;
        const vector<Fact> &invariant_group = mg.getFacts();
        DEBUG_MSG(cout << mg << endl;);
        if (mg.isExactlyOne()) {
            BDD bddInvariant = zeroBDD();
            int var = numeric_limits<int>::max();
            int val = 0;
            bool exactlyOneRelevant = true;

            for (auto &fluent : invariant_group) {
                if (!isRelevantVar(fluent.var)) {
                    exactlyOneRelevant = true;
                    break;
                }
                bddInvariant += vars->preBDD(fluent.var, fluent.value);
                if (fluent.var < var) {
                    var = fluent.var;
                    val = fluent.value;
                }
            }

            if (exactlyOneRelevant) {
                num_invariants++;
                if (genMutexBDD) {
                    invariant_bdds_by_fluent[var][val] *= bddInvariant;
                }
                if (genMutexBDDByFluent) {
                    for (auto &fluent : invariant_group) {
                        exactlyOneBDDsByFluent[fluent.var][fluent.value] *= bddInvariant;
                    }
                }
            }
        }


        for (size_t i = 0; i < invariant_group.size(); ++i) {
            int var1 = invariant_group[i].var;
            if (!isRelevantVar(var1))
                continue;
            int val1 = invariant_group[i].value;
            BDD f1 = vars->preBDD(var1, val1);

            for (size_t j = i + 1; j < invariant_group.size(); ++j) {
                int var2 = invariant_group[j].var;
                if (!isRelevantVar(var2))
                    continue;
                int val2 = invariant_group[j].value;
                BDD f2 = vars->preBDD(var2, val2);
                BDD mBDD = !(f1 * f2);
                if (genMutexBDD) {
                    num_mutex++;
                    mBDDByVar[min(var1, var2)] *= mBDD;
                    if (mBDDByVar[min(var1, var2)].nodeCount() > p.max_mutex_size) {
                        notMutexBDDs.push_back(mBDDByVar[min(var1, var2)]);
                        mBDDByVar[min(var1, var2)] = vars->oneBDD();
                    }
                }
                if (genMutexBDDByFluent) {
                    notMutexBDDsByFluent[var1][val1] *= mBDD;
                    notMutexBDDsByFluent[var2][val2] *= mBDD;
                }
            }
        }
    }

    if (genMutexBDD) {
        for (size_t var = 0; var < g_variable_domain.size(); ++var) {
            if (!mBDDByVar[var].IsOne()) {
                notMutexBDDs.push_back(mBDDByVar[var]);
            }
            for (const BDD &bdd_inv : invariant_bdds_by_fluent[var]) {
                if (!bdd_inv.IsOne()) {
                    notMutexBDDs.push_back(bdd_inv);
                }
            }
        }

        DEBUG_MSG(dumpMutexBDDs(fw););
        merge(vars, notMutexBDDs, mergeAndBDD,
              p.max_mutex_time, p.max_mutex_size);
        std::reverse(notMutexBDDs.begin(), notMutexBDDs.end());
        DEBUG_MSG(cout << "Mutex initialized " << (fw ? "fw" : "bw") << ". Total mutex added: " << num_mutex << " Invariant groups: " << num_invariants << endl;);
        DEBUG_MSG(dumpMutexBDDs(fw););
    }
}

void SymStateSpaceManager::addDeadEndStates(bool fw, BDD bdd) {
    //There are several options here, we could follow with edeletion
    //and modify the TRs, so that the new spurious states are never
    //generated. However, the TRs are already merged and the may get
    //too large. Therefore we just keep this states in another vectors
    //and spurious states are always removed. TODO: this could be
    //improved.
    if (fw || isAbstracted()) {
        if (isAbstracted())
            bdd = shrinkForall(bdd);
        notDeadEndFw.push_back(!bdd);
        mergeBucketAnd(notDeadEndFw);
    } else {
        notDeadEndBw.push_back(!bdd);
        mergeBucketAnd(notDeadEndBw);
    }
}


void SymStateSpaceManager::addDeadEndStates(const std::vector<BDD> &fw_dead_ends,
                                            const std::vector<BDD> &bw_dead_ends) {
    for (BDD bdd : fw_dead_ends) {
        bdd = shrinkForall(bdd);
        if (!bdd.IsZero()) {
            notDeadEndFw.push_back(!bdd);
        }
    }

    for (BDD bdd : bw_dead_ends) {
        bdd = shrinkForall(bdd);
        if (!bdd.IsZero()) {
            notDeadEndFw.push_back(!bdd);
        }
    }
    mergeBucketAnd(notDeadEndFw);
}


void SymStateSpaceManager::dumpMutexBDDs(bool fw) const {
    if (fw) {
        cout << "Mutex BDD FW Size(" << p.max_mutex_size << "):";
        for (const auto &bdd : notMutexBDDsFw) {
            cout << " " << bdd.nodeCount();
        }
        cout << endl;
    } else {
        cout << "Mutex BDD BW Size(" << p.max_mutex_size << "):";
        for (const auto &bdd : notMutexBDDsBw) {
            cout << " " << bdd.nodeCount();
        }
        cout << endl;
    }
}

void SymStateSpaceManager::zero_preimage(const BDD &bdd, vector <BDD> &res, int nodeLimit) const {
    for (const auto &tr : transitions.at(0)) {
        res.push_back(tr.preimage(bdd, nodeLimit));
    }
}

void SymStateSpaceManager::zero_image(const BDD &bdd, vector <BDD> &res, int nodeLimit) const {
    for (const auto &tr : transitions.at(0)) {
        res.push_back(tr.image(bdd, nodeLimit));
    }
}

void SymStateSpaceManager::cost_preimage(const BDD &bdd, map<int, vector<BDD>> &res,
                                         int nodeLimit) const {
    for (auto trs : transitions) {
        int cost = trs.first;
        if (cost == 0)
            continue;
        for (size_t i = res[cost].size(); i < trs.second.size(); i++) {
            BDD result = trs.second[i].preimage(bdd, nodeLimit);
            res[cost].push_back(result);
        }
    }
}

void SymStateSpaceManager::cost_image(const BDD &bdd,
                                      map<int, vector<BDD>> &res,
                                      int nodeLimit) const {
    for (auto trs : transitions) {
        int cost = trs.first;
        if (cost == 0)
            continue;
        for (size_t i = res[cost].size(); i < trs.second.size(); i++) {
            //cout << "Img: " << trs.second[i].nodeCount() << " with bdd " << bdd.nodeCount() << " node limit: " << nodeLimit << endl;
            BDD result = trs.second[i].image(bdd, nodeLimit);
            //cout << "Res: " << result.nodeCount() << endl;
            res[cost].push_back(result);
        }
    }
}

BDD SymStateSpaceManager::filter_mutex(const BDD &bdd, bool fw,
                                       int nodeLimit, bool initialization) {
    BDD res = bdd;
    const vector<BDD> &notDeadEndBDDs = ((fw || isAbstracted()) ? notDeadEndFw : notDeadEndBw);
    for (const BDD &notDeadEnd : notDeadEndBDDs) {
        DEBUG_MSG(cout << "Filter: " << res.nodeCount() << " and dead end " << notDeadEnd.nodeCount() << flush;);
        res = res.And(notDeadEnd, nodeLimit);
        DEBUG_MSG(cout << ": " << res.nodeCount() << endl;);
    }

    const vector<BDD> &notMutexBDDs = (fw ? notMutexBDDsFw : notMutexBDDsBw);


    switch (p.mutex_type) {
    case MutexType::MUTEX_NOT:
        break;
    case MutexType::MUTEX_EDELETION:
        if (initialization) {
            for (const BDD &notMutexBDD : notMutexBDDs) {
                DEBUG_MSG(cout << res.nodeCount() << " and " << notMutexBDD.nodeCount() << flush;);
                res = res.And(notMutexBDD, nodeLimit);
                DEBUG_MSG(cout << ": " << res.nodeCount() << endl;);
            }
        }
        break;
    case MutexType::MUTEX_AND:
        for (const BDD &notMutexBDD : notMutexBDDs) {
            DEBUG_MSG(cout << "Filter: " << res.nodeCount() << " and " << notMutexBDD.nodeCount() << flush;);
            res = res.And(notMutexBDD, nodeLimit);
            DEBUG_MSG(cout << ": " << res.nodeCount() << endl;);
        }
        break;
    case MutexType::MUTEX_RESTRICT:
        for (const BDD &notMutexBDD : notMutexBDDs)
            res = res.Restrict(notMutexBDD);
        break;
    case MutexType::MUTEX_NPAND:
        for (const BDD &notMutexBDD : notMutexBDDs)
            res = res.NPAnd(notMutexBDD);
        break;
    case MutexType::MUTEX_CONSTRAIN:
        for (const BDD &notMutexBDD : notMutexBDDs)
            res = res.Constrain(notMutexBDD);
        break;
    case MutexType::MUTEX_LICOMP:
        for (const BDD &notMutexBDD : notMutexBDDs)
            res = res.LICompaction(notMutexBDD);
        break;
    }
    return res;
}

int SymStateSpaceManager::filterMutexBucket(vector<BDD> &bucket, bool fw,
                                            bool initialization, int maxTime, int maxNodes) {
    int numFiltered = 0;
    setTimeLimit(maxTime);
    try{
        for (size_t i = 0; i < bucket.size(); ++i) {
            DEBUG_MSG(cout << "Filter spurious " << (fw ? "fw" : "bw") << ": " << *this
                           << " from: " << bucket[i].nodeCount() <<
                      " maxTime: " << maxTime << " and maxNodes: " << maxNodes;);

            bucket[i] = filter_mutex(bucket[i], fw, maxNodes, initialization);
            DEBUG_MSG(cout << " => " << bucket[i].nodeCount() << endl;);
            numFiltered++;
        }
    }catch (BDDError e) {
        DEBUG_MSG(cout << " truncated." << endl;);
    }
    unsetTimeLimit();

    return numFiltered;
}

void SymStateSpaceManager::filterMutex(Bucket &bucket, bool fw, bool initialization) {
    filterMutexBucket(bucket, fw, initialization,
                      p.max_aux_time, p.max_aux_nodes);
}

void SymStateSpaceManager::mergeBucket(Bucket &bucket) const {
    mergeBucket(bucket, p.max_aux_time, p.max_aux_nodes);
}

void SymStateSpaceManager::mergeBucketAnd(Bucket &bucket) const {
    mergeBucketAnd(bucket, p.max_aux_time, p.max_aux_nodes);
}

void SymStateSpaceManager::shrinkBucket(Bucket &bucket, int maxNodes) {
    for (size_t i = 0; i < bucket.size(); ++i) {
        bucket[i] = shrinkExists(bucket[i], maxNodes);
    }
}


void SymStateSpaceManager::init_mutex(const std::vector<MutexGroup> &mutex_groups) {
    //If (a) is initialized OR not using mutex OR edeletion does not need mutex
    if (mutexInitialized || p.mutex_type == MutexType::MUTEX_NOT)
        return;     //Skip mutex initialization

    if (p.mutex_type == MutexType::MUTEX_EDELETION) {
        SymStateSpaceManager::init_mutex(mutex_groups, true, true);
    } else {
        SymStateSpaceManager::init_mutex(mutex_groups, true, false);
    }
}


void SymStateSpaceManager::init_transitions() {
    if (!transitions.empty())
        return;                          //Already initialized!
    if (parent_mgr.expired() || abs_trs_strategy == AbsTRsStrategy::REBUILD_TRS) {
        init_individual_trs();
        init_transitions_from_individual_trs();
    } else {
        auto parent = parent_mgr.lock();
        assert(!parent->transitions.empty());

        map<int, vector <TransitionRelation>> failedToShrink;
        switch (abs_trs_strategy) {
        case AbsTRsStrategy::TR_SHRINK:
            for (const auto &trsParent : parent->transitions) {
                int cost = trsParent.first;     //For all the TRs of cost cost
                DEBUG_MSG(cout << "Init trs: " << cost << endl;);
                set <const GlobalOperator *> failed_ops;

                double num_trs = parent->transitions.size() * trsParent.second.size();
                for (const auto &trParent : trsParent.second) {
                    TransitionRelation absTransition = TransitionRelation(trParent);
                    DEBUG_MSG(cout << "SHRINK: " << absTransition << " time_out: "
                                   << 1 + p.max_aux_time / num_trs << " max nodes: "
                                   << 1 + p.max_aux_nodes << endl;);
                    try{
                        vars->setTimeLimit(1 + p.max_aux_time / num_trs);
                        absTransition.shrink(*this, 1 + p.max_aux_nodes);
                        transitions[cost].push_back(move(absTransition));
                        vars->unsetTimeLimit();
                    }catch (BDDError e) {
                        vars->unsetTimeLimit();
                        DEBUG_MSG(cout << "Failed shrinking TR" << endl;);
                        //Failed some
                        const set <const GlobalOperator *> &tr_ops = trParent.getOps();
                        failed_ops.insert(begin(tr_ops), end(tr_ops));
                    }
                }

                if (!failed_ops.empty()) {   //Add all the TRs related with it.
                    cout << "Failed ops" << endl;
                    for (const auto &trParent : indTRs.at(cost)) {
                        if (trParent.hasOp(failed_ops)) {
                            TransitionRelation absTransition = TransitionRelation(trParent);
                            vars->setTimeLimit(p.max_aux_time);
                            try{
                                absTransition.shrink(*this, p.max_aux_nodes);
                                transitions[cost].push_back(absTransition);
                                vars->unsetTimeLimit();
                            }catch (BDDError e) {
                                failedToShrink[cost].push_back(absTransition);
                                vars->unsetTimeLimit();
                            }
                        }
                    }
                }
                merge(vars, transitions[cost], mergeTR, p.max_aux_time / parent->transitions.size(), p.max_tr_size);
            }
            break;
        case AbsTRsStrategy::IND_TR_SHRINK:
            for (const auto &indTRsCost : parent->indTRs) {
		for (const auto &trParent : indTRsCost.second) {
                    TransitionRelation absTransition = TransitionRelation(trParent);
		    assert (absTransition.getOps().size() == 1);
		    if(!is_relevant_op(**(absTransition.getOps().begin()), fullVars)) continue;
		    int cost = cost_type->get_adjusted_cost(*(absTransition.getOps().begin()));
		    if(cost != absTransition.getCost()) absTransition.set_cost(cost);
                    try{
                        vars->setTimeLimit(p.max_aux_time);
                        absTransition.shrink(*this, p.max_aux_nodes);
                        vars->unsetTimeLimit();
                        transitions[cost].push_back(absTransition);
                    }catch (BDDError e) {
                        vars->unsetTimeLimit();
                        failedToShrink[cost].push_back(absTransition);
                    }
                }
            }

	    for (auto & trs : transitions) {
                merge(vars, trs.second, mergeTR, p.max_aux_time, p.max_tr_size);
	    }

            break;
        case AbsTRsStrategy::SHRINK_AFTER_IMG:
            //SetAbsAfterImage
            for (const auto &t : parent->transitions) {
                int cost = t.first;

                for (const auto &tr : t.second) {
                    TransitionRelation newTR = tr;
                    newTR.setAbsAfterImage(this);
                    transitions[cost].push_back(newTR);
                }
            }
            break;

        case AbsTRsStrategy::REBUILD_TRS:
            assert(false);
            break;
        }

        //Use Shrink after img in all the transitions that failedToShrink
        DEBUG_MSG(cout << "Failed to shrink: " << (failedToShrink.empty() ? "no" : "yes") << endl;);
        for (auto &failedTRs : failedToShrink) {
            merge(vars, failedTRs.second, mergeTR, p.max_aux_time, p.max_tr_size);
            for (auto &tr : failedTRs.second) {
                tr.setAbsAfterImage(this);
                transitions[failedTRs.first].push_back(tr);
            }
        }
    }

    DEBUG_MSG(cout << "Finished init trs: " << transitions.size() << endl;);
}


SymParamsMgr::SymParamsMgr(const options::Options &opts) :
    max_tr_size(opts.get<int>("max_tr_size")),
    max_tr_time(opts.get<int>("max_tr_time")),
    mutex_type(MutexType(opts.get_enum("mutex_type"))),
    max_mutex_size(opts.get<int>("max_mutex_size")),
    max_mutex_time(opts.get<int>("max_mutex_time")),
    max_aux_nodes(opts.get<int>("max_aux_nodes")),
    max_aux_time(opts.get<int>("max_aux_time")) {
    //Don't use edeletion with conditional effects
    if (mutex_type == MutexType::MUTEX_EDELETION && has_conditional_effects()) {
        cout << "Mutex type changed to mutex_and because the domain has conditional effects" << endl;
        mutex_type = MutexType::MUTEX_AND;
    }
}

SymParamsMgr::SymParamsMgr() :
    max_tr_size(100000),
    max_tr_time(60000),
    mutex_type(MutexType::MUTEX_EDELETION),
    max_mutex_size(100000),
    max_mutex_time(60000),
    max_aux_nodes(1000000), max_aux_time(2000) {
    //Don't use edeletion with conditional effects
    if (mutex_type == MutexType::MUTEX_EDELETION && has_conditional_effects()) {
        cout << "Mutex type changed to mutex_and because the domain has conditional effects" << endl;
        mutex_type = MutexType::MUTEX_AND;
    }
}

void SymParamsMgr::print_options() const {
    cout << "TR(time=" << max_tr_time << ", nodes=" << max_tr_size << ")" << endl;
    cout << "Mutex(time=" << max_mutex_time << ", nodes=" << max_mutex_size << ", type=" << mutex_type << ")" << endl;
    cout << "Aux(time=" << max_aux_time << ", nodes=" << max_aux_nodes << ")" << endl;
}

void SymParamsMgr::add_options_to_parser(options::OptionParser &parser) {
    parser.add_option<int> ("max_tr_size", "maximum size of TR BDDs", "100000");

    parser.add_option<int> ("max_tr_time",
                            "maximum time (ms) to generate TR BDDs", "60000");

    parser.add_enum_option("mutex_type", MutexTypeValues,
                           "mutex type", "MUTEX_EDELETION");

    parser.add_option<int> ("max_mutex_size",
                            "maximum size of mutex BDDs", "100000");

    parser.add_option<int> ("max_mutex_time",
                            "maximum time (ms) to generate mutex BDDs", "60000");

    parser.add_option<int> ("max_aux_nodes", "maximum size in pop operations", "1000000");
    parser.add_option<int> ("max_aux_time", "maximum time (ms) in pop operations", "2000");
}

std::ostream &operator<<(std::ostream &os, const SymStateSpaceManager &abs) {
    abs.print(os, false);
    return os;
}
}
