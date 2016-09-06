#include "sym_pdb.h"

#include "sym_util.h"
#include "sym_transition.h"

#include "../utils/system.h"


using namespace std;

namespace symbolic {
SymPDB::SymPDB(SymVariables *bdd_vars, const SymParamsMgr &params, OperatorCost cost_type_) :
    SymStateSpaceManager(bdd_vars, params, cost_type_) {
    for (size_t i = 0; i < g_variable_name.size(); i++) {
        fullVars.insert(i);
    }

    nonRelVarsCube = bdd_vars->oneBDD();
    nonRelVarsCubeWithPrimes = bdd_vars->oneBDD();
    if (!nonRelVarsCube.IsCube()) {
        cout << "Error in sym_pdb: nonRelVars should be a cube";
        nonRelVarsCube.print(0, 1);
        cout << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
}

SymPDB::SymPDB(SymVariables *bdd_vars, const SymParamsMgr &params,
               OperatorCost cost_type_, AbsTRsStrategy,
               const set<int> &relevantVars) :
    SymStateSpaceManager(bdd_vars, params, cost_type_) {
    fullVars = relevantVars;
    for (size_t i = 0; i < g_variable_name.size(); i++) {
        if (!fullVars.count(i)) {
            nonRelVars.insert(i);
        }
    }

    nonRelVarsCube = vars->getCubePre(nonRelVars);    // * vars->getCubep(nonRelVars);
    nonRelVarsCubeWithPrimes = nonRelVarsCube * vars->getCubeEff(nonRelVars);
    if (!nonRelVarsCube.IsCube()) {
        cout << "Error in sym_pdb: nonRelVars should be a cube";
        nonRelVarsCube.print(0, 1);
        cout << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
}

SymPDB::SymPDB(shared_ptr<SymStateSpaceManager> parent,
               AbsTRsStrategy absTRsStrategy,
               const std::set<int> &relevantVars) :
    SymStateSpaceManager(parent, absTRsStrategy, relevantVars) {
    nonRelVarsCube = vars->getCubePre(nonRelVars);    // * vars->getCubep(nonRelVars);
    nonRelVarsCubeWithPrimes = nonRelVarsCube * vars->getCubeEff(nonRelVars);
    if (!nonRelVarsCube.IsCube()) {
        cout << "Error in sym_pdb: nonRelVars should be a cube";
        nonRelVarsCube.print(0, 1);
        cout << endl;
        utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
    }
}


BDD SymPDB::shrinkExists(const BDD &bdd, int maxNodes) const {
    return bdd.ExistAbstract(nonRelVarsCube, maxNodes);
}

BDD SymPDB::shrinkTBDD(const BDD &bdd, int maxNodes) const {
    return bdd.ExistAbstract(nonRelVarsCubeWithPrimes, maxNodes);
}

BDD SymPDB::shrinkForall(const BDD &bdd, int maxNodes) const {
    return bdd.UnivAbstract(nonRelVarsCube, maxNodes);
}

void SymPDB::init_individual_trs() {
    cout << "Not implemented SymPDB::init_individual_trs ()" << endl;
    utils::exit_with(utils::ExitCode::UNSUPPORTED);
}

// void SymPDB::init_mutex(const std::vector<MutexGroup> & mutex_groups) {
//      // if(/*p.init_mutex_from_parent &&*/ parentMgr){
//      //     setTimeLimit(p.max_mutex_time);
//      //     DEBUG_MSG(cout << "Init mutex from parent" << endl;);
//      //     mutexInitialized = true;
//      //     //Initialize mutexes from other manager
//      //     try{
//      //      for(auto & bdd : parentMgr->notMutexBDDsFw){
//      //          BDD shrinked = abstraction->shrinkExists(bdd, p.max_mutex_size);
//      //          notMutexBDDsFw.push_back(shrinked);
//      //      }
//      //      for(auto & bdd : parentMgr->notMutexBDDsBw){
//      //          BDD shrinked = abstraction->shrinkExists(bdd, p.max_mutex_size);
//      //          notMutexBDDsBw.push_back(shrinked);
//      //      }
//      //      unsetTimeLimit();
//      //     }catch(BDDError e){
//      //      unsetTimeLimit();
//      //      //Forget about it
//      //      vector<BDD>().swap(notMutexBDDsFw);
//      //      vector<BDD>().swap(notMutexBDDsBw);
//      //      init_mutex(mutex_groups, true, false);
//      //     }
// }


/*void SymPDB::getTransitions(map<int, std::vector <SymTransition> > & trs) const{
  cout << "Initialize trs "<< *this << endl;
  for(int i = 0; i < g_operators.size(); i++){
  const Operator * op = &(g_operators[i]);
  // Skip spurious operators
  if (op->spurious){
  continue;
  }
  int cost = op->get_cost();

  if(cost == 0){
  trs[0].push_back(SymTransition(vars, op, cost, *this));
  }else{
  trs[cost].push_back(SymTransition(vars, op, cost, *this));
  }
  }
  }*/


void SymPDB::init_initial_state() {
    vector<pair<int, int>> abstract_ini;
    for (int var : fullVars) {
        abstract_ini.push_back(std::pair<int, int> (var, g_initial_state_data[var]));
    }
    initialState = vars->getPartialStateBDD(abstract_ini);
}


void SymPDB::init_goal() {
    vector<pair<int, int>> abstract_goal;
    for (auto goal_var : g_goal) {
        if (isRelevantVar(goal_var.first)) {
            abstract_goal.push_back(goal_var);
        }
    }
    goal = vars->getPartialStateBDD(abstract_goal);
}

std::string SymPDB::tag() const {
    return "PDB";
}

void SymPDB::print(std::ostream &os, bool fullInfo) const {
    os << "PDB (" << fullVars.size() << "/" << (nonRelVars.size() + fullVars.size()) << "): ";
    for (int v : fullVars) {
        os << v << " ";
    }
    if (fullInfo && !nonRelVars.empty()) {
        os << " [";
        for (int v : fullVars)
            os << v << " ";
        os << "]";
        os << endl << "Abstracted propositions: ";
        for (int v : nonRelVars) {
            os << v << ": ";
            for (auto &prop : g_fact_names[v])
                cout << prop << ", ";
            os << endl;
        }
        os << endl << "Considered propositions: ";
        for (int v : fullVars) {
            os << v << ": ";
            for (auto &prop : g_fact_names[v])
                os << prop << ", ";
            os << endl;
        }
        os << endl;
    }
}

ADD SymPDB::getExplicitHeuristicADD(bool /*fw*/) {
    return vars->getADD(0);
}
void SymPDB::getExplicitHeuristicBDD(bool /*fw*/, map<int, BDD> &res) {
    res[0] = vars->oneBDD();
}
}
