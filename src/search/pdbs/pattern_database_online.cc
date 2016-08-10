#include "pattern_database_online.h"

#include "../priority_queue.h"
#include "../task_tools.h"

#include "../utils/collections.h"
#include "../utils/logging.h"
#include "../utils/math.h"
#include "../utils/timer.h"

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <boost/lexical_cast.hpp>

using namespace std;

namespace pdbs {
  //NOTE: WE ARE DOING FORWARD SEARCH, SO PRE AND EFF ARE REVERSED FROM NORMAL PDB BACKWARD SEARCH FROM GOAL
  //THIS IS WHY WE USE AbstractOperatorOnline instead of the regular
AbstractOperatorOnline::AbstractOperatorOnline(const vector<pair<int, int>> &prev_pairs,
                                   const vector<pair<int, int>> &pre_pairs,
                                   const vector<pair<int, int>> &eff_pairs,
                                   int cost,
                                   const vector<size_t> &hash_multipliers)
    : cost(cost),
      regression_preconditions(prev_pairs) {
	cout<<"\tAbstract_operator_cost:"<<cost<<endl;
    //regression_preconditions.insert(regression_preconditions.end(), eff_pairs.begin(), eff_pairs.end());
    regression_preconditions.insert(regression_preconditions.end(),
                                    pre_pairs.begin(),
                                    pre_pairs.end());
    // Sort preconditions for MatchTreeOnline construction.
    sort(regression_preconditions.begin(), regression_preconditions.end());
    for (size_t i = 1; i < regression_preconditions.size(); ++i) {
        assert(regression_preconditions[i].first !=
               regression_preconditions[i - 1].first);
    }
    hash_effect = 0;
    assert(pre_pairs.size() == eff_pairs.size());
    for (size_t i = 0; i < pre_pairs.size(); ++i) {
        int var = pre_pairs[i].first;
        assert(var == eff_pairs[i].first);
        //int old_val = eff_pairs[i].second;
        //int new_val = pre_pairs[i].second;
//We are doing forward search when doing online pdb generation!
        int old_val = pre_pairs[i].second;
        int new_val = eff_pairs[i].second;
        assert(new_val != -1);
        size_t effect = (new_val - old_val) * hash_multipliers[var];
        hash_effect += effect;
    }
}

AbstractOperatorOnline::~AbstractOperatorOnline() {
}

void AbstractOperatorOnline::dump(const Pattern &pattern,
                            const TaskProxy &task_proxy) const {
    cout << "AbstractOperatorOnline:" << endl;
    cout << "Regression preconditions:" << endl;
    for (size_t i = 0; i < regression_preconditions.size(); ++i) {
        int var_id = regression_preconditions[i].first;
        int val = regression_preconditions[i].second;
        cout << "Variable: " << var_id << " (True name: "
             << task_proxy.get_variables()[pattern[var_id]].get_name()
             << ", Index: " << i << ") Value: " << val << endl;
    }
    cout << "Hash effect:" << hash_effect << endl;
}

PatternDatabaseOnline::PatternDatabaseOnline(
    const TaskProxy &task_proxy,
    const Pattern &pattern,
    bool dump,
    const vector<int> &operator_costs)
    : task_proxy(task_proxy),
      pattern(pattern),
      match_tree(task_proxy,pattern,hash_multipliers){
    verify_no_axioms(task_proxy);
    verify_no_conditional_effects(task_proxy);
    assert(operator_costs.empty() ||
           operator_costs.size() == task_proxy.get_operators().size());
    assert(utils::is_sorted_unique(pattern));

    utils::Timer timer;
    hash_multipliers.reserve(pattern.size());
    num_states = 1;
    for (int pattern_var_id : pattern) {
        hash_multipliers.push_back(num_states);
        VariableProxy var = task_proxy.get_variables()[pattern_var_id];
        if (utils::is_product_within_limit(num_states, var.get_domain_size(),
                                           numeric_limits<int>::max())) {
            num_states *= var.get_domain_size();
        } else {
            cerr << "Given pattern is too large! (Overflow occured): " << endl;
            cerr << pattern << endl;
            utils::exit_with(utils::ExitCode::CRITICAL_ERROR);
        }
    }
    create_pdb(operator_costs);
    if (dump)
        cout << "PDB construction time: " << timer << endl;
}

void PatternDatabaseOnline::multiply_out(
    int pos, int cost, vector<pair<int, int>> &prev_pairs,
    vector<pair<int, int>> &pre_pairs,
    vector<pair<int, int>> &eff_pairs,
    const vector<pair<int, int>> &effects_without_pre,
    vector<AbstractOperatorOnline> &operators) {
    if (pos == static_cast<int>(effects_without_pre.size())) {
        // All effects without precondition have been checked: insert op.
        if (!eff_pairs.empty()) {
            operators.push_back(
                AbstractOperatorOnline(prev_pairs, pre_pairs, eff_pairs, cost,
                                 hash_multipliers));
        }
    } else {
        // For each possible value for the current variable, build an
        // abstract operator.
        int var_id = effects_without_pre[pos].first;
        int eff = effects_without_pre[pos].second;
        VariableProxy var = task_proxy.get_variables()[pattern[var_id]];
        for (int i = 0; i < var.get_domain_size(); ++i) {
            if (i != eff) {
                pre_pairs.push_back(make_pair(var_id, i));
                eff_pairs.push_back(make_pair(var_id, eff));
            } else {
                prev_pairs.push_back(make_pair(var_id, i));
            }
            multiply_out(pos + 1, cost, prev_pairs, pre_pairs, eff_pairs,
                         effects_without_pre, operators);
            if (i != eff) {
                pre_pairs.pop_back();
                eff_pairs.pop_back();
            } else {
                prev_pairs.pop_back();
            }
        }
    }
}

void PatternDatabaseOnline::build_abstract_operators(
    const OperatorProxy &op, int cost,
    const std::vector<int> &variable_to_index,
    vector<AbstractOperatorOnline> &operators) {
    // All variable value pairs that are a prevail condition
    vector<pair<int, int>> prev_pairs;
    // All variable value pairs that are a precondition (value != -1)
    vector<pair<int, int>> pre_pairs;
    // All variable value pairs that are an effect
    vector<pair<int, int>> eff_pairs;
    // All variable value pairs that are a precondition (value = -1)
    vector<pair<int, int>> effects_without_pre;

    size_t num_vars = task_proxy.get_variables().size();
    vector<bool> has_precond_and_effect_on_var(num_vars, false);
    vector<bool> has_precondition_on_var(num_vars, false);

    for (FactProxy pre : op.get_preconditions())
        has_precondition_on_var[pre.get_variable().get_id()] = true;

    for (EffectProxy eff : op.get_effects()) {
        int var_id = eff.get_fact().get_variable().get_id();
        int pattern_var_id = variable_to_index[var_id];
        int val = eff.get_fact().get_value();
        if (pattern_var_id != -1) {
            if (has_precondition_on_var[var_id]) {
                has_precond_and_effect_on_var[var_id] = true;
                eff_pairs.push_back(make_pair(pattern_var_id, val));
            } else {
                effects_without_pre.push_back(make_pair(pattern_var_id, val));
            }
        }
    }
    for (FactProxy pre : op.get_preconditions()) {
        int var_id = pre.get_variable().get_id();
        int pattern_var_id = variable_to_index[var_id];
        int val = pre.get_value();
        if (pattern_var_id != -1) { // variable occurs in pattern
            if (has_precond_and_effect_on_var[var_id]) {
                pre_pairs.push_back(make_pair(pattern_var_id, val));
            } else {
                prev_pairs.push_back(make_pair(pattern_var_id, val));
            }
        }
    }
    multiply_out(0, cost, prev_pairs, pre_pairs, eff_pairs, effects_without_pre,
                 operators);
}

//Modified, only want to create abstract operators and match_tree
void PatternDatabaseOnline::create_pdb(const std::vector<int> &operator_costs) {
  //needed by the online operator search function OnlineDistCalculator
  operator_costs_copy=operator_costs;

    VariablesProxy vars = task_proxy.get_variables();
    vector<int> variable_to_index(vars.size(), -1);
    for (size_t i = 0; i < pattern.size(); ++i) {
        variable_to_index[pattern[i]] = i;
    }

    // compute all abstract operators
    vector<AbstractOperatorOnline> operators;
    for (OperatorProxy op : task_proxy.get_operators()) {
        int op_cost;
        if (operator_costs.empty()) {
            op_cost = op.get_cost();
        } else {
            op_cost = operator_costs[op.get_id()];
        }
        build_abstract_operators(op, op_cost, variable_to_index, operators);
    }

    // build the match tree
    //MatchTreeOnline match_tree(task_proxy, pattern, hash_multipliers);
    match_tree.update(pattern, hash_multipliers);
    cout<<"operators for match tree size:"<<operators.size()<<endl;
    for (const AbstractOperatorOnline &op : operators) {
        match_tree.insert(op);
    }
    cout<<"match_tree dump:"<<endl;
    match_tree.dump();
    //exit(0);
    
    for (FactProxy goal : task_proxy.get_goals()) {
        int var_id = goal.get_variable().get_id();
        int val = goal.get_value();
        if (variable_to_index[var_id] != -1) {
            abstract_goals.push_back(make_pair(variable_to_index[var_id], val));
        }
    }
    return;

}
bool PatternDatabaseOnline::is_goal_state(
    const size_t state_index){
    for (pair<int, int> abstract_goal : abstract_goals) {
        int pattern_var_id = abstract_goal.first;
        int var_id = pattern[pattern_var_id];
        VariableProxy var = task_proxy.get_variables()[var_id];
        int temp = state_index / hash_multipliers[pattern_var_id];
        int val = temp % var.get_domain_size();
        if (val != abstract_goal.second) {
            return false;
        }
    }
    return true;
}

size_t PatternDatabaseOnline::hash_index(const State &state) const {
    size_t index = 0;
    for (size_t i = 0; i < pattern.size(); ++i) {
        index += hash_multipliers[i] * state[pattern[i]].get_value();
    }
    return index;
}

int PatternDatabaseOnline::get_value(const State &state) const {
    return distances[hash_index(state)];
}

double PatternDatabaseOnline::compute_mean_finite_h() const {
    double sum = 0;
    int size = 0;
    for (size_t i = 0; i < distances.size(); ++i) {
        if (distances[i] != numeric_limits<int>::max()) {
            sum += distances[i];
            ++size;
        }
    }
    if (size == 0) { // All states are dead ends.
        return numeric_limits<double>::infinity();
    } else {
        return sum / size;
    }
}

bool PatternDatabaseOnline::is_operator_relevant(const OperatorProxy &op) const {
    for (EffectProxy effect : op.get_effects()) {
        int var_id = effect.get_fact().get_variable().get_id();
        if (binary_search(pattern.begin(), pattern.end(), var_id)) {
            return true;
        }
    }
    return false;
}
}
