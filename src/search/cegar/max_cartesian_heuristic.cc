#include "max_cartesian_heuristic.h"

#include "abstraction.h"
#include "cost_saturation.h"
#include "utils.h"

#include "../utils/rng.h"

#include <cassert>

using namespace std;

namespace cegar {
MaxCartesianHeuristic::MaxCartesianHeuristic(
    const options::Options &opts,
    std::vector<std::unique_ptr<Abstraction>> &&abstractions,
    int num_orders)
    : Heuristic(opts) {
    subtasks.reserve(abstractions.size());
    refinement_hierarchies.reserve(abstractions.size());
    for (auto &abstraction : abstractions) {
        subtasks.push_back(abstraction->get_task());
        refinement_hierarchies.push_back(abstraction->get_refinement_hierarchy());
    }

    vector<int> indices(abstractions.size());
    iota(indices.begin(), indices.end(), 0);

    for (int order = 0; order < num_orders; ++order) {
        if (order != 0) {
            // Always keep original order in the set of orders.
            g_rng()->shuffle(indices);
        }
        h_maps.push_back(create_additive_h_maps(abstractions, indices));
    }
}

vector<MaxCartesianHeuristic::HMap> MaxCartesianHeuristic::create_additive_h_maps(
    const vector<unique_ptr<Abstraction>> &abstractions,
    const vector<int> &order) const {
    assert(abstractions.size() == order.size());

    vector<int> remaining_costs;
    remaining_costs.reserve(task_proxy.get_operators().size());
    for (OperatorProxy op : task_proxy.get_operators()) {
        remaining_costs.push_back(op.get_cost());
    }

    vector<HMap> h_maps(abstractions.size());
    for (int pos : order) {
        const unique_ptr<Abstraction> &abstraction = abstractions[pos];
        abstraction->set_operator_costs(remaining_costs);
        h_maps[pos] = abstraction->compute_h_map();
        reduce_costs(remaining_costs, abstraction->get_saturated_costs());
    }
    return h_maps;
}

int MaxCartesianHeuristic::compute_sum(
    const vector<State> &local_states, const vector<HMap> &order_h_maps) const {
    int sum_h = 0;
    assert(local_states.size() == order_h_maps.size());
    assert(local_states.size() == refinement_hierarchies.size());
    for (size_t i = 0; i < local_states.size(); ++i) {
        const State &local_state = local_states[i];
        const HMap &h_map = order_h_maps[i];
        const Node *node = refinement_hierarchies[i]->get_node(local_state);
        int value = h_map.at(node);
        assert(value >= 0);
        if (value == INF)
            return DEAD_END;
        sum_h += value;
    }
    assert(sum_h >= 0);
    return sum_h;
}

int MaxCartesianHeuristic::compute_heuristic(const State &state) {
    vector<State> local_states;
    local_states.reserve(subtasks.size());
    for (const shared_ptr<AbstractTask> subtask : subtasks) {
        TaskProxy subtask_proxy(*subtask);
        local_states.push_back(subtask_proxy.convert_ancestor_state(state));
    }
    int max_h = 0;
    for (const vector<HMap> &order_h_maps : h_maps) {
        int sum_h = compute_sum(local_states, order_h_maps);
        max_h = max(max_h, sum_h);
    }
    return max_h;
}

int MaxCartesianHeuristic::compute_heuristic(const GlobalState &global_state) {
    State state = convert_global_state(global_state);
    return compute_heuristic(state);
}
}
