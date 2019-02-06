#ifndef SEARCH_ENGINES_EAGER_SEARCH_INTERLEAVED_H
#define SEARCH_ENGINES_EAGER_SEARCH_INTERLEAVED_H

#include "../search_engine.h"

#include "../open_lists/open_list.h"

#include <memory>
#include <vector>

class GlobalOperator;
class Heuristic;
class PruningMethod;
class ScalarEvaluator;

namespace options {
class Options;
}

namespace eager_interleaved_search {
class EagerSearchInterleaved : public SearchEngine {
    const bool reopen_closed_nodes;
    const bool use_multi_path_dependence;

    std::shared_ptr<StateOpenList> open_list;
    ScalarEvaluator *f_evaluator;

    std::vector<Heuristic *> heuristics;
    std::vector<Heuristic *> preferred_operator_heuristics;

    std::shared_ptr<PruningMethod> pruning_method;

    std::pair<SearchNode, bool> fetch_next_node();
    void start_f_value_statistics(EvaluationContext &eval_context);
    void update_f_value_statistics(const SearchNode &node);
    void reward_progress();
    void print_checkpoint_line(int g) const;
    int max_search_time=5;
    int start_time=0;
    int improvement_time=0;
    bool improvement_found=false;
    std::vector<int> last_key_removed;

protected:
    virtual void initialize() override;
    virtual SearchStatus step() override;

public:
    explicit EagerSearchInterleaved(const options::Options &opts);
    virtual ~EagerSearchInterleaved() = default;

    virtual void print_statistics() const override;

    void dump_search_space() const;
    void restart_search();
};
}

#endif
