#include "tiebreaking_open_list.h"

#include "open_list.h"

#include "../option_parser.h"
#include "../plugin.h"

#include "../utils/memory.h"

#include <cassert>
#include <deque>
#include <map>
#include <utility>
#include <vector>

using namespace std;


template<class Entry>
class TieBreakingOpenList : public OpenList<Entry> {
    using Bucket = deque<Entry>;

    map<const vector<int>, Bucket> buckets;
    int size;

    vector<ScalarEvaluator *> evaluators;
    /*
      If allow_unsafe_pruning is true, we ignore (don't insert) states
      which the first evaluator considers a dead end, even if it is
      not a safe heuristic.
    */
    bool allow_unsafe_pruning;

    int dimension() const;

protected:
    virtual void do_insertion(EvaluationContext &eval_context,
                              const Entry &entry) override;

public:
    explicit TieBreakingOpenList(const Options &opts);
    virtual ~TieBreakingOpenList() override = default;

    virtual Entry remove_min(vector<int> *key = nullptr) override;
    virtual bool empty() const override;
    virtual void clear() override;
    virtual void get_involved_heuristics(set<Heuristic *> &hset) override;
    virtual bool is_dead_end(
        EvaluationContext &eval_context) const override;
    virtual bool is_reliable_dead_end(
        EvaluationContext &eval_context) const override;
    virtual void load_states(int /*max_number*/, std::shared_ptr<std::vector<Entry> > /*state_list*/ ) override;
};


template<class Entry>
TieBreakingOpenList<Entry>::TieBreakingOpenList(const Options &opts)
    : OpenList<Entry>(opts.get<bool>("pref_only")),
      size(0), evaluators(opts.get_list<ScalarEvaluator *>("evals")),
      allow_unsafe_pruning(opts.get<bool>("unsafe_pruning")) {
}

template<class Entry>
void TieBreakingOpenList<Entry>::do_insertion(
    EvaluationContext &eval_context, const Entry &entry) {
    vector<int> key;
    key.reserve(evaluators.size());
    for (ScalarEvaluator *evaluator : evaluators)
        key.push_back(eval_context.get_heuristic_value_or_infinity(evaluator));

    buckets[key].push_back(entry);
    ++size;
}

template<class Entry>
Entry TieBreakingOpenList<Entry>::remove_min(vector<int> *key) {
    assert(size > 0);
    typename map<const vector<int>, Bucket>::iterator it;
    it = buckets.begin();
    assert(it != buckets.end());
    assert(!it->second.empty());
    --size;
    if (key) {
        assert(key->empty());
        *key = it->first;
    }
    Entry result = it->second.front();
    it->second.pop_front();
    if (it->second.empty())
        buckets.erase(it);
    return result;
}
template<class Entry>
void TieBreakingOpenList<Entry>::load_states(int max_number, std::shared_ptr<vector<Entry> > state_list ){
  size_t states_to_collect=min(max_number,size);
  int counter=0;
  if(states_to_collect<1){
    cerr<<"can't call load_states from open_list when there are no states to collect!!"<<endl;exit(1);
  }
  states_loaded_from_open_list.reset();
  cout<<"startes_to_collect:,"<<states_to_collect<<",states_in_open_list:,"<<size<<endl;
  typename map<const vector<int>, Bucket>::iterator it1;
  for(auto it1=buckets.begin();it1!=buckets.end();it1++){
    for(auto it2=it1->second.begin();it2!=it1->second.end();it2++){
      //cout<<"\tcollected state,f=,"<<counter<<flush<<",state:"<<*it2<<flush<<endl;
      if(state_list->size()>states_to_collect){
	return;
      }
      state_list->push_back(*it2);
    }
    counter++;
  }
}

template<class Entry>
bool TieBreakingOpenList<Entry>::empty() const {
    return size == 0;
}

template<class Entry>
void TieBreakingOpenList<Entry>::clear() {
    buckets.clear();
    size = 0;
}

template<class Entry>
int TieBreakingOpenList<Entry>::dimension() const {
    return evaluators.size();
}

template<class Entry>
void TieBreakingOpenList<Entry>::get_involved_heuristics(
    set<Heuristic *> &hset) {
    for (ScalarEvaluator *evaluator : evaluators)
        evaluator->get_involved_heuristics(hset);
}

template<class Entry>
bool TieBreakingOpenList<Entry>::is_dead_end(
    EvaluationContext &eval_context) const {
    // TODO: Properly document this behaviour.
    // If one safe heuristic detects a dead end, return true.
    if (is_reliable_dead_end(eval_context))
        return true;
    // If the first heuristic detects a dead-end and we allow "unsafe
    // pruning", return true.
    if (allow_unsafe_pruning &&
        eval_context.is_heuristic_infinite(evaluators[0]))
        return true;
    // Otherwise, return true if all heuristics agree this is a dead-end.
    for (ScalarEvaluator *evaluator : evaluators)
        if (!eval_context.is_heuristic_infinite(evaluator))
            return false;
    return true;
}

template<class Entry>
bool TieBreakingOpenList<Entry>::is_reliable_dead_end(
    EvaluationContext &eval_context) const {
    for (ScalarEvaluator *evaluator : evaluators)
        if (eval_context.is_heuristic_infinite(evaluator) &&
            evaluator->dead_ends_are_reliable())
            return true;
    return false;
}

TieBreakingOpenListFactory::TieBreakingOpenListFactory(const Options &options)
    : options(options) {
}

unique_ptr<StateOpenList>
TieBreakingOpenListFactory::create_state_open_list() {
    return utils::make_unique_ptr<TieBreakingOpenList<StateOpenListEntry>>(options);
}

unique_ptr<EdgeOpenList>
TieBreakingOpenListFactory::create_edge_open_list() {
    return utils::make_unique_ptr<TieBreakingOpenList<EdgeOpenListEntry>>(options);
}

static shared_ptr<OpenListFactory> _parse(OptionParser &parser) {
    parser.document_synopsis("Tie-breaking open list", "");
    parser.add_list_option<ScalarEvaluator *>("evals", "scalar evaluators");
    parser.add_option<bool>(
        "pref_only",
        "insert only nodes generated by preferred operators", "false");
    parser.add_option<bool>(
        "unsafe_pruning",
        "allow unsafe pruning when the main evaluator regards a state a dead end",
        "true");
    Options opts = parser.parse();
    opts.verify_list_non_empty<ScalarEvaluator *>("evals");
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<TieBreakingOpenListFactory>(opts);
}

static PluginShared<OpenListFactory> _plugin("tiebreaking", _parse);
