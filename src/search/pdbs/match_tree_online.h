#ifndef PDBS_MATCH_TREE_ONLINE_H
#define PDBS_MATCH_TREE_ONLINE_H

#include "types.h"

#include "../task_proxy.h"

#include <cstddef>
#include <vector>

namespace pdbs {
class AbstractOperatorOnline;

// Successor Generator for abstract operators.
class MatchTreeOnline {
    TaskProxy task_proxy;
    struct Node;
    // See PatternDatabase for documentation on pattern and hash_multipliers.
    Pattern pattern;
    std::vector<size_t> hash_multipliers;
    Node *root;
    void insert_recursive(const AbstractOperatorOnline &op,
                          int pre_index,
                          Node **edge_from_parent);
    void get_applicable_operators_recursive(
        Node *node, const size_t state_index,
        std::vector<const AbstractOperatorOnline *> &applicable_operators) const;
    void dump_recursive(Node *node) const;
public:
    // Initialize an empty match tree.
    //MatchTreeOnline();//empty, needed for first declaration
    void update(const std::vector<int> &pattern_, const std::vector<size_t> &hash_multipliers_);//initialize Matchtree
    MatchTreeOnline(const TaskProxy &task_proxy,
              const Pattern &pattern,
              const std::vector<size_t> &hash_multipliers);
    ~MatchTreeOnline();
    /* Insert an abstract operator into the match tree, creating or
       enlarging it. */
    void insert(const AbstractOperatorOnline &op);

    /*
      Extracts all applicable abstract operators for the abstract state given
      by state_index (the index is converted back to variable/values pairs).
    */
    void get_applicable_operators(
        size_t state_index,
        std::vector<const AbstractOperatorOnline *> &applicable_operators) const;
    void dump() const;
};
}

#endif
