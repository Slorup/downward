#ifndef PDBS_PDB_HEURISTIC_ONLINE_H
#define PDBS_PDB_HEURISTIC_ONLINE_H

#include "pattern_database_online.h"

#include "../heuristic.h"
#include <map>

class GlobalState;
class State;

namespace options {
class Options;
}

namespace pdbs {
// Implements a heuristic for a single PDB.
class PDBHeuristicOnline : public Heuristic {
    PatternDatabaseOnline pdb_online;
    // multipliers for each variable for perfect hash function
    std::vector<std::size_t> hash_multipliers;
    int helper_max_size=0;
    std::vector<std::vector<int> > subset_patterns;
    std::vector<std::vector<size_t> >subset_hash_multipliers;
    double overall_extra_helper_gen_time=0;
protected:
    virtual int compute_heuristic(const GlobalState &global_state) override;
    /* TODO: we want to get rid of compute_heuristic(const GlobalState &state)
       and change the interface to only use State objects. While we are doing
       this, the following method already allows to get the heuristic value
       for a State object. */
    int compute_heuristic(const State &state);
    int OnlineDistanceCalculator(const State current_state,std::vector<PDBHeuristic*> &candidate_pdbs_offline,int h_value_to_beat);
    bool backward_search_fully_finished=false;
    std::map<size_t,int> stored_abstract_distance;
    Pattern pattern;
    std::map<size_t,size_t> state_vars_values;
    std::vector<int> operator_costs_copy;
    bool solving_heur=false;
    const CausalGraph &causal_graph;
public:
    /*
      Important: It is assumed that the pattern (passed via Options) is
      sorted, contains no duplicates and is small enough so that the
      number of abstract states is below numeric_limits<int>::max()
      Parameters:
       operator_costs: Can specify individual operator costs for each
       operator. This is useful for action cost partitioning. If left
       empty, default operator costs are used.
    */
    PDBHeuristicOnline(const options::Options &opts);
    virtual ~PDBHeuristicOnline() override = default;
    void get_var_values(size_t set_id);
    size_t get_subset_hash_unoptimized(size_t pdb_helper_index);
    int get_pattern_size(std::vector<int> pattern);
    void set_transformer_subset(std::vector<int> subset_pat);
    std::vector<std::vector<int> >subsets_missing_vars;
    std::vector<std::vector<int> >subsets_missing_index;
    void remove_irrelevant_variables_util(std::vector<int> &pattern);
};
}

#endif
