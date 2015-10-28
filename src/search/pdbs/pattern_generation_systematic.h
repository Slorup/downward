#ifndef PDBS_PATTERN_GENERATION_SYSTEMATIC_H
#define PDBS_PATTERN_GENERATION_SYSTEMATIC_H

#include "../utilities.h"
#include "../task_tools.h"

#include <vector>
#include <unordered_set>

class CanonicalPDBsHeuristic;
class CausalGraph;
class Options;


// Invariant: patterns are always sorted.
typedef std::vector<int> Pattern;

class PatternGenerationSystematicNaive {
    const std::shared_ptr<AbstractTask> task;
    TaskProxy task_proxy;
    std::vector<Pattern> patterns;
public:
    explicit PatternGenerationSystematicNaive(const Options &opts);
    ~PatternGenerationSystematicNaive();
    const std::vector<Pattern> &get_patterns() const {return patterns; }
    CanonicalPDBsHeuristic *get_pattern_collection_heuristic(const Options &opts) const;
};

class PatternGenerationSystematic {
    typedef std::unordered_set<Pattern> PatternSet;

    const std::shared_ptr<AbstractTask> task;
    TaskProxy task_proxy;
    const size_t max_pattern_size;
    std::vector<Pattern> patterns;
    PatternSet pattern_set;  // Cleared after pattern computation.

    void enqueue_pattern_if_new(const Pattern &pattern);
    void compute_eff_pre_neighbors(const CausalGraph &cg,
                                   const Pattern &pattern,
                                   std::vector<int> &result) const;
    void compute_connection_points(const CausalGraph &cg,
                                   const Pattern &pattern,
                                   std::vector<int> &result) const;

    void build_sga_patterns(const CausalGraph &cg);
    void build_patterns();
public:
    explicit PatternGenerationSystematic(const Options &opts);
    ~PatternGenerationSystematic();

    // TODO: Something is wrong with the interface if we need both of
    // the following methods:
    const std::vector<Pattern> &get_patterns() const;
    CanonicalPDBsHeuristic *get_pattern_collection_heuristic(const Options &opts) const;
};

#endif
