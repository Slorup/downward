#ifndef PDBS_PATTERN_COLLECTION_GENERATOR_GENETIC_ONLINE_SS_H
#define PDBS_PATTERN_COLLECTION_GENERATOR_GENETIC_ONLINE_SS_H

#include "pattern_generator.h"
#include "types.h"
//#include "../state_registry.h"

#include <memory>
#include <vector>
#include "../ss/ss_search.h"
#include "../global_state.h"
#include "group_zero_one_pdbs.h"


class AbstractTask;

namespace options {
class Options;
}
namespace utils {
class CountdownTimer;
}

namespace pdbs {
class PDBFactory;
class ZeroOnePDBs;

struct SS_state
{
  size_t id;
  int g;
  double weight;
};
/*
  Implementation of the pattern generation algorithm by Edelkamp. See:
  Stefan Edelkamp, Automated Creation of Pattern Database Search
  Heuristics. Proceedings of the 4th Workshop on Model Checking and
  Artificial Intelligence (MoChArt 2006), pp. 35-50, 2007.
*/
class PatternCollectionGeneratorGeneticSS : public PatternCollectionGenerator {
    std::shared_ptr<PDBFactory> pdb_factory_candidate;
    std::shared_ptr<PDBFactory> pdb_factory_selected;
    utils::CountdownTimer *genetic_SS_timer;
    vector<SS_state> SS_states_vector;
    map<size_t,pair<int,double> > SS_states;
    // Maximum number of states for each pdb
    int modifier=1;
    int pdb_max_size;
    int num_collections;
    int num_episodes;
    double mutation_probability;
  
    bool bin_packed_episode=false;

    int initial_h=0;
    int sampling_threshold=0;
    double prev_current_collector=0;
    double max_collector=0;
    double overall_pdb_gen_time=0;
    double pdb_gen_time_limit=600;
    double overall_sample_generation_timer=0;
    double overall_online_samp_time=0;
    double overall_probe_time=0;
    int total_online_samples=0;
    int overall_sampled_states=0;
    int current_episode=0;
    double avg_sampled_states=0;
    
    //SS data
    std::set<SSQueue, classcomp> L;
    std::set<SSNode, classcomp2> check;
    TypeSystem * sampler;
    int threshold;

    double sampler_time=0;
    double last_pdb_max_size=0;
    double last_pdb_min_size=0;
    bool last_sampler_too_big=false;
    float min_improvement_ratio=0.1;

    std::shared_ptr<AbstractTask> task;
    /* Specifies whether patterns in each pattern collection need to be disjoint
       or not. */
    bool disjoint_patterns;
    

    // Store best pattern collection over all episodes and its fitness value.
    std::shared_ptr<PatternCollection> best_patterns;
    std::vector<std::vector<std::vector<bool>>> best_pattern_collection;
    double best_fitness;
    // pointer to the heuristic in evaluate from the episode before, used to free memory.
    GroupZeroOnePDBs best_heuristic;
    ZeroOnePDBs *current_heuristic;
    //ZeroOnePDBsHeuristic *current_heuristic;
    double average_operator_cost;
    int num_samples;
    vector<State> samples;
    map<size_t,State> unique_samples;

    
    vector<int> best_heuristic_values;
    std::vector<std::vector<std::vector<bool> > > pattern_collections; // all current pattern collections
    bool best_fitness_was_duplicate;
    std::shared_ptr<std::vector<std::vector<int> > > chosen_pattern_collections;
    set<vector<int> > chosen_patterns;
    //PDBHeuristicOnline *current_heuristic,

    /*
      The fitness values (from evaluate) are used as probabilities. Then
      num_collections many pattern collections are chosen from the vector of all
      pattern collections according to their probabilities. If all fitness
      values are 0, we select uniformly randomly. Note that the set of pattern
      collection where we select from is only changed by mutate, NOT by
      evaluate. All changes there (i.e. transformation and removal of irrelevant
      variables) are just temporary for improved PDB computation.
    */
    void select(const std::vector<double> &fitness_values);

    /*
      Iterate over all patterns and flip every variable (set 0 if 1 or 1 if 0)
      with the given probability from options. This method does not check for
      pdb_max_size or disjoint patterns.
    */
    void mutate();
    //Mutate2 checks for irrelevant vars, otherwise can end up with useless vars in the patter
    int mutate2();
    void transform_to_pattern_bitvector_form(vector<bool> &bitvector,
	const vector<int> &pattern) const ;

    /*
      Transforms a bit vector (internal pattern representation in this class,
      mainly for easy mutation) to the "normal" pattern form vector<int>, which
      we need for ZeroOnePDBsHeuristic.
    */
    Pattern transform_to_pattern_normal_form(
        const std::vector<bool> &bitvector) const;

    /*
      Calculates the mean h-value (fitness value) for each pattern collection.
      For each pattern collection, we iterate over all patterns, first checking
      whether they respect the size limit, then modifying them in a way that
      only causally relevant variables remain in the patterns. Then the zero one
      partitioning pattern collection heuristic is constructed and its fitness
      ( = summed up mean h-values (dead ends are ignored) of all PDBs in the
      collection) computed. The overall best heuristic is eventually updated and
      saved for further episodes.
    */
    void evaluate(std::vector<double> &fitness_values);
    bool is_pattern_too_large(const Pattern &pattern) const;

    /*
      Mark used variables in variables_used and return true iff
      anything was already used (in which case we do not mark the
      remaining variables).
    */
    bool mark_used_variables(const Pattern &pattern,
                             std::vector<bool> &variables_used) const;
    void remove_irrelevant_variables(Pattern &pattern) const;

    /*
      Calculates the initial pattern collections with a next-fit bin packing
      algorithm. Variables are shuffled to obtain a different random ordering
      for every pattern collection. Then, variables are inserted into a pattern
      as long as pdb_max_sizes allows. If the "bin" is full, a new one is
      opened. This is repeated until all variables have been treated. Every
      initial pattern collection contains every variable exactly once and all
      initial patterns of each pattern collection are disjoint (regardless of
      the disjoint_patterns flag).
    */
    void bin_packing();

    /*
      Main genetic algorithm loop. All pattern collections are initialized with
      bin packing. In each iteration (or episode), all patterns are first
      mutated, then evaluated and finally some patterns in each collection are
      selected to be part of the next episode. Note that we do not do any kind
      of recombination.
    */
    void genetic_algorithm(std::shared_ptr<AbstractTask> task);
    double probe_best_only(int threshold);
    int get_pattern_size(Pattern pattern);
public:
    PatternCollectionGeneratorGeneticSS(const options::Options &opts);
    virtual ~PatternCollectionGeneratorGeneticSS() = default;

    virtual PatternCollectionInformation generate(
        std::shared_ptr<AbstractTask> task) override;
    virtual std::shared_ptr<PDBFactory> get_factory () override {
	return pdb_factory_selected;
    }
    static bool compare_SS_states(SS_state i, SS_state j) { return (i.weight>j.weight); }
    static bool compare_pattern_sizes_sort (pair<int ,int > i,pair<int ,int > j) { return (i.second<j.second); }
};
ostream& operator<<(ostream& os, const vector<bool>& v);
template<class T>
ostream& operator<<(ostream& stream, const std::vector<T>& values)
{
    stream << "[ ";
    copy( begin(values), end(values), ostream_iterator<T>(stream, " ") );
    stream << ']';
    return stream;
}
struct delete_empty_vector: public std::unary_function<std::vector<int>, bool>
{
    bool operator() (const vector<int> &a) const {
        return a.size() == 0;
    }
};

//struct delete_empty_vector : public std::unary_function<const std:vector<int>&, bool> 
//{
//    bool operator()(const std::vector<int>& vctPath) const
//    {
 //     if(vctPath.size()>0){
//        return true;
//      }
//      else{
//	return false;
//      }
//    }
//}

}

#endif
