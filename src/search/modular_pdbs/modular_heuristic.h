#ifndef MODULAR_PDBS_HEURISTIC_H
#define MODULAR_PDBS_HEURISTIC_H

//#include "pattern_database.h"

#include "../heuristic.h"
#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_evaluator.h"
#include "pattern_collection_local_search.h"
#include "canonical_symbolic_pdbs.h"
//#include "pdb_factory.h"

class GlobalState;
class State;

namespace options {
class Options;
}

namespace pdbs3 {
  class PDBFactory;
// Implements a heuristic for a single PDB.
class ModularHeuristic : public Heuristic {
    //PatternDatabase pdb;
    PatternCollectionContainer best_collection;
    std::shared_ptr<PatternCollectionGeneratorComplementary> pattern_generator;
    std::shared_ptr<PatternCollectionEvaluator> pattern_evaluator;
    std::shared_ptr<PatternCollectionLocalSearch> pattern_local_search;
    int modular_time_limit;
    int always_CBP_or_RBP_or_UCB;
    bool terminate_creation;
    bool create_perimeter;
    bool only_gamer;
    bool gamer_classic;
    bool gamer_excluded;
    int CBP_counter=0;
    int RBP_counter=0;
    int Disj_counter=0;
    int Not_Disj_counter=0;
    bool doing_local_search=false;
    bool doing_dominated_sampling=false;
    std::shared_ptr<PDBFactory> pdb_factory;
    std::shared_ptr<PatternCollectionInformation> result;
    utils::CountdownTimer *modular_heuristic_timer;
    enum{GAMER_LOCAL_SEARCH=0,GA_LOCAL_SEARCH=1,IPDB_LOCAL_SEARCH=2};
    enum{GAMER_GENERATOR=0,RBP_CBP_GENERATOR=1};
    enum{ALWAYS_CBP=0,ALWAYS_RBP=1,ALWAYS_UCB=2,ALWAYS_50_50=3};
    std::vector<std::shared_ptr<MaxAdditivePDBSubsets> > cleaned_best_pdb_collections; //Store PDB Collections for clean check
    float clear_dominated_in_situ_spent_time=0;
    std::unique_ptr<CanonicalSymbolicPDBs> canonical_pdbs;

protected:
    virtual int compute_heuristic(const GlobalState &global_state) override;
    //virtual void initialize() override;
    /* TODO: we want to get rid of compute_heuristic(const GlobalState &state)
       and change the interface to only use State objects. While we are doing
       this, the following method already allows to get the heuristic value
       for a State object. */
//    int compute_heuristic(const State &state) const;
public:
//    int compute_heuristic_id(size_t state_id);
    /*
      Important: It is assumed that the pattern (passed via Options) is
      sorted, contains no duplicates and is small enough so that the
      number of abstract states is below numeric_limits<int>::max()
      Parameters:
       operator_costs: Can specify individual operator costs for each
       operator. This is useful for action cost partitioning. If left
       empty, default operator costs are used.
    */
    ModularHeuristic(const options::Options &opts);
    void clear_dominated_heuristics_in_situ_sampling(std::shared_ptr<ModularZeroOnePDBs> candidate_ptr);
    virtual ~ModularHeuristic() override = default;
    bool do_local_search (PatternCollectionContainer old_PC);
};
}

#endif
