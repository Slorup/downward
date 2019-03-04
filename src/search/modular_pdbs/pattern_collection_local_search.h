#ifndef MODULAR_PDBS_PATTERN_COLLECTION_LOCAL_SEARCH_H
#define MODULAR_PDBS_PATTERN_COLLECTION_LOCAL_SEARCH_H

#include "types.h"
#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_information.h"
#include "pattern_collection_evaluator.h"
#include "zero_one_pdbs.h"
#include <memory>
#include <vector>
#include <set>
#include <random>
#include "../globals.h"
#include "../task_tools.h"


class AbstractTask;

namespace options {
class Options;
}
namespace utils {
class CountdownTimer;
}

namespace pdbs3 {
//class PDBFactory;
class PatternCollectionContainer;
class PatternCollectionLocalSearch {
  std::set<int> initial_remaining_goal_vars;
  std::shared_ptr<TaskProxy> task_proxy;
  protected:
  PatternCollectionContainer new_patterns;
  int num_vars;
  public:
    void initialize(std::shared_ptr<AbstractTask> task);
    //virtual void initialize(std::shared_ptr<AbstractTask> task) = 0;
    virtual PatternCollectionContainer generate_next_candidate(PatternCollectionContainer candidate_pattern) = 0;
    virtual void forbid_last_var(){};
    virtual void reset_forbidden_vars(){};
    virtual void print_last_var(){};
    virtual int get_time_limit(){ return 0;};
    virtual int get_episodes(){ return 0;};
    void remove_irrelevant_variables(Pattern &pattern,Pattern &removed_vars);
    virtual void print_name() = 0;
    virtual std::string get_name() = 0;
    virtual bool do_local_search(std::shared_ptr<PatternCollectionInformation> current_result, 
	std::shared_ptr<PatternCollectionEvaluator> evaluation_method,
	std::shared_ptr<PDBFactory> pdb_factory) = 0;

};

}

#endif
