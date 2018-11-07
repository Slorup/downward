#ifndef MODULAR_PDBS_PATTERN_COLLECTION_LOCAL_SEARCH_H
#define MODULAR_PDBS_PATTERN_COLLECTION_LOCAL_SEARCH_H

//#include "types.h"
#include "pattern_collection_generator_complementary.h"
//#include "pattern_collection_information.h"
//#include "zero_one_pdbs.h"
#include <memory>
#include <vector>
#include <random>


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
  protected:
  PatternCollectionContainer new_patterns;
  public:
    //virtual void initialize(std::shared_ptr<AbstractTask> task) = 0;
    virtual PatternCollectionContainer generate_next_candidate(PatternCollectionContainer candidate_pattern) = 0;
    virtual void forbid_last_var() = 0;
    virtual void reset_forbidden_vars() = 0;
    virtual void print_last_var() = 0;
    virtual int get_time_limit(){ return 0;};
};

}

#endif
