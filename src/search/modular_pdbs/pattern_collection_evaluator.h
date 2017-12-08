#ifndef PDBS_PATTERN_COLLECTION_EVALUATOR_H
#define PDBS_PATTERN_COLLECTION_EVALUATOR_H

//#include "types.h"
#include "pattern_collection_generator_complementary.h"
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
//class PatternCollectionContainer;
class PatternCollectionEvaluator {
    //std::shared_ptr<PatternCollection> patterns;
  unsigned threshold=20000;
  public:
    virtual void initialize(std::shared_ptr<AbstractTask> task) = 0;
    virtual bool evaluate(const PatternCollectionContainer & best_pc,const PatternCollectionContainer & candidate_pc)=0;
    virtual void sample_states(const PatternCollectionContainer & best_pc,const PatternCollectionContainer & candidate_pc)=0;
    void set_threshold(const unsigned thres){
      threshold=thres;
    }
    unsigned get_threshold(){
      return threshold;
    }
    //virtual void set_reward(const PatternCollectionContainer & pc, double reward) = 0;

};

}

#endif
