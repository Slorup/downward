#ifndef MODULAR_PDBS_PATTERN_COLLECTION_GENERATOR_COMPLEMENTARY_H
#define MODULAR_PDBS_PATTERN_COLLECTION_GENERATOR_COMPLEMENTARY_H

//#include "pattern_generator.h"
#include "types.h"

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

class PatternCollectionGeneratorComplementary {
    //std::shared_ptr<PatternCollection> patterns;
  public:
    virtual void initialize(std::shared_ptr<AbstractTask> task) = 0;
    virtual PatternCollectionContainer generate() = 0;
    //virtual void set_reward(const PatternCollectionContainer & pc, double reward) = 0;

};
class PatternCollectionContainer {
    PatternCollection pc;
  public:
    void add_pc(Pattern input){
      pc.push_back(input);
    }
    int get_size() const{
      return pc.size();
    }
};

}

#endif
