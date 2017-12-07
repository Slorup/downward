#ifndef PDBS_PATTERN_COLLECTION_GENERATOR_RBP_H
#define PDBS_PATTERN_COLLECTION_GENERATOR_RBP_H

#include "pattern_collection_generator_complementary.h"
#include "types.h"

#include <memory>
#include <vector>
#include <random>


class AbstractTask;

namespace options {
class OptionParser;
class Options;
}
namespace utils {
class CountdownTimer;
}

namespace pdbs3 {
//class PDBFactory;
class PatternCollectionGeneratorRBP : public PatternCollectionGeneratorComplementary {
  unsigned num_vars=0;
	int time_limit=100;
  std::shared_ptr<PatternCollection> patterns;
  
  public:
  
  virtual void initialize(std::shared_ptr<AbstractTask> task) override;
	explicit PatternCollectionGeneratorRBP(const options::Options &options);
  virtual PatternCollectionContainer generate() override;
  //virtual void set_reward(const PatternCollectionContainer & pc, double reward) = 0;
};

}

#endif
