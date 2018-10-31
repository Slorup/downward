#ifndef MODULAR_PDBS_PATTERN_COLLECTION_GENERATOR_GAMER_H
#define MODULAR_PDBS_PATTERN_COLLECTION_GENERATOR_GAMER_H

#include "pattern_collection_generator_complementary.h"
#include "types.h"

#include <memory>
#include <vector>
#include <random>
#include <set>


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
class PatternCollectionGeneratorGamer : public PatternCollectionGeneratorComplementary {
  int time_limit=100;
  Pattern pattern;
  
  public:
  
  explicit PatternCollectionGeneratorGamer(const options::Options &options);
  virtual PatternCollectionContainer generate() override;
};

}

#endif
