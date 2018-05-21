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
  //std::shared_ptr<PatternCollection> patterns;
  std::set<int> pattern;
  std::set<int> child_pattern;
  std::vector<int> candidate_pattern;
  double avg_h_val=0;
  //GamerPDBsHeuristic * spdbheuristic;
  //std::shared_ptr <SymStateSpaceManager> state_space;
  //std::unique_ptr <UniformCostSearch> uc_search;
  
  public:
  
  //virtual void initialize(std::shared_ptr<AbstractTask> task) override;
  void initialize();
	explicit PatternCollectionGeneratorGamer(const options::Options &options);
  virtual PatternCollectionContainer generate() override;
  std::vector<int> candidate_vars() const ;
  void confirm(std::set<int> input_pattern);
  void check_improv(double input_avg_h);
  PatternCollectionContainer get_PC();
  //virtual void set_reward(const PatternCollectionContainer & pc, double reward) = 0;
};

}

#endif
