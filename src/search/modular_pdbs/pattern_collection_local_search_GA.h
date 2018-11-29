#ifndef MODULAR_PDBS_PATTERN_COLLECTION_LOCAL_SEARCH_GA__H
#define MODULAR_PDBS_PATTERN_COLLECTION_LOCAL_SEARCH_GA__H

#include "pattern_collection_local_search.h"
#include "types.h"

#include <memory>
#include <vector>
#include <set>
#include <map>
#include <random>
#include "../global_state.h"
#include "../task_proxy.h"
#include "../successor_generator.h"
#include "pattern_collection_information.h"


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
//class PatternCollectionContainer;
class PatternCollectionLocalSearchGA : public PatternCollectionLocalSearch {
  int last_var=0;
  std::shared_ptr<TaskProxy> task_proxy;
  //So we avoid re-trying vars which have already been fully tested,
  //forbiddent_vars is cleared whenever we add a variable because then
  //previously discarded vars might actually help the new pattern
  std::set<int> forbidden_vars;
  float mutation_rate=0.2;
  int time_limit=40;
  int episodes=60;
  public:
  PatternCollectionLocalSearchGA(const options::Options & opts);
  virtual PatternCollectionContainer generate_next_candidate(PatternCollectionContainer candidate_pattern) override;
  int mutate(PatternCollectionContainer& candidate_collection);
  void transform_to_pattern_bit_vector_form(std::vector<bool> &bitvector,const std::vector<int> pattern);
  void transform_to_pattern_int_vector_form(const std::vector<bool> &bitvector, std::vector<int> &pattern);
  virtual void print_name() override {std::cout<<"Doing LocalSearchGA"<<std::endl;}
  virtual std::string get_name() override {std::string output="LocalSearchGA";return output;};
  virtual int get_time_limit() override { return time_limit;};
  virtual int get_episodes() override { return episodes;};
};

}

#endif
