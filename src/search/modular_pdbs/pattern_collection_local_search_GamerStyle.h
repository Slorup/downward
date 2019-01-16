#ifndef MODULAR_PDBS_PATTERN_COLLECTION_LOCAL_SEARCH_GAMER_STYLE_H
#define MODULAR_PDBS_PATTERN_COLLECTION_LOCAL_SEARCH_GAMER_STYLE_H

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
class PatternCollectionLocalSearchGamerStyle : public PatternCollectionLocalSearch {
  Pattern current_pattern;
  int time_limit=100;
  int last_var=0;
  std::shared_ptr<TaskProxy> task_proxy;
  //So we avoid re-trying vars which have already been fully tested,
  //forbiddent_vars is cleared whenever we add a variable because then
  //previously discarded vars might actually help the new pattern
  std::set<int> forbidden_vars;
  //utils::CountdownTimer *local_search_timer;
  public:
  explicit PatternCollectionLocalSearchGamerStyle(const options::Options &options);
  //explicit PatternCollectionLocalSearchGamerStyle();
  virtual PatternCollectionContainer generate_next_candidate(PatternCollectionContainer candidate_pattern) override;
  virtual void forbid_last_var() override;
  virtual void reset_forbidden_vars() override;
  virtual void print_last_var() override{std::cout<<"last_var:,"<<last_var<<std::endl;}
  virtual int get_time_limit() override {return time_limit;};
  virtual void print_name() override {std::cout<<"Doing LocalSearchGamerStyle"<<std::endl;}
  virtual std::string get_name() override {std::string output="LocalSearchGamerStyle";return output;};
  std::shared_ptr<ModularZeroOnePDBs> compound_local_search(PatternCollectionContainer old_PC,std::shared_ptr<PatternCollectionInformation> current_result,std::shared_ptr<PDBFactory> pdb_factory);
};

}

#endif
