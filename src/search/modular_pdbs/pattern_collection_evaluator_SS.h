#ifndef MODULAR_PDBS_PATTERN_COLLECTION_EVALUATOR_SS_H
#define MODULAR_PDBS_PATTERN_COLLECTION_EVALUATOR_SS_H

#include "pattern_collection_evaluator.h"
#include "types.h"

#include <memory>
#include <vector>
#include <map>
#include <random>
#include "../global_state.h"
#include "../task_proxy.h"
#include "../successor_generator.h"
#include "pattern_collection_information.h"
#include "../ss/ss_search.h"


class AbstractTask;

namespace options {
class OptionParser;
class Options;
}
namespace utils {
class CountdownTimer;
}

namespace pdbs3 {

  struct SS_state
  {
    size_t id;
    int g;
    double weight;
  };
//class PDBFactory;
//class PatternCollectionContainer;
class PatterCollectionEvaluatorSS : public PatternCollectionEvaluator {
  int time_limit=20;
  //unsigned increased_states=0;//We keep this value so we can pass overall progress 
  //std::vector<State> samples;
  std::vector<SS_state> SS_states_vector;
  std::map<size_t,std::pair<int,double> > SS_states;
  //std::shared_ptr<AbstractTask> task;
  std::shared_ptr<TaskProxy> task_proxy;
  std::unique_ptr<SuccessorGenerator> successor_generator;
  std::shared_ptr<PatternCollectionInformation> result;
  utils::CountdownTimer *evaluator_timer;
  std::map<size_t,std::pair<State,int> > unique_samples;
  int increased_states=0;
  //SS data
  std::set<SSQueue, classcomp> L;
  std::set<SSNode, classcomp2> check;
  TypeSystem * sampler;
  unsigned best_initial_value=0;
  std::map<size_t,pair<int,double> > SS_states_copy;
  double overall_probe_time=0;
  //State *initial_state;
  int best_heur_dead_ends=0;
  double total_SS_gen_nodes=0;
  double pruned_states=0;
  int sampled_states=0;
  long total_online_samples=0;
  int sampling_threshold=0;
  double prev_current_collector=0;
  double max_collector=0;
  double node_gen_and_exp_cost=0;
  double current_heur_time_cost=0;
  
  double probe_best_only(std::shared_ptr<PatternCollectionInformation> current_result);
  
  public:
  
  virtual void initialize(std::shared_ptr<AbstractTask> task) override;
  explicit PatterCollectionEvaluatorSS(const options::Options &options);
  virtual bool evaluate(std::shared_ptr<ModularZeroOnePDBs> candidate_PC) override;
  virtual void sample_states(std::shared_ptr<PatternCollectionInformation> current_result) override;
  virtual void clear_dominated_heuristics(std::shared_ptr<PatternCollectionInformation> current_result,std::shared_ptr<PatternCollectionInformation> &new_result) override;
  virtual int get_reward() override {return increased_states;};//So metric is improved_states
  //virtual void set_reward(const PatternCollectionContainer & pc, double reward) = 0;
};

}

#endif
