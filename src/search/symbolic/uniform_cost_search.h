s
#ifndef SYMBOLIC_UNIFORM_COST_SEARCH_H
#define SYMBOLIC_UNIFORM_COST_SEARCH_H
 
#include "unidirectional_search.h"
#include "sym_heuristic.h"
#include "sym_state_space_manager.h"
#include "sym_bucket.h"
#include "sym_astar_closed.h"
#include "sym_estimate.h"
#include "sym_util.h"
#include <vector>
#include <map>
#include <memory>

namespace symbolic {
/*
 * This class allows to perform a BDD search.  It is designed to
 * mantain the current state in the search.  We consider four
 * different points at which we may truncate the search:
 * pop(), filter_mutex(), expand_zero(), expand_cost()
 * We mantain 3 BDDs to know the current state: Stmp, S and Szero.
 * Briefly:
 * 1) if Sfilter, Szero and S are empty => pop() => Szero.
 * 2) else if Stfilter => filter_mutex() => Szero
 * 3) else if Szero => expand_zero => S (passing by Sfilter)
 * 4) else (S must have something) => expand_cost()
 * 
 * Zero cost operators have been expanded iff !S.IsZero() && Szero.IsZero()
 */
class SymController;
class ClosedList;

class UniformCostSearch : public UnidirectionalSearch  {  
    UnidirectionalSearch * parent; //Parent of the search
    
  //Current state of the search:
  std::map<int, Bucket> open_list;

  Bucket Sfilter;   //current g-bucket without duplicates and h-classified (still not filtered mutexes)
  Bucket Smerge;    // bucket before applying merge
  Bucket Szero;     // bucket to expand 0-cost transitions
  Bucket S;         // bucket to expand cost transitions
  int g;            // f and g value of current bucket (S, Szero, Sfilter and Smerge)


  //bucket to store temporary image results in expand_zero() and expand_cost()
  //For each BDD in Szero or S, stores a map with pairs <cost, resImage>
  std::vector<std::map<int, Bucket>> Simg;

  std::unique_ptr<ClosedList> closed;    // Closed list 
  
  SymStepCostEstimation estimationCost, estimationZero;//Time/nodes estimated
  // NOTE: This was used to estimate the time and nodes needed to
  //perform a step in case that the next bucket is still not prepared.
  //Now, we always prepare the next bucket and when that fails no
  //estimation is needed (the exploration is deemed as not searchable
  //and is worse than any other exploration which has its next bucket
  //to expand ready)
  //SymStepCostEstimation estimationDisjCost, estimationDisjZero;
  bool lastStepCost; //If the last step was a cost step (to know if we are in estimationDisjCost or Zero)

  SymController * engine; //Access to the bound and notification of new solutions

  SymExpStatistics stats;

  bool bucketReady() const {
    return !(Szero.empty() && S.empty() && Sfilter.empty() && Smerge.empty());
  }

  inline bool expansionReady() const {
    return Sfilter.empty() && Smerge.empty() && 
      !(Szero.empty() && S.empty());
  }
  virtual bool initialization() const{
    return g==0 && lastStepCost;
  }


  /*
   * Check generated or closed states with other frontiers.  In the
   * original state space we obtain a solution (maybe suboptimal if
   * the states are not closed). 
   */
  void checkCutOriginal(Bucket & bucket, int g);

  void closeStates(Bucket & bucket, int g);

  /*Get the next set from the open list and update g and f.
    Remove duplicate and spurious states. */
  void pop();

  bool prepareBucket(/*int maxTime, int maxNodes, bool afterPop*/);

  
  /* Apply 0-cost operators over Szero. */
  /* Puts the result on Szero. */
  /* Includes the result on S or open (depends on the heuristic value).*/  
  bool expand_zero(int maxTime, int maxNodes);
  
  /* Apply cost-operators over S. */
  /* Insert S on closed. */
  /* Insert successors on open.  */
  bool expand_cost(int maxTime, int maxNodes);

  // Returns the subset with h_value h
  BDD compute_heuristic(const BDD & from, int fVal, int hVal, bool store_eval); 

  void computeEstimation(bool prepare);

  //void debug_pop();
  
  //////////////////////////////////////////////////////////////////////////////
 public:
  UniformCostSearch(SymController * eng, const SymParamsSearch & params);
  UniformCostSearch(const UniformCostSearch & ) = delete;
  UniformCostSearch(UniformCostSearch &&) = default;
  UniformCostSearch& operator=(const UniformCostSearch& ) = delete;
  UniformCostSearch& operator=(UniformCostSearch &&) = default;
  ~UniformCostSearch() {}


  virtual bool finished() const {
      return open_list.empty() && !bucketReady(); 
  }

  virtual bool stepImage(int maxTime, int maxNodes);


  bool init(SymStateSpaceManager * manager, bool fw); //Init forward or backward search



  //Then, relaxFrontier only relaxes the first bucket to expand. 
  //The caller should check if expansion is feasible and useful
  //Finally, all the open list is relaxed to the new abstract state space
  bool relaxFrontier(SymStateSpaceManager * manager, int maxTime, int maxNodes);

  bool relax_open(int , int ){
      std::cerr << "Not implemented relax_open in sym_ucs.h" << std::endl;
      utils::exit_with(utils::ExitCode::UNSUPPORTED);
      return false;
  }
  void relaxClosed();

  virtual void getHeuristic(std::vector<ADD> & heuristics,
			    std::vector <int> & maxHeuristicValues) const ;

  void notifyPrunedBy(int fVal, int gVal);
  void notify(const Bucket & bucket, int fNotClosed = 0); //May prune  
  void notifyNotClosed(int fValue, int hValue);

  void getPossiblyUsefulExplorations(std::vector <UniformCostSearch *> & potentialExps);

  bool isBetter(const UniformCostSearch & other) const;

  UniformCostSearch * getOpposite() const;

  virtual bool isSearchableWithNodes(int maxNodes) const;

  bool isUseful(const std::vector<BDD> & evalStates, 
		       const std::vector<BDD> & newFrontier, 
		       double ratio) const;

  virtual bool isUseful(double ratio) const ;

  /* inline bool isOtherUseful(Bucket & closedAbstract,  */
  /* 			    double ratio) const { */
  /*     assert (!isAbstracted()); // We should only call this method */
  /* 				//over original state space searches */
  /*     double rUseful = ratioUseful(closedAbstract); */
  /*     return rUseful > 0 && rUseful >= ratio  ; */
  /* } */

  double ratioUseful(Bucket & bucket) const;

  // Pointer to the closed list Used to set as heuristic of other explorations.
  inline ClosedList * getClosed() const{
    return closed.get();
  }

  inline SymController * getEngine() const {
      return engine;
  }

  inline const SymStateSpaceManager * get_mgr() const{
    return mgr;
  }

  inline Bucket getSfilter() const{
    return Sfilter;
  }

  inline Bucket getSmerge() const{
    return Smerge;
  }

  inline Bucket getSzero() const{
    return Szero;
  }

  inline Bucket getS() const{
    return S;
  }

  inline void getFrontier(Bucket & res) const {
      for(const BDD & bdd : Sfilter){
	  res.push_back(bdd);
      }
      for(const BDD & bdd : Smerge){
	  res.push_back(bdd);
      }
      for(const BDD & bdd : S){
	  res.push_back(bdd);
      }
      for(const BDD & bdd : Szero){
	  res.push_back(bdd);
      }
  }

  inline int getG() const{
    return g;
  }

  BDD getClosedTotal();

  BDD notClosed();

  void desactivate(); 

  void filterDuplicates(Bucket & bucket);

  virtual long nextStepTime() const;
  virtual long nextStepNodes() const;
  virtual long nextStepNodesResult() const;

  //Returns the nodes that have been expanded by the algorithm (closed without the current frontier)
  BDD getExpanded() const;
  void getNotExpanded(Bucket & res) const;

  //void write(const std::string & file) const;

  void filterMutex (Bucket & bucket) {
      mgr->filterMutex(bucket, fw, initialization());
  }
 private: 
  
  virtual void printFrontier() const;

  int frontierNodes() const{
    if(!Szero.empty()){
      return nodeCount(Szero);
    }else if (!S.empty()){
      return nodeCount(S);
    }else{
      return nodeCount(Sfilter) + nodeCount(Smerge);
    }
  }

  double frontierStates() const{
    if(!Szero.empty()){
      return mgr->stateCount(Szero);
    }else if (!S.empty()){
      return mgr->stateCount(S);
    }else{
      return mgr->stateCount(Sfilter) + mgr->stateCount(Smerge);
    }
  }


  int frontierBuckets() const{
    if(!Szero.empty()){
      return Szero.size();
    }else if (!S.empty()){
      return S.size();
    }else{
      return Sfilter.size() + Smerge.size();
    }
  }

  void violated(TruncatedReason reason , double time, 
		int maxTime, int maxNodes);

  friend std::ostream & operator<<(std::ostream &os, const UniformCostSearch & bdexp);
  
};

}
#endif

