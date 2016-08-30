#ifndef SYMBOLIC_SYM_ENGINE_H
#define SYMBOLIC_SYM_ENGINE_H

#include <vector>
#include <memory>

#include "sym_controller.h"
#include "sym_enums.h"
#include "../search_engine.h"

namespace options{ 
    class Options;
}

namespace symbolic {
class SymHeuristicGenerator;
class SymHeuristic;
class SymExploration;
class SymSolution;
class HNode;
class HTree;

class SymPH;
class SymBDExp;
class SymPruneHeuristic;

class SymEngine : public SearchEngine, public SymController { 
 protected:
  Dir searchDir; //Direction of search in the original state space
  
  //Tree of state spaces
  std::unique_ptr<HTree> tree;

  // List of hierarchy policies to derive new abstractions
  std::vector <SymPH *> phs;
  
  // List of heuristics to use in fw or bw original search. Allows to
  // use precomputed heuristics such as M&S or PDBs. 
  std::vector <SymHeuristicGenerator *> heuristic_generators;
  
  // List of heuristics generated by the heuristic generators
  std::vector <std::shared_ptr<SymHeuristic> > heuristics;

  //Prune heuristic
  std::unique_ptr<SymPruneHeuristic> prune_heuristic;
    
  // List of abstract state spaces. We store a list with the unique
  // pointer so that if we want ever to delete some hNode we just
  // remove it from this list. TODO: maybe we could use
  // shared_pointers instead....
  //std::vector <std::unique_ptr<HNode>> nodes;

  //Variable to keep the current lower bound. Used to know when we have proven an optimal solution.
  int lower_bound;

  //Inherited methods
  virtual void initialize();
  virtual SearchStatus step() = 0;

  //Auxiliar method to get the return code for step()
  int stepReturn() const;
 public:
  SymEngine(const Options &opts);
  virtual ~SymEngine();

  void statistics() const;
  void dump_search_space();
  virtual void new_solution(const SymSolution & sol);

  virtual void setLowerBound(int lower){
    if(lower > lower_bound){
      lower_bound = lower;
      std::cout << "BOUND: " << lower_bound << " < " << bound
		<< ", total time: " << utils::g_timer << std::endl;
    }
  }

  virtual int getUpperBound() const{
    return bound;
  }

   virtual int getLowerBound() const{
    return lower_bound;
  }
  virtual bool solved () const {
    return lower_bound >= bound;
  }

  //virtual SymBDExp * relax(SymBDExp * exp) const;

  static void add_options_to_parser(OptionParser &parser);
  static void set_default_options(Options & opts);
};

}

#endif
