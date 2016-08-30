#ifndef SYMBOLIC_TRANSITION_H
#define SYMBOLIC_TRANSITION_H

#include "../global_operator.h"
#include "sym_variables.h"
//#include "../cudd-2.5.0/cudd/cudd.h"
//#include "../cudd-2.5.0/obj/cuddObj.hh"

#include <set>
#include <vector>

namespace symbolic {

class SymStateSpaceManager;
class SymPDB;
class SymSMAS;
class SimulationRelation;
class DominanceRelation;

/*
 * Represents a symbolic transition.
 * It has two differentiated parts: label and abstract state transitions
 * Label refers to variables not considered in the merge-and-shrink
 * Each label has one or more abstract state transitions
 */
class SymTransition{
  SymVariables * sV; //To call basic BDD creation methods
  int cost; // transition cost
  BDD tBDD; // bdd for making the relprod
  
  std::vector<int> effVars; //FD Index of eff variables. Must be sorted!!
  BDD existsVars, existsBwVars;     // Cube with variables to existentialize
  std::vector<BDD> swapVarsS, swapVarsSp; // Swap variables s to sp and viceversa
  std::vector<BDD> swapVarsA, swapVarsAp; // Swap abstraction variables

  std::set<const GlobalOperator *> ops; //List of operators represented by the TR

  const SymStateSpaceManager * absAfterImage;
 public:
  //Constructor for abstraction transitions
 SymTransition(SymStateSpaceManager * mgr, 
	       const DominanceRelation & sim_relations);

  //Constructor for transitions irrelevant for the abstraction
  SymTransition(SymVariables * sVars, 
		const GlobalOperator * op, int cost_);

  //Copy constructor
  SymTransition(const SymTransition &) = default;

  BDD image(const BDD & from) const;
  BDD preimage(const BDD & from) const;
  BDD image(const BDD & from, int maxNodes) const;
  BDD preimage(const BDD & from, int maxNodes) const;

  void edeletion(SymStateSpaceManager & mgr); //Includes mutex into the transition

  void merge(const SymTransition & t2,
	     int maxNodes);

  //shrinks the transition to another abstract state space (useful to preserve edeletion)
  void shrink(const SymStateSpaceManager & abs, int maxNodes);

  bool setMaSAbstraction(const SymStateSpaceManager & abs,
			 const BDD & bddSrc, const BDD &  bddTarget);

  inline void setAbsAfterImage(const SymStateSpaceManager * abs){
    absAfterImage = abs;
  }

  inline int getCost() const{
    return cost;
  }
  inline int nodeCount() const{
    return tBDD.nodeCount();
  }
  inline const std::set<const GlobalOperator *> &getOps() const {
    return ops;
  }

  inline bool hasOp(std::set<const GlobalOperator *> ops) const {
    for(const auto & op : ops){
      if(ops.count(op)){
	return true;
      }
    }
    return false;
  }

  inline const BDD & getBDD() const {
    return tBDD;
  }

  friend std::ostream & operator<<(std::ostream &os, const SymTransition & tr);
};

}
#endif
