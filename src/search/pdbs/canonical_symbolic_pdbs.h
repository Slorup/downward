#ifndef PDBS_CANONICAL_SYMBOLIC_PDBS_H
#define PDBS_CANONICAL_SYMBOLIC_PDBS_H

#include "types.h"
#include "cuddObj.hh"

#include <memory>
#include "../utils/dynamic_bitset.h"
#include <iostream>

class State;

namespace symbolic {
    class SymVariables;
}

namespace pdbs {
class CanonicalSymbolicPDBs {

    std::vector<ADD> pdbs;
    std::vector<std::vector<int> > max_additive_subsets;
    std::vector<ADD> singlePDBs;

    std::vector<BDD> dead_end_detection;    


    mutable utils::DynamicBitset<> valid_cache;
    mutable std::vector<int> cache;


    std::shared_ptr <symbolic::SymVariables> symbolic_vars;

public:
    CanonicalSymbolicPDBs(std::shared_ptr<PDBCollection> pattern_databases,
			  std::shared_ptr<MaxAdditivePDBSubsets> max_additive_subsets,
			  bool dominance_pruning, int compress_nodes, int compress_time);
    ~CanonicalSymbolicPDBs() = default;

    int get_value(const State &state) const;
    virtual int count_pdbs(){
      int count=0;
      for(auto pattern_collection : max_additive_subsets){
	count+=pattern_collection.size();
      }
      return count;
    }

};
}

#endif
