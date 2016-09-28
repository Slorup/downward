#include "canonical_symbolic_pdbs.h"

#include "dominance_pruning.h"
#include "pattern_database.h"

#include "../symbolic/sym_variables.h"
#include "../symbolic/sym_util.h"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>

using namespace std;

namespace pdbs {
CanonicalSymbolicPDBs::CanonicalSymbolicPDBs(
    shared_ptr<PDBCollection> pattern_databases,
    shared_ptr<MaxAdditivePDBSubsets> max_additive_subsets_,
    bool dominance_pruning, int compress_nodes, int compress_time) : valid_cache(pattern_databases->size()), 
								     cache(pattern_databases->size())  {
    assert(max_additive_subsets);

    for (const auto &pdb : *pattern_databases) {
	if (pdb->get_symbolic_variables()) {
	    symbolic_vars = pdb->get_symbolic_variables();
	    break;
	}
    }

    if (dominance_pruning) {
        max_additive_subsets_ = prune_dominated_subsets(*pattern_databases, 
							*max_additive_subsets_);
    }
   
    for (const auto & subset : *max_additive_subsets_) {
	vector <ADD> subsetADDs;
	for (const auto & pdb : subset) {
	    subsetADDs.push_back(pdb->get_ADD());
	}

	if(compress_nodes) {
	    merge(symbolic_vars.get(), subsetADDs, symbolic::mergeSumADD, compress_time, compress_nodes);
	}

	assert(!subsetADDs.empty());
	if(subsetADDs.size() == 1) {
	    singlePDBs.push_back(subsetADDs[0]);
	} else {
	    vector<int> indexes (subsetADDs.size());
	    int i = 0;
	    for (const auto & pdb : subsetADDs) {
		int pos = find(pdbs.begin(), pdbs.end(), pdb) - pdbs.begin();
		if (pos < (int) pdbs.size()) { 
		    indexes [i++] = pos;
		}else {
		    indexes [i++] = pdbs.size();
		    pdbs.push_back(pdb);
		}
	    }
	    max_additive_subsets.push_back(indexes);
	}
    }

    if (compress_nodes) {
	merge(symbolic_vars.get(), singlePDBs, symbolic::mergeMaxADD, compress_time, compress_nodes);
    }

    cout << "Single PDBs: ";
    for (const auto & pdb : singlePDBs) cout << pdb.nodeCount() << " ";
    cout << endl << "Shared PDBs: ";
    for (const auto & pdb : pdbs) cout << pdb.nodeCount() << " ";
    cout << endl;
    cout << "Max additive subsets: " << max_additive_subsets.size() << endl;
   
}

int CanonicalSymbolicPDBs::get_value(const State &state) const {
    // If we have an empty collection, then max_additive_subsets = { \emptyset }.
    assert(!max_additive_subsets->empty());
    valid_cache.reset();
    int * inputs = symbolic_vars->getBinaryDescription(state.get_values());

    int max_h = 0;
    for (const auto & pdb : singlePDBs) {
	int val = Cudd_V(pdb.Eval(inputs).getRegularNode());
	if (val == -1) return numeric_limits<int>::max();

        max_h = max(max_h, val);
    }

    for (const auto &subset : max_additive_subsets) {
        int subset_h = 0;
        for (int pdb : subset) {
	    if (!valid_cache[pdb]) {
		valid_cache.set(pdb);
		cache[pdb] = Cudd_V(pdbs[pdb].Eval(inputs).getRegularNode());
	    }
	    
	    if (cache[pdb] == -1) return numeric_limits<int>::max();

            subset_h += cache[pdb];
        }
        max_h = max(max_h, subset_h);
    }

    return max_h;
}

}
