#ifndef PDBS_PDB_FACTORY_H
#define PDBS_PDB_FACTORY_H

#include <memory>
#include <string>
#include <vector>
#include <map>
#include <limits>

#include "../utils/system.h"

#include "types.h"

class TaskProxy;


namespace pdbs {

class PatternDatabaseInterface;

class PDBKey {
    Pattern pattern;
    std::vector<int> operator_costs;

public:
    PDBKey (Pattern pattern_, 
	    std::vector<int> operator_costs_) : 
    pattern(pattern_), operator_costs(operator_costs_) {

    }

    bool operator<(const PDBKey& other) const {
	if(pattern.size() < other.pattern.size()) return true;
	else if (pattern.size() > other.pattern.size()) return false;
	for(size_t i = 0; i < pattern.size(); i++) {
	    if (pattern[i] < other.pattern[i]) return true;
	    else if (pattern[i] > other.pattern[i]) return false;
	}

	for(size_t i = 0; i < operator_costs.size(); i++) {
	    if (operator_costs[i] < other.operator_costs[i]) return true;
	    else if (operator_costs[i] > other.operator_costs[i]) return false;
	}

	return false;
    }

    bool operator==(const PDBKey& other) const {
	if(pattern.size() != other.pattern.size()) return false;
	for(size_t i = 0; i < pattern.size(); i++) {
	    if (pattern[i] != other.pattern[i]) return false;
	}
	for(size_t i = 0; i < operator_costs.size(); i++) {
	    if (operator_costs[i] != other.operator_costs[i]) return false;
	}
	return true;
    }
};
class PDBFactory {
    //std::map <PDBKey, std::weak_ptr<PatternDatabaseInterface>> stored_pdbs;
    std::map <PDBKey, std::shared_ptr<PatternDatabaseInterface>> stored_pdbs;

    int num_patterns_created;
    int num_patterns_requested; 
    int num_patterns_regenerated;
protected:
    virtual void dump_strategy_specific_options() const = 0;

    // Type is shared because, in certain configurations, the factories
    // might want to store a copy of the result. 
    virtual std::shared_ptr<PatternDatabaseInterface> 
	create_pdb(const TaskProxy & task, 
		    const Pattern &pattern, 
		   const std::vector<int> &operator_costs = std::vector<int>(),
		   int time_limit = std::numeric_limits<int>::max()
	    ) = 0;    
public:
PDBFactory() : num_patterns_created(0), num_patterns_requested(0), num_patterns_regenerated(0) {}
    virtual ~PDBFactory() = default;
    void dump_options() const;
    
    std::shared_ptr<PatternDatabaseInterface> 
	compute_pdb(const TaskProxy & task, 
		    const Pattern &pattern, 
		    const std::vector<int> &operator_costs = std::vector<int>(), 
		    int time_limit = std::numeric_limits<int>::max()
	    );

    virtual std::string name() const = 0;
    void statistics() const;

    virtual bool is_solved () const {
	return false;
    }

    std::shared_ptr<PDBCollection> terminate_creation (const PDBCollection & pdb_collection) {
	//By default we just make a copy
	return std::make_shared<PDBCollection>(pdb_collection);
    }
};
}

#endif
