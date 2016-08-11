#ifndef PDBS_PDB_FACTORY_H
#define PDBS_PDB_FACTORY_H

#include <memory>
#include <string>
#include <vector>

#include "types.h"

class TaskProxy;

namespace pdbs {

class PatternDatabaseInterface;

class PDBFactory {
protected:
    virtual void dump_strategy_specific_options() const = 0;

public:
    PDBFactory() = default;
    virtual ~PDBFactory() = default;
    void dump_options() const;
    
    // Type is shared because, in certain configurations, the factories
    // might want to store a copy of the result. 
    virtual std::shared_ptr<PatternDatabaseInterface> 
	compute_pdb(std::shared_ptr<TaskProxy> task, 
		    const Pattern &pattern, 
		    const std::vector<int> &operator_costs = std::vector<int>()
	    ) = 0;

    virtual std::string name() const = 0;
};
}

#endif
