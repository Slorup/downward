#ifndef PDBS_PDB_FACTORY_ONLINE_PLUS_H
#define PDBS_PDB_FACTORY_ONLINE_PLUS_H

#include "pdb_factory.h"

#include "../symbolic/sym_enums.h"
#include "../symbolic/original_state_space.h"
#include "../symbolic/sym_controller.h"

#include <vector>

class Heuristic;

namespace options {
class OptionParser;
class Options;
}

namespace symbolic {
    class SymSolution;
}

namespace utils {
class RandomNumberGenerator;
}

namespace pdbs {

    class PatternDatabaseOnlinePlus;

class PDBFactoryOnlinePlus : public PDBFactory, public symbolic::SymController {
    const double precomputationTime, precomputationNodes;
    const double terminationTime, terminationNodes;
    const double onlineTime, onlineExpansions;

    const bool use_pdbs_in_online_search;
    const bool online_use_canonical_pdbs;
    const bool online_prune_dominated_pdbs; 

    const bool use_online_during_search;
    const bool dump;
	

        std::shared_ptr<symbolic::OriginalStateSpace> manager;

    protected:
	virtual void dump_strategy_specific_options() const override;
    public:
        explicit PDBFactoryOnlinePlus(const options::Options &options);
	virtual ~PDBFactoryOnlinePlus() override = default;

    // Type is shared because, in certain configurations, the factories
    // might want to store a copy of the result. 
	virtual std::shared_ptr<pdbs::PatternDatabaseInterface> 
	create_pdb(const TaskProxy & task, 
		   const Pattern &pattern, 
		   const std::vector<int> &operator_costs = std::vector<int>(), 
		   double time_limit = std::numeric_limits<int>::max(),
		   double memory_limit=2000
	    );

    void get_heuristics_for (const PatternDatabaseOnlinePlus & pdb, 
			     std::vector<std::shared_ptr<Heuristic>> & heuristics);

    virtual std::string name() const override;

    virtual bool is_solved () const override {
	return solved();
    }
    
    virtual symbolic::Bucket get_dead_ends() const override;

    virtual std::shared_ptr<PDBCollection> terminate_creation (PDBCollection & pdb_collection); 
};
}

#endif

