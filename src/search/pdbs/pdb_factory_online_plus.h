#ifndef PDBS_PDB_FACTORY_ONLINE_PLUS_H
#define PDBS_PDB_FACTORY_ONLINE_PLUS_H

#include "pdb_factory.h"

#include "../symbolic/sym_enums.h"
#include "../symbolic/original_state_space.h"
#include "../symbolic/sym_controller.h"

#include <vector>

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


class PDBFactoryOnlinePlus : public PDBFactory, public symbolic::SymController {
	const double generationTime;
	const double generationMemoryGB; 
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

    virtual std::string name() const override;

    virtual bool is_solved () const override {
	return solved();
    }
    
    virtual symbolic::Bucket get_dead_ends() const override;


};
}

#endif

