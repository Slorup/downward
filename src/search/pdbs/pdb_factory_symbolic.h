#ifndef PDBS_PDB_FACTORY_SYMBOLIC_H
#define PDBS_PDB_FACTORY_SYMBOLIC_H

#include "pdb_factory.h"

#include "../symbolic/sym_enums.h"
#include "../symbolic/original_state_space.h"
#include "../symbolic/sym_controller.h"

#include <vector>

namespace options {
class OptionParser;
class Options;
}

namespace utils {
class RandomNumberGenerator;
}

namespace pdbs {

class PDBFactorySymbolic : public PDBFactory, public symbolic::SymController {
	const int generationTime;
	const double generationMemoryGB; 
	const symbolic::AbsTRsStrategy absTRsStrategy;
        const bool dump;

        std::shared_ptr<symbolic::OriginalStateSpace> manager;

    protected:
	virtual void dump_strategy_specific_options() const override;
    public:
        explicit PDBFactorySymbolic(const options::Options &options);
	virtual ~PDBFactorySymbolic() override = default;

    // Type is shared because, in certain configurations, the factories
    // might want to store a copy of the result. 
	virtual std::shared_ptr<pdbs::PatternDatabaseInterface> 
	compute_pdb(const TaskProxy & task, 
		    const Pattern &pattern, 
		    const std::vector<int> &operator_costs = std::vector<int>()
	    );

    virtual std::string name() const override;
};
}

#endif
