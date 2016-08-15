#ifndef PDBS_PDB_FACTORY_EXPLICIT_H
#define PDBS_PDB_FACTORY_EXPLICIT_H

#include "pdb_factory.h"

#include <vector>

namespace options {
class OptionParser;
class Options;
}

namespace utils {
class RandomNumberGenerator;
}

namespace pdbs {

    class PDBFactoryOnline : public PDBFactory {
	const bool dump;
    protected:
	virtual void dump_strategy_specific_options() const override;
    public:
	explicit PDBFactoryOnline(const options::Options &options);
	virtual ~PDBFactoryOnline() override = default;

    // Type is shared because, in certain configurations, the factories
    // might want to store a copy of the result. 
	virtual std::shared_ptr<pdbs::PatternDatabaseInterface> 
	compute_pdb(const TaskProxy & task, 
		    const Pattern &pattern, 
		    const std::vector<int> &operator_costs = std::vector<int>()
	    ) override;

    virtual std::string name() const override;
};
}

#endif
