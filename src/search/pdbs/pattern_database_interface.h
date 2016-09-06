#ifndef PDBS_PATTERN_DATABASE_INTERFACE_H
#define PDBS_PATTERN_DATABASE_INTERFACE_H

#include "types.h"

#include "../task_proxy.h"

#include <utility>
#include <vector>

namespace pdbs {

// Implements a single pattern database
class PatternDatabaseInterface {
protected: 
    TaskProxy task_proxy;

    Pattern pattern;

    std::vector<int> operator_costs;

public:
    /*
      Important: It is assumed that the pattern (passed via Options) is
      sorted, contains no duplicates and is small enough so that the
      number of abstract states is below numeric_limits<int>::max()
      Parameters:
       operator_costs: Can specify individual operator costs for each
       operator. This is useful for action cost partitioning. If left
       empty, default operator costs are used.
    */
    PatternDatabaseInterface(
        const TaskProxy &task_proxy,
        const Pattern &pattern,
        const std::vector<int> &operator_costs = std::vector<int>());

    virtual ~PatternDatabaseInterface() = default;

    // Returns the pattern (i.e. all variables used) of the PDB
    const Pattern &get_pattern() const {
        return pattern;
    }

    // Returns true iff op has an effect on a variable in the pattern.
    bool is_operator_relevant(const OperatorProxy &op) const;

    virtual int get_value(const State &state) const = 0;

    virtual int get_value(const std::vector<int> & state) const = 0;


    // Returns the size (number of abstract states) of the PDB
    virtual std::size_t get_size() const {
	return 0;
	//std::cerr << "Error: method not implemented in this type of PDBs" << std::endl;
	// utils::exit_with(utils::ExitCode::UNSUPPORTED);
    }

    /*
      Returns the average h-value over all states, where dead-ends are
      ignored (they neither increase the sum of all h-values nor the
      number of entries for the mean value calculation). If all states
      are dead-ends, infinity is returned.
      Note: This is only calculated when called; avoid repeated calls to
      this method!
    */
    virtual double compute_mean_finite_h() const = 0;
};
}

#endif
