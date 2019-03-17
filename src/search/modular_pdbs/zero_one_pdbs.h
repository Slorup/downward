#ifndef MODULAR_PDBS_ZERO_ONE_PDBS_H
#define MODULAR_PDBS_ZERO_ONE_PDBS_H

#include "types.h"

class State;
class TaskProxy;

namespace pdbs3 {
class PDBFactory;

class ModularZeroOnePDBs {
    PDBCollection pattern_databases;
    bool dead_end=false;
public:
    ModularZeroOnePDBs(TaskProxy task_proxy, const PatternCollection &patterns, PDBFactory & factory);
    ModularZeroOnePDBs(const std::shared_ptr<PatternDatabaseInterface> Start_PDB);//Mostly to pass to pattern_collection_evaluator::evaluate
    ~ModularZeroOnePDBs() = default;

    int get_value(const State &state) const;
    float fitness=0;
    /*
      Returns the sum of all mean finite h-values of every PDB.
      This is an approximation of the real mean finite h-value of the Heuristic,
      because dead-ends are ignored for the computation of the mean finite
      h-values for a PDB. As a consequence, if different PDBs have different
      states which are dead-end, we do not calculate the real mean h-value for
      these states.
    */
    double compute_approx_mean_finite_h() const;
    void print() const;

    const PDBCollection & get_pattern_databases () const {
	return pattern_databases;
    } 
    void set_fitness(float fitness_) {fitness=fitness_;};
    float get_fitness() const {return fitness;};
		size_t new_pdbs=0;
		size_t get_new_pdbs(){return new_pdbs;};
    size_t get_size(){return pattern_databases.size();};
    bool is_dead_end(const State &state) const;
    bool is_finished() const;
};
}

#endif
