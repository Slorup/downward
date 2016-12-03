#ifndef PDBS_PATTERN_DATABASE_ONLINE_PLUS_H
#define PDBS_PATTERN_DATABASE_ONLINE_PLUS_H

#include "types.h"
#include "pattern_database_interface.h"
#include "pattern_database_symbolic.h"

#include "../task_proxy.h"

#include <utility>
#include <vector>
#include <map>
#include "match_tree_online.h"
#include "pdb_heuristic.h"
#include "../successor_generator.h"

namespace options {
class Options;
}

namespace pdbs {

typedef int LocalStateID;

class SearchStateInfo {
public:
    bool closed;
    int h;
    int g;

   SearchStateInfo() : closed(false), 
	h(-1), g(-1) {
    }

    int f() {
	return g + h;
    }
};


class SearchInfo {

    struct LocalStateIDSemanticHash {
        const std::vector<State> &state_data_pool;
        LocalStateIDSemanticHash(
            const std::vector<State> &state_data_pool)
            : state_data_pool(state_data_pool)
	    {
        }

        size_t operator()(LocalStateID id) const {
            return state_data_pool[id].hash();
        }
    };

    struct LocalStateIDSemanticEqual {
        const std::vector<State> &state_data_pool;
    LocalStateIDSemanticEqual(const std::vector<State> &state_data_pool) :
	state_data_pool(state_data_pool) {
        }

        bool operator()(LocalStateID lhs, LocalStateID rhs) const {
            return state_data_pool[lhs] == state_data_pool[rhs];
        }
    };

    std::vector<SearchStateInfo> state_info;
    std::vector<State> data_pool;
    std::unordered_set<LocalStateID,  LocalStateIDSemanticHash, LocalStateIDSemanticEqual> idSet;

public:
    SearchInfo() : idSet (0,
			  LocalStateIDSemanticHash(data_pool), 
			  LocalStateIDSemanticEqual(data_pool)) { 
    }

    LocalStateID get_id(const State & state){ 
	LocalStateID id(data_pool.size());
	data_pool.push_back(state);
	auto result = idSet.insert(id);
	bool is_new_entry = result.second;
	if (!is_new_entry) {
	    data_pool.pop_back();
	}
	assert(idSet.size() == data_pool.size());
	return id;	
    }

    SearchStateInfo & get_state_info(LocalStateID id)  {
	if((int)(state_info.size()) <= id) {
	    state_info.resize(id + 1);
	}
	assert(id < (int)(state_info.size()));

	return state_info[id];
    }

    State & get_state(LocalStateID id) {
	assert(id < (int)(data_pool.size()));
	return data_pool[id];
    }

};

class PatternDatabaseOnlinePlus : public PatternDatabaseInterface {
    
    std::shared_ptr<AbstractTask> pdb_task;
    std::vector <Heuristic *> heuristics;
    TaskProxy task_proxy;
    SuccessorGenerator successor_generator;

    std::unique_ptr<PatternDatabaseSymbolic> symbolic_pdb;

    int compute_heuristic(const State & state) const;
    int get_goal_cost(const State & state) const;

public:
    PatternDatabaseOnlinePlus(const TaskProxy &task_proxy, 
			      const Pattern &pattern,
			      const std::vector<int> &operator_costs,
			      std::shared_ptr<AbstractTask> pdb_task,
			      std::vector<Heuristic *> heuristics_, 
			      symbolic::SymController * engine, 
			      std::shared_ptr<symbolic::SymVariables> vars, 
			      std::shared_ptr<symbolic::SymStateSpaceManager> manager, 
			      const symbolic::SymParamsSearch & params, 
			      int generationTime, double generationMemoryGB);

    virtual ~PatternDatabaseOnlinePlus() = default;

    virtual int get_value(const State &state) const override;
    virtual int get_value(const std::vector<int> &state) const override {
	return get_value (State(*pdb_task, std::vector<int>(state)));
    }

    virtual double compute_mean_finite_h() const override;

    virtual const BDD & get_dead_ends() const override{
	assert(symbolic_pdb);
	return symbolic_pdb->get_dead_ends();
    }


};
}

#endif
