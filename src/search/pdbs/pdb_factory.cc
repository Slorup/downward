#include "pdb_factory.h"

#include "../options/plugin.h"
#include "pattern_database_interface.h"

#include <iostream>
#include "../utils/timer.h"

using namespace std;

namespace pdbs {

void PDBFactory::dump_options() const {
    cout << "PDB factory: " << name () << endl;
    dump_strategy_specific_options();
}

static options::PluginTypePlugin<PDBFactory> _type_plugin(
    "PDBFactory",
    "This page describes the various pattern database factories supported"
    "by the planner.");



    std::shared_ptr<PatternDatabaseInterface> 
    PDBFactory::compute_pdb(const TaskProxy & task, 
			   const Pattern &pattern, 
			    const std::vector<int> &operator_costs, 
			    double time_limit,
			    double memory_limit) {
	assert(!pattern.empty ());
	num_patterns_requested ++;
	auto item = stored_pdbs.find(PDBKey(pattern, operator_costs));
	if (item != stored_pdbs.end()) {
	    return item->second;
	    // if(!item->second.expired()){
	    // 	return item->second.lock();
	    // } else {
	    // 	num_patterns_regenerated ++;
	    // }
	}

	num_patterns_created ++;
	shared_ptr<PatternDatabaseInterface> result = create_pdb(task, pattern, operator_costs, time_limit,memory_limit); 
        //cout<<"compute_after create,g_timer:"<<utils::g_timer()<<endl;
	stored_pdbs[PDBKey(pattern, operator_costs)] = result;    
        //cout<<"compute_after stored,g_timer"<<utils::g_timer()<<endl;
	
	return result;
    }



    void PDBFactory::statistics() const {
	cout << num_patterns_created << " patterns were generated from which " << num_patterns_regenerated << " were regenerated. " << num_patterns_requested << " patterns provided" << endl;   
    }


}
