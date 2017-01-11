#include "pdb_factory_online_plus.h"

#include "pattern_database_online_plus.h"
#include "canonical_pdbs_heuristic.h"

#include "../symbolic/sym_params_search.h"
#include "../symbolic/sym_pdb.h"
#include "../symbolic/sym_solution.h"

#include "../tasks/pdb_task.h"

#include "../options/option_parser.h"
#include "../options/plugin.h"
#include "../utils/debug_macros.h"
#include "../utils/system.h"
#include "../operator_cost_function.h"

using namespace std;
using namespace symbolic;

namespace pdbs {

    PDBFactoryOnlinePlus::PDBFactoryOnlinePlus(const options::Options & opts) : 
	SymController(opts), 
	precomputationTime(opts.get<int>("precomputation_time")), 
	precomputationNodes(opts.get<int>("precomputation_nodes")), 
	terminationTime(opts.get<int>("termination_time")), 
	terminationNodes(opts.get<int>("termination_nodes")), 
	onlineTime(opts.get<int>("online_time")), 
	onlineExpansions(opts.get<int>("online_expansions")), 
	use_pdbs_in_online_search (opts.get<bool>("use_pdbs_in_online_search")), 
	online_use_canonical_pdbs (opts.get<bool>("online_use_canonical_pdbs")), 
	online_prune_dominated_pdbs(opts.get<bool>("online_prune_dominated_pdbs")), 
	use_online_during_search(opts.get<bool>("use_online_during_search")), 
	dump (opts.get<bool>("dump")) {
	manager = make_shared<OriginalStateSpace>(vars.get(), mgrParams,
						  OperatorCostFunction::get_cost_function());
    }

    std::shared_ptr<PatternDatabaseInterface> 
    PDBFactoryOnlinePlus::create_pdb(const TaskProxy & task, 
				     const Pattern &pattern, 
				     const std::vector<int> &operator_costs, 
				     double time_limit,
				     double /*memory_lmit*/) {
	
	assert(!pattern.empty());
	assert(!solved());
	DEBUG_MSG(cout << "COMPUTE SYMBOLIC PDB" << endl;);
	std::set<int> pattern_set (pattern.begin(), pattern.end()); 
	DEBUG_MSG(cout << "Pattern: "; for (int v : pattern_set) { cout << " " << v; }cout << endl;);
	
	assert(manager);

	shared_ptr<SymStateSpaceManager> state_space_mgr;
	if (pattern.size() == g_variable_domain.size()) {
	    state_space_mgr = manager;
	} else  {
	    state_space_mgr = make_shared<SymPDB> (*manager, pattern_set, 
						   OperatorCostFunction::get_cost_function(operator_costs).get());
	}
	DEBUG_MSG(cout << "INIT PatternDatabaseSymbolic" << endl;);

	auto pdb_task = make_shared<extra_tasks::PDBTask> (g_root_task(), pattern, operator_costs); 

	DEBUG_MSG(cout << "PDB Task created" << endl;);
	auto new_pdb = 
	    make_shared<PatternDatabaseOnlinePlus> (this, task, pattern, operator_costs, 
						    pdb_task, vars, state_space_mgr, 
						    searchParams, 
						    std::min(precomputationTime, time_limit), 
						    precomputationNodes);

	if(new_pdb->is_finished()) {
	    DEBUG_MSG(cout << "Dead end states discovered: " << new_pdb->get_dead_ends().nodeCount() << endl;);

	    if(!(new_pdb->get_dead_ends()*manager->getInitialState()).IsZero()) {
		cout << "Problem proved unsolvable by: " << *new_pdb << endl;
		utils::exit_with(utils::ExitCode::UNSOLVABLE);
	    }
	    manager->addDeadEndStates(true, new_pdb->get_dead_ends());

	}
	
	return new_pdb;
    }

    std::shared_ptr<PDBCollection> PDBFactoryOnlinePlus::terminate_creation (PDBCollection & pdb_collection) {
	auto result = std::make_shared<PDBCollection> ();
	for(auto & pdb : pdb_collection ) {
	    pdb->terminate_creation(terminationTime, terminationNodes);
	    result->push_back(pdb);
	}
	return result;
    }


    void PDBFactoryOnlinePlus::dump_strategy_specific_options() const {
	cout << " dump: " << (dump ? "true" : "false") << endl;
    }

    string PDBFactoryOnlinePlus::name() const {
	return "symbolic";
    }


    static shared_ptr<PDBFactory>_parse(options::OptionParser &parser) {
	symbolic::SymController::add_options_to_parser(parser, 30e3, 1e7);

	parser.add_option<bool> ("dump", "If set to true, prints the construction time.", "false");
	parser.add_option<int> ("max_time", "Maximum construction time for each PDB.", "1800");
	parser.add_option<double> ("max_memory_gb", "Maximum memory in GB.", "4.0");
	parser.add_option<bool> ("use_pdbs_in_online_search", "Whether smaller PDBs are used by the heuristic.", "false");
	parser.add_option<bool> ("online_use_canonical_pdbs", "Use canonical PDBs or just max PDBs.", "false");
	parser.add_option<bool> ("online_prune_dominated_pdbs", "Prune dominated PDBs in those selected for symbolic online.", "false");
	parser.add_option<int> ("termination_time", "Maximum construction time for the termination phase of each PDB.", "1800");
	parser.add_option<double> ("termination_memory_gb", "Maximum memory in GB for the termination phase.", "4.0");


	options::Options options = parser.parse();
	parser.document_synopsis(
	    "PDB Factory Symbolic",
	    "Symbolic-search PDBS");
	if (parser.dry_run())
	    return nullptr;
	else
	    return make_shared<PDBFactoryOnlinePlus>(options);
    }

    symbolic::Bucket PDBFactoryOnlinePlus::get_dead_ends() const {

	const auto & non_mutex = manager->getNotMutexBDDs(true);
	const auto & non_dead_ends = manager->getNotDeadEnds(true);

	//cout << "Obtaining dead ends from factory: " << non_dead_ends.size() << endl;
	symbolic::Bucket dead_ends;
	dead_ends.reserve(non_dead_ends.size() + non_mutex.size());

	for (const auto & bdd  : non_dead_ends) {
	    dead_ends.push_back(!bdd);
	}

	for (const auto & bdd  : non_mutex) {
	    dead_ends.push_back(!bdd);
	}
    
	return dead_ends;
    }

    void PDBFactoryOnlinePlus::get_heuristics_for (const PatternDatabaseOnlinePlus & pdb, 
						   std::vector<shared_ptr<Heuristic>> & heuristics) {
	if(!use_pdbs_in_online_search) {
	    return;
	}

	vector<shared_ptr<PatternDatabaseInterface>> pdb_collection;
	for(const auto & stored_pdb : stored_pdbs) {
	    if(stored_pdb.first.is_pdb_for(pdb.get_pattern(), pdb.get_operator_costs())) {
		pdb_collection.push_back(stored_pdb.second->get_offline_pdb());
	    }
	}
	if(!pdb_collection.empty()) { 
	    heuristics.push_back(make_shared<CanonicalPDBsHeuristic> (pdb_collection, online_use_canonical_pdbs, online_prune_dominated_pdbs));
	}
    }


    static options::PluginShared<PDBFactory> _plugin("online_symbolic", _parse);

}
