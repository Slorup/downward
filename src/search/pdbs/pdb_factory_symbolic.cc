#include "pdb_factory_symbolic.h"

#include "pattern_database_symbolic.h"

#include "../symbolic/sym_params_search.h"
#include "../symbolic/sym_pdb.h"

#include "../options/option_parser.h"
#include "../options/plugin.h"
#include "../utils/debug_macros.h"
#include "../operator_cost_function.h"


using namespace std;
using namespace symbolic;

namespace pdbs {

    PDBFactorySymbolic::PDBFactorySymbolic(const options::Options & opts) : 
	SymController(opts), 
	generationTime(opts.get<int>("max_time")), 
	generationMemoryGB(opts.get<double>("max_memory_gb")), 
	absTRsStrategy(AbsTRsStrategy(opts.get_enum("tr_st"))),
	dump (opts.get<bool>("dump")) {
	manager = make_shared<OriginalStateSpace>(vars.get(), mgrParams,
						  OperatorCostFunction::get_cost_function());
	manager->init();
    }

    std::shared_ptr<PatternDatabaseInterface> 
PDBFactorySymbolic::create_pdb(const TaskProxy & task, 
		    const Pattern &pattern, 
		    const std::vector<int> &operator_costs){
	

	DEBUG_MSG(cout << "COMPUTE SYMBOLIC PDB" << endl;);
	std::set<int> pattern_set (pattern.begin(), pattern.end()); 
	DEBUG_MSG(cout << "Make copy" << endl;);
	
	assert(manager);
	auto state_space_mgr = make_shared<SymPDB> (manager, absTRsStrategy, pattern_set, 
						    OperatorCostFunction::get_cost_function(operator_costs));
	DEBUG_MSG(cout << "INIT PatternDatabaseSymbolic" << endl;);

	return make_shared<PatternDatabaseSymbolic> (task, pattern, operator_costs,
						     this, vars, state_space_mgr, 
						     searchParams, generationTime, generationMemoryGB);
    }


void PDBFactorySymbolic::dump_strategy_specific_options() const {
    cout << " dump: " << (dump ? "true" : "false") << endl;
}

string PDBFactorySymbolic::name() const {
    return "symbolic";
}


static shared_ptr<PDBFactory>_parse(options::OptionParser &parser) {
    symbolic::SymController::add_options_to_parser(parser, 30e3, 1e7);

    parser.add_option<bool> ("dump", "If set to true, prints the construction time.", "false");
    parser.add_option<int> ("max_time", "Maximum construction time for each PDB.", "1800");
    parser.add_option<double> ("max_memory_gb", "Maximum memory in GB.", "4.0");

    parser.add_enum_option("tr_st", AbsTRsStrategyValues,
                           "Strategy to initialize the transition relation for the abstract state spaces",
			   "IND_TR_SHRINK");


    options::Options options = parser.parse();
    parser.document_synopsis(
        "PDB Factory Symbolic",
        "Symbolic-search PDBS");
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<PDBFactorySymbolic>(options);
}

static options::PluginShared<PDBFactory> _plugin("symbolic", _parse);

}
