 #include "pdb_factory_explicit.h"

#include "pattern_database.h"


#include "../options/option_parser.h"
#include "../options/plugin.h"


using namespace std;

namespace pdbs {


PDBFactoryExplicit::PDBFactoryExplicit(const options::Options & opts) : 
    dump (opts.get<bool>("dump")) {
}
PDBFactoryExplicit::PDBFactoryExplicit() : 
    dump (false) {
}

    std::shared_ptr<PatternDatabaseInterface> 
PDBFactoryExplicit::create_pdb(const TaskProxy & task, 
		    const Pattern &pattern, 
		    const std::vector<int> &operator_costs){
	return make_shared<PatternDatabase> (task, pattern, dump, operator_costs);
}


void PDBFactoryExplicit::dump_strategy_specific_options() const {
    cout << " dump: " << (dump ? "true" : "false") << endl;
}

string PDBFactoryExplicit::name() const {
    return "explicit";
}



static shared_ptr<PDBFactory>_parse(options::OptionParser &parser) {
    parser.add_option<bool> ("dump", "If set to true, prints the construction time.", "false");
    options::Options options = parser.parse();
    parser.document_synopsis(
        "PDB Factory Explicit",
        "Explicit-search PDBS");
    if (parser.dry_run())
        return nullptr;
    else
        return make_shared<PDBFactoryExplicit>(options);
}

static options::PluginShared<PDBFactory> _plugin("explicit", _parse);

}
