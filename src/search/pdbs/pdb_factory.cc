#include "pdb_factory.h"

#include "../options/plugin.h"

#include <iostream>

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
}

