#include "pattern_collection_local_search.h"

#include "../option_parser.h"
#include "../plugin.h"
//#include "../task_proxy.h"

#include "../task_tools.h"


    
using namespace std;
namespace pdbs3 {
//    virtual void initialize(std::shared_ptr<AbstractTask> task) {
//    cout << "Manual pattern collection: " << *patterns << endl;
//    return PatternCollectionInformation(task, patterns);
//
//   }

static options::PluginTypePlugin<PatternCollectionLocalSearch> _type_plugin(
"PatterCollectionLocalSearch",
"The various local search algorithms for pattern(s) improvement usable by the complementary_modular_pdbs heurisitic.");
}
