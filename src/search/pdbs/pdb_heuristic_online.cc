#include "pdb_heuristic_online.h"

#include "pattern_generator.h"

#include "../option_parser.h"
#include "../plugin.h"
#include "../task_proxy.h"
#include "../priority_queue.h"
#include <climits>
#include "../utils/timer.h"
#include <unordered_set>

#include <limits>
#include <memory>

#include "../causal_graph.h"


using namespace std;

namespace pdbs {
  //Name already taken in pdb_heuristc.cc, hence t he numeral
PatternDatabaseOnline get_pdb_from_options2(const shared_ptr<AbstractTask> task,
                                     const Options &opts) {
    shared_ptr<PatternGenerator> pattern_generator =
        opts.get<shared_ptr<PatternGenerator>>("pattern");
    Pattern pattern = pattern_generator->generate(task);
    TaskProxy task_proxy(*task);
    return PatternDatabaseOnline(task_proxy, pattern, true);
}

PDBHeuristicOnline::PDBHeuristicOnline(const Options &opts)
    : Heuristic(opts),
      pdb_online(get_pdb_from_options2(task, opts)), 
      causal_graph(task_proxy.get_causal_graph())
{
}
int PDBHeuristicOnline::compute_heuristic(const GlobalState &global_state) {
    State state = convert_global_state(global_state);
    return compute_heuristic(state);
}

int PDBHeuristicOnline::compute_heuristic(const State &state) {
    //int h = pdb_online.get_value(state);
    vector<PatternDatabase*> candidate_pdbs;//empty for now, until online sampling and selection of patterns is added
    //int h = OnlineDistanceCalculator(state,candidate_pdbs,0);
    int h = pdb_online.OnlineDistanceCalculator2(state,candidate_pdbs,0);
    //cout<<"h:"<<h<<endl;debug
    if (h == numeric_limits<int>::max())
        return DEAD_END;
    return h;
}

int PDBHeuristicOnline::get_pattern_size(vector<int> pattern){
    // test if the pattern respects the memory limit
    int mem = 1;
    for (size_t j = 0; j < pattern.size(); ++j) {
        int domain_size = g_variable_domain[pattern[j]];
        mem *= domain_size;
    }
    return mem;
}
static Heuristic *_parse(OptionParser &parser) {
  cout<<"hello!"<<endl;
    parser.document_synopsis("Pattern database heuristic", "TODO");
    parser.document_language_support("action costs", "supported");
    parser.document_language_support("conditional effects", "not supported");
    parser.document_language_support("axioms", "not supported");
    parser.document_property("admissible", "yes");
    parser.document_property("consistent", "yes");
    parser.document_property("safe", "yes");
    parser.document_property("preferred operators", "no");

    parser.add_option<shared_ptr<PatternGenerator>>(
        "pattern",
        "pattern generation method",
        "greedy()");
    Heuristic::add_options_to_parser(parser);

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;

    return new PDBHeuristicOnline(opts);
}

static Plugin<Heuristic> _plugin("pdb_online", _parse);
}
