#include "modular_heuristic.h"
#include "types.h"

//#include "pattern_generator.h"
#include "pattern_collection_generator_RBP.h"
#include "pattern_collection_evaluator_RandWalk.h"
#include "pdb_factory_symbolic.h"


#include "../option_parser.h"
#include "../plugin.h"
#include "../task_proxy.h"

#include <limits>
#include <memory>
#include <climits>
#include "zero_one_pdbs.h"
//#include "pdb_factory_symbolic.h"
//#include "pdb_factory.h"

using namespace std;

//PatternCollectionGeneratorComplementary is to be the driver for 8 PDB-based options
//RandomCollectionGeneration: CBP, RBP, CGamer, the one used by *Pommerening et al. 
//Local Search: iPDB, gaPDB, CGamer, VPN, *none, *changing order of patterns in gaPDB mutation. (this one matters to 0-1 greedy cost partitioning).
//GenPDB: Symbolic, Explicit, Online, expressed on pdb_factory class
//PDBEval: AvgH, Random sampling, Stratified sampling, *original iPDB method
//CombPDBs->Canonical, hPO, Max
//CostPartition->*None, Saturated, 0-1 greedy 
//Learning: UCB1 to choose bin packing, pdb size. 
//Re-evaluate: None, RemovedPDBsDominated, Run GHS (note: I think if we run GHS here, it may help to generate PDBs that are complementary to LM-cut)
// And the pseudo-code logic for it:
//While (time < 900 seconds)
//     PC<-RandomCollectionGeneration(MaxSize,CostPartition)
//     setInterestingPCs<-LocalSearch(P,PDBEval,MaxSize,CostPartition)
//     selectedPCs<-SubsetSelection(setInterestingPCs,PDBEval,ComPDBs)
//     generatedPDBs <-GenPDB(selectedPCs)
//     H<-Re-evaluate (H¿generatedPCs)
//     Learning(BinPackingRewards)
//ComPDBs (Hl)
namespace pdbs3 {
//Pattern get_pattern_from_options(const shared_ptr<AbstractTask> task,
//				 const Options &opts) {
//    shared_ptr<PatternGenerator> pattern_generator =
//        opts.get<shared_ptr<PatternGenerator>>("pattern");
//    return pattern_generator->generate(task);
//}

ModularHeuristic::ModularHeuristic(const Options &opts)
    : Heuristic(opts),
    pattern_generator(opts.get<shared_ptr<PatternCollectionGeneratorComplementary>>("patterns")),
    pattern_evaluator(opts.get<shared_ptr<PatternCollectionEvaluator>>("evaluator")),
    pdb_factory (opts.get<shared_ptr<PDBFactory>>("pdb_factory")) {
    TaskProxy task_proxy(*task);
	  //need result here to store final PDB collection
    result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
    //std::shared_ptr<TaskProxy> task_proxy=make_shared<TaskProxy>(task_proxy_temp);
	
    cout<<"initial pdb type:"<<pdb_factory->name()<<endl;
//    shared_ptr<PatternCollectionGeneratorComplementary> pattern_generator =
//        opts.get<shared_ptr<PatternCollectionGeneratorComplementary>>("patterns");
    
//    shared_ptr<PatternCollectionEvaluator> pattern_evaluator =
//        opts.get<shared_ptr<PatternCollectionEvaluator>>("evaluator");
    
    pattern_generator->initialize(task);
    PatternCollectionContainer Initial_collection=pattern_generator->generate();
    Initial_collection.print();
    PatternCollection temp_pc=Initial_collection.get_PC();
    cout<<"temp_pc.size:"<<temp_pc.size()<<endl;
		const State &initial_state = task_proxy.get_initial_state();
    //ModularZeroOnePDBs candidate(task_proxy, Initial_collection.get_PC(), *pdb_factory);
    shared_ptr<ModularZeroOnePDBs> candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, Initial_collection.get_PC(), *pdb_factory);
    cout<<"Initial collection zero-one h value:"<<candidate_ptr->get_value(initial_state)<<endl;
    best_collection=Initial_collection;
    //generate sample states:
    pattern_evaluator->initialize(task);
    //result included in sample_states call because needed for dead_end detection
    //set_dead_ends add dead_ends for symbolic, NEED TO ASK ALVARO ABOUT THIS
    result->set_dead_ends(pdb_factory->get_dead_ends());
    //NOTE: Missing dead_ends for latest heuristic,change sampling so it is aware of this
    pattern_evaluator->sample_states(candidate_ptr,result);
    
    ///DISCUSS WITH ALVARO:adding pdbs to current set if evaluator says new collection is helpful
    //WHEN ADDING THE PDB, terminate_creation makes comparisons biased
    //because new candidate has less time to generate pdb, should we wait for 
    //terminate_pdb till the end of subset selection???
    //result->include_additive_pdbs(candidate_ptr->get_pattern_databases());
    result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));

    }

//void ModularHeuristic::initialize(){
//  PatternCollectionContainer Initial_collection;
  //Initial_collection=
//    shared_ptr<PatternCollectionGeneratorComplementary> pattern_generator =
//        opts.get<shared_ptr<PatternCollectionGeneratorComplementary>>("patterns");
//        pattern_generator->generate(task);
//}


int ModularHeuristic::compute_heuristic(const GlobalState &global_state) {
    State state = convert_global_state(global_state);
    cout<<"Testing modular_heuristic constructor finished"<<endl;
    exit(1);
    return 1;
//   return compute_heuristic(state);
}

//int ModularHeuristic::compute_heuristic(const State &state) const {
//    int h = pdb.get_value(state);
//    if (h == numeric_limits<int>::max())
//        return DEAD_END;
//    return h;
//}

//int ModularHeuristic::compute_heuristic_id(size_t state_id) {
//  //cout<<"calling offline_compute_heuristic_id"<<endl;fflush(stdout);
//  //cout<<"state_id="<<state_id<<",entries:"<<num_states<<endl;fflush(stdout);
//    int h = pdb.distances[state_id];
//    //cout<<"h_offline:"<<h<<endl;fflush(stdout);
//    if (h == numeric_limits<int>::max())
//        return INT_MAX/2;//Better when doing maxes
//    return h;
//}

static Heuristic *_parse(OptionParser &parser) {
    parser.document_synopsis("Pattern database heuristic", "TODO");
    parser.document_language_support("action costs", "supported");
    parser.document_language_support("conditional effects", "not supported");
    parser.document_language_support("axioms", "not supported");
    parser.document_property("admissible", "yes");
    parser.document_property("consistent", "yes");
    parser.document_property("safe", "yes");
    parser.document_property("preferred operators", "no");


    parser.add_option<shared_ptr<PatternCollectionGeneratorComplementary>>(
        "patterns",
        "pattern Collection generation method",
        "modular_rbp");
    parser.add_option<shared_ptr<PatternCollectionEvaluator>>(
        "evaluator",
        "pattern Collection evaluation method",
        "rand_walk");
    parser.add_option<shared_ptr<PDBFactory>>(
        "pdb_factory",
        "See detailed documentation for pdb factories. ",
	      "modular_symbolic");
    
    Heuristic::add_options_to_parser(parser);

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;

    return new ModularHeuristic(opts);
}

static Plugin<Heuristic> _plugin("modular_pdb", _parse);
}
