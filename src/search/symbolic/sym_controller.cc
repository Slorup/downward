#include "sym_controller.h"

#include "opt_order.h"

#include "sym_manager.h"
#include "sym_bdexp.h"
#include "sym_hnode.h"
#include "../debug.h" 
#include "../globals.h"
// #include "../merge_and_shrink/ld_simulation.h"
// #include "../merge_and_shrink/abstraction_builder.h"


#include "../option_parser.h"

using namespace std;

namespace symbolic {

SymController::SymController(const Options &opts)
    : vars(new SymVariables(opts)), mgrParams(opts),  searchParams(opts) {
	
    mgrParams.print_options();
    searchParams.print_options();

    //TODO: This should be done before computing the var order and
    //initializing vars. Done here to avoid memory errors
    // if(abstractionBuilder) {
    // 	unique_ptr<LDSimulation> ldSim;
    // 	std::vector<std::unique_ptr<Abstraction> > abstractions;
    // 	// TODO: This irrelevance pruning is only safe for detecting unsolvability
    // 	abstractionBuilder->build_abstraction
    // 	    (true, OperatorCost::ZERO, ldSim, abstractions);
    // 	cout << "LDSimulation finished" << endl;
	
    // 	ldSim->release_memory();
    // }

    vars->init();
}

SymHNode * SymController::createHNode(SymHNode * node,
				      unique_ptr <SymAbstraction> && abs, 
				      SymPH * ph){
    SymHNode * newNode = new SymHNode(node, ph, move (abs));
    nodes.push_back(unique_ptr<SymHNode> (newNode));
    node->addChildren(newNode);
    newNode->addParent(node);
    return newNode;
}


void SymController::add_options_to_parser(OptionParser &parser, int maxStepTime, int maxStepNodes) {
    SymVariables::add_options_to_parser(parser);
    SymParamsMgr::add_options_to_parser(parser);
    SymParamsSearch::add_options_to_parser(parser, maxStepTime, maxStepNodes);
}
}
