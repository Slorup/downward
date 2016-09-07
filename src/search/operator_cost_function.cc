#include "operator_cost_function.h"

using namespace std;


std::shared_ptr<OperatorCostFunction> OperatorCostFunction::default_cost_function;

std::shared_ptr<OperatorCostFunction> OperatorCostFunction::get_cost_function(const std::vector<int> & costs) {
    if(!costs.empty()) {
	return make_shared<OperatorCostPredefined> (costs);
    } else if (!OperatorCostFunction::default_cost_function) {
	OperatorCostFunction::default_cost_function = make_shared<OperatorCostConstant> (OperatorCost::NORMAL);	
    }
    return OperatorCostFunction::default_cost_function;
}
