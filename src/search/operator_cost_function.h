#ifndef OPERATOR_COST_FUNCTION_H
#define OPERATOR_COST_FUNCTION_H

#include "globals.h"
#include "operator_cost.h"
#include "global_operator.h"

class OperatorCostFunction {
public:
    virtual int get_adjusted_cost(int op_id) const = 0;
};

class OperatorCostConstant : public OperatorCostFunction {
    OperatorCost cost_type;
public:
    
OperatorCostConstant(OperatorCost c = OperatorCost::NORMAL) : cost_type(c) {}
    virtual int get_adjusted_cost(int op_id) const override {
	return get_adjusted_action_cost(g_operators[op_id].get_cost(), cost_type);
    }	
};

class OperatorCostPredefined : public OperatorCostFunction{
    std::vector<int> costs;
public:
    
    OperatorCostPredefined(const std::vector<int> & c) : costs(c) {}

    virtual int get_adjusted_cost(int op_id) const override {
	return costs[op_id];
    }	
};

#endif
