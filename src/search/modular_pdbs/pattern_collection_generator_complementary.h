#ifndef MODULAR_PDBS_PATTERN_COLLECTION_GENERATOR_COMPLEMENTARY_H
#define MODULAR_PDBS_PATTERN_COLLECTION_GENERATOR_COMPLEMENTARY_H

//#include "pattern_generator.h"
#include "types.h"
#include "../globals.h"
#include "../task_tools.h"

#include <memory>
#include <vector>
#include <random>
#include <iostream>
#include <iterator>



class AbstractTask;

namespace options {
class Options;
}
namespace utils {
class CountdownTimer;
}

namespace pdbs3 {
class PatternCollectionContainer {
    PatternCollection pc;
  public:
    PatternCollectionContainer() = default;
    bool operator==(const PatternCollectionContainer input){ 
      if(this->get_size()!=input.get_size())
	return false;

    PatternCollection pc1=this->get_PC();
    PatternCollection pc2=input.get_PC();
      for(int i=0;i<this->get_size(); i++){
	if(pc1[i]!=pc2[i]){
	  return false;
	}
      }
      return true;//All patterns equal and in the same position, which matters for cost partitioning, so 
      //same patterns in different position make a different PC
    }
    PatternCollectionContainer(PatternCollection Input_pc){
      pc=Input_pc;
    }
    void add_pc(Pattern input){
      pc.push_back(input);
    }
    void restart_pc(Pattern input){//For generators like Gamer, which keep a single PC
      pc.clear();
      pc.push_back(input);
    }
    void clear(){
      pc.clear();
    }
    bool empty(){
      if(pc.size()==0){
        return true;
      }
      std::cout<<"PatternCollectionContainer returns not-empty,NumberOfPatterns="<<pc.size()<<",1stPatternSize:"<<pc.begin()->size()<<std::endl;
      return false;
    }

    int get_size() const{
      return pc.size();
    }
    PatternCollection get_PC() const{
      return pc;
    }
    void pop_back() {
      std::cout<<"complementary_pop_back of pc,inital_size:"<<pc.size()<<",";
      pc.pop_back();
      std::cout<<"final_size:"<<pc.size()<<std::endl;
    }
    void print(){
      std::cout<<"PC";
      for (auto pattern : pc){
        std::cout<<"[";
        for (auto var : pattern){
          std::cout<<var<<",";
        }
        std::cout<<"]";
      }
      std::cout<<std::endl;
    }

    double get_overall_size(){
      if(pc.size()==0){
        return 0;
      }
      
      double overall_size=0;
      for (auto pattern : pc){
        double mem = 1;
        for (size_t j = 0; j < pattern.size(); ++j) {
          //cout<<"g_variable_domain[pattern["<<j<<"]]:"<<g_variable_domain[pattern[j]]<<",mem:"<<mem<<endl;
          double domain_size = g_variable_domain[pattern[j]];
          mem *= domain_size;
        }
        overall_size+=mem;
      }
      return overall_size;
    }
    
    //Moved from pattern_collection_generator_genetic_online so they are available
    //to all pattern generation methods
    void transform_to_pattern_bitvector_form(std::vector<bool> &bitvector,const std::vector<int> &pattern) const {
      bitvector.assign(g_variable_name.size(), false);
      for (size_t i = 0; i < pattern.size(); ++i) {
        bitvector[pattern[i]]=true;
      }
    }

    Pattern get_top_pattern(){
      assert(pc.size()>0);
      return pc[0];
    }
};
//class PatternCollectionContainer;
class PatternCollectionGeneratorRBP;
class PatternCollectionGeneratorBinPackingV1;

class PatternCollectionGeneratorComplementary {
  double max_single_PDB_size=8;
  bool disjunctive_patterns=true;
  unsigned num_vars=0;
  unsigned goals_to_add;
  bool InSituCausalCheck=true;
    //std::shared_ptr<PatternCollection> patterns;
  
  public:
    


  PatternCollectionContainer generate_perimeter(){
      PatternCollectionContainer PC;
      Pattern temp_pattern;
      for (unsigned i=0;i< num_vars;i++){
        temp_pattern.push_back(i);
      }
      std::cout<<"temp_pattern.size:"<<temp_pattern.size()<<",num_vars:"<<num_vars<<std::endl;
      PC.add_pc(temp_pattern);
      return PC;
  }
  void set_pdb_max_size(double pdb_max_size){
    max_single_PDB_size=pdb_max_size;
  }
  double get_pdb_max_size(){
    return max_single_PDB_size;
  }
  void set_disjunctive_patterns(){
    disjunctive_patterns=true;
  }
  void set_goals_to_add(unsigned number_goals){
    goals_to_add=number_goals;
  }
  unsigned get_goals_to_add(){
    return goals_to_add;
  }
  void set_disjunctive_patterns(bool condition){
    disjunctive_patterns=condition;
    //std::cout<<"disjunctive patterns for pdb_size:,"<<max_single_PDB_size<<",set to:,"<<disjunctive_patterns<<std::endl;
  }
  bool get_disjunctive_patterns(){
    return disjunctive_patterns;
  }
  void set_InSituCausalCheck(bool status){
    InSituCausalCheck=status;
  }

    virtual void initialize(std::shared_ptr<AbstractTask> task) {
      num_vars= task->get_num_variables();
      std::cout<<"num_vars:"<<num_vars<<std::endl;
    }

    virtual PatternCollectionContainer generate() = 0;
    //virtual void set_reward(const PatternCollectionContainer & pc, double reward) = 0;
    
    double get_pattern_size(Pattern pattern){
      if(pattern.size()==0){
        return 0;
      }

    double mem = 1;
    for (size_t j = 0; j < pattern.size(); ++j) {
      //cout<<"g_variable_domain[pattern["<<j<<"]]:"<<g_variable_domain[pattern[j]]<<",mem:"<<mem<<endl;
      double domain_size = g_variable_domain[pattern[j]];
      mem *= domain_size;
    }   
    return mem;
    }
    
    //Moved from pattern_collection_generator_genetic_online so they are available
    //to all pattern generation methods
    void transform_to_pattern_bitvector_form(std::vector<bool> &bitvector,const std::vector<int> &pattern) const {
      bitvector.assign(g_variable_name.size(), false);
      for (size_t i = 0; i < pattern.size(); ++i) {
        bitvector[pattern[i]]=true;
      }
    }

    Pattern transform_to_pattern_normal_form(const std::vector<bool> &bitvector) const {
      Pattern pattern;
      for (size_t i = 0; i < bitvector.size(); ++i) {
          if (bitvector[i]){
        pattern.push_back(i);
          }
      }
      return pattern;
    }

    /*
    void PatternCollectionGeneratorGeneticSS::remove_irrelevant_variables( Pattern &pattern) const {
      TaskProxy task_proxy(*task);

      unordered_set<int> in_original_pattern(pattern.begin(), pattern.end());
      unordered_set<int> in_pruned_pattern;

      vector<int> vars_to_check;
      for (FactProxy goal : task_proxy.get_goals()) {
        int var_id = goal.get_variable().get_id();
        if (in_original_pattern.count(var_id)) {
      // Goals are causally relevant.
      vars_to_check.push_back(var_id);
      in_pruned_pattern.insert(var_id);
        }
      }

	while (!vars_to_check.empty()) {
	    int var = vars_to_check.back();
	    vars_to_check.pop_back();
	    *//*
	      A variable is relevant to the pattern if it is a goal variable or if
	      there is a pre->eff arc from the variable to a relevant variable.
	      Note that there is no point in considering eff->eff arcs here.
	    *//*
	    const CausalGraph &cg = task_proxy.get_causal_graph();

	    const vector<int> &rel = cg.get_eff_to_pre(var);
	    for (size_t i = 0; i < rel.size(); ++i) {
		int var_no = rel[i];
		if (in_original_pattern.count(var_no) &&
		    !in_pruned_pattern.count(var_no)) {
		    // Parents of relevant variables are causally relevant.
		    vars_to_check.push_back(var_no);
		    in_pruned_pattern.insert(var_no);
		}
	    }
	}

	pattern.assign(in_pruned_pattern.begin(), in_pruned_pattern.end());
	sort(pattern.begin(), pattern.end());
    }

    bool PatternCollectionGeneratorGeneticSS::is_pattern_too_large(const Pattern &pattern) const {
	// Test if the pattern respects the memory limit.
	TaskProxy task_proxy(*task);
	VariablesProxy variables = task_proxy.get_variables();
	double mem = 1;
	for (size_t i = 0; i < pattern.size(); ++i) {
	    VariableProxy var = variables[pattern[i]];
	    double domain_size = var.get_domain_size();
	    if (!utils::is_product_within_limit(mem, domain_size, pdb_max_size))
		return true;
	    mem *= domain_size;
	}
	return false;
    }

    bool PatternCollectionGeneratorGeneticSS::mark_used_variables(
	const Pattern &pattern, vector<bool> &variables_used) const {
	for (size_t i = 0; i < pattern.size(); ++i) {
	    int var_id = pattern[i];
	    if (variables_used[var_id])
		return true;
	    variables_used[var_id] = true;
	}
	return false;
    }*/
void print_vect_int(const std::vector<int>& values)
{
  std::cout << "[ ";
    for(auto element : values)
      std::cout<<element<<",";
    std::cout << ']';
}
friend PatternCollectionGeneratorRBP;//mostly for max_single_PDB_size
friend PatternCollectionGeneratorBinPackingV1;
virtual unsigned get_RBP_calls(){ return 0;};
virtual unsigned get_CBP_calls(){ return 0;};
};

}
#endif
