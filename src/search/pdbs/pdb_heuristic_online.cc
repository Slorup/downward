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
  hash_multipliers=pdb_online.get_hash_multipliers();
  pattern=pdb_online.get_pattern();

 //use smaller PDB as a heuristic if pattern size cross product is big enough to 
 //require it (likely, not for sure, based on size of pattern) 
  if(pdb_online.get_size()<100000){
    helper_max_size=get_pattern_size(pattern)/1;
  }
  else{
    helper_max_size=get_pattern_size(pattern)/1000;
  }
}

int PDBHeuristicOnline::compute_heuristic(const GlobalState &global_state) {
    State state = convert_global_state(global_state);
    return compute_heuristic(state);
}

int PDBHeuristicOnline::compute_heuristic(const State &state) {
    //int h = pdb_online.get_value(state);
    vector<PDBHeuristic*> candidate_pdbs;//empty for now, until online sampling and selection of patterns is added
    int h = OnlineDistanceCalculator(state,candidate_pdbs,0);
    cout<<"h:"<<h<<endl;exit(0);
    if (h == numeric_limits<int>::max())
        return DEAD_END;
    return h;
}

void PDBHeuristicOnline::get_var_values(size_t set_id){
  int temp=0;
  //cout<<endl;
  for (size_t var = 0; var < pattern.size(); ++var) {
        temp = set_id/hash_multipliers[var];
        state_vars_values[pattern[var]] = temp % g_variable_domain[pattern[var]];
	//cout<<"pattern_var:"<<pattern[var]<<",state_vars_values["<<pattern[var]<<"]="<<state_vars_values[pattern[var]]<<endl;
  }
}

int PDBHeuristicOnline::OnlineDistanceCalculator(const State current_state, 
						 vector<PDBHeuristic*> & candidate_pdbs_offline,
						 int h_value_to_beat) {  
  if(pattern.size()==0){
    return 0;
  }
    
  size_t state_index=pdb_online.hash_index(current_state);
    
  if(stored_abstract_distance.find(state_index)!=stored_abstract_distance.end()){//We already know this state ;-)
      return stored_abstract_distance[state_index];
  }
  else if(backward_search_fully_finished){//Not storing all dead_ends, simply returning when state not found and PDB has been fully generated
    //cout<<"initial_state is not stored"<<endl;fflush(stdout);
    return DEAD_END;
  }

    size_t initial_state_index=state_index;

    AdaptiveQueue<pair<size_t,size_t> > pq; // (first implicit entry: priority,) second entry: index for an abstract state
    int initial_h=0;
    if(candidate_pdbs_offline.size()>0){
      get_var_values(state_index);
      //cout<<"candidate_pdbs_ofline.size:"<<candidate_pdbs_offline.size()<<endl;fflush(stdout);
      for(size_t i=0;i<candidate_pdbs_offline.size();i++){
	initial_h=max(candidate_pdbs_offline[i]->compute_heuristic_id(get_subset_hash_unoptimized(i)),initial_h);
      }
    }

    if(initial_h==0){//heuristic has to be admissible!
      if (pdb_online.is_goal_state(state_index)) {//no need to search, this state is abstract goal already!
	//cout<<"initial_state is goal"<<endl;fflush(stdout);
	return 0;
      }
    }

    if(initial_h==INT_MAX/2){
      cout<<"initial_abstract state:,"<<state_index<<", is DEAD_END according to helper_pdb"<<endl;
      return DEAD_END;
    }
    
    pq.push(initial_h,make_pair(state_index,0) );
    map<size_t,pair<int,int> > backtracking_ops;//storing hash_effect and op_cost
    backtracking_ops[state_index]=make_pair(0, 0);

    // Dijkstra loop
    int stored_alternative_cost;
    int alternative_cost;
    int best_stored_goal_distance=INT_MAX;
    size_t goal_state_index=-1;
    int current_boundary=0;
    int counter=0;
    double start_timer=utils::g_timer();
    //START SEARCH FOR GOAL
    while (!pq.empty()) {
      counter++;
	cout<<"counter:,"<<counter<<", memory:,"<<utils::get_current_memory_in_kb()<<endl;fflush(stdout);
      if(counter%1000==0){
	if(utils::g_timer()-start_timer>0.01){
	  //cout<<"\thelper_max_size:"<<helper_max_size<<endl;
	    if(subset_patterns.size()<10&&helper_max_size<=get_pattern_size(pattern)&&helper_max_size<500000){
	      //we are going to add a subset of variables as a pdb_to hopefully speed up things
	      double start_extra_helper_gen_time=utils::g_timer();
	      vector<int> candidate_subset_pattern=pattern;
		//skip, already gen biggest PDB, we should do instead code to generate the full PDB when online search is worse
	      helper_max_size=helper_max_size*5;
	      while(get_pattern_size(candidate_subset_pattern)>helper_max_size&&candidate_subset_pattern.size()>1){
		size_t var_to_remove=rand()%candidate_subset_pattern.size();
		candidate_subset_pattern.erase(candidate_subset_pattern.begin()+var_to_remove); //pattern_collection_helper.erase(pattern_collection_helper.begin());
		remove_irrelevant_variables_util(candidate_subset_pattern);
	      }
	      if(candidate_subset_pattern.size()>1){
		//cout<<"generating extra pdb_helper["<<subset_patterns.size()<<"],subset_par:"<<candidate_subset_pattern<<",mem size:"<<get_pattern_size(candidate_subset_pattern)<<endl;
		Options opts2;
		//opts2.set<TaskProxy *>("task", task);
		opts2.set<int>("cost_type", cost_type);
		opts2.set<vector<int> >("pattern", candidate_subset_pattern);
		opts2.set<vector<int> >("operator_costs", pdb_online.operator_costs);
		PDBHeuristic *pdb_heuristic_helper1=new PDBHeuristic(opts2 );
		candidate_pdbs_offline.push_back(pdb_heuristic_helper1);
		set_transformer_subset(candidate_subset_pattern);
		overall_extra_helper_gen_time+=utils::g_timer()-start_extra_helper_gen_time;
		cout<<"extra,helper["<<subset_patterns.size()<<",helper_max_size:"<<helper_max_size<<",pattern_size:"<<get_pattern_size(pattern)<<",overall_extra_helper_gen_time:"<<overall_extra_helper_gen_time<<endl;
		//if(g_timer()-start_extra_helper_gen_time>0.2){
		//  create_pdb_time_limit(operator_costs_copy,0.2);
		//}
	      }
	      //So we need to at least run another 0.01 secs of search before we build another PDB
	      start_timer=utils::g_timer();
	    }
	  if(solving_heur){
	    //IF USING HEURISTIC AS SOLVER, BETTER TO GENERATE PDB  UNTIL MEMORY USAGE IS TOO GREAT
	    //CURRENTLY NEEDS RE_IMPLEMENTATION, ADD TO TO DO LIST
	    //if((double(utils::get_current_memory_in_kb(false))/1024)<2500){
	    //  cout<<"call_create_pdb_time_limit, memory:"<<get_current_memory_in_kb(false)<<endl;
	    //  vector<int> void_operator_costs;
	    //  create_pdb_time_limit(void_operator_costs,2.0);
	    //}
	    if(subset_patterns.size()<30&&(double(utils::get_current_memory_in_kb())/1024)<3000){
	      double start_extra_helper_gen_time=utils::g_timer();
	      //we are going to add a new pdb_helper to hopefully speed up things
	      vector<int> candidate_subset_pattern=pattern;
	      helper_max_size=helper_max_size*10;
	      cout<<"call_extra_pdb_helper, memory:"<<utils::get_current_memory_in_kb()<<",helper_max_size:"<<helper_max_size<<endl;
	      while(get_pattern_size(candidate_subset_pattern)>=helper_max_size&&candidate_subset_pattern.size()>1){
		size_t var_to_remove=rand()%candidate_subset_pattern.size();
		candidate_subset_pattern.erase(candidate_subset_pattern.begin()+var_to_remove); //pattern_collection_helper.erase(pattern_collection_helper.begin());
		remove_irrelevant_variables_util(candidate_subset_pattern);
	      }
	      if(candidate_subset_pattern.size()>1){
		//cout<<"generating extra pdb_helper["<<subset_patterns.size()<<"],subset_par:"<<candidate_subset_pattern<<",mem size:"<<get_pattern_size(candidate_subset_pattern)<<endl;
		Options opts2;
		//opts2.set<TaskProxy *>("task", task);
		opts2.set<int>("cost_type", cost_type);
		opts2.set<vector<int> >("pattern", candidate_subset_pattern);
		opts2.set<vector<int> >("operator_costs", pdb_online.operator_costs);
		PDBHeuristic *pdb_heuristic_helper1=new PDBHeuristic(opts2);
		candidate_pdbs_offline.push_back(pdb_heuristic_helper1);
		set_transformer_subset(candidate_subset_pattern);
		overall_extra_helper_gen_time+=utils::g_timer()-start_extra_helper_gen_time;
		cout<<"extra,pdb_helper_patterns:"<<subset_patterns.size()<<",overall_extra_helper_gen_time:"<<overall_extra_helper_gen_time<<endl;
	      }
	    }
	    return current_boundary;
	  }
	}
      }
        pair<int, pair<size_t,size_t> > node = pq.pop();
        int current_f = node.first;
	if(h_value_to_beat>0){
	  if(current_f>h_value_to_beat){//past h_value_to_beat so we are finished
	    return current_f;
	  }
	}
        size_t state_index = node.second.first;
        int current_g = node.second.second;
	if(current_f>current_boundary){
	  //cout<<"current_boundary:"<<current_f<<endl;fflush(stdout);
	  current_boundary=current_f;
	}
	cout<<"\t next state:"<<state_index<<",current_g:"<<current_g<<",current_f:"<<current_f<<endl;
	cout<<"\t best_stored_goal_distance:"<<best_stored_goal_distance<<endl;fflush(stdout);
	if(current_f>=best_stored_goal_distance){//if pulling states past or equal to best stored optimal path, we already found optimal path
	  //cout<<"found cached goal is optimal path, returning "<<best_stored_goal_distance<<",expansion_counter:"<<expansion_counter<<",current_f:"<<current_f<<endl;
	  //if(expansion_counter>best_stored_goal_distance){
	    //cout<<"goal distance<expansion_counter!,goal_distance:"<<best_stored_goal_distance<<",expansion_counter:"<<expansion_counter<<endl;
	  //}
	  //cout<<"found cached goal is optimal path, returning "<<best_stored_goal_distance<<endl;
	  while(initial_state_index!=goal_state_index){
	      goal_state_index-= backtracking_ops[goal_state_index].first;
	      if(stored_abstract_distance.find(goal_state_index)!=stored_abstract_distance.end()){//we are finished backtracking
		//cout<<"breaking, goal_state stored seen before"<<endl;
		break;
	      }
	      stored_alternative_cost=best_stored_goal_distance-backtracking_ops[goal_state_index].second;
	      stored_abstract_distance[goal_state_index]=stored_alternative_cost;
	      //cout<<"\tbacktracked_successor:"<<goal_state_index<<",stored_abstract_distance:"<<stored_alternative_cost<<endl;
	  }
	  //cout<<"returning best_stored_goal_distance"<<endl;
	  return best_stored_goal_distance;
	}

	//Checking if state already in hashed states   
	//If it is, check the recored goal is smaller than actual goal, if it is, update
	if(stored_abstract_distance.find(state_index)!=stored_abstract_distance.end()){
	  cout<<"\t state is stored with goal distance:"<<stored_abstract_distance[state_index]<<endl;fflush(stdout);
	  if(stored_abstract_distance[state_index]==DEAD_END){//Skip DEAD_ENDs!
	    //cout<<"\tskipping stored DEAD_END"<<endl;
	    continue;
	  }
	   alternative_cost=current_g+stored_abstract_distance[state_index];
	   if(alternative_cost<best_stored_goal_distance){
	     //cout<<"state_index:"<<state_index<<",maximum goal distance is stored_abstract_distance:"<<stored_abstract_distance[state_index]<<"+current_g:"<<current_g<<"="<<alternative_cost<<endl;
	      best_stored_goal_distance=alternative_cost;
	      goal_state_index=state_index;
	      //NO POINT KEEPING DONIG SEARCH IF GUARANTEED BIGGEST POSSIBLE GOAL DIST IS LESS THAN 
	      //BEST HEURISTIC VALUE ALREADY
	      if(h_value_to_beat>best_stored_goal_distance){
		//cout<<"\t\th_value_to_beat:"<<h_value_to_beat<<",best_stored_goal_distance:"<<best_stored_goal_distance<<",g_timer:"<<g_timer()<<",finished"<<endl;
		return best_stored_goal_distance;
	      }
	   }
	   continue;
	}
	cout<<"\t state_index:"<<state_index<<",g:"<<current_g<<endl;fflush(stdout);
	if(backtracking_ops.find(state_index)!=backtracking_ops.end()){
	  if (current_g > backtracking_ops[state_index].second) {
	    //cout<<"\tskipping state_index at g:"<<current_g<<"it is dup, prev dist:"<<backtracking_ops[state_index].second<<endl;
	      continue;
	  }
	  /*  else if(current_g == backtracking_ops[state_index].second) {
	    if(expanded_depth.find(state_index)!=expanded_depth.end()){
	      if(expanded_depth[state_index]==current_g){
		//cout<<"\tskipping state_index at g:"<<current_g<<"it is dup, already expanded "<<backtracking_ops[state_index].second<<endl;
		continue;//already expanded at current depth
	      }
	    }
	  }*/
	}

        // regress abstract_state
        vector<const AbstractOperatorOnline *> applicable_operators;
        pdb_online.match_tree.get_applicable_operators(state_index, applicable_operators);
	cout<<"\t applicable operators:"<<applicable_operators.size()<<endl;fflush(stdout);
	//expansion_counter++;
	//expanded_depth[state_index]=current_g;//node first expanded at this depth
	//cout<<"\tstate_index added to expanded_depth:"<<current_g<<endl;
	//cout<<"Pulled state_index:"<<state_index<<",g:"<<current_g<<endl;
	//cout<<"\t Applicable operators:"<<applicable_operators.size()<<endl;fflush(stdout);
        for (size_t i = 0; i < applicable_operators.size(); i++) {
	  cout<<"\t\ti:"<<i<<" out of "<<applicable_operators.size()<<endl;fflush(stdout);
            size_t successor = state_index + applicable_operators[i]->get_hash_effect();
	    if(stored_abstract_distance.find(successor)!=stored_abstract_distance.end()){
	      if(stored_abstract_distance[successor]==DEAD_END){//Skip DEAD_ENDs!
		cout<<"Skipping adding stored DEAD_END!"<<endl;fflush(stdout);
		continue;
	      }
	    }
	    applicable_operators[i]->dump(pdb_online.pattern, pdb_online.task_proxy);
	    cout<<"\tSuccessor["<<i<<"]:";fflush(stdout);cout<<successor<<endl;fflush(stdout);
	    cout<<"\thash_effect2:"<<applicable_operators[i]->get_hash_effect()<<endl;fflush(stdout);
	    cout<<",op cost:"<<applicable_operators[i]->get_cost()<<endl;fflush(stdout);
	    int alternative_cost = current_g + applicable_operators[i]->get_cost();
	    cout<<"\talternative_cost:"<<alternative_cost<<",parent_g:"<<current_g;
		
	    int h=0;
	    if(backtracking_ops.find(successor)!=backtracking_ops.end()){
	      //cout<<"found prev state_index:"<<successor<<"with distance:";fflush(stdout);cout<<backtracking_ops[successor].second<<endl;
	      if (alternative_cost < backtracking_ops[successor].second) {
		if(alternative_cost<best_stored_goal_distance){//Found new goal distance
		  if (pdb_online.is_goal_state(successor)) {
		     //cout<<"\tfound shortest goal distance, new shortest goal_state:"<<successor<<",g:"<<alternative_cost<<endl;fflush(stdout);
		     best_stored_goal_distance=alternative_cost;
		     goal_state_index=successor;
		     if(h_value_to_beat>best_stored_goal_distance){
		       //cout<<"\t\th_value_to_beat:"<<h_value_to_beat<<",best_stored_goal_distance:"<<best_stored_goal_distance<<",g_timer:"<<g_timer()<<",finished"<<endl;
		       return best_stored_goal_distance;
		     }
		   }
		}
		//cout<<"successor:"<<successor<<",parent:"<<state_index<<",g:"<<alternative_cost<<endl;
		//cout<<"\top:"<<i<<",alternative_cost:"<<alternative_cost<<"is smaller than stored distance:"<<backtracking_ops[successor].second<<",so updating"<<endl;
		if(candidate_pdbs_offline.size()>0){
		  get_var_values(successor);
		  for(size_t j=0;j<candidate_pdbs_offline.size();j++){
		    h=max(candidate_pdbs_offline[j]->compute_heuristic_id(get_subset_hash_unoptimized(j)),h);
		    //cout<<",pdb_helper["<<j<<"]:"<<candidate_pdbs_offline[j]->compute_heuristic_id(get_subset_hash(successor,j))<<endl;//fflush(stdout);
		  }
		  //cout<<"h:"<<h<<endl;fflush(stdout);
		  if(h==INT_MAX/2){
		    stored_abstract_distance[successor]=DEAD_END;
		    //cout<<"Skipping adding helper_pdb DEAD_END!"<<endl;fflush(stdout);
		    continue;//skipping dead_end
		  }
		}
		else{
		  h=0;
		}
		  pq.push(alternative_cost+h, make_pair(successor,alternative_cost));
		  backtracking_ops[successor]=make_pair(applicable_operators[i]->get_hash_effect(), alternative_cost);//storing hash_effect and op_cost
		  //expanded_depth.erase(successor);
		  //cout<<"\talternative successor:"<<successor<<",new g:"<<alternative_cost<<endl;
		/*if(successor==2||successor==3||successor==0){
		  cout<<"new alternative successor:"<<successor<<",cost:"<<applicable_operators[i]->get_cost()<<endl;;
		}*/
	      }
	      else{
		//cout<<"\top:"<<i<<",dup,prev state_index:"<<successor<<"with distance:"<<backtracking_ops[successor].second<<",alternative_cost:"<<alternative_cost<<endl;
		continue;
	      }
	    }
	    else{
	      //cout<<"\top:"<<i<<",successor_state:"<<successor<<" is new,distance:"<<alternative_cost<<endl;
	      backtracking_ops[successor]=make_pair(applicable_operators[i]->get_hash_effect(), alternative_cost);//storing hash_effect and op_cost
	      //distances_map[successor] = alternative_cost;
	      if(candidate_pdbs_offline.size()>0){
		get_var_values(successor);
		for(size_t j=0;j<candidate_pdbs_offline.size();j++){
		  h=max(candidate_pdbs_offline[j]->compute_heuristic_id(get_subset_hash_unoptimized(j)),h);
		  //cout<<"\th["<<j<<"]:"<<candidate_pdbs_offline[j]->compute_heuristic_id(get_subset_hash(successor,j))<<endl;//fflush(stdout);
		}
		if(h==INT_MAX/2){
		  stored_abstract_distance[successor]=DEAD_END;
		  continue;//Not adding a dead end to queue!
		}
	      }
	      else{
		h=0;
	      }
	    
	    //cout<<"\th:"<<h<<",succesor_f:"<<h+alternative_cost<<endl;fflush(stdout);
	      pq.push(alternative_cost+h, make_pair(successor,alternative_cost));
	      /*  if(successor==2||successor==3||successor==0){
		  cout<<"new successor:"<<successor<<",hash_effect:"<<successor-state_index<<",cost:"<<applicable_operators[i]->get_cost()<<endl;;
	      }*/
	      //cout<<"\tnew successor:"<<successor<<",g:"<<alternative_cost<<endl;
	      if (pdb_online.is_goal_state(successor)) {
	       if(alternative_cost<best_stored_goal_distance){
		 //cout<<"\tgenerated new shortest goal_state:"<<successor<<",g:"<<alternative_cost<<endl;
		 best_stored_goal_distance=alternative_cost;
		 goal_state_index=successor;
		 if(h_value_to_beat>best_stored_goal_distance){
		   //cout<<"\t\th_value_to_beat:"<<h_value_to_beat<<",best_stored_goal_distance:"<<best_stored_goal_distance<<",g_timer:"<<g_timer()<<",finished"<<endl;
		   return best_stored_goal_distance;
		 }
	       }
	      }
	    }
        }
    }
    if(best_stored_goal_distance<INT_MAX){//if we found a stored goal distance, return the optimal one
      cout<<"best_stored_goal_distance:"<<best_stored_goal_distance<<", but queue is finished, so it should be DEAD_END, pls DEBUG ME!!!"<<endl;
      exit(0);
    }
    cout<<"DEAD_END FOUND"<<endl;
    return DEAD_END;
}

size_t PDBHeuristicOnline::get_subset_hash_unoptimized(size_t pdb_helper_index){
  //Make sure get_var_values was called, otherwise not doing correct state
  //only needs to be done once per state
  size_t subset_index=0;
  for (size_t var = 0; var < subset_patterns[pdb_helper_index].size(); ++var) {
    subset_index+=state_vars_values[subset_patterns[pdb_helper_index][var]]*subset_hash_multipliers[pdb_helper_index][var];
    //cout<<"subset_index:"<<subset_index<<",pdb_helper_index:"<<pdb_helper_index;fflush(stdout);
    //cout<<",var:"<<subset_patterns[pdb_helper_index][var]<<",hash_multiplier:"<<subset_hash_multipliers[pdb_helper_index][var]<<",value:"<<subset_patterns[pdb_helper_index][var]<<endl;
  }
  //cout<<",final set_id:"<<subset_index<<endl;fflush(stdout);
  return subset_index;
}

int get_pattern_size(vector<int> pattern){
    // test if the pattern respects the memory limit
    int mem = 1;
    for (size_t j = 0; j < pattern.size(); ++j) {
        int domain_size = g_variable_domain[pattern[j]];
        mem *= domain_size;
    }   
    return mem;
}

void PDBHeuristicOnline::set_transformer_subset(vector<int> subset_pat){
  std::vector<int> subset_missing_vars;
  std::vector<int> subset_missing_index;
  subset_patterns.push_back(subset_pat);
  
    
  vector<size_t> temp_hash_multipliers;
  int temp_num_states=1;
  for (size_t i = 0; i < subset_pat.size(); ++i) {
      temp_hash_multipliers.push_back(temp_num_states);
      temp_num_states *= g_variable_domain[subset_pat[i]];
  }
  subset_hash_multipliers.push_back(temp_hash_multipliers);
  //sort(subset_pat.begin(),subset_pat.end());//just in case it is not sorted already
  std::set_difference(
    pattern.begin(), pattern.end(),
    subset_pat.begin(), subset_pat.end(),
    std::back_inserter( subset_missing_vars )
    );
  subsets_missing_vars.push_back(subset_missing_vars);

  //cout<<"subset pattern:"<<subset_pat<<endl;fflush(stdout);
  //cout<<"subset_missing_vars:"<<subset_missing_vars<<endl;
  //cout<<"transformer_subset online pattern:"<<subset_pat<<",subset_hash_multiplier:"<<subset_hash_multipliers.back()<<endl;
  for(size_t i=0;i<subset_missing_vars.size();i++){
    for(size_t j=0;j<pattern.size();j++){
      if(pattern[j]==subset_missing_vars[i]){
	subset_missing_index.push_back(j);
	break;
      }
    }
  }
  subsets_missing_index.push_back(subset_missing_index);
  //cout<<"subset_missing_index:"<<subset_missing_index<<endl;
}

void PDBHeuristicOnline::remove_irrelevant_variables_util(
    vector<int> &pattern) {
  
  unordered_set<int> in_original_pattern(pattern.begin(), pattern.end());
  unordered_set<int> in_pruned_pattern;

    vector<int> vars_to_check;
    for (size_t i = 0; i < g_goal.size(); ++i) {
        int var_no = g_goal[i].first;
        if (in_original_pattern.count(var_no)) {
            // Goals are causally relevant.
            vars_to_check.push_back(var_no);
            in_pruned_pattern.insert(var_no);
        }
    }

    while (!vars_to_check.empty()) {
        int var = vars_to_check.back();
        vars_to_check.pop_back();
        // A variable is relevant to the pattern if it is a goal variable or if
        // there is a pre->eff arc from the variable to a relevant variable.
        // Note that there is no point in considering eff->eff arcs here.
        const vector<int> &rel = causal_graph.get_eff_to_pre(var);
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

    //if(in_pruned_pattern.size()!=in_original_pattern.size()){
      //cout<<"util_version,found "<<in_original_pattern.size()-in_pruned_pattern.size()<<" irrelevant pattern vars"<<endl;
      //cout<<"util,in_original_pattern.size() "<<in_original_pattern.size()<<",prunned size:"<<in_pruned_pattern.size()<<endl;
    //}
    pattern.assign(in_pruned_pattern.begin(), in_pruned_pattern.end());
    sort(pattern.begin(), pattern.end());
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
