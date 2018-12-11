//#include "pattern_collection_generator_complementary.h"
#include "pattern_collection_evaluator_SS.h"
#include "../sampling.h"
#include "../utils/timer.h"
#include "../task_tools.h"
#include "../utils/countdown_timer.h"
#include "../options/option_parser.h"
#include "../options/plugin.h"
#include "pattern_database_interface.h"


//#include "../causal_graph.h"
//#include "../globals.h"
//#include "../task_proxy.h"

//#include "../utils/markup.h"
#include "../utils/math.h"
#include "../utils/rng.h"
//#include "../utils/timer.h"
//#include "../heuristic.h"

//#include <algorithm>
//#include <cassert>
//#include <iostream>
//#include <unordered_set>
//#include <vector>
//#include <math.h>
//Hack to use SS get_type, it needs heuristic object in constructor
#include "../heuristics/blind_search_heuristic.h"
//#include "../heuristics/lm_cut_heuristic.h"
#include "pdb_factory.h"
//#include "pattern_database_interface.h"
#include "../utils/debug_macros.h"
//#include <random>
//#include "../sampling.h"

    
using namespace std;
//PatternCollectionEvaluator is to be the driver for 8 PDB-based options
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
PatterCollectionEvaluatorSS::PatterCollectionEvaluatorSS(const options::Options & opts) :
	time_limit (opts.get<int>("time_limit")),
	time_or_size_selection(opts.get<int>("time_or_size_selection")){
    cout<<"hello EvaluatorSS"<<endl;fflush(stdout);
  //num_vars=task->get_num_variables();
}
  void PatterCollectionEvaluatorSS::initialize(std::shared_ptr<AbstractTask> task) {
    set_threshold(0);
    int num_vars= task->get_num_variables();
    cout<<"num_vars:"<<num_vars<<flush<<endl;
    TaskProxy task_proxy_temp(*task);
    task_proxy=make_shared<TaskProxy>(task_proxy_temp);
    successor_generator=utils::make_unique_ptr<SuccessorGenerator>(task);
    cout<<"succesor generator created"<<flush<<endl;
   //Need a heuristic for SS h-basedtype system,using blind as a place holder
    const State &initial_state = task_proxy->get_initial_state();
    //Do A* related timings for predictions
    if(time_or_size_selection==TIME_SELECTION){
      cout<<"SS_based TIME_SELECTION,starting timings"<<flush<<endl;
      utils::Timer node_gen_and_exp_timings;
      int node_counter=0;
      State succ_state = initial_state;
      vector<State> succ_states;
      set<State> visited_states;
      vector<OperatorProxy> applicable_ops; 
      while(node_gen_and_exp_timings()<1.0){
	node_counter++;
	applicable_ops.clear();
	successor_generator->generate_applicable_ops(succ_state, applicable_ops);
	for(auto op : applicable_ops){
	  succ_states.push_back(succ_state.get_successor(op));
	}
	//const OperatorProxy &op = applicable_ops[rand()%applicable_ops.size()];//choose operator at random
	if(applicable_ops.size()==0){//dead_end,go back to initial state
	  succ_state = task_proxy->get_initial_state();
	  node_counter--;
	  continue;
	}
	else{
	  //cout<<"\t\tapplicable_ops:"<<applicable_ops.size()<<endl;
	  succ_state = succ_state.get_successor(applicable_ops[rand()%applicable_ops.size()]);
	  visited_states.insert(succ_state);
	}
      }
      visited_states.clear();
      node_gen_and_exp_cost=node_gen_and_exp_timings.stop()/double(node_counter);
      cout<<"node gen_and_exp_cost:"<<node_gen_and_exp_cost<<",node_counter:"<<node_counter<<endl;
    }
    else{
      cout<<"SS_based SIZE_SELECTION"<<endl;
    }
  }
  bool PatterCollectionEvaluatorSS::evaluate(std::shared_ptr<ModularZeroOnePDBs> candidate_PC){
    //cout<<"ss, evaluate,candidate_PC.size:"<<candidate_PC->get_size()<<endl;
    //Update SS_states if necessary and get how many states are improved by the candidate_PC
    pruned_states=0;
    //fitness=0;
    sampled_states=0;
    increased_states=0;
    g_rng()->shuffle(SS_states_vector);
    double start_sampler_time=utils::g_timer();
    vector<SS_state>::iterator SS_iter;

    DEBUG_MSG(cout<<"time;,"<<utils::g_timer()<<",SS_states_vector.size:"<<SS_states_vector.size()<<flush<<endl;);
	
    //current_heur_initial_value=candidate_PC->get_value(initial_state);
    total_online_samples++;
    total_SS_gen_nodes=0;
    for(SS_iter=SS_states_vector.begin();SS_iter!=SS_states_vector.end();){
	//cout<<"time:,"<<utils::g_timer<<",working on state:"<<sampled_states<<endl;
	if(unique_samples.find(SS_iter->id)==unique_samples.end()){
	    cerr<<"state not in unique_samples!!!"<<endl;exit(22);
	}

	if(sampled_states%1000==0){
	    if(pruned_states==0){
		if(utils::g_timer()-start_sampler_time>0.2){
		    //cout<<"\tcurrent_episode:"<<current_episode<<",exiting candidate vs best_heuristic SS_states comparison, 0.2 secs iterating without a single better h value"<<endl;
		    break;
		}
	    }
	    else if(pruned_states>0){
		if(utils::g_timer()-start_sampler_time>0.5){
		  cout<<"\t exiting candidate vs best_heuristic SS_states comparison, spent max 0.5 secs,sampled states:"<<sampled_states<<endl;
		    //cout<<"\t exiting candidate vs best_heuristic SS_states comparison, spent max 0.5 secs"<<endl;
		    break;
		}
	    }
	}
	//int best_h=get_best_value(unique_samples.at(SS_iter->id).first);
	int best_h=unique_samples.at(SS_iter->id).second;
	if(best_h==numeric_limits<int>::max()){
	    best_heur_dead_ends++;
	    SS_states.erase(SS_iter->id);
	    SS_iter=SS_states_vector.erase(SS_iter);
	    //cout<<"\tstate is already known dead_end, removing"<<endl;
	    continue;
	}
	else if(best_h+SS_iter->g>sampling_threshold){
	    SS_states.erase(SS_iter->id);
	    SS_iter=SS_states_vector.erase(SS_iter);
	    //cout<<"\tstored best_h is above prunning threshold,f:"<<best_h+SS_iter->g<< ",threshold:"<<sampling_threshold<<endl;
	    continue;
	}
	sampled_states++;
	total_SS_gen_nodes+=SS_iter->weight;
	//cout<<"sampled_state:"<<sampled_states<<",new_f="<<current_heuristic->get_heuristic()+SS_iter->g<<",old f:"<<best_heuristic->get_heuristic()+SS_iter->g<<",g:"<<SS_iter->g<<",weight:"<<SS_iter->weight<<",sampling_threshold:"<<sampling_threshold<<endl;
	int candidate_h=candidate_PC->get_value(unique_samples.at(SS_iter->id).first);
	//cout<<"candidate_h:"<<candidate_h<<",best_h:"<<best_h<<endl;
	/*if(candidate_h<candidate_explicit.get_value(unique_samples.at(SS_iter->id))){
	  cout<<"candidate_h:"<<candidate_h<<",candidate_explicit:"<<candidate_explicit.get_value(unique_samples.at(SS_iter->id))<<endl;
	  exit(0);
	  }*/
	if(candidate_h==numeric_limits<int>::max()){
	    increased_states++;
	    pruned_states+=SS_iter->weight;
	    //cout<<"\t\tsampled_state:,"<<sampled_states<<",out of "<<SS_states.size()<<"is now pruned by dead_end, weight:"<<SS_iter->weight<<",current_total:"<<total_SS_gen_nodes<<",best_h:"<<best_h<<endl;
	    SS_iter++;
	    //if(float(increased_states)/float(sampled_states)>min_improvement_ratio){
	    //  break;
	    //}
	    continue;
	}

	//fitness+=candidate_h;

	/*  if(candidate_h+SS_iter->g>sampling_threshold){
	    pruned_states+=SS_iter->weight;
	    increased_states++;
	    //cout<<"id:,"<<SS_iter->id<<",sampled_state:,"<<sampled_states<<",out of "<<SS_states.size()<<"is now pruned by higher F, weight:"<<SS_iter->weight<<",current_total:"<<total_SS_gen_nodes<<endl;
	    //cout<<"h1="<<current_heuristic->get_heuristic()<<"+g="<<SS_iter->g<<",f:"<<current_heuristic->get_heuristic()+SS_iter->g<<",sampling_threshold:"<<sampling_threshold<<endl;
	    //cout<<"h2="<<best_heuristic->get_heuristic()<<"+g="<<SS_iter->g<<",f:"<<best_heuristic->get_heuristic()+SS_iter->g<<",sampling_threshold:"<<sampling_threshold<<endl;
	}
	else {*/
	    if(candidate_h>best_h){
		pruned_states+=SS_iter->weight;
		increased_states++;
		//if(float(increased_states)/float(sampled_states)>min_improvement_ratio){
		 // break;
		//}
		//cout<<"id:,"<<SS_iter->id<<",candidate_h:"<<candidate_h<<",best_h:"<<best_h<<endl;
		//cout<<"sampled_state:,"<<sampled_states<<",out of "<<SS_states.size()<<"is now pruned by higher h, weight:"<<SS_iter->weight<<",current_total:"<<total_SS_gen_nodes<<endl;
		//cout<<"\tsampled_state:,"<<sampled_states<<",out of "<<SS_states.size()<<"is now pruned by candidate_h:"<<candidate_h<<",best_h:"<<best_h<<",current_total:"<<total_SS_gen_nodes<<endl;
	    }
	    //else{
	      //cout<<"\tstored best_h is above candidate_h:"<<candidate_h<<",best_h:"<<best_h<<endl;
	    //}
	//}
	SS_iter++;
    }
    
    //overall_sampled_states+=sampled_states;
    DEBUG_MSG(cout<<"finished sampling,sampled_states:"<<sampled_states<<",increased_states:"<<increased_states<<",pruned_states:"<<pruned_states<<endl;fflush(stdout););
    //double sampler_time=utils::g_timer()-start_sampler_time;
    
    double saved_time=0;
    if(increased_states>0){

      if(time_or_size_selection==TIME_SELECTION){
	utils::Timer heur_timer;
	//New, we calculate the current time cost for the current result
	int counter_temp=0;
	for (auto sample=unique_samples.begin();sample!=unique_samples.end();sample++){
	  counter_temp++;
	  candidate_PC->get_value(sample->second.first);
	  if(counter_temp>20000){
	  break;
	  }
	}
	double candidate_heur_time_cost=heur_timer.stop()/counter_temp;
	saved_time=total_SS_gen_nodes*(node_gen_and_exp_cost+current_heur_time_cost)-(total_SS_gen_nodes-pruned_states)*(node_gen_and_exp_cost+candidate_heur_time_cost+current_heur_time_cost);
      }
      else{
	//We want to use the same ratio of increased states as if we were doing random_walk, this is 100(threshold) over 20K
	//but since we are doing a SS sample, it applies to the SS prediction of generated and pruned states.
       if((pruned_states/total_SS_gen_nodes)>double(get_threshold())/double(get_num_samples()))
	return increased_states;//True if >0, so choosing any heuristic if it lowers search space in the slightest
       else
	 return false;
      }
     //cout<<"ratio:"<<float(increased_states)/float(sampled_states)<<",node_gen_and_exp_cost:"<<node_gen_and_exp_cost<<",sampling_time:"<<sampler_time<<",candidate_heur_time_cost:"<<candidate_heur_time_cost<<",current_heur_time_cost:"<<current_heur_time_cost<<",best_prev_time:"<<total_SS_gen_nodes*(node_gen_and_exp_cost+best_pdb_collections.size()*heur_time_cost)<<",new_time:"<<(total_SS_gen_nodes-pruned_states)*(node_gen_and_exp_cost+heur_time_cost*(best_pdb_collections.size()+1))<<",saved_time:"<<saved_time<<"total_SS_gen_nodes:"<<total_SS_gen_nodes<<",new_nodes:"<<total_SS_gen_nodes-pruned_states<<",initial_h:"<<candidate.get_value(initial_state)<<endl;
    }
    if(saved_time>0){
      return true;
    }
    else{
      return false;
    }
}

void PatterCollectionEvaluatorSS::sample_states(std::shared_ptr<PatternCollectionInformation> current_result){
  //cout<<"ss,sample_states being called"<<flush<<endl;
  best_heur_dead_ends=0;
  SS_states.clear();
  SS_states_vector.clear();
  const State &initial_state = task_proxy->get_initial_state();
  unsigned candidate_threshold=max(unsigned(current_result->get_value(initial_state)),max(unsigned(1),get_threshold()));
  if(get_threshold()<candidate_threshold){
    set_threshold(candidate_threshold);
    cout<<"ss, sample_states raised SS threhsold to:,"<<get_threshold()<<",current initial h value:"<<current_result->get_value(initial_state)<<endl;
  }
  else
    cout<<"ss, sample_states kept SS threhsold to:,"<<get_threshold()<<",current initial h value:"<<current_result->get_value(initial_state)<<endl;

  double start_probe_time=utils::g_timer();
  int repetition=0;
  vector<double> probe_avgs;
  //threshold=44;
  //Now sampling
  //cout<<"ss, sample_states 1), setting threshold"<<flush<<endl;
  //1) Set initial Threshold
  if(get_threshold()==unsigned(current_result->get_value(initial_state))){
    //set_threshold(200);
    //cout<<"ss,initial threshold:"<<get_threshold()<<endl;
    for (int prob_index=0;prob_index<2000;prob_index++){
      SS_states.clear();SS_states_vector.clear();
      for(size_t j=1;j<20;j++){
	probe_best_only(current_result);
	//cout<<"ss,SS_states.size:"<<SS_states.size()<<endl;
	if(SS_states.size()>500){
	  //cout<<"ss,breaking,SS_states.size:"<<SS_states.size()<<endl;
	  break;
	}
      }
      if(SS_states.size()>500){
	cout<<"ss,best initial threshold raised to:"<<get_threshold()<<",probe states:"<<SS_states.size()<<endl;
	SS_states.clear();SS_states_vector.clear();
	break;
      }
      set_threshold(max(get_threshold()+1,unsigned(get_threshold()*1.2)));
      cout<<"ss,next_threshold:"<<get_threshold()<<",SS_states.size:"<<SS_states.size()<<endl;
    }
  }

  cout<<"ss, sample_states 2), Run probes"<<flush<<endl;
    //2) Run Probes till time limit or size limit
    for (repetition=0;repetition<10;repetition++){
	vector<double> probe_data;
	for (int prob_index=0;prob_index<1000;prob_index++){
	  DEBUG_MSG(cout<<"prob_index:"<<prob_index<<",starting probe_best_only"<<flush<<endl;;);
	    probe_data.push_back(probe_best_only(current_result));
	    if(utils::g_timer()-start_probe_time>10.0){
		cout<<"exceeded 10 seconds limit for probes, number of repetitions completed:"<<repetition<<endl;
		break;
	    } else if(SS_states.size()>20000){
		cout<<"exceeded 20K max SS_states sampled,current size:"<<SS_states.size()<<",probe:"<<prob_index<<",threshold:"<<get_threshold()<<endl;
		break;
	    }
	}
	cout<<"\tss,repetition:"<<repetition<<",threshold:"<<get_threshold()<<",SS_states:"<<SS_states.size()<<endl;
	pair<double,double> avg_and_dev=utils::avg_and_standard_deviation(probe_data);
	//cout<<scientific<<"repetition:"<<repetition<<",probe data average:"<<avg_and_dev.first<<endl;
	//cout<<scientific<<"repetition:"<<repetition<<",probe data standard deviation:"<<avg_and_dev.second<<endl;
	probe_avgs.push_back(avg_and_dev.first);
  
	if(utils::g_timer()-start_probe_time>10.0){
	    cout<<"ss,exceeded 10 seconds limit for probes, number of repetions completed:"<<repetition<<",threshold:"<<get_threshold()<<endl;
	    break;
	} else if(SS_states.size()>20000){
	    cout<<"ss,exceeded 20K max SS_states sampled,current size:"<<SS_states.size()<<",threshold:"<<get_threshold()<<endl;
	    break;
	} else if(avg_and_dev.first>pow(10,100)){
	    cout<<"ss,avg probe:"<<avg_and_dev.first<<" past 10^100 nodes, not going further, gets very unprecise"<<",threshold:"<<get_threshold()<<endl;
	    break;
	}

	SS_states_copy=SS_states;
	//update threshold, minimum raise of one unit
	set_threshold(max(get_threshold()+1,unsigned(get_threshold()*1.2)));
	cout<<"ss, sample_states2,next_threshold:"<<get_threshold()<<endl;
	SS_states.clear();
	SS_states_vector.clear();
    }
  
    overall_probe_time+=utils::g_timer()-start_probe_time;
  
    cout<<"ss, populate SS_iter_map 3)"<<flush<<endl;
    //3)Now populate SS_iter_map to speed up evaluations
			
    //cout<<"time:,"<<utils::g_timer()<<",starting sorting SS states,size:"<<SS_states.size()<<endl;
    map<size_t,pair<int,double> >::iterator SS_iter_map;
			
    for(SS_iter_map=SS_states.begin();SS_iter_map!=SS_states.end();){
      //cout<<",SS_iter_map->first:"<<SS_iter_map->first<<endl;
      if(unique_samples.find(SS_iter_map->first)==unique_samples.end()) {
	cout<<"ss,state:"<<SS_iter_map->first<<" not in unique_samples!!!"<<endl;exit(0);
      }
			    
      State current_state(unique_samples.at(SS_iter_map->first).first);
      int current_h=current_result->get_value(current_state);
      if(current_h==numeric_limits<int>::max()){
	best_heur_dead_ends++;
	SS_iter_map=SS_states.erase(SS_iter_map++);
	//cout<<"eliminating best_heur dead_end, state_id:"<<SS_iter_map->first<<endl;
	continue;
      }

      SS_state temp;
      temp.id=SS_iter_map->first;
      //cout<<"temp.id:"<<temp.id<<",SS_iter_map->first:"<<SS_iter_map->first<<endl;
      temp.g=SS_iter_map->second.first;
      temp.weight=SS_iter_map->second.second;
      SS_states_vector.push_back(temp);
      SS_iter_map++;
    }
    g_rng()->shuffle(SS_states_vector);
    //4) Measure the time costs of current heuristic combination, needed for evaluate when doing time predictions
    if(time_or_size_selection==TIME_SELECTION){
      utils::Timer heur_timer;
      int counter_temp=0;
      for (auto sample=unique_samples.begin();sample!=unique_samples.end();sample++){
	counter_temp++;
	current_result->get_value(sample->second.first);
	if(counter_temp>20000){
	  break;
	}
      }
      current_heur_time_cost=heur_timer.stop()/counter_temp;
    }
    DEBUG_MSG(cout<<"time:,"<<utils::g_timer()<<",finished randomzing SS states vector,size:"<<SS_states.size()<<",best_heur_dead_ends"<<best_heur_dead_ends<<endl;);
    cout<<"ss,time:,"<<utils::g_timer()<<",finished randomzing SS states vector,size:"<<SS_states.size()<<",best_heur_dead_ends"<<best_heur_dead_ends<<endl;
  }

double PatterCollectionEvaluatorSS::probe_best_only(std::shared_ptr<PatternCollectionInformation> current_result){
  DEBUG_MSG(cout<<"hello probe_best_only"<<flush<<endl;);
  //Creating fake heuristic (blind) to store h-value results needed by SS sampler class interface
  options::Options temp_options2;
  temp_options2.set<int>(
      "cost_type", NORMAL);
  temp_options2.set<bool>(
      "cache_estimates", false);
  blind_search_heuristic::BlindSearchHeuristic temp_blind_heuristic(temp_options2);
  DEBUG_MSG(cout<<"blind_search_heuristic created"<<flush<<endl;);
  sampler = new TypeSystem(&temp_blind_heuristic);
  DEBUG_MSG(cout<<"sampler created,"<<flush<<endl;);
  set<vector<int> > F_culprits;
  map<Type, SSNode> queue;
  const State &initial_state = task_proxy->get_initial_state();
  //cout<<"after initial_state"<<flush<<endl;
  set_threshold(max(unsigned(current_result->get_value(initial_state)),max(unsigned(1),get_threshold())));
  DEBUG_MSG(cout<<"ss, calling probe_best_only,threshold:"<<get_threshold()<<endl;fflush(stdout););
  
  if(current_result->pdb_counts()>0){
      int initial_h=current_result->get_value(initial_state);
      DEBUG_MSG(cout<<initial_h<<flush<<endl;);
      //initial_h needs to be at least 1
      initial_h=max(1,initial_h);
  }
  else{
      cerr<<"cant call probe_best_only if there are no pdbs already added! current_result collections are empty,DEBUG ME!!!"<<endl;
      exit(1);
  }	  
  sampling_threshold=get_threshold();
  
  vector<OperatorProxy> applicable_ops; 
  successor_generator->generate_applicable_ops(initial_state,applicable_ops);
  //cout<<"before amount_initial"<<endl;fflush(stdout);
  double amount_initial = (double)applicable_ops.size();
  
  //probe_best 2) Need to initialize child state before using	
  State child=initial_state; 
  prev_current_collector=1 + amount_initial;
  max_collector=1 + amount_initial;
  size_t initial_state_id = initial_state.hash();
  unique_samples.insert(make_pair(initial_state_id,make_pair(child,0)));
	
  //3) Cycle check, this is an advantage over random_walk, we are going further this way
  //It could also be implemented in random_walk, TO_DO
  map<size_t,int> cycle_check;
  cycle_check[initial_state.hash()]=0;
  //cout<<"initial_state.hash:"<<initial_state_id<<endl;
	
  //4) Initialize SSNode root
  SSNode node;
  node.setId(initial_state_id);
  node.setWeight(1.0);
  node.setGreal(0);  //g_real value of the node
  int initial_value=current_result->get_value(initial_state);
  node.setH(initial_value);
  
  SS_states[initial_state_id].first=0;
  SS_states[initial_state_id].second=amount_initial;
  // Seeding the prediction with the children of the start state
  Type type = sampler->getType(initial_value);
  type.setLevel( 0 ); // level where the node is located
  queue.insert( pair<Type, SSNode>( type, node ) );
  int nraiz = 0;
  long queue_counter=0;
  
  while( !queue.empty() )
  {
    queue_counter++;
    if(queue_counter%1000==0){
      if(SS_states.size()>22000){
	cout<<"SS_states past limit,size:"<<SS_states.size()<<", no more SS generation, getting out"<<endl;
	return 0;
      }
      if(utils::g_timer()>time_limit){
	cout<<"Search_timer past maximum sampling_time"<<endl;fflush(stdout);
	//cout<<"selecting best heuristic after search_time: "<<search_time()<<", seconds,g_timer:"<<g_timer()<<endl;
	return(-1);
      }
    }
	    
    //cout<<"queue.size:"<<queue.size()<<endl;//",search_time:"<<search_time<<endl;
	    
    Type out = queue.begin()->first;
    SSNode s = queue.begin()->second;
	    
    int g_real =  s.getGreal();
    int level = out.getLevel();
    double w = s.getWeight(); 
	    
    std::map<Type, SSNode>::iterator rt;
    rt = queue.find(out);
	    
    queue.erase( rt );
    
    nraiz++;                
	    
    /*cout<<"\tRaiz: "<<nraiz<<" h  = "<<max_h;fflush(stdout);
    cout<<", g_real = "<<g_real<<", f = ";fflush(stdout);
    cout<<max_h + g_real;fflush(stdout);
    cout<<", level = "<<level;fflush(stdout);*/
    
		
  vector<OperatorProxy> applicable_ops;
  map<size_t,pair<State,int> >::iterator it=unique_samples.find(s.get_id());
  if (it == unique_samples.end()){
      cerr<<"any retrieved state should be on the queue, FIX ME!!!!"<<endl;exit(1);
  }

    State* current_state=&(it->second.first);
    successor_generator->generate_applicable_ops(*current_state,applicable_ops); //count nodes generated

    //cout<<"\t_____________________begin Childs________________________\n";fflush(stdout);
    int h =  INT_MAX/2;
    L.clear();
    check.clear();
    
    for (size_t i = 0; i < applicable_ops.size(); ++i)
    {
	const OperatorProxy &op = applicable_ops[i];
	child = current_state->get_successor(op);

	size_t child_hash=child.hash();
	//check if we already have state
	std::map<size_t, int>::iterator cycle_check_iterator=cycle_check.find(child_hash);
		
	//Cycle check
	if(cycle_check.find(child_hash)!=cycle_check.end()){
	  //check if dup state is at lower depth than recorded, otherwise, is cycle so continue
	  if(cycle_check_iterator->second>g_real+op.get_cost()){
	    cycle_check_iterator->second=g_real+op.get_cost();
	  }
	  else{
	    continue;
	  }
		
	}
	else{//add state
	  cycle_check[child_hash]=g_real+op.get_cost();
	}
	//SS_states is global to all runs, used here to avoid opening paths where g is known to be suboptimal
	//cycle_check is just to avoid cycles, specially zero-op cost
	if(SS_states.find(child_hash)!=SS_states.end()){
	    if(g_real+op.get_cost()>SS_states[child_hash].first){
		continue;
	    }
	}
			    
	h=current_result->get_value(child);
	if(h==numeric_limits<int>::max()){
	  pair<map<size_t,pair<State,int> >::iterator,bool> ret;
	  ret=unique_samples.insert(make_pair(child_hash,make_pair(child,h)));
	  if(!ret.second){//keep max h value stored, so we do not need to recalculate
	    ret.first->second.second=max(h,ret.first->second.second);
	  }
	  continue;
	  //cout<<"prev_h_deade_end_found"<<endl;
	}
//#ifdef _SS_DEBUG
	//cout<<"\tg_real = "<<g_real + op.get_cost()<<" f_min = "<< h + g_real + op.get_cost()<<endl;//",b_child_v.count:"<<b_child_v.count()<<endl;
	//cout<<"\tget_adjusted_cost(*op) = "<<get_adjusted_action_cost(*op,cost_type)<<"\n";
	//cout<<"\tChild_"<<(i+1)<<" : h = "<<h<<",b_child_v:"<<b_child_v<<endl; 
		
	//cout<<h + g_real + get_adjusted_action_cost(*op,cost_type)<<endl;
	//cout<<", level = "<<(level + 1);
	//cout<<", w = "<<w<<"\n";
//#endif

	vector<OperatorProxy> applicable_ops2; 
	//cout<<"S:"<<endl;global_state_2.dump_inline();fflush(stdout);
	successor_generator->generate_applicable_ops(child,applicable_ops2); //count nodes generated
             
	int amount = applicable_ops2.size();
                
	//add to the collector
	pair<map<size_t,pair<State,int> >::iterator,bool> ret;
	ret=unique_samples.insert(make_pair(child_hash,make_pair(child,h)));
	//ret=unique_samples.insert(make_pair(child_hash,make_pair(child,5)));
	if(!ret.second){//keep max h value stored, so we do not need to recalculate
	  ret.first->second.second=max(h,ret.first->second.second);
	}
	//visited_states.insert(child.hash());
	if ( unsigned(h + g_real + op.get_cost())  <= get_threshold()) {
	    //Keep a record of all sampled states and their maximum weights and minimum depths
	    if(SS_states.find(child_hash)!=SS_states.end()){
	      if(g_real+op.get_cost()<SS_states[child_hash].first){
		  max_collector+=amount*w-SS_states[child_hash].second;
		  //cout<<"Reviewed Stored id:,"<<child_hash<<",new_g:"<<g_real + op.get_cost()<<",new F:"<<h + g_real + op.get_cost()<<endl;
		  SS_states[child_hash].first=g_real + op.get_cost();
		  //If state is found in different probes with same depth, we want to keep record of maximum impact
		  SS_states[child_hash].second=max(SS_states[child_hash].second,amount*w);
	      }
	      else if(g_real+op.get_cost()<SS_states[child_hash].first){
		  max_collector+=SS_states[child_hash].second;
	      }
		    
	    }
	    else{
		max_collector+=amount*w;
		SS_states[child_hash].first=g_real + op.get_cost();
		SS_states[child_hash].second=amount*w;
		//cout<<"New Stored id:,"<<child_hash<<",new_g:"<<g_real + op.get_cost()<<",new F:"<<h + g_real + op.get_cost()<<endl;
	    }

	    //Make pruning
	    Type object = sampler->getType(h);
			   
		    object.setLevel( level + 1 );
                           
		    SSNode child_node;
                           
		    child_node.setId(child_hash);
		    child_node.setWeight(w);
		    child_node.setGreal(g_real + op.get_cost()); 
		    //child_node.setHC(h_child_v);
		    child_node.setH(h);

				
#ifdef _SS_DEBUG
		    cout<<"\t\tChild f<=threshold: h = "<<h; 
		    cout<<", g_real = "<<g_real + get_adjusted_action_cost(*op,cost_type)<<" f = ";
		    cout<<h + g_real  +  get_adjusted_action_cost(*op,cost_type);
		    cout<<", level = "<<level + 1<<"\n";
#endif
		    //ZERO COST OPERATORS NEED EXTRA BFS, until all descendant nodes found for this F-level  
		    //if (op.get_cost() == 0) {//ZERO COST OPERATORS
		    //  cout<<"TO-Do:Update BFS to current version so we can add SS_States to zero-cost operations!"<<endl;exit(0);
//#ifdef _SS_DEBUG
//				cout<<"\t\tget_adjusted_cost(*op) == 0\n";fflush(NULL);
//#endif
//				//TO-DO:UPDATE BFS function to current FD version
//			   	BFS(child_node, object,best_heuristic);
//#ifdef _SS_DEBUG
//				cout<<"after BFS, L size:"<<L.size()<<endl;
//#endif
//				std::set<SSQueue>::iterator it;
//				for (it = L.begin(); it != L.end(); it++) {	
//					SSQueue sst = *it;
//					SSNode node = sst.getNode();
//					Type t = sst.getT();
//					//double new_g_real = node.getGreal();
//					int new_state_id = node.get_id();
//					pair<set<StateID>::iterator, bool> ret2;
//					//std::pair<std::map<boost::dynamic_bitset<>, double>::iterator, bool> ret2;
//					State global_state3 = g_state_registry->lookup_state(new_state_id);
//					unique_samples.insert(global_state3);
//					ret2=visited_states.insert(new_state_id);
//
//					if(SS_states.find(new_state_id)!=SS_states.end()){
//					  if(g_real<SS_states[new_state_id].first){
//					    //cout<<"Reviewed Stored id:,"<<child.get_id()<<",new_g:"<<g_real + get_adjusted_action_cost(*op,cost_type)<<",new F:"<<h + g_real + get_adjusted_action_cost(*op,cost_type)<<endl;
//					    SS_states[new_state_id].first=g_real+get_adjusted_action_cost(*op,cost_type);
//					    SS_states[new_state_id].second=amount*w;
//					  }
//					} 
//					else{
//					  SS_states[new_state_id].first=g_real+get_adjusted_action_cost(*op,cost_type);
//					  SS_states[new_state_id].second=amount*w;
//					}
//					
//					//Check for state in cyclical path
//					if (!ret2.second) {
//					  continue;//the state was already visited
//					}
//					double w2 = node.getWeight();
//
//					//cout<<"\n\t\tNode restored: h = "<<t.getH()<<", g_real = "<<node.getGreal()<<", f = "<<t.getH() + node.getGreal()<<", level = "<<t.getLevel()<<", w = "<<w2<<"\n";
//                           			
//					map<Type, SSNode>::iterator queueIt = queue.find( t );
//#ifdef _SS_DEBUG
//					cout<<"state_id:"<<new_state_id<<endl;
//#endif
//			   		if( queueIt != queue.end() )
//			   		{
//                                	SSNode snode = queueIt->second;
//
//#ifdef _SS_DEBUG
//                                		cout<<"\t\t\tzc: The duplicate node is: h = "<<queueIt->first.getH()<<", g = "<<snode.getGreal()<<", f = "<< queueIt->first.getH() + snode.getGreal()<<", w = "<<snode.getWeight()<<", level = "<<queueIt->first.getLevel()<<"\n";
//#endif
//                                
//						double wa = (double)snode.getWeight();
//						//snode.setWeight( wa + w);
//                                		queueIt->second.setWeight(wa + w2); // set w the node that already exists
//                                		//cout<<"\t\t\tzc: before ss process starts, the w of the duplicate node is updated to: "<<queueIt->second.getWeight()<<endl; 
//                                		//std::pair<std::map<Type, SSNode>::iterator, bool> ret0;
//                                		//ret0 = queue.insert(pair<Type, SSNode>(object, snode));
//                                		//cout<<"\tsnode.getWeight() = "<<snode.getWeight()<<endl;
//                                		//queueIt->second.setWeight(snode.getWeight());
//						double prob = ( double )w2 / (double)( wa + w2);
//						//int rand_100 =  RanGen2->IRandom(0, 99);  //(int)g_rng.next(100);
//						int rand_100 =  rand()%100;  //(int)g_rng.next(100);
//                          	 
//                                		double a = (( double )rand_100) / 100;
//                                		//cout<<"a = "<<a<<" prob = "<<prob<<endl;
//                                
//						if (a < prob) 
//						{
//#ifdef _SS_DEBUG
//                                        		cout<<"\t\t\tzc: Added even though is duplicate.\n";                               
//#endif
//				        		node.setWeight( wa + w2);
//#ifdef _SS_DEBUG
//                                        		cout<<"\t\t\tzc: the w is updated to = "<<node.getWeight()<<endl;
//#endif
//                                        		std::pair<std::map<Type, SSNode>::iterator, bool> ret3;
//                                     			queue.erase(t); 
//                                        
//                                        		ret3 = queue.insert( pair<Type, SSNode>( t, node ));      
//                                        		queueIt = ret3.first;
//                                        		queueIt->second.setWeight(node.getWeight());
//						} else {
//#ifdef _SS_DEBUG
//                                        		cout<<"\t\t\tzc: Not added.\n";
//                                        		cout<<"\t\t\tbut the w is updated for the node that already exists to: "<<queueIt->second.getWeight()<<endl;
//#endif
//                                		}
//			   		} 
//			   		else
//			   		{
//#ifdef _SS_DEBUG
//                                		cout<<"\t\t\tzc: New L node added.\n";
//#endif
//						queue.insert( pair<Type, SSNode>( t, node ) );
//                                		//cout<<"\t\tsucc_node2.getWeight() = "<<succ_node2.getWeight()<<"\n";
//                                
//                                		//cout<<"\t\t\tzc: Child: h = "<< t.getH() <<", g_real = "<< new_g_real <<", f = "<< t.getH() + new_g_real << " threshold: " << threshold <<" w = "<<node.getWeight()<<endl;
//                           		}// End queueIt != queue.end()
//				}   //End for set lopp
//#ifdef _SS_DEBUG
//				cout<<"Finished with L"<<endl;
//#endif
		    //} else {

						       
		    map<Type, SSNode>::iterator queueIt = queue.find( object );
		    if( queueIt != queue.end() )
		    {
			SSNode snode = queueIt->second;

#ifdef _SS_DEBUG
			cout<<"\t\tThe duplicate node is: h = "<<h;
				      
			cout<<", g_real = "<<g_real + get_adjusted_action_cost(*op,cost_type)<<" f = ";
			cout<<h + g_real  +  get_adjusted_action_cost(*op,cost_type);
			cout<<", w = "<<snode.getWeight();
			cout<<", level = "<<level + 1<<"\n";
#endif

			double wa = (double)snode.getWeight();
			//snode.setWeight( wa + w);
			queueIt->second.setWeight(wa + w);
#ifdef _SS_DEBUG
			cout<<"\t\tbefore ss process starts, the w of the duplicate node is updated to: "<<queueIt->second.getWeight()<<endl;
#endif
			//std::pair<std::map<Type, SSNode>::iterator, bool> ret0;

			//ret0 = queue.insert(pair<Type, SSNode>(object, snode));
			//cout<<"\tsnode.getWeight() = "<<snode.getWeight()<<endl;
			//queueIt->second.setWeight(snode.getWeight());
   
   
			double prob = ( double )w / (double)( wa + w );
			//int rand_100 =  RanGen2->IRandom(0, 99);  //(int)g_rng.next(100);
			int rand_100 =  rand()%100;  //(int)g_rng.next(100);
				   
			double a = (( double )rand_100) / 100;
			//cout<<"a = "<<a<<" prob = "<<prob<<endl; 
				  
			if( (a < prob))
			{
			    //unique_samples.insert(make_pair(child.hash(),child));
			    //cout<<"a<prob,hash:"<<child.hash()<<endl;
#ifdef _SS_DEBUG
			    cout<<"\t\tAdded even though is duplicate.\n";
#endif
					  
			    child_node.setWeight( wa + w);
#ifdef _SS_DEBUG
			    cout<<"\t\tthe w is updated to = "<<child_node.getWeight()<<endl;
#endif
			    std::pair<std::map<Type, SSNode>::iterator, bool> ret;
			    queue.erase(object); 

			    ret = queue.insert( pair<Type, SSNode>( object, child_node ));      

			    queueIt = ret.first;
			    queueIt->second.setWeight(child_node.getWeight());
					  
					  
			} else {
#ifdef _SS_DEBUG
			    cout<<"\t\tNot added.\n";
			    cout<<"\t\tbut the w is updated for the node that already exists to: "<<queueIt->second.getWeight()<<endl;
#endif
			}
		    } 
		    else
		    {
			//unique_samples.insert(make_pair(child.hash(),child));
			//cout<<"new_SS_state,hash:"<<child.hash()<<endl;
#ifdef _SS_DEBUG
			cout<<"\t\tNew node added\n";
			//Now update the non-prunning set of heuristics for the node
#endif
			queue.insert( pair<Type, SSNode>( object, child_node ) );
		    }
		    //}
		}
		else 
		{
#ifdef _SS_DEBUG
		    cout << "\tNode was pruned!" << endl;
#endif
		}
#ifdef _SS_DEBUG
		cout<<"\tend Child_"<<(i+1)<<"\n";
#endif
	    }
	}
  delete sampler;
  return max_collector;
}


  void PatterCollectionEvaluatorSS::clear_dominated_heuristics(std::shared_ptr<PatternCollectionInformation> current_result,std::shared_ptr<PatternCollectionInformation> &new_result,
      shared_ptr<ModularZeroOnePDBs> candidate_ptr){
	double start_time=utils::g_timer();
	vector<int> current_best_h_values;
	current_best_h_values.reserve(unique_samples.size());

	//cout<<"calling clear_dominated_heuristics with "<<current_result->get_max_additive_subsets()->size()+1<<" best heuristics and unique_samples:"<<unique_samples.size()<<endl;
	
	//First we get all the values for sampled states with the latest PC
	int i=0;
	for(map<size_t,pair<State,int> >::iterator it=unique_samples.begin(); it!=unique_samples.end();it++){
	  current_best_h_values[i++]=candidate_ptr->get_value(it->second.first);
	}
	cout<<"\ttime to populate current_best_h_values:"<<utils::g_timer()-start_time<<",unique_samples:"<<i<<endl;

	//Now go through each additive subset and check if they are dominated for the sampled set of states
	shared_ptr<MaxAdditivePDBSubsets> current_max_additive_subsets=current_result->get_max_additive_subsets();


	int h=0;
	i=0;
	for(auto additive_subset : *current_max_additive_subsets){
	  bool dominated_heur=true;
	  int j=0;
	 for(map<size_t,pair<State,int> >::iterator it=unique_samples.begin(); it!=unique_samples.end();it++){
	    //If dead_end skip
	    if(current_best_h_values[j]==INT_MAX){
	      j++;
	      continue;
	    }
	    h=calculate_max_additive_subset(additive_subset,it->second.first);
	    //NO BREAKS BECAUSE WE WANT TO CALCULATE ALL THE NEW HIGHER H VALUES
	    //IF HEUR IS NOT DOMINATED
	    if (h == numeric_limits<int>::max()){
	      dominated_heur=false;
	      //cout<<"\tcolleciton ["<<i<<" is undominated because of dead_end, prev_val:"<<current_best_h_values[j]<<",h:"<<h<<endl;
	      current_best_h_values[j]=INT_MAX;
	    }
	    else if(h>current_best_h_values[j]){
	      dominated_heur=false;
	      //cout<<"\tcolleciton ["<<i<<" is undominated because of higher_h, prev_val:"<<current_best_h_values[j]<<",h:"<<h<<endl;
	      current_best_h_values[j]=h;
	    }
	    j++;
	  }
    
	  if(!dominated_heur){
	    //cout<<"adding heur["<<i<<"] to list of heurs"<<endl;

	    shared_ptr<vector<shared_ptr<pdbs3::PatternDatabaseInterface> > > additive_subset_ptr;
	    additive_subset_ptr=make_shared<vector<shared_ptr<pdbs3::PatternDatabaseInterface> > >(additive_subset);
	    new_result->include_additive_pdbs(additive_subset_ptr);
	    /*if(cleaned_best_pdb_collections.size()>15){
	      cout<<"max of 15 pdb_collections, otherwise timewise takes too long"<<endl;
	     break;
	    } */
	  }
	  else{
	    //cout<<"collection["<<i<<"] is dominated,eliminating "<<endl;
	  }
	  i++;
	}

	
	//cout<<"clear_dominated_heuristics finished:"<<new_result->get_max_additive_subsets()->size()<<","<<",time:"<<utils::g_timer()-start_time<<endl;
  }
  int PatterCollectionEvaluatorSS::calculate_max_additive_subset(PDBCollection max_subset,State current_state){
    int h=0;
    int h_temp=0;
      for (auto pdb : max_subset){
	h_temp=pdb->get_value(current_state);
	if (h_temp == numeric_limits<int>::max()){
	  h=numeric_limits<int>::max();
	  break;
	}
	else{
	  h+=h_temp;
	}
      }
      return h;
  }




static shared_ptr<PatternCollectionEvaluator>_parse(options::OptionParser &parser) {
  parser.add_option<int> ("time_limit", "If populated,stop construction on first node past boundary and time limit", "100");
  parser.add_option<int> ("time_or_size_selection", "If TIME_SELECTION do time predictions, if SIZE_SELECTION do size prediction", "1");
  options::Options options = parser.parse();
  parser.document_synopsis(
      "Pattern Generator RBP",
      "RBP-stype selection of variables to generate Pattern Collection");
  options::Options opts = parser.parse();
  if (parser.dry_run())
      return 0;

  return make_shared<PatterCollectionEvaluatorSS>(opts);
}

static options::PluginShared<PatternCollectionEvaluator> _plugin("ss_walk", _parse);
}
