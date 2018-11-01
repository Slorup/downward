#include "modular_heuristic.h"
#include "types.h"

//#include "pattern_generator.h"
#include "pattern_collection_generator_RBP.h"
//#include "pattern_collection_generator_bin_packing_v1.h"
#include "pattern_collection_evaluator_RandWalk.h"
#include "pdb_factory_symbolic.h"


#include "../option_parser.h"
#include "../plugin.h"
#include "../task_proxy.h"
#include "../globals.h"

#include <limits>
#include <memory>
#include <climits>
#include "zero_one_pdbs.h"
#include "../utils/countdown_timer.h"
#include "learning.h"
#include "pattern_database_symbolic.h"
#include "pattern_collection_generator_GamerStyle.h"
#include <boost/range/adaptor/reversed.hpp>
//#include "pdb_factory_symbolic.h"
//#include "pdb_factory.h"

using namespace std;
std::ostream & operator<<(std::ostream &os, set<int> pattern){
    for (int v : pattern) os << "," << v;  
    return os;
}
std::ostream & operator<<(std::ostream &os, vector<int> pattern){
    for (int v : pattern) os << "," << v;  
    return os;
}

bool recompute_additive_sets=false;

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
    pattern_local_search(opts.get<shared_ptr<PatternCollectionLocalSearch>>("local_search")),
	  modular_time_limit(opts.get<int>("modular_time_limit")),
	  always_CBP_or_RBP_or_UCB(opts.get<int>("always_cbp_or_rbp_or_ucb")), 
	  terminate_creation(opts.get<bool>("terminate_creation")),
	  create_perimeter(opts.get<bool>("create_perimeter")), 
	  only_gamer(opts.get<bool>("only_gamer")), 
	  gamer_classic(opts.get<bool>("gamer_classic")), 
	  only_CBP(opts.get<bool>("only_cbp")), 
    pdb_factory (opts.get<shared_ptr<PDBFactory>>("pdb_factory")) {
      cout<<"Hi nonagnostic_v2"<<endl;
      bool unterminated_pdbs=false;
      cout<<"modular_time_limit:"<<modular_time_limit<<endl;
      cout<<"terminate_creation:"<<terminate_creation<<endl;
      cout<<"pdb_type:"<<pdb_factory->name()<<endl;
      cout<<"only_gamer:"<<only_gamer<<",gamer_classic:"<<gamer_classic<<endl;
      cout<<"only_CBP:"<<only_CBP<<endl;
      unsigned num_goals_to_group=0;
      bool disjunctive_choice=true;
      modular_heuristic_timer = new utils::CountdownTimer(modular_time_limit);
      TaskProxy task_proxy(*task);
      int initial_h=0;
      int new_initial_h=0;
      int num_episodes=0;
      int PC_counter=0;
      double pdb_max_size=0;
      vector<shared_ptr<ModularZeroOnePDBs> > pdb_ptr_collection;
      Options opts2=opts;
      Options opts3=opts;
      Options opts4=opts;
      Options opts5=opts;
      opts2.set<int>("time_limit",40);
      opts3.set<int>("packer_selection",0);
      opts4.set<int>("packer_selection",1);
      opts5.set<int>("packer_selection",2);
      PatternCollectionGeneratorGamer alternative_pattern_generator(opts2);
      


      vector<pair<double,Pattern > > improving_patterns;
      //double best_average_h_value=0;//For Gamer-Style selection
      //1 means Random split into two patterns, 0 means CBP
      Learning UCB_generator; UCB_generator.insert_choice(GAMER_GENERATOR); UCB_generator.insert_choice(RBP_CBP_GENERATOR);
      Learning UCB_local_search; UCB_local_search.insert_choice(GAMER_LOCAL_SEARCH);
      //Initialize reward to 2 so it does not go too fast for initial selection
      UCB_generator.increase_reward(RBP_CBP_GENERATOR,4);UCB_generator.increase_reward(GAMER_GENERATOR,4);
      Learning UCB_RBP_vs_CBP; UCB_RBP_vs_CBP.insert_choice(1); UCB_RBP_vs_CBP.insert_choice(0);
      UCB_RBP_vs_CBP.increase_reward(1,10);UCB_RBP_vs_CBP.increase_reward(0,10);

      Learning UCB_sizes;//PDB size selector
      map<double,Learning> UCB_Disjunctive_patterns;//Disjunctive (or not) pattern selector,one per size;
      //UCB_sizes.insert_choice(pow(10,4));
      //UCB_sizes.insert_choice(8);
      Learning binary_choice;binary_choice.insert_choice(1);binary_choice.insert_choice(0);
      //Initialize reward to 2 so it does not go too fast for initial selection
      binary_choice.increase_reward(1);binary_choice.increase_reward(0);
      size_t overall_num_goals=task_proxy.get_goals().size();
      cout<<"overall_num_goals:"<<overall_num_goals<<endl;
      //bool partial_gamer_run=false;//In case we did not process all available variabes due to taking too long
      //we need to know so that if Gamer come back with a partial list and none was used then we do not block
      //more Gamer turns unless none of those partial runs resulted on adding a variable
      bool partial_gamer_improvement_found = false;
      
      Learning goals_choice;
      for (size_t i=0;i<overall_num_goals;i++){
        goals_choice.insert_choice(i+1);
      //Initialize reward to 2 so it does not go too fast for initial selection
        goals_choice.increase_reward(i+1);
      }
      map<double,Learning> UCB_goals_to_group;//Disjunctive (or not) pattern selector,one per size;
      
      Learning terminate_choice;terminate_choice.insert_choice(1);terminate_choice.insert_choice(0);
      //terminate_choice.increase_reward(0,10);//biasing towards not terminating in the begining
      
      UCB_sizes.insert_choice(9);
      UCB_sizes.insert_choice(10);
      UCB_sizes.insert_choice(11);
      UCB_sizes.insert_choice(12);
      UCB_sizes.insert_choice(13);
      UCB_sizes.insert_choice(14);
      UCB_sizes.insert_choice(15);
      UCB_sizes.insert_choice(16);
      UCB_sizes.insert_choice(17);
      UCB_sizes.insert_choice(18);
      UCB_sizes.insert_choice(19);
      UCB_sizes.insert_choice(20);
      //UCB_sizes.insert_choice(25);
      //UCB_sizes.insert_choice(30);
      //UCB_sizes.insert_choice(35);
      //UCB_sizes.insert_choice(50);
      //UCB_sizes.insert_choice(60);

      /*UCB_Disjunctive_patterns[pow(10,8)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,9)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,10)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,11)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,12)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,13)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,14)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,15)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,16)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,17)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,18)]=binary_choice;
      UCB_Disjunctive_patterns[pow(10,19)]=binary_choice;
      //UCB_Disjunctive_patterns[pow(10,20)]=binary_choice;
      
      UCB_goals_to_group[pow(10,8)]=goals_choice;
      UCB_goals_to_group[pow(10,9)]=goals_choice;
      UCB_goals_to_group[pow(10,10)]=goals_choice;
      UCB_goals_to_group[pow(10,11)]=goals_choice;
      UCB_goals_to_group[pow(10,12)]=goals_choice;
      UCB_goals_to_group[pow(10,13)]=goals_choice;
      UCB_goals_to_group[pow(10,14)]=goals_choice;
      UCB_goals_to_group[pow(10,15)]=goals_choice;
      UCB_goals_to_group[pow(10,16)]=goals_choice;
      UCB_goals_to_group[pow(10,17)]=goals_choice;
      UCB_goals_to_group[pow(10,18)]=goals_choice;
      UCB_goals_to_group[pow(10,19)]=goals_choice;*/



      //double initial_pdb_size=pow(10.0,UCB_sizes.make_choice());
      double initial_pdb_size=pow(10.0,9);//start safely!
      //double initial_disjunctive=UCB_Disjunctive_patterns[initial_pdb_size].make_choice();
      double initial_disjunctive=binary_choice.make_choice();
      //cout<<"after disjunctive_choice"<<endl;
      double initial_goals_to_group=goals_choice.make_choice();
      cout<<"after goals_choice"<<endl;
      pdb_max_size=initial_pdb_size;
      cout<<"initial_pdb_size:"<<initial_pdb_size<<endl;
      //need result here to store final PDB collection
      result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
    
      cout<<"initial pdb type:"<<pdb_factory->name()<<endl;

      const State &initial_state = task_proxy.get_initial_state();
      //ModularZeroOnePDBs candidate(task_proxy, Initial_collection.get_PC(), *pdb_factory);
      //best_collection=Initial_collection;
      //generate sample states:
      pattern_evaluator->initialize(task);
      pattern_generator->initialize(task);
      pattern_generator->set_pdb_max_size(initial_pdb_size);
      pattern_generator->set_disjunctive_patterns(initial_disjunctive);
      pattern_generator->set_goals_to_add(initial_goals_to_group);
      
      shared_ptr<ModularZeroOnePDBs> candidate_ptr;
      
      //PatternCollectionContainer initial_Gamer_Collection=alternative_pattern_generator.generate();
      //cout<<"Initial Gamer PDB:";initial_Gamer_Collection.print();
      //candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, initial_Gamer_Collection.get_PC(), *pdb_factory);
      if(pdb_factory->is_solved()){
	cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
	result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
	return;
      }
      //cout<<"initial avg_h for Gamer-Style:"<<candidate_ptr->compute_approx_mean_finite_h();
      //cout<<"initial h value for Gamer-Style:"<<candidate_ptr->get_value(initial_state)<<endl;
      //alternative_pattern_generator.check_improv(candidate_ptr->compute_approx_mean_finite_h());//This way we populate initial h value

      
      //result included in sample_states call because needed for dead_end detection
      //set_dead_ends add dead_ends for symbolic, NEED TO ASK ALVARO ABOUT THIS
      //OK, so this was a call to collapse de dead_ends so adding new symbolic PDBs
      //Would not take ages or something like that, needs to check on it and see what 
      //we do now.
      //result->set_dead_ends(pdb_factory->get_dead_ends());
      
      ///DISCUSS WITH ALVARO:adding pdbs to current set if evaluator says new collection is helpful
      //WHEN ADDING THE PDB, terminate_creation makes comparisons biased
      //because new candidate has less time to generate pdb, should we wait for 
      //terminate_pdb till the end of subset selection???
      //result->include_additive_pdbs(candidate_ptr->get_pattern_databases());
      
      //Always adding first collection
      double overall_problem_size=0;
      PatternCollectionContainer perimeter_collection=pattern_generator->generate_perimeter();
      overall_problem_size=perimeter_collection.get_overall_size();
      cout<<"overall_problem_size:"<<overall_problem_size<<endl;
      if(create_perimeter){
        create_perimeter=false;
        cout<<"seeding with creating_perimeter,time:"<<utils::g_timer<<endl;
        cout<<"perimeter PC:";perimeter_collection.print();
        candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, perimeter_collection.get_PC(), *pdb_factory);
        result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases(), 250000, 50000, 10000000));
	result->set_dead_ends(pdb_factory->get_dead_ends());
        //result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases(), 85000, 50000, 10000000));
        if(pdb_factory->is_solved()){
          cout<<"Perimeter solved the problem!"<<endl;
        }
        cout<<"seeding with creating_perimeter finished,time:"<<utils::g_timer<<",h_val:"<<candidate_ptr->get_value(initial_state)<<endl;
        //cout<<"second terminate finished,time:"<<utils::g_timer<<",h_val:"<<candidate_ptr->get_value(initial_state)<<endl;
        pattern_evaluator->set_threshold(1);//If using perimeter, then we want all heuristics which can see further
      }
      else{
	if(always_CBP_or_RBP_or_UCB==0)
	  pattern_generator->set_InSituCausalCheck(true);
	else if(always_CBP_or_RBP_or_UCB==1)
	  pattern_generator->set_InSituCausalCheck(false);

	if(!only_gamer){
	  PatternCollectionContainer Initial_collection=pattern_generator->generate();
	  cout<<"Initial PC:";
	  Initial_collection.print();
	  PatternCollection temp_pc=Initial_collection.get_PC();
	  cout<<"temp_pc.size:"<<temp_pc.size()<<endl;
	  candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, Initial_collection.get_PC(), *pdb_factory);
	  cout<<"initial avg_h:"<<candidate_ptr->compute_approx_mean_finite_h();
	  //pdb_factory->no_terminate_creation(candidate_ptr->get_pattern_databases());//To update pdb_factory in case problem is solved
	  result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
	  result->set_dead_ends(pdb_factory->get_dead_ends());//NOT SURE WHAT IT DOES, CHECK WITH ALVARO
	  //NOTE:Skipped on purpose this for perimeter, because not want to give option to terminate for UBC
	  //I think perimeter PDB would be too big for this to be effective
	  //pdb_ptr_collection.push_back(candidate_ptr);
	}
      }
        
      if(pdb_factory->is_solved()){
          cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
          result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
          return;
      }
      else{
        cout<<"Sadly,solution not found yet after first PDB ;-) so we keep going"<<endl;
      }
      
      
      initial_h=result->get_value(initial_state);
      cout<<"Initial collection zero-one h value:"<<initial_h<<endl;
      
      PC_counter++;
      
      pattern_evaluator->sample_states(result);
      cout<<"first pdb_max_size:"<<pattern_generator->get_pdb_max_size()<<endl;
      
      bool check_to_terminate=false;
      //bool terminate_or_not=true;
      //For Gamer-like local search when doing gamer only we store last successful pattern 
      //in gamer_current_pattern, new_candidate_gamer is the local_search generated
      //pattern which needs to be evaluated
      PatternCollectionContainer gamer_current_pattern, new_candidate_Gamer;
      //separate so clearing it for other bin packing algorithms does not clear partial gamer collections which we might want to try
      PatternCollectionContainer gamer_current_pattern_Gamer;
      PatternCollectionContainer selected_collection_Gamer;
      bool generator_type=true;
      while(!modular_heuristic_timer->is_expired()){
	if(always_CBP_or_RBP_or_UCB==2){
	  generator_type=UCB_RBP_vs_CBP.make_choice();
	  if(generator_type==1){
	    CBP_counter++;
	    pattern_generator->set_InSituCausalCheck(true);
	  }
	  else{
	    RBP_counter++;
	    pattern_generator->set_InSituCausalCheck(false);
	  }
	}
	else if(always_CBP_or_RBP_or_UCB==3){
	  if(rand()%2==1){
	    CBP_counter++;
	    pattern_generator->set_InSituCausalCheck(true);
	  }
	  else{
	    RBP_counter++;
	    pattern_generator->set_InSituCausalCheck(false);
	  }
	}
	 if(double(utils::get_current_memory_in_kb())/1024.0>memory_limit){
	    cout<<"break-3,memory limit breached,current_memory(MB):"<<utils::get_current_memory_in_kb()/1024.0<<",memory_limit:"<<memory_limit<<endl;
	    return;
	 }
        //First we decide whether to try improving existing pdbs or keep generating new ones
        //0 means no terminate, 1 means to try terminating the pdbs
          
        if(unterminated_pdbs){
          int initial_terminate_time=utils::g_timer();
          bool terminate_or_not=terminate_choice.make_choice();
          terminate_or_not=false;
          DEBUG_COMP(cout<<"TERMINATE WHILE SEARCHING FOR PATTERNS DISABLED,DEBUGGING GAMER-STYLE ON ITS OWN"<<endl;);
          DEBUG_COMP(cout<<"terminate_or_not:"<<terminate_or_not<<endl;);
          if(terminate_or_not){
            bool success=false;
            int pdb_counter=0;
            int finished_counter=0;
              auto pdb_collection=result->get_pdbs(); 
              for (auto pdb : *pdb_collection){
                pdb_counter++;
                if(pdb->is_finished()){//looking for first unfinished pdb
                  finished_counter++;
                  continue;
                }
                else{
                  cout<<"pdb is unfinished, trying to terminate it now"<<endl;
                }
                double initial_mean_finite_h=pdb->compute_mean_finite_h();
                cout<<"time:"<<utils::g_timer()<<",average h value before Terminate:"<<initial_mean_finite_h<<endl;
                pdb->terminate_creation(20000,10000,2000000,4000);
                cout<<"time:"<<utils::g_timer()<<",after Terminate average h value:"<<pdb->compute_mean_finite_h()<<endl;
                if(pdb->compute_mean_finite_h()>initial_mean_finite_h){//has terminate improved the pdb
                  if(pdb->is_finished()){
                    result->set_dead_ends(pdb_factory->get_dead_ends());//NOT SURE WHAT IT DOES, CHECK WITH ALVARO
                    cout<<"pdb_finished after terminate,time:"<<utils::g_timer()<<",after setting dead_ends "<<endl;
                    success=true;
                  }
                  //terminate_or_not=false;//already done with the terminate action
                  check_to_terminate=true;//need to check if it was the last pdb to terminate
                  success=true;
                  break;
                }
                //Check if we need to update initial_h_value after terminate
                int new_h_value=result->get_value(initial_state);
                if(result->get_value(initial_state)>initial_h){
                  cout<<"new inital_h_value raised to:"<<new_h_value<<endl;
                  initial_h=new_h_value;
                }
              }
              if(success){
                break;
              }
            if(pdb_counter==finished_counter){
              cout<<"Debug me, trying to terminate but pdb_counter:"<<pdb_counter<<"==finished_counter!!!"<<endl;
            }
            else{
              cout<<"pdb_counter:"<<pdb_counter<<",finished:"<<finished_counter<<endl;
            }
            //Always add cost for the choice,no matter the result!
            terminate_choice.increase_cost(terminate_or_not,utils::g_timer()-initial_terminate_time);
            if(success){//terminate helped
              terminate_choice.increase_reward(terminate_or_not,utils::g_timer()-initial_terminate_time);
              cout<<"terminate helped!"<<endl;
              //Now check if there are more unterminated pdbs
              unterminated_pdbs=false;
            for(auto i : pdb_ptr_collection){
              const PDBCollection & pdb_collection=i->get_pattern_databases(); 
              for (auto pdb : pdb_collection){
                if(pdb->is_finished()){
                  unterminated_pdbs=true;
                  check_to_terminate=false;
                  cout<<"pdb is not terminated,UCB to choose whether to terminate"<<endl;
                  break;
                }
                if(unterminated_pdbs)
                  break;
              }
            }
            continue;//Going back to decide whether to do more terminates
          }
        }
        }

        num_episodes++;

        if(pdb_factory->is_solved()){
          cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
          //best_pdb_collections.push_back(pdb_factory->terminate_creation(candidate.get_pattern_databases()));
          result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
          return;
        }
        //Now choosing between RandomSplit and CBP pattern generation
        int generator_choice=UCB_generator.make_choice();
        if(only_CBP){
          generator_choice=1; cout<<"Forcing generator_choice to non-CGamer-style always, for debugging!"<<endl;
        }
        else if(only_gamer){
          generator_choice=0; cout<<"Forcing generator_choice to Gamer-style always, for debugging!"<<endl;
        }
        DEBUG_COMP(cout<<"time:"<<utils::g_timer()<<",generator_choice:"<<generator_choice<<flush<<endl;);
        float start_time=utils::g_timer();
        float pdb_time=0;
        if(generator_choice==0){//Gamer
	    if(gamer_current_pattern_Gamer.get_size()==0){
	    gamer_current_pattern_Gamer=alternative_pattern_generator.generate();
	    cout<<"time:"<<utils::g_timer<<"Gamer,starting new_gamer_run with num_patterns:"<<gamer_current_pattern_Gamer.get_size()<<endl;
	    //partial_gamer_run=false;
	    improving_patterns.clear();
	    //partial_gamer_improvement_found=false;
	    DEBUG_COMP(cout<<"modular_heuristic, new set of candidate_patterns from Gamer:";gamer_current_pattern_Gamer.print();cout<<endl;);
	  }
	  else if(gamer_current_pattern_Gamer.empty()&&!partial_gamer_improvement_found){
            cout<<"No more Gamer-style pattern selection is possible"<<endl;
            UCB_generator.increase_cost(generator_choice,INT_MAX);//So we do not pick this option again
            if(only_gamer){
              cout<<"Finished, we are doing only_gamer as pattern generator"<<endl;
              return;
            }
            continue;
          }
	  else{
	    //partial_gamer_run=true;
	    cout<<"time:"<<utils::g_timer<<"Gamer,starting partial_gamer_run with num_patterns:"<<gamer_current_pattern_Gamer.get_size()<<endl;
	  }
          DEBUG_COMP(cout<<"modular_heuristic, continuing with remaining set of candidate_patterns from Gamer:";gamer_current_pattern_Gamer.print();cout<<endl;);
        }
        else{
	  gamer_current_pattern.clear();
          //Now getting next pdb size
          //cout<<"previos pdb_max_size:"<<pdb_max_size;
          pdb_max_size=pow(10.0,UCB_sizes.make_choice());
          pattern_generator->set_pdb_max_size(pdb_max_size);
          DEBUG_COMP(cout<<",new pdb_max_size:,"<<pattern_generator->get_pdb_max_size(););
          //bool disjunctive_choice=UCB_Disjunctive_patterns[pdb_max_size].make_choice();
          bool disjunctive_choice=binary_choice.make_choice();
          DEBUG_COMP(cout<<",new disjunctive_choice:"<<disjunctive_choice;);
          pattern_generator->set_disjunctive_patterns(disjunctive_choice);
          if(disjunctive_choice){//1 goal per pattern
            //Note that we may end up with more than one goal per disjunctive patern,
            //All this does is to disable greedy goal grouping, which makes sense to try when
            //variables are reusable among patterns but not if grouping all the goal variables 
            //in one patern means discarding most of variables not in first pattern due to inability
            //to reuse.  However, if goals happen to be dividide between a few patterns, that is not a 
            //problem either, can happen if goal variables get chosen randomly into patterns with existing goals
            num_goals_to_group=1;
          }
          else{
            num_goals_to_group=goals_choice.make_choice();
          }
          DEBUG_COMP(cout<<",new goals_to_group:"<<num_goals_to_group;);
          pattern_generator->set_goals_to_add(num_goals_to_group);

          //Now generating next set of patterns
	  pattern_generator->set_InSituCausalCheck(true);
          gamer_current_pattern=pattern_generator->generate();
          candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, gamer_current_pattern.get_PC(), *pdb_factory);
          pdb_time=utils::g_timer()-start_time;
          DEBUG_COMP(cout<<"pdb_max_size:"<<pdb_max_size<<",pdb_time:"<<pdb_time<<endl;);
          if(pdb_factory->is_solved()){
		  cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
		  result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
		  return;
	  }
        }
          
	if(generator_choice==0){
	  //Pattern selection in traditional Gamer uses Avg-heuristic value
	  //but any of the sampling techniques are applicable as well
	  //For my implementation of the traditional Gamer selection techqnique
	  //use the flag gamer_classic
//	  if(gamer_classic){
//	    cout<<"gamer_classic"<<endl;
//            double new_best_average_h_value=-1;
//            //Raising best_average_h_value by a thousandth so that it knows when 
//            DEBUG_COMP(cout<<"best_average_h_value raised by a thousandth to:"<<best_average_h_value<<",to ensure next set of variables added actually raises the avg_heuristic_value"<<endl;);
//            double candidate_average_h_value=0;
//	    vector<Pattern> patterns=gamer_current_pattern_Gamer.get_PC();
//       std::reverse(patterns.begin(),patterns.end());//So we can look at the same order as Gamer for comparison purposes
//	    float start_gamer_iteration_time=utils::g_timer();
//	    for(auto imprv : improving_patterns){
//	      best_average_h_value=max(best_average_h_value,imprv.first);
//	    }
//	    cout<<"time:"<<utils::g_timer()<<",Gamer_iteration,partial_gamer_run="<<partial_gamer_run<<",partial_gamer_improvement_found="<<partial_gamer_improvement_found<<",starting best_average_h_value:"<<best_average_h_value<<endl;
//            while(patterns.size()>0){
//	      cout<<"time:"<<utils::g_timer()<<",Gamer,patterns.size:"<<patterns.size()<<endl;
//	      Pattern pattern=patterns.back();patterns.pop_back();
//	      cout<<"Gamer,working on pattern:"<<pattern<<endl;
//	      gamer_current_pattern_Gamer.pop_back();
//              if(float(utils::g_timer())/float(modular_time_limit)>0.80&&new_best_average_h_value>best_average_h_value){
//		//PDB generation in Gamer selection can take too long
//                //We want to leave some time to add any improving patterns if found
//                cout<<"modular_time_limit:"<<modular_time_limit<<",time:"<<utils::g_timer<<",break-1,no more Gamer-style candidate PDB generation (or any other) because we are too close to time limit and we want remaining time to add any improving patterns from last iteration"<<endl;
//                break;
//              }
//              PatternCollectionContainer GamerCandidateContainer;GamerCandidateContainer.add_pc(pattern);
//                
//              float indiv_start_time=utils::g_timer();
//              unique_ptr<ModularZeroOnePDBs> Gamer_candidate_ptr=make_unique<ModularZeroOnePDBs>(task_proxy, GamerCandidateContainer.get_PC(),*pdb_factory);
//              //Need to terminate creation to properly evaluate avg_h_value, limit is 30 secs and one million nodes, if bigger it is "deemed" too expensive to consider
//              //Just for selection purposes for me
//              pdb_factory->terminate_creation(Gamer_candidate_ptr->get_pattern_databases(), 20000, 5000, 1000000);
//              float indiv_pdb_time=utils::g_timer()-indiv_start_time;
//              auto pdb=candidate_ptr->get_pattern_databases()[0];
//              DEBUG_COMP(cout<<"\tGamer-style, indiv pdb gen time:"<<indiv_pdb_time<<",finished after terminate:"<<pdb->is_finished()<<",for pattern:"<<pattern<<endl;);
//              if(modular_heuristic_timer->is_expired()){//PDB generation in Gamer selection can take too long
//		cout<<"timer_is_expired,break-2"<<endl;
//                break;
//              }
//              if(indiv_pdb_time>150.0){//so if the time to generate all options is large, bettern to stop now with Gamer candidates
//                cout<<"Stopping gamer-style selection, indiv pdb gen time:,"<<indiv_pdb_time<<",limit:150"<<",pattern:"<<pattern<<endl;
//                UCB_generator.increase_cost(generator_choice,INT_MAX);//And make sure we do not try it again!
//                break;
//              }
//              candidate_average_h_value=Gamer_candidate_ptr->compute_approx_mean_finite_h();
//              DEBUG_COMP(cout<<"time:"<<utils::g_timer()<<",evaluating pattern:"<<pattern<<endl;);
//              if(candidate_average_h_value>best_average_h_value){
//                cout<<"\ttime:,"<<utils::g_timer<<",Gamer-Style,pattern:";for(auto i :pattern) cout<<i<<",";cout<<",raised average_h_value to:,"<<candidate_average_h_value<<"vs prevous iteration best:,"<<best_average_h_value<<endl;
//                new_best_average_h_value=max(candidate_average_h_value,new_best_average_h_value);
//                improving_patterns.push_back(make_pair(candidate_average_h_value,pattern));
//                DEBUG_COMP(cout<<"modular_heuristic,Gamer candidate avg_h_val:"<<best_average_h_value<<",initial h:"<<Gamer_candidate_ptr->get_value(initial_state)<<endl;);
//              }
//              else{
//                DEBUG_COMP(cout<<"candidate_pattern:";for(auto i :pattern) cout<<i<<",";cout<<",avg_h:"<<candidate_average_h_value<<",is smaller than current_avg_h:"<<best_average_h_value<<endl;);
//              }
//
//	      if(utils::g_timer()-start_gamer_iteration_time>100){
//		cout<<"break-3,interrupting gamer candidate generation, taking too long, invested time:"<<utils::g_timer()-start_gamer_iteration_time<<",limit:100"<<endl;
//		break;
//	      }
//            }
//            cout<<"Gamer,removing all added patterns smaller than :"<<0.999*new_best_average_h_value<<",improving_patterns_size:"<<improving_patterns.size()<<endl;
//              
//	    for(auto imprv : improving_patterns){
//	      new_best_average_h_value=max(new_best_average_h_value,imprv.first);
//	    }
//	
//            improving_patterns.erase(std::remove_if(
//                  improving_patterns.begin(), 
//                  improving_patterns.end(),
//                  [new_best_average_h_value](pair<double, Pattern> & x){
//                  return x.first < 0.999*new_best_average_h_value;
//                  }), improving_patterns.end());
//            set<int> Gamer_Combined;
//            for (auto i : improving_patterns){
//              cout<<"\t\tpattern:";
//              for(int j : i.second){ 
//                cout<<j<<",";
//                Gamer_Combined.insert(j);
//              }
//              cout<<"is improving,avg_h:"<<i.first<<endl;
//            }
//            vector<int> Gamer_Combined_vect;for (auto i : Gamer_Combined) Gamer_Combined_vect.push_back(i);
//            if(Gamer_Combined_vect.size()==0){
//              cout<<"Gamer_Combined_vect:"<<Gamer_Combined_vect<<",is empty,so we are finished with Gamer-style for this problem"<<endl;
//              UCB_generator.increase_cost(generator_choice,INT_MAX);
//	      if(only_gamer){
//		cout<<"only-gamer=true, so we are finished building the heuristic"<<endl;
//		return;
//	      }
//              continue;
//            }
//	    if(Gamer_Combined_vect.size()>0&&partial_gamer_run){
//	      partial_gamer_improvement_found=true;
//	      cout<<"Gamer,Marking partial_gamer_improvement_found as true because we just did a partial run where improvements on avg_h value were found"<<endl;
//	    }
//              
//            selected_collection_Gamer.restart_pc(Gamer_Combined_vect);
//	    if(patterns.size()==0){
//	      cout<<"Gamer,confirming pattern:"<<Gamer_Combined<<"as there are no more patterns to test for this iteration"<<endl;
//	      alternative_pattern_generator.confirm(Gamer_Combined);
//	    }
//            candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, selected_collection_Gamer.get_PC(), *pdb_factory);
//	    //To get decent h values when combining all the variables need to terminate the combined vector as well
//	    cout<<"time:,"<<utils::g_timer()<<",finished:"<<pdb_factory->is_finished()<<endl;
//	    auto pdb=candidate_ptr->get_pattern_databases()[0];
//	    pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases(), 50000, 50000, 5000000);
//            best_average_h_value=candidate_ptr->compute_approx_mean_finite_h();
//            cout<<"Gamer_Combined_vect:,"<<Gamer_Combined_vect<<",initial_h:,"<<candidate_ptr->get_value(initial_state)<<",previous initial_h:,"<<initial_h<<",finished after terminate:,"<<pdb->is_finished()<<",new best_average_h_value:,"<<best_average_h_value<<endl;
//
//            pdb_time=utils::g_timer()-start_time;//So includes all the pdb generations for this Gamer step
//	  }
//	  else{//use the evaluator for Gamer patterns as specified by the user
//	    cout<<"gamer_classic is off"<<endl;
	    if(gamer_current_pattern_Gamer.get_size()==0){//starting new gamer pattern
	      gamer_current_pattern_Gamer=alternative_pattern_generator.generate();
	    }
	    //vector<Pattern> patterns=gamer_current_pattern_Gamer.get_PC();
	    cout<<"Starting gamer collection:";gamer_current_pattern_Gamer.print();
            new_candidate_Gamer=pattern_local_search->generate_next_candidate(gamer_current_pattern_Gamer);
	    //If we can not find any new better improvement
	    //We add 1 secs as a penalty everytime this happens if mixing local_search_methods, to avoid loops
	    //otherwise UCB would keep choosing local_search on the same pattern repeatedly until the time spent 
	    //is enough to prefer an alternative.
	    //If only gamer, we simply return, no improvements possible by adding a single var.
	    //Note: we could keep looking for better pattenrs by adding pairs, n-tupples, etc.  That is future work.
	    if(new_candidate_Gamer.get_top_pattern()==gamer_current_pattern_Gamer.get_top_pattern()){
	      if(only_gamer){
		cout<<"No improvement possible, and only doing gamer, so we are finished."<<endl;
		break;//breaking because we want to do the terminate of the selected PDB 
		//in case it is not finished
	      }
	      else{
		UCB_local_search.increase_cost(GAMER_LOCAL_SEARCH,1.0);
	      }
	    }
            candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, new_candidate_Gamer.get_PC(), *pdb_factory);
 //         }
	}
	else{
            UCB_sizes.increase_cost(log10(pattern_generator->get_pdb_max_size()),pdb_time);
            binary_choice.increase_cost(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
            goals_choice.increase_cost(double(pattern_generator->get_goals_to_add()),pdb_time);
	    UCB_RBP_vs_CBP.increase_cost(generator_type,pdb_time);
          }
          UCB_generator.increase_cost(generator_choice,pdb_time);
          //UCB_Disjunctive_patterns[pdb_max_size].increase_cost(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
          //UCB_goals_to_group[pdb_max_size].increase_cost(double(pattern_generator->get_goals_to_add()),pdb_time);
          PC_counter++;
          if(pdb_factory->is_solved()){
		  cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
		  result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
		  return;
	  }

          new_initial_h=candidate_ptr->get_value(initial_state);
	  initial_h=result->get_value(initial_state);
	  //cout<<"Initial value before evaluations:,"<<initial_h<<",new_initial_h:"<<new_initial_h<<endl;
          
          if(generator_choice==0){
            if(new_initial_h>initial_h){
              cout<<"Gamer,increased initial_h from:,"<<initial_h<<",to,"<<new_initial_h<<endl;
            }
          }
          
          if (initial_h == numeric_limits<int>::max()) {
            cout<<"initial state is dead_end according to PDB, problem unsolvable!!!"<<endl;
            exit(1);
          }
          else if(new_initial_h>initial_h){//we always add the collection and re-sample if initial_h increased
	    if(generator_choice==0){
	      cout<<"Gamer improved initial h value from:,"<<initial_h<<",to,"<<new_initial_h<<endl;
	      gamer_current_pattern_Gamer=new_candidate_Gamer;
	      //We just added a new variable, so all previously 
	      //tested vars had to be given a new chance
	      pattern_local_search->reset_forbidden_vars();
	      if(only_gamer){//recompute PDBs, very likely new collection fully dominates previous ones
		cout<<"time:"<<utils::g_timer<<",calling recompute_max_additive_subsets"<<endl;
		result->recompute_max_additive_subsets();
		cout<<"time:"<<utils::g_timer<<",after recompute_max_additive_subsets"<<endl;
	      }
	    }
            //int temp_initial_h=candidate_ptr->get_value(initial_state);
            //cout<<"time:"<<utils::g_timer()<<",Initial h value before terminate:"<<temp_initial_h<<endl;
            result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
            cout<<"time:,"<<utils::g_timer()<<"pdb_max_size:,"<<pdb_max_size<<",generator_choice:,"<<generator_choice<<",goals_choice:"<<num_goals_to_group<<",disjoint:"<<disjunctive_choice<<",Selecting PC and resampling because initial_h has been raised from "<<initial_h<<"to "<<new_initial_h<<endl;
            initial_h=result->get_value(initial_state);//Might be higher than simply new cadidate collection value due to max additive pattern combinations beyond current collection, unlikely but possible!
	      
	    //Clean dominated PDBs after addition of improving patterns
	    if(recompute_additive_sets){
	      cout<<"time:"<<utils::g_timer<<",calling recompute_max_additive_subsets"<<endl;
	      result->recompute_max_additive_subsets();
	      cout<<"time:"<<utils::g_timer<<",after recompute_max_additive_subsets"<<endl;
	    }
            
	    check_to_terminate=true;
            pdb_ptr_collection.push_back(candidate_ptr);

            //UCB_Disjunctive_patterns[pdb_max_size].increase_reward(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
            //UCB_goals_to_group[pdb_max_size].increase_reward(double(pattern_generator->get_goals_to_add()),pdb_time);
            UCB_generator.increase_reward(generator_choice,pdb_time);
            if(generator_choice!=0){//Not Gamer options
              binary_choice.increase_reward(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
              goals_choice.increase_reward(double(pattern_generator->get_goals_to_add()),pdb_time);
              UCB_sizes.increase_reward(log10(pdb_max_size),pdb_time);
	      UCB_RBP_vs_CBP.increase_reward(generator_type,pdb_time);
            }
            result->set_dead_ends(pdb_factory->get_dead_ends());//NOT SURE WHAT IT DOES, CHECK WITH ALVARO
            pattern_evaluator->sample_states(result);
          }
          else{//OK,so lets check if candidate_PC is good enough to add to current collection
            if(pattern_evaluator->evaluate(candidate_ptr)){
	      if(generator_choice==0){
		gamer_current_pattern_Gamer=new_candidate_Gamer;
		pattern_local_search->reset_forbidden_vars();
		cout<<"Gamer improved prunning but initial_h still,"<<initial_h<<endl;
		if(only_gamer){//recompute PDBs, very likely new collection fully dominates previous ones
		  cout<<"time:"<<utils::g_timer<<",calling recompute_max_additive_subsets"<<endl;
		  result->recompute_max_additive_subsets();
		  cout<<"time:"<<utils::g_timer<<",after recompute_max_additive_subsets"<<endl;
		}
	      }
              //NEED TO CHECK WITHOUT TERMINATING CREATION UNTIL ALL PDBs ARE SELECTED
              cout<<"time:"<<utils::g_timer()<<"pdb_max_size:"<<pdb_max_size<<",generator_choice:"<<generator_choice<<",disjoint:"<<disjunctive_choice<<",goals_choice:"<<num_goals_to_group<<",modular_heuristic_selecting PC"<<endl;
              result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
	 
	      //Clean dominated PDBs after addition of improving patterns

	    if(recompute_additive_sets){
	      cout<<"time:"<<utils::g_timer<<",calling recompute_max_additive_subsets"<<endl;
	      result->recompute_max_additive_subsets();
	      cout<<"time:"<<utils::g_timer<<",after recompute_max_additive_subsets"<<endl;
	    }
              check_to_terminate=true;
              if(generator_choice!=0){
                UCB_sizes.increase_reward(log10(pdb_max_size),pdb_time);
                //UCB_Disjunctive_patterns[pdb_max_size].increase_reward(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
                //UCB_goals_to_group[pdb_max_size].increase_reward(double(pattern_generator->get_goals_to_add()),pdb_time);
                binary_choice.increase_reward(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
                goals_choice.increase_reward(double(pattern_generator->get_goals_to_add()),pdb_time);
		UCB_RBP_vs_CBP.increase_reward(generator_type,pdb_time);
              }
              UCB_generator.increase_reward(generator_choice,pdb_time);
              result->set_dead_ends(pdb_factory->get_dead_ends());//NOT SURE WHAT IT DOES, CHECK WITH ALVARO
              //should we always re-sample or only if number of improved states is large enough?
              pattern_evaluator->sample_states(result);
            }
            else{
              DEBUG_COMP(cout<<"time:"<<utils::g_timer()<<"pdb_max_size:"<<pdb_max_size<<",modular_heuristic_not_selecting_PC"<<endl;);
	      if(generator_choice==0){
		cout<<"time:"<<utils::g_timer()<<"pdb_max_size:"<<pdb_max_size<<",modular_heuristic_not_selecting_PC"<<endl;
	      }
	      if(only_gamer){
		//Need to check if latest PDB was fully grown
		//Then we remove last var from set of testable vars
		//until an improving gaming pattern is actually found.
		if(candidate_ptr->is_finished()){
		  pattern_local_search->forbid_last_var();
		  cout<<"\tforbidding last_var, PDB was finished,last_var:,";pattern_local_search->print_last_var();cout<<endl;
		}
		else{
		  cout<<"\tnot forbidding last_var, PDB was not finished,last_var:,";pattern_local_search->print_last_var();cout<<endl;
		  candidate_ptr->print();
		}
	      }
            }
          }
          //Checking if any of the pdbs were not terminated
          //If they are, allow UCB to choose whether to terminate them instead of generating new pdbs
          if(check_to_terminate&&unterminated_pdbs==false){
            check_to_terminate=false;
            unterminated_pdbs=false;
              auto pdb_collection=result->get_pdbs(); 
              for (auto pdb : *pdb_collection){
                if(!pdb->is_finished()){
                  unterminated_pdbs=true;
                  check_to_terminate=false;
                  cout<<"pdb is not terminated,UCB to choose whether to terminate"<<endl;
                  break;
                }
              }
          }
      }
      cout<<"CBP calls:,"<<pattern_generator->get_CBP_calls()<<",RBP calls:,"<<pattern_generator->get_RBP_calls()<<endl;
      //Now terminate creation of unfinished selected PDBs
      float start_time=utils::g_timer();
      auto pdb_collection=result->get_pdbs(); 
      for(auto i : boost::adaptors::reverse(*pdb_collection)){//better to start terminating the best PDBs
        //usually the last added, in case we run out of time
        if(i){
          int temp_h=i->get_value(initial_state);
          cout<<"time:"<<utils::g_timer()<<",Initial pdb h value before Terminate:"<<temp_h<<endl;
          i->terminate_creation(60000,20000,2000000,memory_limit);
          temp_h=i->get_value(initial_state);
          cout<<"time:"<<utils::g_timer()<<",Initial pdb h value after Terminate:"<<temp_h<<endl;
        }
        if((utils::g_timer()-start_time)>200.0){
          cout<<"time:"<<utils::g_timer()<<",terminate_time:"<<utils::g_timer()-start_time<<",Interrupting terminate because it is taking too long"<<endl;
          break;
        }
      }
      //Check if any collection is useless because it is dominated
      //pattern_evaluator->clear_dominated_heuristics(result,&);
	    
      float terminate_time=utils::g_timer()-start_time;
      cout<<"time:"<<utils::g_timer()<<",before_recompute_max_additive_subset,Testing modular_heuristic constructor finished,time:"<<utils::g_timer()<<",episodes:"<<num_episodes<<",PC created:"<<PC_counter<<",final_pdbs:"<<result->get_patterns()->size()<<",terminate_time:"<<terminate_time<<endl;
      cout<<"CBP_counter:,"<<CBP_counter<<",RBP_counter:,"<<RBP_counter<<endl;
      //result->recompute_max_additive_subsets();
      //cout<<"time:"<<utils::g_timer()<<",after recompute_max_additive_subset,Testing modular_heuristic constructor finished,time:"<<utils::g_timer()<<",episodes:"<<num_episodes<<",PC created:"<<PC_counter<<",final_pdbs:"<<result->get_patterns()->size()<<",terminate_time:"<<terminate_time<<endl;

    }

int ModularHeuristic::compute_heuristic(const GlobalState &global_state) {
    State state = convert_global_state(global_state);
    int h = result->get_value(state);
    if (h == numeric_limits<int>::max()) {
        return DEAD_END;
    } else {
        return h;
    }
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

void ModularHeuristic::clear_dominated_heuristics(){
  //no point checking if we have only one collection
  if(result->get_pdbs()->size()<2){
    return;
  }
  std::shared_ptr<PatternCollectionInformation> new_result;
	      
  new_result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
  pattern_evaluator->clear_dominated_heuristics(result,new_result);
}
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
    parser.add_option<int>(
        "modular_time_limit",
        "time limit in seconds for modular_pdb_heuristic initialization cut off",
        "900");
    parser.add_option<int> (
	"always_cbp_or_rbp_or_ucb", 
	"If 0,ensure any added variable is causally connected(CBP), if 1, only check causal connection after all patterns are selected(RBP), if 2, use UCB to learn which is better, if 3 do 50/50", 
	"0");
    parser.add_option<shared_ptr<PDBFactory>>(
        "pdb_factory",
        "See detailed documentation for pdb factories. ",
	      "modular_symbolic");
    parser.add_option<bool>(
        "terminate_creation",
        "give extra generation time to selected PDBs but not to candidate PDBs",
	      "false");
    parser.add_option<bool>(
        "create_perimeter",
        "do initial perimeter",
	      "false");
    parser.add_option<bool>(
        "only_gamer",
        "only gamer-style w/wo Perimeter",
	      "false");
    parser.add_option<bool>(
        "only_cbp",
        "only CBP-style w/wo Perimeter",
	      "false");
    parser.add_option<bool>(
        "gamer_classic",
        "only gamer-style w avg_h value as selector",
	      "false");
    parser.add_option<shared_ptr<PatternCollectionLocalSearch>>(
        "local_search",
        "pattern Collection local search method",
        "local_search_gamer");
    
    Heuristic::add_options_to_parser(parser);

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;

    return new ModularHeuristic(opts);
}

static Plugin<Heuristic> _plugin("modular_pdb", _parse);
}
