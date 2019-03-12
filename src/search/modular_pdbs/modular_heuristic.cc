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
#include "pattern_database_symbolic.h"
#include "pattern_collection_generator_GamerStyle.h"
#include "pattern_collection_local_search.h"
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
    gamer_classic(opts.get<bool>("gamer_classic")), 
    gamer_excluded(opts.get<bool>("gamer_excluded")), 
    doing_local_search(opts.get<bool>("doing_local_search")), 
    doing_canonical_search(opts.get<bool>("doing_canonical_search")), 
    pdb_factory (opts.get<shared_ptr<PDBFactory>>("pdb_factory")) {
      cout<<"Hi nonagnostic_v2"<<endl;
      cout<<"modular_time_limit:"<<modular_time_limit<<endl;
      cout<<"terminate_creation:"<<terminate_creation<<endl;
      cout<<"pdb_type:"<<pdb_factory->name()<<endl;
      cout<<"gamer_excluded:"<<gamer_excluded<<",always_cbp_or_rbp_or_ucb:"<<always_CBP_or_RBP_or_UCB<<endl;
      cout<<"doing_local_search:"<<doing_local_search<<endl;
      cout<<"doing_canonical_search:"<<doing_canonical_search<<endl;
      

      TaskProxy task_proxy(*task);
      const State &initial_state = task_proxy.get_initial_state();
      //need result here to store final PDB collection
      result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());

      //DEBUG BLOCK, testing recompute//
//      PatternCollectionContainer PC;
//      Pattern temp_pattern1{2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,20,21,22,23,25,26,27,28,29,30,31,32,33,34,35,36,37,38,40,41,42,44,45,46,48,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,72,73,74,75,76,79,80,81,82,83,84,85,86,88,89,90,91,92,93,94,95,96,98,100,101,102,104,105,106,107,109,110,111,112,113,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,138,139,140,141,142,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,190,191,192,193,194,195,196,197,198,199,200,201,203,204,205,206,207,208,209,210,211,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,231,233,234,235,236,237,238,239,240,242,243,244,245,246,247,248,249,250,251,252,253,255,256,258,261,262,263,264,265,266,268};
//      Pattern temp_pattern2{1,16,19,24,39,43,47,49,71,77,78,87,97,99,103,108,114,137,143,165,166,189,202,212,229,230,232,241,254,257,259,260,267};
//      Pattern temp_pattern3{1,2,3,4,5,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,24,26,27,29,30,31,32,33,34,35,36,37,38,40,43,44,45,46,47,48,49,51,52,54,55,56,58,59,60,61,62,63,64,65,66,68,70,71,73,74,75,76,77,79,80,81,82,83,84,87,88,90,91,92,93,94,95,96,97,99,100,101,102,103,104,105,106,107,108,109,110,111,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,150,151,152,153,154,155,156,158,159,162,163,164,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,218,219,220,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,239,240,241,242,243,244,246,247,249,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,268};
//      Pattern temp_pattern4{6,7,23,25,28,39,41,42,50,53,57,67,69,72,78,85,86,89,98,112,149,157,160,161,165,200,217,221,237,238,245,248,250};
//      Pattern temp_pattern5{1,2,4,5,8,9,11,12,13,14,15,17,18,19,20,22,24,26,27,29,30,31,32,33,34,35,36,37,38,40,42,43,44,45,46,47,48,49,51,52,54,55,58,59,60,61,62,63,64,66,68,70,73,74,75,76,77,79,80,81,82,83,84,87,88,89,90,91,92,93,94,95,96,97,99,101,102,103,104,106,107,108,109,110,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,131,132,133,134,135,136,137,139,141,143,144,145,146,147,148,150,151,152,153,154,155,156,158,159,162,164,166,167,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,185,186,187,189,190,191,193,194,195,196,197,198,199,201,202,203,204,205,206,208,209,210,211,212,213,214,215,216,218,219,220,222,224,225,226,227,228,229,230,231,232,233,234,236,239,240,241,242,243,244,246,247,249,251,252,254,255,256,257,258,259,260,261,262,263,265,266,267};
//      Pattern temp_pattern6{6,7,16,23,25,27,41,42,50,52,53,57,67,69,73,78,83,85,86,88,90,93,98,107,112,118,124,139,141,149,157,160,161,165,200,204,207,213,217,220,221,222,237,238,246,248,249,258,261,264,266};
//
//      PC.clear();PC.add_pc(temp_pattern1);PC.add_pc(temp_pattern2);
//      candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, PC.get_PC(), *pdb_factory);
//      cout<<"PC1:,initial h_value before terminate:"<<candidate_ptr->get_value(initial_state)<<flush<<endl;
//      result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
//      cout<<"PC1:,initial h_value after terminate:"<<candidate_ptr->get_value(initial_state)<<",result:"<<result->get_value(initial_state)<<flush<<endl;
//
//      PC.clear();PC.add_pc(temp_pattern3);PC.add_pc(temp_pattern4);
//      candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, PC.get_PC(), *pdb_factory);
//      cout<<"PC2:,initial h_value before terminate:"<<candidate_ptr->get_value(initial_state)<<flush<<endl;
//      result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
//      cout<<"PC2:,initial h_value after terminate:"<<candidate_ptr->get_value(initial_state)<<",result:"<<result->get_value(initial_state)<<flush<<endl;
//      
//      PC.clear();PC.add_pc(temp_pattern5);PC.add_pc(temp_pattern6);
//      candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, PC.get_PC(), *pdb_factory);
//      cout<<"PC3:,initial h_value before terminate:"<<candidate_ptr->get_value(initial_state)<<flush<<endl;
//      result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
//      cout<<"PC3:,initial h_value after terminate:"<<candidate_ptr->get_value(initial_state)<<",result:"<<result->get_value(initial_state)<<flush<<endl;
//
//      cout<<"time:"<<utils::g_timer()<<",initial_h before recompute:,"<<result->get_value(initial_state)<<endl;
//      result->recompute_max_additive_subsets();
//      cout<<"time:"<<utils::g_timer()<<",initial_h after recompute:,"<<result->get_value(initial_state)<<endl;
//      exit(1);

      //DEBUG BLOCK FINISHED
      PatternCollectionContainer PC;
      //Pattern temp_pattern1{2,5,11,14,19,23,27,40,42,44,45,46,50,51,52,53,54,56,58,59,60,61};
      //Pattern temp_pattern1{7,8,9,10,11,12,13,14,15,16,17,18};
      //Pattern temp_pattern1{0,1,2,3,4,5,6,7,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79};
      //Pattern temp_pattern1{0,1,2,3,4,5,6,7,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79};
	//Pattern temp1a{0,3,6,7,28,32,33,37,38,42,43,47,58,60,63,64,65,66,67,69,71,73,74,76,77,78,79};
	//Pattern temp2a{0,3,6,7,28,32,33,37,38,42,43,47,58,60,63,64,65,67,68,69,71,73,74,76,77,78,79};
	//Pattern temp3a{,3,6,7,28,32,33,37,38,42,43,47,58,60,63,64,65,67,69,70,71,73,74,76,77,78,79};
	//Pattern temp_pattern1{0,3,6,7,28,32,33,37,38,42,43,47,58,60,63,64,65,66,67,68,69,70,71,72,73,74,76,77,78,79};
        //Pattern temp_pattern1{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50};
        /*PC.clear();PC.add_pc(temp_pattern1);
      cndidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, PC.get_PC(), *pdb_factory);
      cout<<"initial h value for Gamer-Style before_terminate:"<<candidate_ptr->get_value(initial_state)<<endl;
      result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
      cout<<"initial h value for Gamer-Style after:"<<candidate_ptr->get_value(initial_state)<<endl;
      result->set_dead_ends(pdb_factory->get_dead_ends());
      if(doing_canonical_search){
	canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,false, 0, 0);//no need to prune anything, solution found
      }
      return;*/
      
      
      //unsigned num_goals_to_group=0;
      modular_heuristic_timer =make_unique<utils::CountdownTimer>(modular_time_limit);
      Options opts2=opts;
      Options opts3=opts;
      Options opts4=opts;
      Options opts5=opts;
      opts2.set<int>("time_limit",40);
      opts3.set<int>("packer_selection",0);
      opts4.set<int>("packer_selection",1);
      opts5.set<int>("packer_selection",2);
      //double best_average_h_value=0;//For Gamer-Style selection
      //1 means Random split into two patterns, 0 means CBP
      //Initialize reward to 2 so it does not go too fast for initial selection
      //Learning UCB_RBP_vs_CBP; UCB_RBP_vs_CBP.insert_choice(1); UCB_RBP_vs_CBP.insert_choice(0);
      //UCB_RBP_vs_CBP.increase_reward(1,10);UCB_RBP_vs_CBP.increase_reward(0,10);

      map<double,Learning> UCB_Disjunctive_patterns;//Disjunctive (or not) pattern selector,one per size;
      //UCB_sizes.insert_choice(pow(10,4));
      //UCB_sizes.insert_choice(8);
      UCB_generator.insert_choice(GAMER_GENERATOR); UCB_generator.insert_choice(RBP_CBP_GENERATOR);
      UCB_generator.increase_reward(RBP_CBP_GENERATOR,4);UCB_generator.increase_reward(GAMER_GENERATOR,4);
      binary_choice.insert_choice(1);binary_choice.insert_choice(0);
      UCB_local_search.insert_choice(GAMER_LOCAL_SEARCH);
      //Initialize reward to 2 so it does not go too fast for initial selection
      binary_choice.increase_reward(1);binary_choice.increase_reward(0);
      size_t overall_num_goals=task_proxy.get_goals().size();
      cout<<"overall_num_goals:"<<overall_num_goals<<endl;
      
      for (size_t i=0;i<overall_num_goals;i++){
        goals_choice.insert_choice(i+1);
      //Initialize reward to 2 so it does not go too fast for initial selection
        goals_choice.increase_reward(i+1);
      }
      map<double,Learning> UCB_goals_to_group;//Disjunctive (or not) pattern selector,one per size;
      
      terminate_choice.insert_choice(1);terminate_choice.insert_choice(0);
      //terminate_choice.increase_reward(0,10);//biasing towards not terminating in the begining
      
      //We divide the pdb_max_size limit in 10 groups, each a 10% bigger relative to 
      //the overall problem_size, and we round up.  if problem size<10, we add as many categories
      //until we reach the max size limit
      double overall_problem_size=0;
      pattern_generator->initialize(task);
      
      cout<<"initializing local_search"<<flush<<endl;
      if(doing_local_search)
	pattern_local_search->initialize(task);
      cout<<"finished initializing local_search"<<flush<<endl;
     
      PatternCollectionContainer perimeter_collection=pattern_generator->generate_perimeter();
      PatternCollectionContainer candidate_collection;
      overall_problem_size=perimeter_collection.get_overall_size();
      max_size_step=log10(overall_problem_size)/10.0;
      cout<<"overall_problem_size:"<<overall_problem_size<<",max_size_step:"<<max_size_step<<endl;
      for(int i=1;i<11;i++){
	if(i*max_size_step<4)//skipping tiny patterns!
	  continue;
	UCB_sizes.insert_choice(int(i*max_size_step)+1);
	cout<<"added to UCB_sizes size_limit:"<<int(i*max_size_step)+1<<",overall_problem_size:"<<overall_problem_size<<endl;
	if((i*max_size_step)>log10(overall_problem_size)){
	  cout<<"no more max_size_steps, overall_problem_size too small"<<endl;
	  break;
	}
      }
      //Need to at least have one option!
      if(max_size_step<4||UCB_sizes.size()<1){
	UCB_sizes.insert_choice(log10(overall_problem_size)+1);
      }

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
      double initial_goals_to_group=goals_choice.make_choice();
      cout<<"after goals_choice"<<endl;
      pdb_max_size=initial_pdb_size;
      cout<<"initial_pdb_size:"<<initial_pdb_size<<endl;
    
      cout<<"initial pdb type:"<<pdb_factory->name()<<endl;

      //ModularZeroOnePDBs candidate(task_proxy, Initial_collection.get_PC(), *pdb_factory);
      //best_collection=Initial_collection;
      //generate sample states:
      pattern_evaluator->initialize(task);
      pattern_generator->set_pdb_max_size(initial_pdb_size);
      pattern_generator->set_disjunctive_patterns(initial_disjunctive);
      pattern_generator->set_goals_to_add(initial_goals_to_group);
      
      //shared_ptr<ModularZeroOnePDBs> candidate_ptr;
      
      //PatternCollectionContainer initial_Gamer_Collection=alternative_pattern_generator.generate();
      //cout<<"Initial Gamer PDB:";initial_Gamer_Collection.print();
      //candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, initial_Gamer_Collection.get_PC(), *pdb_factory);
      if(pdb_factory->is_solved()){
	cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
	result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
	result->set_dead_ends(pdb_factory->get_dead_ends());
	if(doing_canonical_search){
	  canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,false, 0, 0);//no need to prune anything, solution found
	}
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
      
      //For Gamer-like local search when doing gamer only we store last successful pattern 
      //in gamer_current_pattern, new_candidate_gamer is the local_search generated
      //pattern which needs to be evaluated
      //separate so clearing it for other bin packing algorithms does not clear partial gamer collections which we might want to try
      
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
	if(always_CBP_or_RBP_or_UCB==ALWAYS_CBP){
	  pattern_generator->set_InSituCausalCheck(true);
	}
	else if(always_CBP_or_RBP_or_UCB==ALWAYS_RBP){
	  pattern_generator->set_InSituCausalCheck(false);
	}

	PatternCollectionContainer Initial_collection=pattern_generator->generate();
	cout<<"Initial PC:";
	Initial_collection.print();
	PatternCollection temp_pc=Initial_collection.get_PC();
	cout<<"temp_pc.size:"<<temp_pc.size()<<endl;
	candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, Initial_collection.get_PC(), *pdb_factory);
	cout<<"initial avg_h:"<<candidate_ptr->compute_approx_mean_finite_h();

	int initial_h_pre_terminate=candidate_ptr->get_value(initial_state);
	result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
	result->set_dead_ends(pdb_factory->get_dead_ends());
	int initial_h_post_terminate=candidate_ptr->get_value(initial_state);
	if(initial_h_post_terminate>initial_h_pre_terminate){
	  cout<<"terminate when adding PDB has raised initial h value from:"<<initial_h_pre_terminate<<",to:,"<<initial_h_post_terminate<<endl;
	}
      
      if(pdb_factory->is_solved()){
	    cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
	    result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
	    result->set_dead_ends(pdb_factory->get_dead_ends());
	    if(doing_canonical_search){
	      canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,false, 0, 0);//no need to prune anything, solution found
	    }
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
      }
      
      if(modular_heuristic_timer->is_expired()){
	if(doing_canonical_search){
	  canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,true, 0, 0);
	}
	return;
      }

      
      
      //bool terminate_or_not=true;
      generator_type=true;

      cout<<"FINDING IMPROVEMENTES"<<flush<<endl;
      if(!pdb_factory->is_solved()){
	bool improvement_found=find_improvements(modular_time_limit);
	cout<<"time:,"<<utils::g_timer()<<",improvement_found:"<<improvement_found<<endl;
	cout<<"CBP calls:,"<<pattern_generator->get_CBP_calls()<<",RBP calls:,"<<pattern_generator->get_RBP_calls()<<endl;
	cout<<"CBP_counter:,"<<CBP_counter<<",RBP_counter:,"<<RBP_counter<<",Disj_counter:,"<<Disj_counter<<",Not_Disj_counter:"<<Not_Disj_counter<<endl;
      }
      else{
	cout<<"skipping looking for improvement because valid optimal solution already found in pdb_factory"<<endl;
      }
    }

int ModularHeuristic::compute_heuristic(const GlobalState &global_state) {
      //cout<<"using compute_heuristic"<<flush<<endl;
    State state = convert_global_state(global_state);
    int h=0;
    if(doing_canonical_search){
      //cout<<"using canonical_search instead"<<flush<<endl;
      h = canonical_pdbs->get_value(state);
    }
    else {
      //cout<<"using result instead"<<flush<<endl;
      h = result->get_value(state);
    }

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

//This does a cannonical check, alternatively we
//can clear quasi-dominated heuristics with in situ
//sampling test
void ModularHeuristic::clear_dominated_heuristics_in_situ_sampling(shared_ptr<ModularZeroOnePDBs> candidate_ptr){
  float start_time=utils::g_timer();
  //no point checking if we have only one collection
  if(result->get_max_additive_subsets()->size()<2){
    return;
  }
    std::shared_ptr<PatternCollectionInformation> new_result;
    new_result=make_shared<PatternCollectionInformation>(task, make_shared<PatternCollection>());
    pattern_evaluator->clear_dominated_heuristics(result,new_result,candidate_ptr);
    cout<<"clear_dominated_heuristics,time_spent,"<<utils::g_timer()-start_time<<",Reduced number of collections from:,"<<result->get_max_additive_subsets()->size()<<",to:,"<<new_result->get_max_additive_subsets()->size()<<endl;
    result=new_result;
}
    
int ModularHeuristic::recompute_heuristic(const GlobalState& current_state){
  return compute_heuristic(current_state);
}
//This is to continue search for better patterns
//Usually called when we have spent too long without finding solution 
//To start with we pass search_time as limit to improvement_found
bool ModularHeuristic::find_improvements(int time_limit) {
  cout<<"helo find_impovements,time_limit:,"<<time_limit<<flush<<endl;
  modular_heuristic_timer =make_unique<utils::CountdownTimer>(time_limit);
      const State &initial_state = task_proxy.get_initial_state();
      bool improvement_found=false;
      while(!modular_heuristic_timer->is_expired()){
	cout<<"find_improvements,remaining_time:"<<modular_heuristic_timer->get_remaining_time()<<",elapsed_time:,"<<modular_heuristic_timer->get_elapsed_time()<<endl;
	//cout<<"\tcurrent_h_modular:"<<result->get_value(initial_state)<<endl;
	if(always_CBP_or_RBP_or_UCB==ALWAYS_UCB){
	  //generator_type=UCB_RBP_vs_CBP.make_choice();
	  if(generator_type==1){
	    CBP_counter++;
	    pattern_generator->set_InSituCausalCheck(true);
	  }
	  else{
	    RBP_counter++;
	    pattern_generator->set_InSituCausalCheck(false);
	  }
	}
	else if(always_CBP_or_RBP_or_UCB==ALWAYS_50_50){
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
	    cout<<"time:"<<utils::g_timer()<<",break-3,memory limit breached,current_memory(MB):"<<utils::get_current_memory_in_kb()/1024.0<<",memory_limit:"<<memory_limit<<endl;
	    if(doing_canonical_search&&improvement_found){//No point callling to recompute if no improvement has been found
	      if(recompute_additive_sets){
		cout<<"time:"<<utils::g_timer()<<",initial_h before recompute:,"<<result->get_value(initial_state)<<endl;
		result->recompute_max_additive_subsets();
		cout<<"time:"<<utils::g_timer()<<",initial_h after recompute:,"<<result->get_value(initial_state)<<endl;
	      }
		    
	    //Now make canonical_symbolic_pdbs
	    //CanonicalSymbolicPDBs canonical_pdbs (PatternCollectionInformation & info,bool dominance_pruning, int compress_nodes, int compress_time);
	      cout<<"time:"<<utils::g_timer()<<",initial_h before canonical_pdbs:,"<<result->get_value(initial_state)<<endl;
	      canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,true, 0, 0);
	      cout<<"time:"<<utils::g_timer()<<",initial_h after canonical_pdbs:,"<<canonical_pdbs->get_value(initial_state)<<endl;
	      cout<<"canonical_pdbs_count:"<<canonical_pdbs->count_pdbs()<<endl;
	    }
	    result->set_dead_ends(pdb_factory->get_dead_ends());//NOT SURE WHAT IT DOES, CHECK WITH ALVARO
	    return improvement_found;
	 }
        //First we decide whether to try improving existing pdbs or keep generating new ones
        //0 means no terminate, 1 means to try terminating the pdbs
        if(unterminated_pdbs){
	  //cout<<"unterminated_pdbs"<<endl;
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
                    result->set_dead_ends(pdb_factory->get_dead_ends());
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
	  result->set_dead_ends(pdb_factory->get_dead_ends());
	  if(doing_canonical_search){
	    canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,false, 0, 0);//no need to prune anything, solution found
	  }
          return true;
        }
	  
	//Now choosing between RandomSplit and CBP pattern generation
	int generator_choice=0;
	if(gamer_excluded){
	  generator_choice=RBP_CBP_GENERATOR;
	}
	else{
	  generator_choice=UCB_generator.make_choice();
	}

	float pdb_time=0;
	bool disjunctive_choice=false;
	
      DEBUG_COMP(cout<<"time:"<<utils::g_timer()<<",generator_choice:"<<generator_choice<<flush<<endl;);
      float start_time=utils::g_timer();
    
      //Now getting next pdb size
      //cout<<"previos pdb_max_size:"<<pdb_max_size;
      pdb_max_size=pow(10.0,UCB_sizes.make_choice());
      cout<<"pdb_max_size_chosen;"<<pdb_max_size<<endl;
      //We add a random amount to pdb_max_size if the max_step_size>1 so that 
      //we increase the diversity of patterns generated in big problems
      if(max_size_step>1){
	increase_level=floor((float(rand()%100)/100.0)*max_size_step);
	//cout<<"increase_level:"<<increase_level<<endl;
      }
      pattern_generator->set_pdb_max_size(pdb_max_size+pow(10.0,increase_level));
      DEBUG_COMP(cout<<",new pdb_max_size:,"<<pattern_generator->get_pdb_max_size(););
      //bool disjunctive_choice=UCB_Disjunctive_patterns[pdb_max_size].make_choice();
      disjunctive_choice=binary_choice.make_choice();
      DEBUG_COMP(cout<<",new disjunctive_choice:"<<disjunctive_choice;);
      pattern_generator->set_disjunctive_patterns(disjunctive_choice);
      if(disjunctive_choice)
       Disj_counter++;
      else
	Not_Disj_counter++;
      /*if(disjunctive_choice){//1 goal per pattern
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
      }*/
      //DEBUG_COMP(cout<<",new goals_to_group:"<<num_goals_to_group;);
      //pattern_generator->set_goals_to_add(num_goals_to_group);

      //Now generating next set of patterns and PDB
      if(pattern_generator->get_name()=="GamerStyle"){
	if(doing_local_search){
	  pattern_local_search->do_local_search(result,pattern_evaluator,pdb_factory,modular_heuristic_timer->get_remaining_time());
	  if(check_for_solution()){
	    return true;
	  }
	}
      }
      else{
	candidate_collection=pattern_generator->generate(); 
	candidate_ptr=make_shared<ModularZeroOnePDBs>(task_proxy, candidate_collection.get_PC(), *pdb_factory);
	pdb_time=utils::g_timer()-start_time;
	
	DEBUG_COMP(cout<<"pdb_max_size:"<<pdb_max_size<<",pdb_time:"<<pdb_time<<endl;);
	UCB_sizes.increase_cost(log10(pdb_max_size),pdb_time);
	binary_choice.increase_cost(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
	goals_choice.increase_cost(double(pattern_generator->get_goals_to_add()),pdb_time);
	//UCB_RBP_vs_CBP.increase_cost(generator_type,pdb_time);

	//Now evaluate and update UCB parameters as necessary
	//Adding a minmum of 0.5 in case the PDB size limit is so low
	//that we keep generating almost empty patterns but still calling it
	//because the generation time is almost nill or also because the pattern
	//is already stored
	UCB_generator.increase_cost(generator_choice,max(pdb_time,float(0.5)));

	  //UCB_Disjunctive_patterns[pdb_max_size].increase_cost(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
	  //UCB_goals_to_group[pdb_max_size].increase_cost(double(pattern_generator->get_goals_to_add()),pdb_time);
	  PC_counter++;
	  if(pdb_factory->is_solved()){
		  cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
		  result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
		  result->set_dead_ends(pdb_factory->get_dead_ends());
		  if(doing_canonical_search){
		    canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,false, 0, 0);//no need to prune anything, solution found
		  }
		  return true;
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
	      cerr<<"initial state is dead_end according to PDB, problem unsolvable!!!"<<endl;
	      exit(1);
	    }
	    else if(new_initial_h>initial_h){//we always add the collection and re-sample if initial_h increased
	      improvement_found=true;
	      cout<<"time:,"<<utils::g_timer()<<",improvement found"<<endl;
	      
	      cout<<"time:,"<<utils::g_timer()<<"pdb_max_size:,"<<pdb_max_size<<",generator_choice:,"<<generator_choice<<",disjoint:,"<<disjunctive_choice<<",Selecting PC and resampling because initial_h has been raised from,"<<initial_h<<"to:,";

	      result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
	      result->set_dead_ends(pdb_factory->get_dead_ends());
	  
	      if(recompute_additive_sets){//recompute PDBs, very likely new collection fully dominates previous ones
		cout<<"time:"<<utils::g_timer<<",calling recompute_max_additive_subsets"<<endl;
		result->recompute_max_additive_subsets();
		cout<<"time:"<<utils::g_timer<<",after recompute_max_additive_subsets"<<endl;
	      }
	      //int temp_initial_h=candidate_ptr->get_value(initial_state);
	      //cout<<"time:"<<utils::g_timer()<<",Initial h value before terminate:"<<temp_initial_h<<endl;
	      if(doing_dominated_sampling)
		clear_dominated_heuristics_in_situ_sampling(candidate_ptr);
	      
	      initial_h=result->get_value(initial_state);//Might be higher than simply new cadidate collection value due to max additive pattern combinations beyond current collection, unlikely but possible!
	      
	      cout<<initial_h;
	      
		
	      //Clean dominated PDBs after addition of improving patterns
		//RECOMPUTE_MAX_ADDITIVE_SUBSETS IS BUGGY
		//SO USE clear_dominated_heuristics_in_situ_sampling instead
		//if(only_gamer){//recompute PDBs, very likely new collection fully dominates previous ones
	      //if(recompute_additive_sets){
	      //  cout<<"time:"<<utils::g_timer<<",calling recompute_max_additive_subsets"<<endl;
	       // result->recompute_max_additive_subsets();
		//cout<<"time:"<<utils::g_timer<<",after recompute_max_additive_subsets"<<endl;
	      //}
	      
	      check_to_terminate=true;
	      //Need to reorganize pdb_ptr_collection if choosing to terminate
	      //whenever we clear_dominate heuristics
	      //pdb_ptr_collection.push_back(candidate_ptr);

	      //UCB_Disjunctive_patterns[pdb_max_size].increase_reward(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
	      //UCB_goals_to_group[pdb_max_size].increase_reward(double(pattern_generator->get_goals_to_add()),pdb_time);
	      UCB_generator.increase_reward(generator_choice,pdb_time);
	      if(generator_choice!=0){//Not Gamer options
		binary_choice.increase_reward(double(pattern_generator->get_disjunctive_patterns()),pdb_time);
		goals_choice.increase_reward(double(pattern_generator->get_goals_to_add()),pdb_time);
		UCB_sizes.increase_reward(log10(pdb_max_size),pdb_time);
		UCB_RBP_vs_CBP.increase_reward(generator_type,pdb_time);
	      }
	      //Always trying to further improve any already improving
	      //pattern collection,set by doing_local_search
	      if(doing_local_search){
		pattern_local_search->do_local_search(result,pattern_evaluator,pdb_factory,modular_heuristic_timer->get_remaining_time());
		if(check_for_solution()){
		  return true;
		}
	      }
	      else{//do_local_search does sample_states
		pattern_evaluator->sample_states(result);
	      }
	      if(pdb_factory->is_solved()){
		cout<<"local_search found Solution, exiting"<<endl;
		if(doing_canonical_search){
		  canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,false, 0, 0);//no need to prune anything, solution found
		}
		result->set_dead_ends(pdb_factory->get_dead_ends());
		return true;
	      }

	    }
	    else{//OK,so lets check if candidate_PC is good enough to add to current collection
	      if(pattern_evaluator->evaluate(candidate_ptr)){
		improvement_found=true;
		cout<<"time:,"<<utils::g_timer()<<",improvement found"<<endl;
		
		result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
		result->set_dead_ends(pdb_factory->get_dead_ends());
		
		if(recompute_additive_sets){//recompute PDBs, very likely new collection fully dominates previous ones
		  cout<<"time:"<<utils::g_timer<<",calling recompute_max_additive_subsets"<<endl;
		  result->recompute_max_additive_subsets();
		  cout<<"time:"<<utils::g_timer<<",after recompute_max_additive_subsets"<<endl;
		}
		//NEED TO CHECK WITHOUT TERMINATING CREATION UNTIL ALL PDBs ARE SELECTED
		//cout<<"time:"<<utils::g_timer()<<"pdb_max_size:"<<pdb_max_size<<",generator_choice:"<<generator_choice<<",disjoint:"<<disjunctive_choice<<",goals_choice:"<<num_goals_to_group<<",modular_heuristic_selecting PC"<<endl;
		cout<<"time:"<<utils::g_timer()<<"pdb_max_size:"<<pdb_max_size<<",generator_choice:"<<generator_choice<<",disjoint:"<<disjunctive_choice<<",modular_heuristic_selecting PC"<<endl;
		
		//IN SITU SAMPLING DOMINATION CHECK
		if(doing_dominated_sampling)
		  clear_dominated_heuristics_in_situ_sampling(candidate_ptr);
		
		check_to_terminate=true;
		result->set_dead_ends(pdb_factory->get_dead_ends());
		//should we always re-sample or only if number of improved states is large enough?
		pattern_evaluator->sample_states(result);
		
		//Always trying to further improve any already improving
		//pattern collection,set by doing_local_search
		if(doing_local_search){
		  pattern_local_search->do_local_search(result,pattern_evaluator,pdb_factory,modular_heuristic_timer->get_remaining_time());
		  if(check_for_solution()){
		    return true;
		  }
		}
		else{//do_local_search does sample_states
		  pattern_evaluator->sample_states(result);
		}
	      }
	      else{
		DEBUG_COMP(cout<<"time:"<<utils::g_timer()<<",pdb_max_size:"<<pdb_max_size<<",modular_heuristic_not_selecting_PC"<<endl;);
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
      }
      
	 
      //Now terminate creation of unfinished selected PDBs
      float start_time=utils::g_timer();
      auto pdb_collection=result->get_pdbs(); 
      for(auto i : boost::adaptors::reverse(*pdb_collection)){//better to start terminating the best PDBs
        //usually the last added, in case we run out of time
        if(i){
          int temp_h=i->get_value(initial_state);
	  if(!i->is_finished()){//looking for first unfinished pdb
	    cout<<"pdb unterminated,time:"<<utils::g_timer()<<",Initial pdb h value before Terminate:"<<temp_h<<endl;
	    i->terminate_creation(60000,20000,2000000,memory_limit);
	    temp_h=i->get_value(initial_state);
	    cout<<"time:"<<utils::g_timer()<<",Initial pdb h value after Terminate:"<<temp_h<<endl;
	  }
        }
        if((utils::g_timer()-start_time)>200.0){
          cout<<"time:"<<utils::g_timer()<<",terminate_time:"<<utils::g_timer()-start_time<<",Interrupting terminate because it is taking too long"<<endl;
          break;
        }
      }
	    
      float terminate_time=utils::g_timer()-start_time;
      cout<<"time:"<<utils::g_timer()<<",before_recompute_max_additive_subset,Testing modular_heuristic constructor finished,episodes:"<<num_episodes<<",PC created:"<<PC_counter<<",final_pdbs:"<<result->get_patterns()->size()<<",terminate_time:"<<terminate_time<<endl;

      if(doing_canonical_search){
	//Now recalculate additive subsets
	if(recompute_additive_sets){
	  cout<<"time:"<<utils::g_timer()<<",initial_h before recompute:,"<<result->get_value(initial_state)<<endl;
	  result->recompute_max_additive_subsets();
	  cout<<"time:"<<utils::g_timer()<<",after recompute_max_additive_subset,Testing modular_heuristic constructor finished,episodes:"<<num_episodes<<",PC created:"<<PC_counter<<",final_pdbs:"<<result->get_patterns()->size()<<endl;
	  cout<<"time:"<<utils::g_timer()<<",initial_h after recompute:,"<<result->get_value(initial_state)<<endl;
	}
	
      //Now make canonical_symbolic_pdbs
      //CanonicalSymbolicPDBs canonical_pdbs (PatternCollectionInformation & info,bool dominance_pruning, int compress_nodes, int compress_time);
	cout<<"time:"<<utils::g_timer()<<",initial_h before canonical_pdbs:,"<<result->get_value(initial_state)<<endl;
	canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,true, 0, 0);
	cout<<"time:"<<utils::g_timer()<<",initial_h after canonical_pdbs:,"<<canonical_pdbs->get_value(initial_state)<<endl;
	cout<<"canonical_pdbs_count:"<<canonical_pdbs->count_pdbs()<<endl;
      }
      result->set_dead_ends(pdb_factory->get_dead_ends());
      return improvement_found;
}
bool ModularHeuristic::check_for_solution(){
  if(pdb_factory->is_solved()){
    cout<<"Solution found while generating PDB candidate of type:"<<pdb_factory->name()<<", adding PDB and exiting generation at time"<<utils::g_timer()<<endl;
    result->include_additive_pdbs(pdb_factory->terminate_creation(candidate_ptr->get_pattern_databases()));
    result->set_dead_ends(pdb_factory->get_dead_ends());
    if(doing_canonical_search){
      canonical_pdbs=make_unique<CanonicalSymbolicPDBs>(result,false, 0, 0);//no need to prune anything, solution found
    }
    return true;
  }
  return false;
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
      parser.add_option<shared_ptr<PatternCollectionLocalSearch>>(
	  "local_search",
	  "pattern Collection local search method",
	  "local_search_ga(time_limit=60)");
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
        "gamer_excluded",
        "No Gamer pattern generation",
	      "false");
    parser.add_option<bool>(
        "gamer_classic",
        "only gamer-style w avg_h value as selector",
	      "false");
    parser.add_option<bool>(
        "doing_local_search",
        "Use local_search to try to furtherly improve any PCs we have selected.",
	"false");
    parser.add_option<bool>(
        "doing_dominated_sampling",
        "Use local sampling to eliminate quasi-dominated PCs.",
	"false");
    parser.add_option<bool>(
        "doing_canonical_search",
        "Combining PDBs in a canonical set",
	"true");
    
    Heuristic::add_options_to_parser(parser);

    Options opts = parser.parse();
    if (parser.dry_run())
        return nullptr;

    return new ModularHeuristic(opts);
}

static Plugin<Heuristic> _plugin("modular_pdb", _parse);
}
