#include "pattern_collection_local_search_GA.h"
#include "../sampling.h"
#include "../utils/timer.h"
#include "../task_tools.h"
#include "../utils/countdown_timer.h"
#include "../options/option_parser.h"
#include "../options/plugin.h"

#include "../causal_graph.h"
//#include "../globals.h"
//#include "../task_proxy.h"

//#include "../utils/markup.h"
//#include "../utils/math.h"
#include "../utils/rng.h"
//#include "../utils/timer.h"
//#include "../heuristic.h"

#include "pdb_factory.h"
#include "../utils/debug_macros.h"
#include "pattern_collection_generator_complementary.h"
#include "../utils/collections.h"

    
using namespace std;
using namespace utils;
namespace pdbs3 {
  PatternCollectionLocalSearchGA::PatternCollectionLocalSearchGA(const options::Options & opts)
  :time_limit (opts.get<int>("time_limit")),
  episodes (opts.get<int>("episodes")){
  //PatternCollectionLocalSearchGA::PatternCollectionLocalSearchGA() 
      cout<<"hello LocalSearchGA,time_limit:"<<time_limit<<",episodes:"<<episodes<<endl;
    //num_vars=task->get_num_variables();
  }

  PatternCollectionContainer PatternCollectionLocalSearchGA::generate_next_candidate(PatternCollectionContainer candidate_collection){
    //cout<<"hello generate_next_candidate_GA"<<flush<<endl;
    //If more than one pattern, we only do local search on the first pattern
    new_patterns=candidate_collection;
    int mutation_count=mutate(new_patterns);
    TaskProxy task_proxy2(*(g_root_task()));
    int num_vars = task_proxy2.get_variables().size();
    cout<<"\t\tmutations:"<<mutation_count<<",patterns:"<<candidate_collection.get_size()<<",num_vars:"<<num_vars<<endl;
    // cout<<"old_patterns:";candidate_collection.print(); cout<<"new_patterns:";new_patterns.print();
    return new_patterns;//Not adding collection
  }

  /*void PatternCollectionLocalSearchGA::forbid_last_var(){
    forbidden_vars.insert(last_var);
  }
  void PatternCollectionLocalSearchGA::reset_forbidden_vars(){
    forbidden_vars.clear();
  }*/
  
  int PatternCollectionLocalSearchGA::mutate(PatternCollectionContainer& candidate_collection) {
    //cout<<"hello mutate"<<flush<<endl;
    //cout<<"Input Collection:";candidate_collection.print();
    int mutations=0;
    vector<int> int_pattern,removed_vars;
    vector<bool> bool_pattern;
    PatternCollection collection=candidate_collection.get_PC();
    PatternCollection new_collection;
    vector<vector<bool> > bit_collection;
    for (auto pattern : collection) {
      transform_to_pattern_bit_vector_form(bool_pattern,pattern);
      bit_collection.push_back(bool_pattern);
      //cout<<"\tbool_pattern:"<<bool_pattern<<endl;
    }
    //cout<<"finished creating bit_collection "<<flush<<endl;
      
    for (vector<bool> &pattern : bit_collection) {
      for (size_t k = 0; k < pattern.size(); ++k) {
	double random = (*g_rng())(); // [0..1)
	if (random < mutation_rate) {
	  mutations++;
	  //cout<<"doing mutation,random:,"<<random<<",mutation_rate:,"<<mutation_rate<<flush<<endl;
	  pattern[k].flip();
	}
      }

      transform_to_pattern_int_vector_form(pattern,int_pattern);
      //cout<<"Before removing irrelevant_variables:"<<int_pattern<<endl;
      remove_irrelevant_variables(int_pattern,removed_vars);
      
      if(int_pattern.size()==0){
	new_collection=collection;
	return 0;
      }
      new_collection.push_back(int_pattern);
      //cout<<"\tadding clean_pattern:";for(auto i : int_pattern) cout<<i<<",";cout<<endl;
    }
    //cout<<"finished processing bit_collection "<<flush<<endl;
    
    if(mutations==0){//No changes were made so we enforce at least one!
      int random = (*g_rng())(bit_collection.size()); // [0..1)
      size_t random_pattern=random;
      vector<bool> pattern=bit_collection[random_pattern];
      size_t random_variable=(*g_rng())(pattern.size());
      pattern[random_variable].flip();
      transform_to_pattern_int_vector_form(pattern,int_pattern);
      remove_irrelevant_variables(int_pattern,removed_vars);
      new_collection[random_pattern]=int_pattern;

      /*cout<<"old_collection:";
      for (auto pat : collection)
	for(auto i : pat) 
	  cout<<i<<",";
      cout<<endl;

      cout<<"new_collection:";
      for (auto collection : new_collection){
	for(auto i : collection) {
	  cout<<i<<",";
	}
	cout<<endl;
      }
      cout<<endl;*/
      //cerr<<"Need to make sure new_collection is updated correctly when mutations==0"<<endl;exit(1);
    }

    /*cout<<"mutations:"<<mutations<<",old_collection:"<<endl;
    candidate_collection.print();
    cout<<"new_collection:"<<endl;
    PatternCollectionContainer tempPCC(new_collection);tempPCC.print();*/

    candidate_collection=new_collection;
    return mutations;
  }
  
  void PatternCollectionLocalSearchGA::transform_to_pattern_bit_vector_form(vector<bool> &bitvector,const vector<int> pattern) {
	bitvector.assign(g_variable_name.size(), false);
	for (size_t i = 0; i < pattern.size(); ++i) {
	    bitvector[pattern[i]]=true;
	}
  }
  
  void PatternCollectionLocalSearchGA::transform_to_pattern_int_vector_form(const vector<bool> &bitvector, vector<int> &pattern){
    pattern.clear();
    for (size_t i = 0; i < bitvector.size(); ++i) {
      if (bitvector[i]){
	pattern.push_back(i);
      }
    }
  }

  static shared_ptr<PatternCollectionLocalSearch>_parse(options::OptionParser &parser) {
    parser.add_option<int> ("time_limit", "Time limit for local search", "100");
    parser.add_option<int> ("episodes", "how many tries to find an improving neighour in pattern space.","60");
    options::Options options = parser.parse();
    parser.document_synopsis(
        "Pattern Generator Local Search module",
        "Selection of variables to generate Pattern Collection");
    options::Options opts = parser.parse();
    if (parser.dry_run())
        return 0;

    return make_shared<PatternCollectionLocalSearchGA>(opts);
  }
  
  static options::PluginShared<PatternCollectionLocalSearch> _plugin("local_search_ga", _parse);
  

}
