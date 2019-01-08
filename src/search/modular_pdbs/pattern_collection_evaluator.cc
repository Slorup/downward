#include "pattern_collection_evaluator.h"

#include "../option_parser.h"
#include "../plugin.h"
//#include "../task_proxy.h"

#include "../task_tools.h"
#include "../utils/timer.h"
#include <climits>
#include "pattern_database_interface.h"


    
using namespace std;
//PatterCollectionEvaluator is to be the driver for 8 PDB-based options
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
//    virtual void initialize(std::shared_ptr<AbstractTask> task) {
//    cout << "Manual pattern collection: " << *patterns << endl;
//    return PatternCollectionInformation(task, patterns);
//
//   }

  void PatternCollectionEvaluator::clear_dominated_heuristics(std::shared_ptr<PatternCollectionInformation> current_result,std::shared_ptr<PatternCollectionInformation> &new_result,
      shared_ptr<ModularZeroOnePDBs> candidate_ptr){
      double start_time=utils::g_timer();
      if(unique_samples.size()<100){//Better to do nothing if unique_samples size is too small
	new_result=current_result;
	return;
      }
      vector<int> current_best_h_values;
      current_best_h_values.reserve(unique_samples.size());


      cout<<"time:"<<utils::g_timer()<<",calling clear_dominated_heuristics with "<<current_result->get_max_additive_subsets()->size()+1<<" best heuristics and unique_samples:"<<unique_samples.size()<<endl;
      
      //First we get all the values for sampled states with the latest PC
      int i=0;
      for(map<size_t,pair<State,int> >::iterator it=unique_samples.begin(); it!=unique_samples.end();it++){
	current_best_h_values[i++]=candidate_ptr->get_value(it->second.first);
      }
      cout<<"\ttime to populate current_best_h_values:"<<utils::g_timer()-start_time<<",unique_samples:"<<i<<endl;

      //Now go through each additive subset and check if they are dominated for the sampled set of states
      shared_ptr<MaxAdditivePDBSubsets> current_max_additive_subsets=new_result->get_max_additive_subsets();

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
	}
	else{
	  //cout<<"collection["<<i<<"] is dominated,eliminating "<<endl;
	}
	i++;
      }
      cout<<"clear_dominated_heuristics finished:"<<new_result->get_max_additive_subsets()->size()<<","<<",time:"<<utils::g_timer()-start_time<<endl;
  }
    
  int PatternCollectionEvaluator::calculate_max_additive_subset(PDBCollection max_subset,State current_state){
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


static options::PluginTypePlugin<PatternCollectionEvaluator> _type_plugin(
"PatterCollectionEvaluator",
"The various pattern evaluation algorithms usable by the complementary_modular_pdbs heurisitic.");
}
