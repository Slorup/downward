#ifndef MODULAR_PDBS_LEARNING_H
#define MODULAR_PDBS_LEARNING_H
#include<map>
#include"../utils/debug_macros.h"

using namespace std;

class data_point{
  double cost;
  double reward;
  double score;

  public:
  data_point(double _cost, double _reward, double _score){
    cost=_cost;
    reward=_reward;
    score=_score;
  }
  data_point(){
    cost=1;
    reward=1;
    score=1;
  }
  void increase_reward(){
    reward++;
  }

  void increase_reward(double _reward){//In case we want to weight rewards,e.g. pdb size weighted by gen time
    reward+=_reward;
  }

  void increase_cost(){
    cost++;
  }
  void increase_cost(double increase){//when using other cost metrics, e.g. time
    cost+=increase;
  }
  void set_score(double _score){
    score=_score;
  }
  double get_reward(){
    return reward;
  }
  double get_cost(){
    return cost;
  }
  double get_score(){
    return score;
  }
};
//for putting data_points on a set
//All we care about for placing ucb data_points is the reference(choice unique identifier)
//, which should be unique!

class Learning {
  //first element is option, e.g. pdb_max_size, second is associated UCB reward
  map<double,data_point> choices;
  double total_cost=0;

  public:
  int size(){
    return choices.size();
  }
  void update_scores(){//Using UCB1 to update scores, need to call after each selection to update scores
    //for next selection
    for (auto& i : choices){
      if(i.second.get_cost()>0){//keep initial value until at least one call
        i.second.set_score((i.second.get_reward()/i.second.get_cost())+sqrt(2*log(i.second.get_cost())/total_cost));
      }
    }
  }
  void update_specific_score(double ref){
    auto search = choices.find(ref);
    if(search != choices.end()) {
      //cout<<"updated score from "<<search->second.get_score();
      search->second.set_score((search->second.get_reward()/search->second.get_cost())+sqrt(2*log(search->second.get_cost())/total_cost));
      //cout<<" to:"<<search->second.get_score()<<endl;
    }
    else{
      cerr<<"Reference "<<ref<<" not in existing choices in learning, pls debug me!!!"<<endl;
      exit(1);
    }
  }
  void increase_cost(double reference,double increase=1.0){
    auto search = choices.find(reference);
    if(search != choices.end()) {
      DEBUG_MSG(cout<<"updated cost from "<<search->second.get_score(););
      search->second.increase_cost(increase);
      //cout<<" to:"<<search->second.get_cost()<<" for reference:"<<reference<<endl;
      total_cost+=increase;
    }
    else {
        std::cerr << "reference value: " << reference << "Does not exist in set of data points, debug me!!!"<<endl;
        exit(1);
    }
  }
void insert_choice(double reference){
  choices[reference];
}

void increase_reward(double reference,double increase=1){
  auto search = choices.find(reference);
  if(search != choices.end()) {
    search->second.increase_reward(increase);
  }
  else{
    std::cerr << "reference value: " << reference << "Does not exist in set of data points, debug me!!!"<<endl;
    exit(1);
  }
}

double make_choice(bool print=false){
  DEBUG_COMP(cout<<"making_choice,avaliable choices:"<<choices.size()<<endl;);
  if(choices.size()==0){
    cerr<<"PLS DEBUG ME, making a choice but none are available!!!"<<endl;
    exit(1);
  }
  if(total_cost>0)
    update_scores();
  double max_score=0;
  vector<double> winners;
  for (auto& i : choices){
    if(print)
      cout<<endl<<"\treference:"<<i.first<<",score:"<<i.second.get_score()<<",cost:"<<i.second.get_cost()<<",reward:"<<i.second.get_reward()<<endl;
    if (i.second.get_score()>max_score){
      max_score=max(max_score,i.second.get_score());
    }
  }
  for (auto& i : choices){
    if (i.second.get_score()==max_score){
      winners.push_back(i.first);
    }
  }
  size_t chosen_one=rand()%winners.size();
  DEBUG_COMP(cout<<"max_score:"<<max_score<<",winners:";);
  DEBUG_COMP(for (auto i : winners) cout<<i<<",";);
  DEBUG_COMP(cout<<",chosen:"<<winners[chosen_one]<<endl;);
  //choices[winners[chosen_one]].increase_cost();//We assume this choice is going to be used!
  //total_cost++;
  //Responsability for calculating reward outside of class
  //Then caller also has to update reward fields after choice run finished
  //And finally call the update_scores function before next call
  return winners[chosen_one];
}





};


#endif
