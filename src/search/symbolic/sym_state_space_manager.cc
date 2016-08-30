#include "sym_state_space_manager.h"

#include "sym_enums.h" 
#include "../utils/debug_macros.h"
#include <queue>
#include <limits>
#include <algorithm>

#include "sym_util.h"
#include "../globals.h"
#include "../global_operator.h"
#include "../mutex_group.h"
#include "../utils/timer.h"
#include "../options/options.h"
#include "../options/option_parser.h"
#include "../abstract_task.h"

using namespace std;

namespace symbolic {

    SymStateSpaceManager::SymStateSpaceManager(shared_ptr<SymStateSpaceManager> & parent,
					       const std::set<int>& relevantVars)
 :
	vars(parent->vars), p(parent->p), cost_type(parent->cost_type), 
	parent_mgr(parent),
	fullVars(relevantVars), 
	min_transition_cost(parent->min_transition_cost), 
	hasTR0(parent->hasTR0), mutexInitialized(false), 
	mutexByFluentInitialized(false) {
    
	for(size_t i = 0; i < g_variable_name.size(); i++){
	    if (!fullVars.count(i)){
		nonRelVars.insert(i);
	    }
	}
    }

    SymStateSpaceManager::SymStateSpaceManager(SymVariables * v,
					       const SymParamsMgr & params, 
					       OperatorCost cost_type_) :
	vars(v), p(params), cost_type(cost_type_),
	initialState(v->zeroBDD()), goal(v->zeroBDD()), 
	min_transition_cost(0), hasTR0(false), 
	mutexInitialized(false), 
	mutexByFluentInitialized(false)  {

	for(const auto & op : g_operators){
	    //if (op.is_dead()) continue;       
      
	    if(min_transition_cost == 0 || min_transition_cost > get_adjusted_action_cost(op, cost_type)){
		min_transition_cost = get_adjusted_action_cost(op, cost_type);      
	    }
	    if(get_adjusted_action_cost(op, cost_type) == 0){
		hasTR0 = true;
	    }
	}
    }


    // SymStateSpaceManager::SymStateSpaceManager(SymVariables * v, const SymParamsMgr & params, 
    // 					       OperatorCost cost_type_, const set<int> & /*relVars*/) :
    // 	vars(v), p(params), cost_type(cost_type_),
    // 	initialState(v->zeroBDD()), goal(v->zeroBDD()), 
    // 	min_transition_cost(0), hasTR0(false), 
    // 	mutexInitialized(false), 
    // 	mutexByFluentInitialized(false)  {

    // 	for(const auto & op : g_operators){
    // 	    //if (op.is_dead()) continue;       
      
    // 	    if(min_transition_cost == 0 || min_transition_cost > get_adjusted_action_cost(op, cost_type)){
    // 		min_transition_cost = get_adjusted_action_cost(op, cost_type);      
    // 	    }
    // 	    if(get_adjusted_action_cost(op, cost_type) == 0){
    // 		hasTR0 = true;
    // 	    }
    // 	}
    // }


    // SymStateSpaceManager::SymStateSpaceManager(SymStateSpaceManager * mgr,
    // 					       const SymParamsMgr & params) : vars(mgr->vars), p(params), cost_type(mgr->cost_type), 
    // 									      parentMgr(mgr),  initialState(mgr->zeroBDD()), 
    // 									      goal(mgr->zeroBDD()), 
    // 									      min_transition_cost(0), hasTR0(false), 

    // 									      mutexInitialized(false), mutexByFluentInitialized(false)   {  
    // 	if(mgr){
    // 	    min_transition_cost = mgr->getMinTransitionCost();
    // 	    hasTR0 = mgr->hasTransitions0();
    // 	}else{
    // 	    for(const auto & op : g_operators){
    // 		//if (op.is_dead()) continue; 

    // 		if(min_transition_cost == 0 || min_transition_cost >  get_adjusted_action_cost(op, cost_type)){
    // 		    min_transition_cost = get_adjusted_action_cost(op, cost_type);      
    // 		}
    // 		if(get_adjusted_action_cost(op, cost_type) == 0){
    // 		    hasTR0 = true;
    // 		}
    // 	    }
    // 	}
    // }

    void SymStateSpaceManager::init_transitions_from_individual_trs () {
	if(!transitions.empty()) return; //Already initialized!
	DEBUG_MSG(cout << "Init transitions" << endl;);

	DEBUG_MSG(cout << "Generate individual TRs" << endl;);
	transitions = map<int, vector <SymTransition>> (indTRs); //Copy
	DEBUG_MSG(cout << "Individual TRs generated" << endl;);
	min_transition_cost = 0;
	hasTR0 = transitions.count(0) > 0;
  
	for(map<int, vector<SymTransition> >::iterator it = transitions.begin(); 
	    it != transitions.end(); ++it){
	    merge(vars, it->second, mergeTR, p.max_tr_time, p.max_tr_size);

	    if(min_transition_cost == 0 || min_transition_cost > it->first){
		min_transition_cost = it->first;      
	    }

	    cout << "TRs cost=" << it->first << " (" << it->second.size() << ")" << endl;
	}
    }


    void SymStateSpaceManager::init_mutex(const std::vector<MutexGroup> & mutex_groups,
					  bool genMutexBDD, bool genMutexBDDByFluent){

	//Check if I should initialize something and return
	if(mutexInitialized) genMutexBDD = false;
	if(mutexByFluentInitialized) genMutexBDDByFluent = false;
	if(!genMutexBDD && !genMutexBDDByFluent)
	    return;
	if(genMutexBDD) mutexInitialized = true;
	if(genMutexBDDByFluent) mutexByFluentInitialized = true;

	if(genMutexBDDByFluent){
	    //Initialize structure for exactlyOneBDDsByFluent (common to both init_mutex calls) 
	    exactlyOneBDDsByFluent.resize(g_variable_domain.size());
	    for (size_t i = 0; i < g_variable_domain.size(); ++i){
		exactlyOneBDDsByFluent[i].resize(g_variable_domain[i]); 
		for(int j = 0; j < g_variable_domain[i]; ++j){
		    exactlyOneBDDsByFluent[i][j] = oneBDD();
		}
	    }
	}
  
	init_mutex(mutex_groups, genMutexBDD, genMutexBDDByFluent, false);
	init_mutex(mutex_groups, genMutexBDD, genMutexBDDByFluent, true);
    }

    void SymStateSpaceManager::init_mutex(const std::vector<MutexGroup> & mutex_groups,
					  bool genMutexBDD, bool genMutexBDDByFluent, bool fw){
	DEBUG_MSG(cout << "Init mutex BDDs " << (fw ? "fw" : "bw") << ": "
		  << genMutexBDD << " " << genMutexBDDByFluent << endl;);

	vector<vector<BDD>> & notMutexBDDsByFluent = 
	    (fw ? notMutexBDDsByFluentFw : notMutexBDDsByFluentBw);
  
	vector<BDD> & notMutexBDDs = 
	    (fw ? notMutexBDDsFw : notMutexBDDsBw);
  
	//BDD validStates = vars->oneBDD();
	int num_mutex = 0;
	int num_invariants = 0;

	if(genMutexBDDByFluent){
	    //Initialize structure for notMutexBDDsByFluent 
	    notMutexBDDsByFluent.resize(g_variable_domain.size());
	    for (size_t i = 0; i < g_variable_domain.size(); ++i){
		notMutexBDDsByFluent[i].resize(g_variable_domain[i]); 
		for(int j = 0; j < g_variable_domain[i]; ++j){
		    notMutexBDDsByFluent[i][j] = oneBDD();
		}
	    }
	}
  
	//Initialize mBDDByVar and invariant_bdds_by_fluent
	vector<BDD>  mBDDByVar;
	mBDDByVar.reserve(g_variable_domain.size());
	vector<vector<BDD>> invariant_bdds_by_fluent (g_variable_domain.size());
	for(size_t i = 0; i < invariant_bdds_by_fluent.size(); i++){
	    mBDDByVar.push_back(oneBDD());
	    invariant_bdds_by_fluent[i].resize(g_variable_domain[i]);
	    for(size_t j = 0; j < invariant_bdds_by_fluent[i].size(); j++){
		invariant_bdds_by_fluent[i][j] = oneBDD();
	    }
	}
  
	for (auto & mg : mutex_groups){
	    if(mg.pruneFW() != fw)
		continue;
	    const vector<Fact> & invariant_group = mg.getFacts();
	    DEBUG_MSG(cout << mg << endl;);
	    if(mg.isExactlyOne()){
		BDD bddInvariant = zeroBDD();
		int var = numeric_limits<int>::max();
		int val = 0;
		bool exactlyOneRelevant = true;

		for(auto & fluent : invariant_group){
		    if(!isRelevantVar(fluent.var)){
			exactlyOneRelevant = true;
			break;
		    }
		    bddInvariant += vars->preBDD(fluent.var, fluent.value);
		    if(fluent.var < var){
			var = fluent.var;
			val = fluent.value;
		    }
		}

		if(exactlyOneRelevant){ 
		    num_invariants++;	
		    if(genMutexBDD){
			invariant_bdds_by_fluent[var][val] *= bddInvariant;
		    }
		    if(genMutexBDDByFluent){
			for(auto & fluent : invariant_group){
			    exactlyOneBDDsByFluent[fluent.var][fluent.value] *= bddInvariant;
			}
		    }
		}
	    }
  
  
	    for (size_t i = 0; i < invariant_group.size(); ++i){
		int var1 = invariant_group[i].var;
		if(!isRelevantVar(var1)) continue;
		int val1 = invariant_group[i].value;
		BDD f1 = vars->preBDD(var1, val1);

		for (size_t j = i+1; j < invariant_group.size(); ++j){
		    int var2 = invariant_group[j].var;
		    if(!isRelevantVar(var2)) continue;
		    int val2 = invariant_group[j].value;
		    BDD f2 = vars->preBDD(var2, val2);
		    BDD mBDD = !(f1*f2);
		    if(genMutexBDD){
			num_mutex++;
			mBDDByVar[min(var1, var2)] *= mBDD;
			if(mBDDByVar[min(var1, var2)].nodeCount() > p.max_mutex_size){
			    notMutexBDDs.push_back(mBDDByVar[min(var1, var2)]);
			    mBDDByVar[min(var1, var2)] = vars->oneBDD();
			}
		    }
		    if(genMutexBDDByFluent){
			notMutexBDDsByFluent[var1][val1] *= mBDD;
			notMutexBDDsByFluent[var2][val2] *= mBDD;
		    }
		}
	    }
	}

	if(genMutexBDD){
	    for(size_t var = 0; var < g_variable_domain.size(); ++var){
		if(!mBDDByVar[var].IsOne()){
		    notMutexBDDs.push_back(mBDDByVar[var]);
		}
		for (const BDD & bdd_inv : invariant_bdds_by_fluent[var]){
		    if(!bdd_inv.IsOne()){
			notMutexBDDs.push_back(bdd_inv);
		    }
		}
	    }

	    DEBUG_MSG(dumpMutexBDDs(fw););
	    merge(vars, notMutexBDDs, mergeAndBDD,
		  p.max_mutex_time,  p.max_mutex_size);
	    std::reverse(notMutexBDDs.begin(), notMutexBDDs.end());
	    DEBUG_MSG(cout << "Mutex initialized "<< (fw ? "fw" : "bw") << ". Total mutex added: " << num_mutex  << " Invariant groups: " << num_invariants  << endl;);
	    dumpMutexBDDs(fw);
	}

	//gst_mutex.check_mutexes(*this);
    }

    void SymStateSpaceManager::addDeadEndStates(bool fw, BDD bdd) {
	//There are several options here, we could follow with edeletion
	//and modify the TRs, so that the new spurious states are never
	//generated. However, the TRs are already merged and the may get
	//too large. Therefore we just keep this states in another vectors
	//and spurious states are always removed. TODO: this could be
	//improved.
	if(fw || isAbstracted()) {
	    if (isAbstracted()) bdd = shrinkForall(bdd);
	    notDeadEndFw.push_back(!bdd);
	    mergeBucketAnd(notDeadEndFw);
	}else{
	    notDeadEndBw.push_back(!bdd);
	    mergeBucketAnd(notDeadEndBw);
	}
    }


    void SymStateSpaceManager::addDeadEndStates(const std::vector<BDD> & fw_dead_ends,
						const std::vector<BDD> & bw_dead_ends) {
	for (BDD bdd : fw_dead_ends){
	    bdd = shrinkForall(bdd);
	    if(!bdd.IsZero()) {
		notDeadEndFw.push_back(!bdd);
	    }
	}

	for (BDD bdd : bw_dead_ends){
	    bdd = shrinkForall(bdd);
	    if(!bdd.IsZero()) {
		notDeadEndFw.push_back(!bdd);
	    }
	}
	mergeBucketAnd(notDeadEndFw);
    }


    void SymStateSpaceManager::dumpMutexBDDs(bool fw) const {
	if(fw){
	    cout << "Mutex BDD FW Size(" <<p.max_mutex_size << "):";
	    for(const auto & bdd : notMutexBDDsFw){
		cout << " " << bdd.nodeCount();
	    }
	    cout << endl;
	}else{

	    cout << "Mutex BDD BW Size(" <<p.max_mutex_size << "):";
	    for(const auto & bdd : notMutexBDDsBw){
		cout << " " << bdd.nodeCount();
	    }
	    cout << endl;
	}

    }


    // void AbstractStateSpace::init_transitions(){
    // 	if(!transitions.empty()) return; //Already initialized!
    // 	init_individual_trs(); 
    // 	init_transitions_from_individual_trs();

    // 	if(parentMgr){
    // 	    DEBUG_MSG(cout << "Init transitions from parent state space" << endl;);
    // 	    SymStateSpaceManager * mgr_parent = parentMgr;
    // 	    while(mgr_parent && mgr_parent->transitions.empty()){
    // 		mgr_parent = mgr_parent->parentMgr;
    // 	    }
    // 	    const auto & trsParent = mgr_parent->getTransitions();
    
    // 	    while(mgr_parent && mgr_parent->indTRs.empty()){
    // 		mgr_parent = mgr_parent->parentMgr;
    // 	    }
    // 	    const auto & indTRsParent = mgr_parent->getIndividualTransitions();

    // 	    shrinkTransitions(trsParent, indTRsParent, transitions,
    // 			      p.max_tr_time, p.max_tr_size);
    // 	    return;
    // 	}else {
    // 	    init_individual_trs(); 
    // 	    init_transitions_from_individual_trs();
    // 	}
    // }

    void SymStateSpaceManager::zero_preimage(const BDD & bdd, vector <BDD> & res, int nodeLimit) const{
	for(const auto & tr : transitions.at(0)){   
	    res.push_back(tr.preimage(bdd, nodeLimit));
	}
    }

    void SymStateSpaceManager::zero_image(const BDD & bdd, vector <BDD> & res, int nodeLimit) const{
	for(const auto & tr : transitions.at(0)){   
	    res.push_back(tr.image(bdd, nodeLimit));
	}
    }
 
    void SymStateSpaceManager::cost_preimage(const BDD & bdd, map<int, vector<BDD> > &res,
					     int nodeLimit) const {
	for(auto trs : transitions){
	    int cost = trs.first;
	    if(cost == 0) continue;
	    for(size_t i = res[cost].size(); i < trs.second.size(); i++){
		BDD result = trs.second[i].preimage(bdd, nodeLimit);
		res[cost].push_back(result);
	    }
	}
    }

    void SymStateSpaceManager::cost_image(const BDD & bdd,
					  map<int, vector<BDD> > &res,
					  int nodeLimit) const{
	for(auto trs : transitions){
	    int cost = trs.first;
	    if(cost == 0) continue;
	    for(size_t i = res[cost].size(); i < trs.second.size(); i++){
		//cout << "Img: " << trs.second[i].nodeCount() << " with bdd " << bdd.nodeCount() << " node limit: " << nodeLimit << endl;
		BDD result = trs.second[i].image(bdd, nodeLimit);
		//cout << "Res: " << result.nodeCount() << endl;
		res[cost].push_back(result);
	    }
	}
    }

    BDD SymStateSpaceManager::filter_mutex(const BDD & bdd, bool fw,
					   int nodeLimit, bool initialization) {  
	BDD res = bdd;
	const vector<BDD> & notDeadEndBDDs = ((fw || isAbstracted()) ? notDeadEndFw : notDeadEndBw);
	for(const BDD & notDeadEnd : notDeadEndBDDs){
	    DEBUG_MSG(cout << "Filter: " << res.nodeCount()  << " and dead end " <<  notDeadEnd.nodeCount()  << flush;);
	    res = res.And(notDeadEnd, nodeLimit);
	    DEBUG_MSG(cout << ": " << res.nodeCount() << endl;);
	}

	const vector<BDD> & notMutexBDDs = (fw ? notMutexBDDsFw : notMutexBDDsBw);


	switch (p.mutex_type){
	case MutexType::MUTEX_NOT:
	    break;
	case MutexType::MUTEX_EDELETION:
	    if(initialization){
		for(const BDD & notMutexBDD : notMutexBDDs){
		    DEBUG_MSG(cout << res.nodeCount()  << " and " <<  notMutexBDD.nodeCount()  << flush;);
		    res = res.And(notMutexBDD, nodeLimit);
		    DEBUG_MSG(cout << ": " << res.nodeCount() << endl;);
		}
	    }
	    break;
	case MutexType::MUTEX_AND:
	    for(const BDD & notMutexBDD : notMutexBDDs){
		DEBUG_MSG(cout << "Filter: " << res.nodeCount()  << " and " <<  notMutexBDD.nodeCount()  << flush;);
		res = res.And(notMutexBDD, nodeLimit);
		DEBUG_MSG(cout << ": " << res.nodeCount() << endl;);
	    }
	    break;
	case MutexType::MUTEX_RESTRICT:
	    for(const BDD & notMutexBDD : notMutexBDDs)
		res = res.Restrict(notMutexBDD);
	    break;
	case MutexType::MUTEX_NPAND:
	    for(const BDD & notMutexBDD : notMutexBDDs)
		res = res.NPAnd(notMutexBDD);
	    break;
	case MutexType::MUTEX_CONSTRAIN:
	    for(const BDD & notMutexBDD : notMutexBDDs)
		res = res.Constrain(notMutexBDD);
	    break;
	case MutexType::MUTEX_LICOMP:
	    for(const BDD & notMutexBDD : notMutexBDDs)
		res = res.LICompaction(notMutexBDD);
	    break;
	}
	return res;
    }

    int SymStateSpaceManager::filterMutexBucket(vector<BDD> & bucket, bool fw,
						bool initialization, int maxTime, int maxNodes){
	int numFiltered = 0;
	setTimeLimit(maxTime);
	try{
	    for (size_t i = 0; i < bucket.size(); ++i){
		DEBUG_MSG(cout <<  "Filter spurious " << (fw ? "fw" : "bw") << ": " << *this
			  << " from: "  << bucket[i].nodeCount() <<
			  " maxTime: " << maxTime << " and maxNodes: " << maxNodes;);
      
		bucket[i] = filter_mutex(bucket[i], fw, maxNodes, initialization);
		DEBUG_MSG(cout << " => " << bucket[i].nodeCount() << endl;);
		numFiltered ++;
	    }  
	}catch(BDDError e){
	    DEBUG_MSG(cout << " truncated." << endl;);
	}
	unsetTimeLimit();

	return numFiltered;
    }

    void SymStateSpaceManager::filterMutex (Bucket & bucket, bool fw, bool initialization) {
	filterMutexBucket(bucket, fw, initialization, 
			  p.max_aux_time, p.max_aux_nodes);
    }

    void SymStateSpaceManager::mergeBucket(Bucket & bucket) const {
	mergeBucket(bucket, p.max_aux_time, p.max_aux_nodes);
    }

    void SymStateSpaceManager::mergeBucketAnd(Bucket & bucket) const {
	mergeBucketAnd(bucket, p.max_aux_time, p.max_aux_nodes);
    }

    void SymStateSpaceManager::shrinkBucket(Bucket & bucket, int maxNodes) {
	for(size_t i = 0; i < bucket.size(); ++i){
	    bucket[i] = shrinkExists(bucket[i], maxNodes);
	}
    }


    void SymStateSpaceManager::init_mutex(const std::vector<MutexGroup> & mutex_groups){
	//If (a) is initialized OR not using mutex OR edeletion does not need mutex
	if(mutexInitialized || p.mutex_type == MutexType::MUTEX_NOT)
	    return; //Skip mutex initialization
 
	if(p.mutex_type == MutexType::MUTEX_EDELETION){
	    SymStateSpaceManager::init_mutex(mutex_groups, true, true);
	}else{
	    SymStateSpaceManager::init_mutex(mutex_groups, true, false);
	}
    }


    void SymStateSpaceManager::init_transitions(){
	if(!transitions.empty()) return; //Already initialized!
	init_individual_trs(); 
	init_transitions_from_individual_trs();
    }





    // void AbstractStateSpace::init_mutex(const std::vector<MutexGroup> & mutex_groups){
    // 	//If (a) is initialized OR not using mutex OR edeletion does not need mutex
    // 	if(mutexInitialized || p.mutex_type == MutexType::MUTEX_NOT)
    // 	    return; //Skip mutex initialization
 
    // 	if(parentMgr){
    // 	    setTimeLimit(p.max_mutex_time);
    // 	    DEBUG_MSG(cout << "Init mutex from parent" << endl;);
    // 	    mutexInitialized = true;
    // 	    //Initialize mutexes from other manager
    // 	    try{
    // 		for(auto & bdd : parentMgr->notMutexBDDsFw){
    // 		    BDD shrinked = shrinkExists(bdd, p.max_mutex_size);
    // 		    notMutexBDDsFw.push_back(shrinked);
    // 		}
    // 		for(auto & bdd : parentMgr->notMutexBDDsBw){
    // 		    BDD shrinked = shrinkExists(bdd, p.max_mutex_size);
    // 		    notMutexBDDsBw.push_back(shrinked);
    // 		}
    // 		unsetTimeLimit();
    // 	    }catch(BDDError e){ 
    // 		unsetTimeLimit();
    // 		//Forget about it
    // 		vector<BDD>().swap(notMutexBDDsFw);
    // 		vector<BDD>().swap(notMutexBDDsBw);
    // 		init_mutex(mutex_groups, true, false);
    // 	    }
    // 	    //We will compute mutex by fluent on demand 
    // 	} else if(p.mutex_type == MutexType::MUTEX_EDELETION){
    // 	    init_mutex(mutex_groups, true, true);
    // 	}else{
    // 	    init_mutex(mutex_groups, true, false);
    // 	}
    // }

    SymParamsMgr::SymParamsMgr(const options::Options & opts) : 
	max_tr_size(opts.get<int>("max_tr_size")),
	max_tr_time(opts.get<int>("max_tr_time")),
	mutex_type(MutexType(opts.get_enum("mutex_type"))),
	max_mutex_size(opts.get<int>("max_mutex_size")),
	max_mutex_time(opts.get<int>("max_mutex_time")),
	max_aux_nodes(opts.get<int>("max_aux_nodes")),
	max_aux_time (opts.get<int>("max_aux_time")) {

	//Don't use edeletion with conditional effects
	if(mutex_type == MutexType::MUTEX_EDELETION && has_conditional_effects()){
	    cout << "Mutex type changed to mutex_and because the domain has conditional effects" << endl;
	    mutex_type = MutexType::MUTEX_AND;
	}
    }

    SymParamsMgr::SymParamsMgr() : 
	max_tr_size(100000),
	max_tr_time(60000),
	mutex_type(MutexType::MUTEX_EDELETION),
	max_mutex_size(100000),
	max_mutex_time(60000), 
	max_aux_nodes(1000000), max_aux_time(2000) {
	//Don't use edeletion with conditional effects
	if(mutex_type == MutexType::MUTEX_EDELETION && has_conditional_effects()){
	    cout << "Mutex type changed to mutex_and because the domain has conditional effects" << endl;
	    mutex_type = MutexType::MUTEX_AND;
	}
    }

    void SymParamsMgr::print_options() const{
	cout << "TR(time=" << max_tr_time << ", nodes=" << max_tr_size << ")" << endl;
	cout << "Mutex(time=" << max_mutex_time << ", nodes=" << max_mutex_size << ", type=" << mutex_type << ")" << endl;
	cout << "Aux(time=" << max_aux_time << ", nodes=" << max_aux_nodes << ")" << endl;
    }

    void SymParamsMgr::add_options_to_parser(options::OptionParser &parser){
  
	parser.add_option<int> ("max_tr_size", "maximum size of TR BDDs", "100000");
  
	parser.add_option<int> ("max_tr_time",
				"maximum time (ms) to generate TR BDDs", "60000");
     
	parser.add_enum_option("mutex_type", MutexTypeValues,
			       "mutex type", "MUTEX_EDELETION");

	parser.add_option<int> ("max_mutex_size",
				"maximum size of mutex BDDs", "100000");

	parser.add_option<int> ("max_mutex_time",
				"maximum time (ms) to generate mutex BDDs", "60000");

	parser.add_option<int> ("max_aux_nodes", "maximum size in pop operations", "1000000");
	parser.add_option<int> ("max_aux_time", "maximum time (ms) in pop operations", "2000");

    }

    std::ostream & operator<<(std::ostream &os, const SymStateSpaceManager & abs){
	abs.print(os, false);
	return os;
    }
}
