#include "uniform_cost_search.h"

#include "closed_list.h"
#include "sym_solution.h"
#include "sym_util.h" 
#include "sym_engine.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>

#include "../utils/debug_macros.h"
#include "../utils/timer.h"
#include "../globals.h"

#include "sym_test.h"

using namespace std;
using utils::g_timer;
using utils::Timer;

namespace symbolic {

UniformCostSearch::UniformCostSearch(SymController * eng,
		   const SymParamsSearch & params): 
    UnidirectionalSearch(params), 
    parent(nullptr), g(0), closed(new ClosedList()),  
    estimationCost(params), estimationZero(params),
    //estimationDisjCost(params), estimationDisjZero(params),
    lastStepCost(true), engine(eng) {}

bool UniformCostSearch::init(SymStateSpaceManager * manager, bool forward){
    mgr = manager;
    fw = forward;
    lastStepCost = true;
    g = 0;
    //Ensure that the mgr of the original state space is initialized (only to get the planner output cleaner) 
    mgr->init();
  
    DEBUG_MSG(cout << "Init exploration: " << dirname(forward) << *this/* << " with mgr: " << manager */<< endl;);

    if(fw){
	BDD init = mgr->getInitialState();
	Sfilter.push_back(init);
    }else{
	Sfilter.push_back(mgr->getGoal());
    }
    DEBUG_MSG (cout << dirname(forward) << " initialized to "; Sfilter[0].print(0,1); cout << "Init closed" << endl;);

    closed->init(this, mgr);
 
    closed->insert(0, Sfilter[0]);
    closed->setHNotClosed(open_list.minNextG(g));
    closed->setFNotClosed(f);

    prepareBucket();

    if(!isAbstracted()) 
	engine->setLowerBound(f);

    return true;
}

void UniformCostSearch::getNotExpanded(Bucket & res) const {
    copyBucket(Sfilter, res);
    copyBucket(Smerge, res);
    copyBucket(Szero, res);
    copyBucket(S, res);
        
    open_list.get_reopen_states(res);
}

void UniformCostSearch::relaxClosed(){
    open_list.init(this, mgr);
    closed->init(this, mgr);  
    
    assert (expansionReady());
    assert (!S.empty());

    for(const BDD & bdd : S){
	closed->insert(g, bdd);
    }
    
    closed->setHNotClosed(open_list.minNextG(g));
    DEBUG_MSG(cout << "HNotClosed initialized to: " << closed->getHNotClosed() << " g=" << g<< endl;);
}

void UniformCostSearch::closeStates(Bucket & bucket, int g_val) {
    if(!isAbstracted() && lastStepCost && g_val == 0){
	return; //Avoid closing init twice
    } 
    DEBUG_MSG (cout <<"Insert g="<< g_val << " states into closed: " << nodeCount(bucket) << " (" << bucket.size() << " bdds)" << endl;);

    for(const BDD & states : bucket){
	DEBUG_MSG (cout <<"Inserting: " << states.nodeCount() << " " << *closed << endl;);

	closed->insert(g_val, states);
    }

    // Check cut with other frontier in the abstract search to set fMainDiagonal
    if(isAbstracted() &&  perfectHeuristic && bdExp->getFMainDiagonal() == -1){ 
	for (const BDD & bdd : bucket){
	    BDD cut = perfectHeuristic->getClosed()*bdd;
	    if(!cut.IsZero()){
		bdExp->setFMainDiagonal(f);
	    }
	}    
    }
}


void UniformCostSearch::checkCutOriginal(Bucket & bucket, int g_val){
    DEBUG_MSG (cout <<"CHECK CUT: " << g_val << ", bucket: " << nodeCount(bucket) << endl;);
    //If it is the original space, maybe we have found a solution, set upper bound  
    if(!perfectHeuristic || p.get_non_stop()){
	return;
    }

    for (BDD & bucketBDD : bucket){
	SymSolution sol = perfectHeuristic->checkCut(bucketBDD, g_val, fw);
	if (sol.solved()){
	    cout << "Solution found with cost " << sol.getCost() << 
		" total time: " << g_timer <<  endl;
	    // Solution found :)
	    engine->new_solution(sol);
	}
	bucketBDD *= perfectHeuristic->notClosed();   //Prune everything closed in opposite direction
    }
}



/* pop() must recompute the current g and f value, 
   extract states with those values from open,
   remove duplicates, classify h and put them on Sfilter. */
void UniformCostSearch::pop(){
    DEBUG_MSG(cout << "POP: bucketReady: " << bucketReady() << endl;);
    if(bucketReady() || expansionReady() || engine->solved()){
	return;
    }

    int prevF = f;
    while(f < numeric_limits<int>::max() && !bucketReady()){ 
	pair<int, int> nextFG = open_list.pop(f, g, perfectHeuristic);
	f = nextFG.first;
	g = nextFG.second;
	if(perfectHeuristic){
	    perfectHeuristic->cleanEvals(isOriginal());
	}
	open_list.extract_states(f, g, Sfilter);
    }
	
    if (open_list.minG () <= open_list.minNextG(g)){
	open_list.closeMinOpen();
    }
    closed->setHNotClosed(open_list.minNextG(g));

    if(f > prevF){
	//States closed, update the f value if is the original state space
	if(isOriginal()) 
	    engine->setLowerBound(f);
	
	closed->setFNotClosed(f);
    }
    
    if(f == numeric_limits<int>::max()) return;

    //Try to prepare next Bucket
    DEBUG_MSG(cout << "I AM COMPUTING ESTIMATION BECAUSE OF POP" << endl;);
    computeEstimation(true);
    DEBUG_MSG(cout << "Estimation computed because of pop" << endl;);
}


bool UniformCostSearch::prepareBucket(){
    int maxTime = p.getAllotedTime(nextStepTime());
    int maxNodes = p.getAllotedNodes(nextStepNodesResult());

    if(!bucketReady()){
	pop();
	DEBUG_MSG(cout << "Pop done for bucket prepared" << endl;);
    }

    if(engine->solved()){
	DEBUG_MSG(cout << "SOLVED!!!: "<< engine->getLowerBound() << " >= " << engine->getUpperBound() << endl;);
	return true; //If it has been solved, return 
    }

    DEBUG_MSG(cout << "Preparing bucket time(" << maxTime << ") nodes(" << maxNodes << ")" << endl;);
    if(initialization()){
	mgr->init_mutex(g_mutex_groups);
    }
    Timer filterTime;
    if(!Sfilter.empty()){
	//First, if possible, attempt to merge the g-Sopen (only uses
	//pop_time). This is only to reuse the most resources possible.
	//mergeBucket(Sfilter, p.max_pop_time, p.max_pop_nodes); it has been merged in pop
	int numFiltered = mgr->filterMutexBucket(Sfilter, fw, initialization(), maxTime, maxNodes);
	if(numFiltered > 0){
	    Smerge.insert(std::end(Smerge), std::begin(Sfilter), std::begin(Sfilter) + numFiltered); 
	}
	if(numFiltered == (int)(Sfilter.size())){
	    Bucket().swap(Sfilter);
	}else{
	    Sfilter.erase(std::begin(Sfilter), std::begin(Sfilter) + numFiltered);
	    violated(TruncatedReason::FILTER_MUTEX, filterTime(), maxTime, maxNodes);
	    return false;
	}
    }
  
    if(!Smerge.empty()){  
	if(Smerge.size() > 1){
	    int remainingTime = maxTime - 1000*filterTime();
	    if(remainingTime < 0 || !mgr->mergeBucket(Smerge, remainingTime, maxNodes)){
		violated(TruncatedReason::MERGE_BUCKET, filterTime(), maxTime, maxNodes);
		return false;
	    }
	}

	DEBUG_MSG(cout << " BUCKET MERGED: " << Smerge.size() << " " << nodeCount(Smerge) << endl;);

	//Successfully merged
	// a) close Smerge

	//b) put result on Szero or S (or both)
	if(mgr->hasTransitions0()){
	    S.insert(std::end(S), std::begin(Smerge), std::end(Smerge));
	    assert(Szero.empty());
	    Szero.swap(Smerge); 
	}else{
	    S.swap(Smerge);      
	}
    }

    //If there are zero cost operators, merge S 
    if(mgr->hasTransitions0() && Szero.empty()){
	if(S.size() > 1){
	    int remainingTime = maxTime - 1000*filterTime();
	    if(remainingTime < 0 || !mgr->mergeBucket(S, remainingTime, maxNodes)){
		violated(TruncatedReason::MERGE_BUCKET_COST, filterTime(), maxTime, maxNodes);
		return false;
	    }
	}
    }

    // if(lastStepCost){
    //   DEBUG_MSG(cout << "EST_DISJ_COST " << filterTime() << "    ";);
    //   estimationDisjCost.stepTaken(1000*filterTime(), max(nodeCount(S), nodeCount(Szero)));
    // }else{
    //   DEBUG_MSG(cout << "EST_DISJ_ZERO " << filterTime() << "    ";);
    //   estimationDisjZero.stepTaken(1000*filterTime(), max(nodeCount(S), nodeCount(Szero)));
    // }
    DEBUG_MSG(cout << " BUCKET PREPARED: " << S.size() << " " << nodeCount(S) << endl;);
    return true;
}


bool UniformCostSearch::expand_zero(int maxTime, int maxNodes){
    DEBUG_MSG(cout<<"expand_zero"<< endl;);
    //Image with respect to 0-cost actions
    assert(expansionReady() && nodeCount(Szero) <= maxNodes);
    // DEBUG_MSG(cout << "    >> Truncated. 0-Frontier size exceeded: " << nodeCount(Szero) << endl;);
    // estimationZero.violated_nodes(nodeCount(Szero));
    // return false;
  
    Timer image_time;

    int nodesStep = nodeCount(Szero);
    double statesStep = mgr->stateCount(Szero);
    mgr->setTimeLimit(maxTime);
    //Compute image, storing the result on Simg
    try{
	int numImagesComputed = 0;
	for(size_t i = 0; i < Szero.size(); i++){	
	    Simg.push_back(map<int, Bucket>());
	    mgr->zero_image(fw, Szero[i], Simg[i][0], maxNodes);
	    DEBUG_MSG(for(auto bdd : Simg[i][0]){ cout << "RESULT OF ZERO_IMG: "; bdd.print(0, 1);});
	    ++ numImagesComputed;
	}
	mgr->unsetTimeLimit();    
    }catch(BDDError e){
	mgr->unsetTimeLimit();
	violated(TruncatedReason::IMAGE_ZERO, image_time(), maxTime, maxNodes);
	// Szero.erase(begin(Szero), begin(Szero) + numImagesComputed);

	stats.add_image_time_failed(image_time());
	return false;
    }
    lastStepCost = false; //Must be set to false before calling checkCut

    DEBUG_MSG(cout << "EXPAND ZERO HAS PUT IN Simg: " << Simg.size() << endl;);
    Bucket().swap(Szero); //Delete Szero because it has been expanded
    long nodesRes = 0;
    //Process Simg, removing duplicates and computing h. Store in Sfilter and reopen.
    for(auto & resimg : Simg){
	if(!isAbstracted()) { // Apply nipping in the original state space
	    checkCutOriginal(resimg[0], g);
	}
	for(auto bdd : resimg[0]){
	    nodesRes = max<long>(nodesRes, bdd.nodeCount());
#ifdef DEBUG_GST
	    gst_plan.checkOpen(bdd, g, this);
#endif
	}
    
	if(perfectHeuristic) 
	    perfectHeuristic->cleanEvals(!isAbstracted());

	open_list.extract_and_insert(resimg[0], f, g, Sfilter);
    }

    vector<map<int, Bucket>>().swap(Simg);
    DEBUG_MSG(cout << "EST_ZERO " << image_time() << "   ";);
    estimationZero.stepTaken(1000*image_time(), max<long>(nodesRes, nodesStep));

    //Try to prepare next Bucket
    computeEstimation(true);

    DEBUG_MSG(cout << "  >> 0-cost step " <<  (fw ? " fw " : " bw ") << ": " << nodesStep << " nodes "
	      << statesStep << " states " << " done in " << image_time << endl;);

    stats.add_image_time(image_time());
    DEBUG_MSG( cout << "Prepared in Sfilter: " << Sfilter.size()  << ": " << nodeCount(Sfilter) << endl;
	       cout << "Prepared in Szero: " << Szero.size()  << ": " << nodeCount(Szero) << endl;
	       cout << "Prepared in S: " << S.size()  << ": " << nodeCount(S) << endl;);

    if (Szero.empty() && !S.empty()){
	// We've finished with Szero, now let's re-evaluate S so that
	// evalOrig is set

	perfectHeuristic->cleanEvals(!isAbstracted());

	for (BDD & bdd : S){ 
	    BDD pruned = perfectHeuristic->evaluate(bdd, f, f-g, this, true);
	    if (!pruned.IsZero()) {
		bdd -= pruned;
		open_list.insert_reopen(pruned, g);	      
	    }	      
	}

	prepareBucket();
    }
    return true;
}


bool UniformCostSearch::expand_cost(int maxTime, int maxNodes){
    DEBUG_MSG(cout<<"expand_cost: " << nodeCount(S) << endl;);
    assert(expansionReady());
    assert(nodeCount(S) <= maxNodes);
    int nodesStep = nodeCount(S);
    double statesStep = mgr->stateCount(S);
    Timer image_time;
    DEBUG_MSG(cout << "Setting maxTime: " << maxTime << endl;);
    mgr->setTimeLimit(maxTime);
    try{
	for(size_t i = 0; i < S.size(); i++){
	    Simg.push_back(map<int, Bucket>());
	    mgr->cost_image(fw, S[i], Simg[i], maxNodes);
	}
	mgr->unsetTimeLimit();
    }catch(BDDError e){
	//Update estimation
	mgr->unsetTimeLimit();
	violated(TruncatedReason::IMAGE_COST, image_time(), maxTime, maxNodes);
	stats.add_image_time_failed(image_time());
	return false;
    }

    //Include new states in the open list 
    int stepNodes = nodesStep;
    for(auto & resImage : Simg){
	for(auto & pairCostBDDs : resImage){
	    int cost = g + pairCostBDDs.first; //Include states of the given cost 
	    mgr->mergeBucket(pairCostBDDs.second);

	    //Check the cut (removing states classified, since they do not need to be included in open)
	    if (!isAbstracted()){
		checkCutOriginal(pairCostBDDs.second, cost); 
	    }

	    for(auto & bdd : pairCostBDDs.second){  
		if(!bdd.IsZero()){
		    //TODO: maybe we can also use the heuristics to prune states
		    //right here. Also, we could prune duplicates. Not sure if it
		    //is worth it.
		    int fVal = cost;
		    if(perfectHeuristic){
			fVal += perfectHeuristic->getHNotClosed();
		    }
		    if (fVal < engine->getUpperBound()){
	
			stepNodes = max(stepNodes, bdd.nodeCount());
#ifdef DEBUG_GST
			gst_plan.checkOpen(bdd, cost, this);
#endif
			open_list.insert(bdd, cost);
		    }
		}
	    }
	}
    }
    vector<map<int, Bucket>>().swap(Simg);
    Bucket().swap(S);
    DEBUG_MSG(cout << "EST_COST " << image_time() << "   " << endl;);
    estimationCost.stepTaken(1000*image_time(), stepNodes);
    lastStepCost = true;

    stats.add_image_time(image_time());

    DEBUG_MSG(cout << "  >> cost step" << (fw ? " fw " : " bw ") << ": " << nodesStep << " nodes " << statesStep << " states " << " done in " << image_time << " total time: " << g_timer << endl;);
    return true; 
}


bool UniformCostSearch::stepImage(int maxTime, int maxNodes){
    if(mgr->getAbstraction())
	cout << ">> Step: " << *(mgr->getAbstraction());
    else
	cout << ">> Step: original";
    cout << (fw ? " fw " : " bw ") << "f=" << f << ", g=" << g;
    cout << " frontierNodes: " << frontierNodes() << " [" << frontierBuckets() << "]"  << " total time: " << g_timer 
	 << " total nodes: " << mgr->totalNodes() << " total memory: " << mgr->totalMemory() << endl;

#ifdef DEBUG_GST
    gst_plan.checkUcs(this );
#endif

    DEBUG_MSG(cout << "Step " << dirname(fw) << " f: " << f
	      << " g: " << g << endl;);
    Timer sTime;
    DEBUG_MSG(cout << "preparing bucket.." << " total time: " << g_timer  << endl;);
    if(!expansionReady() && !prepareBucket(/*maxTime, maxNodes, false*/)){
	cout << "    >> Truncated while preparing bucket" << endl;
	if(sTime()*1000.0 > p.maxStepTime){
	    double ratio = (double)p.maxStepTime/((double)sTime()*1000.0);
	    p.maxStepNodes *= ratio;
	    DEBUG_MSG(cout << "MAX STEP NODES CHANGED TO: " << p.maxStepNodes <<
		      " after truncating with " << sTime() << " seconds" << endl;);
	}
	stats.step_time += sTime();
	return false;    
    } 
    DEBUG_MSG(cout << "... bucket prepared. " << endl;);
    if(engine->solved()) return true; //Skip image if we are done

    bool ok = true;
    mgr->init_transitions(); // Ensure that transitions have been initialized
    int stepNodes = frontierNodes();
    if(!Szero.empty()){
	ok = expand_zero(maxTime, maxNodes);
    }else if(!S.empty()){
	//Image with respect to cost actions
	ok = expand_cost(maxTime, maxNodes);
    }else{
	DEBUG_MSG(cout << "   >> All pop states were filtered: " << open_list.empty() << endl;);
    }

    pop(); //We prepare the next bucket before checking time in doing
    //the step because we consider preparing the bucket as a part
    //of the step.

    if(sTime()*1000.0 > p.maxStepTime){
	double ratio = (double)p.maxStepTime/((double)sTime()*1000.0);
	p.maxStepNodes = stepNodes*ratio;
	DEBUG_MSG(cout << "MAX STEP NODES CHANGED TO: " << p.maxStepNodes << " after taking " << sTime() << " seconds" << endl;);
    }else if (!ok){
	//In case maxAllotedNodes were exceeded we reduce the maximum
	//frontier size by 3/4.  TODO: make this a parameter
	p.maxStepNodes = stepNodes*0.75; 
    }

    DEBUG_MSG(cout << "END STEP: has eval orig: " << perfectHeuristic->hasEvalOrig() << endl;);
    assert(isAbstracted() || perfectHeuristic->hasEvalOrig());
    stats.step_time += sTime();
    return ok;
}



/////////////////////////////////////////////////
///  Methods related to heuristic evaluation  ///
/////////////////////////////////////////////////

void UniformCostSearch::notifyPrunedBy(int fVal, int gVal){
    DEBUG_MSG(cout << "NOTIFIED pruned by: " << fVal << ", " << 
	      gVal << " on " << *this << endl;);
    assert(gVal == g);
    
    acceptedValues[fVal].insert(gVal);
}

void UniformCostSearch::notify(const Bucket & bucket, int fNotClosed /*=0*/){
    DEBUG_MSG(cout << " NOTIFIED BUCKET to: " << *this << endl; );

    //Now, try to reduce current frontier
    Bucket bucketPruned;
    if(extract_states(Sfilter, bucket, bucketPruned) |
       extract_states(Smerge, bucket, bucketPruned) |
       extract_states(Szero, bucket, bucketPruned) |
       extract_states(S, bucket, bucketPruned)){
	if(fNotClosed){
	    DEBUG_MSG(cout << "NOTIFIED in notify pruned by: " << fNotClosed << ", " << g << " on " << *this << endl;);
	    acceptedValues[fNotClosed].insert(g);
	}
	if (!bucketPruned.empty()){
	    open_list.insert_reopen(bucketPruned, g);
	}
    }
  
    if(!expansionReady() && !bucketReady()){
	pop(); // Only will do something if !bucketReady 
    }else{
	computeEstimation(false);
    }
    DEBUG_MSG(cout << "END NOTIFY" << endl;
	      cout << "AFTER: "; printFrontier(); cout << endl;);
}

//In the original state space search, there cannot be any state in the
//frontier closed (because of checkCut). So fNotClosed and hNotClosed
//are the minimum f and h for sure
void UniformCostSearch::notifyNotClosed(int fValue, int hValue){
    assert(!isAbstracted());
    DEBUG_MSG(cout << "NOTIFIED NOT CLOSED to: " << *this << ": f="<< fValue << " h="<< hValue << endl;);
    if(fValue > f || hValue > f-g){
    	// We should increase our f value
	if(open_list.first_position(g)) {
	    DEBUG_MSG(cout << "We were in the first position " << endl;);

	    int newFValue = max(fValue, g + hValue);
	    if (newFValue <= f) return;
	    f = newFValue;
	
	    if(hValuesExplicit.empty() && !perfectHeuristic->hasChildren()) {
		DEBUG_MSG(cout << "SetOrigF: " << f << " hVal: " << (f-g) << endl;);
		//We are starting the diagonal => the g should be maintained
		if(isOriginal()) perfectHeuristic->setOrigF(f, f-g);
		else perfectHeuristic->setOppositeF(f, f-g);
		closed->setFNotClosed(f);	
	    } else {

		DEBUG_MSG(cout << "Need to reevaluate states in open[g]. f =" 
			  << f << " g=" << g << endl;);

		Bucket bucketPruned;
		moveBucket(Sfilter, bucketPruned); 
		moveBucket(Smerge, bucketPruned);
		moveBucket(Szero, bucketPruned);
		moveBucket(S, bucketPruned);
		
		if (!bucketPruned.empty()){
		    open_list.insert_reopen(bucketPruned, g);
		}
		vector<map<int, Bucket>>().swap(Simg);

		perfectHeuristic->cleanEvals(!isAbstracted());

		//Extract from open
		open_list.extract_states(f, g, Sfilter);

		prepareBucket();
	    }
	} else {
	    DEBUG_MSG(cout << "Re-pop " << (fw ? " fw " : " bw ") << endl;);
	    //DEBUG_MSG(cout << "VENGA "  << reopen.at(g).size() << endl; for (auto bdd : reopen[g]) cout << bdd.nodeCount() << endl;);

	    Bucket bucketPruned;
		moveBucket(Sfilter, bucketPruned); 
		moveBucket(Smerge, bucketPruned);
		moveBucket(Szero, bucketPruned);
		moveBucket(S, bucketPruned);
		
		if (!bucketPruned.empty()){
		    open_list.insert_reopen(bucketPruned, g);
		}

	    vector<map<int, Bucket>>().swap(Simg);
	    //f -= 1;
	    //g = numeric_limits<int>::max();
	    pop();
	    // if(!bucketReady() && !expansionReady()){
	    // 	pop(); // Only will do something if !bucketReady 
	    // }else{
	    // 	computeEstimation(false);
	    // }
	}
    } else {
	//cout << "Don't need to do anything" << endl;
    }
}


void UniformCostSearch::addHeuristic(std::shared_ptr<SymHeuristic> newHeuristic){
    mgr->addDeadEndStates(fw, newHeuristic->getDeadEnds());
    if(!perfectHeuristic || 
       newHeuristic->getMaxValue() > perfectHeuristic->getHNotClosed()){
	const auto & hVals = newHeuristic->getHValues();
	hValuesExplicit.insert(begin(hVals), end(hVals));
        heuristics.push_back(newHeuristic);
        Bucket toNotify;
	toNotify.push_back(newHeuristic->prunedStates(f-g));
    
	notify(toNotify);
    }
}

BDD UniformCostSearch::compute_heuristic(const BDD & from, int fVal, int hVal, bool store_eval){
    //assert(isAbstracted() || !perfectHeuristic || hVal <= perfectHeuristic->getHNotClosed());
    assert(hVal >= 0);
    assert(fVal >= 0);
    assert(fVal >= hVal);
    DEBUG_MSG(cout << "Compute heuristic f=" << fVal << " h=" << hVal << endl;);
    // Compute static heuristics
    BDD pruned = mgr->zeroBDD();
    BDD notPruned = from;
    for (auto & heur : heuristics){
	BDD newPruned = from * heur->prunedStates(hVal);
	notPruned -= newPruned; 
	pruned += newPruned;
    }

    // Compute collision with the other perimeter and related abstractions 
    if(perfectHeuristic){
	pruned += perfectHeuristic->evaluate(notPruned, fVal, hVal, this, store_eval);
    }
    DEBUG_MSG(if (!pruned.IsZero()) cout << "Pruned: " << pruned.nodeCount() << endl;);
    return pruned;
}

BDD UniformCostSearch::getExpanded() const{
    BDD res = closed->getClosed();
    for(const BDD & bdd : Sfilter){
	res *= !bdd;
    }
    for(const BDD & bdd : Smerge){
	res *= !bdd;
    }
    for(const BDD & bdd : S){
	res *= !bdd;
    }
    return res;
}

/////////////////////////////////////////////////
///    Methods to decide useful/searchable    ///
/////////////////////////////////////////////////

bool UniformCostSearch::isUseful(const vector<BDD> & evalStates,
			      const std::vector<BDD> & newFrontier, 
			      double ratio) const {
    vector<BDD> bucketwoFrontier; 
    for(const BDD & bdd : evalStates){
	BDD bddwoFrontier  = bdd;
	for(const BDD & fBDD : newFrontier) {
	    bddwoFrontier -= fBDD;
	}
	if(!bddwoFrontier.IsZero()){
	    bucketwoFrontier.push_back(bddwoFrontier);
	}
    }

    //Filters all the states that are no useful anymore
    double rUseful = ratioUseful(bucketwoFrontier);
    if(rUseful > 0 && rUseful >= ratio){
	return true;
    }
    
    return false;    
}

	 
double UniformCostSearch::ratioUseful(Bucket & bucket) const{
    //1) Remove those states from bucket that are no more useful
    for(auto & bucketBDD : bucket){
	BDD newBDD = mgr->zeroBDD();
	for(const BDD & bdd : S) newBDD += bucketBDD*bdd;
	for(const BDD & bdd : Szero) newBDD += bucketBDD*bdd;
	for(const BDD & bdd : Sfilter) newBDD += bucketBDD*bdd;
	for(const BDD & bdd : Smerge) newBDD += bucketBDD*bdd;
	bucketBDD = newBDD;
    }
    double possibleStates = mgr->stateCount(bucket); //This count is not
    //exact so ratioUseful
    //may be greater than 1
    double totalStates = mgr->stateCount(S) +  mgr->stateCount(Szero) +
	mgr->stateCount(Sfilter) +  mgr->stateCount(Smerge);

    return min<double>(1, possibleStates/totalStates);
}

bool UniformCostSearch::isSearchableWithNodes(int maxNodes) const{   
    return expansionReady() && nextStepNodes() <= maxNodes;
}


void UniformCostSearch::computeEstimation(bool prepare){
    if(prepare){
	prepareBucket(/*p.max_pop_time, p.max_pop_nodes, true*/);
	DEBUG_MSG(cout << " bucket prepared for compute estimation" << endl;);
    }

    if(expansionReady()){
	//Succeded, the estimation will be only in image
	if(nodeCount(Szero)){
	    estimationZero.nextStep(nodeCount(Szero));
	}else{
	    estimationCost.nextStep(nodeCount(S));
	}
    }else{
	if (mgr->hasTransitions0()){
	    estimationZero.nextStep(nodeCount(Sfilter) + nodeCount(Smerge) +
				    nodeCount(Szero));
	}else{
	    estimationCost.nextStep(nodeCount(Sfilter) + nodeCount(Smerge) +
				    nodeCount(S));
	}
    } 
    DEBUG_MSG(cout << "estimation computed" << endl;);
}


long UniformCostSearch::nextStepTime() const{
    long estimation = 0;

    if(mgr->hasTransitions0() && (!expansionReady() || !Szero.empty())){
	estimation += estimationZero.time();
    }else {
	estimation += estimationCost.time();
    }
    return estimation;
}

long UniformCostSearch::nextStepNodes() const {
    if(mgr->hasTransitions0() && (!expansionReady() || !Szero.empty())){
	return estimationZero.nextNodes();
    }else {
	return estimationCost.nextNodes();
    }  
}

long UniformCostSearch::nextStepNodesResult() const {
    long estimation = 0;
  
    if(mgr->hasTransitions0() && (!expansionReady() || !Szero.empty())){
	estimation = max(estimation, estimationZero.nodes());
    }else {
	estimation = max(estimation, estimationCost.nodes());
    }
    return estimation;
}

/////////////////////////////////////////////////
////   Auxiliar methods to load/save/print   ////
/////////////////////////////////////////////////

std::ostream & operator<<(std::ostream &os, const UniformCostSearch & exp){
    os << "exp " << dirname(exp.isFW());
    if(exp.mgr){
	os << " in ";
	if(exp.mgr->getAbstraction()){
	    os << "abstract ";
	    os << *(exp.mgr->getAbstraction());
	}else{
	    os << "original state space";
	}
	os << " f=" << exp.getF() << flush;
	os << " g=" << exp.getG() << flush;
	os << exp.open_list << flush;
	os << " est_time: " << exp.nextStepTime() << flush;
	os << " est_nodes: " << exp.nextStepNodes() <<flush;  
    }
    return os;
}

void UniformCostSearch::printFrontier() const{
    if(!Sfilter.empty()) cout << "Sf: " << nodeCount(Sfilter) << " ";
    if(!Smerge.empty()) cout << "Sm: " << nodeCount(Smerge) << " ";
    if(!Szero.empty()) cout <<"Sz: " << nodeCount(Szero) << " ";
    if(!S.empty()) cout << "S: " << nodeCount(S) << " ";
}

bool UniformCostSearch::isBetter(const UniformCostSearch & other) const{
    return f > other.f ||
	(f == other.f && 
	 nextStepTime() < other.nextStepTime());
}

void UniformCostSearch::violated(TruncatedReason reason, double ellapsed_seconds, int maxTime, int maxNodes){
  //DEBUG_MSG(
  cout << "Truncated in " << reason << ", took " << ellapsed_seconds << " s," << 
    " maxtime: " << maxTime << " maxNodes: " << maxNodes<< endl;
  //);
  int time = 1 + ellapsed_seconds*1000;

  if(mgr->hasTransitions0() && 
     (!expansionReady() || !Szero.empty())){
      estimationZero.violated(time, maxTime, maxNodes);
  }else{
    estimationCost.violated(time, maxTime, maxNodes);
  }
}

void UniformCostSearch::filterDuplicates(Bucket & bucket) {
    //For each BDD in the bucket, get states with f
    for(auto & bdd : bucket){
	bdd *= closed->notClosed();
	DEBUG_MSG(std::cout << ", duplicates: " << bdd.nodeCount(););
	if(perfectHeuristic && 
	   perfectHeuristic->getFNotClosed() == std::numeric_limits<int>::max()){
	    bdd *= perfectHeuristic->getClosed();
	    DEBUG_MSG(std::cout << ", dead ends: " << bdd.nodeCount(););
	}
    }
}


UniformCostSearch * UniformCostSearch::getOpposite() const {
    if(perfectHeuristic)
	return perfectHeuristic->getExploration();
    else 
	return nullptr;
}

BDD UniformCostSearch::getClosedTotal() {
    return closed->getClosed();
}

BDD UniformCostSearch::notClosed(){
    return closed->notClosed();
}

std::ostream & operator<<(std::ostream &os, const TruncatedReason & reason){
    switch(reason){
    case TruncatedReason::FILTER_MUTEX:
	return os << "filter_mutex";
    case TruncatedReason::MERGE_BUCKET:
	return os << "merge";
    case TruncatedReason::MERGE_BUCKET_COST:
	return os << "merge_cost";
    case TruncatedReason::IMAGE_ZERO: 
    return os << "0-image";
  case TruncatedReason::IMAGE_COST:
    return os << "cost-image";
  default: 
    cerr << "UniformCostSearch truncated by unkown reason" << endl;
    utils::exit_with(utils::ExitCode::UNSUPPORTED);
  }
}
}
