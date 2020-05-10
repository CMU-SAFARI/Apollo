/** @file HMMDecoder.cpp
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/

#include "HMMDecoder.h"

HMMDecoder::HMMDecoder(){graph = NULL;}
HMMDecoder::HMMDecoder(HMMGraph* graph):graph(graph){}
HMMDecoder::~HMMDecoder(){}

void HMMDecoder::backtrace(unsigned threadCnt, seqan::String<seqan::Dna5String>& decodedOut){
    if(graph == NULL || graph->contigLength <= 0) return;
    
    std::vector<std::thread> threads;
    ind_prec seqIndex = 0;
    ind_prec batchSize = (graph->params.batchSize == 0)?graph->contigLength:graph->params.batchSize;
    resize(decodedOut, ((graph->contigLength-1)/batchSize) + 1);

    for(unsigned i = 0; i < threadCnt; ++i){
        threads.push_back(std::thread(&HMMDecoder::viterbiThreadPool, this, batchSize, std::ref(seqIndex),
                                      std::ref(decodedOut)));
    }
    for(unsigned i = 0; i < threadCnt; ++i) threads[i].join();
}

void HMMDecoder::viterbiThreadPool(ind_prec seqLength, ind_prec& seqIndex,seqan::String<seqan::Dna5String>& decodedOut){
    
    ind_prec curSequence = 0; //0-based. refers to the current piece of sequence
    {//block for obtaining new index for the next contig chunk
        std::lock_guard<std::mutex> lk(indexMutex);
        curSequence = seqIndex++;
    }

    while(curSequence < length(decodedOut)){

        std::vector<std::map<ind_prec, ind_prec> > backtrace; //index
        ind_prec startState = MATCH_OFFSET(curSequence*seqLength, 0, graph->params.maxInsertion);
        ind_prec endState = MATCH_OFFSET(((curSequence+1)* seqLength) + 1, 0, graph->params.maxInsertion);

        if(endState > graph->numberOfStates-1) endState = graph->numberOfStates-1;
        ind_prec totalStates = endState - startState + 1;

        std::vector<bool> allowedParentStates; //index - offset
        allowedParentStates.resize(totalStates);
        std::fill_n(allowedParentStates.begin(), allowedParentStates.size(), false);

        std::set<ind_prec>* curTrSet = new std::set<ind_prec>(); //index
        std::set<ind_prec>* nextTrSet = new std::set<ind_prec>(); //index
        std::set<ind_prec>* tmpTrSet; //index

        std::map<ind_prec, prob_prec>* prevViterbi = new std::map<ind_prec, prob_prec>(); //index
        std::map<ind_prec, prob_prec>* curViterbi = new std::map<ind_prec, prob_prec>(); //index
        std::map<ind_prec, prob_prec>* tmpViterbi; //index
        std::map<ind_prec, prob_prec>::iterator itVit;

        //initialization step t = 1 (curTime = 0)
        std::vector<TransitionInfoNode> curStateTransitions; //pushing real index
        insertNewForwardTransitions(&curStateTransitions, graph->seqGraph[startState], graph->params.maxDeletion,
                                    graph->params.maxInsertion);
        backtrace.push_back(std::map<ind_prec, ind_prec >());

        ind_prec curTime = 0;
        prevViterbi->clear();
        while(!curStateTransitions.empty()){

            ind_prec next = curStateTransitions.begin()->toState;
            ind_prec transitionIndex = (next - startState)/(graph->params.maxInsertion+1); //0:ins., 1: match, 2..:dels

            if(next <= endState){
                (*prevViterbi)[next] = graph->transitionProbs[startState][transitionIndex] +
                graph->maxEmissionProbs[next].first;

                curTrSet->insert(next);
                allowedParentStates[next-startState] = true;
                backtrace[curTime][next] = curStateTransitions.begin()->from; //offset
            }
            curStateTransitions.erase(curStateTransitions.begin());
        }
        curStateTransitions.clear();
        curTime++;

        bool finished = false;
        ind_prec maxFinalStateTime = 0;
        prob_prec maxFinalStateValue = INT_MIN;
        while(!finished && curTime < (graph->params.maxInsertion+1)*seqLength){

            //new storage for the current time
            curViterbi->clear();
            backtrace.push_back(std::map<ind_prec, ind_prec>());

            for(std::set<ind_prec>::iterator itSet = (*curTrSet).begin(); itSet != curTrSet->end(); ++itSet){
                //calculate the viterbi values from this state to the states that it has transitions
                ind_prec fromState = (*itSet);
                //if this state is one of the major ones from the previous time
                if(allowedParentStates[fromState-startState]){

                    insertNewForwardTransitions(&curStateTransitions, graph->seqGraph[fromState],
                                                graph->params.maxDeletion, graph->params.maxInsertion);
                    for(size_t toStateIndex = 0; toStateIndex < curStateTransitions.size(); ++toStateIndex){

                        ind_prec toState = curStateTransitions[toStateIndex].toState;
                        ind_prec matchoff = MATCH_OFFSET(graph->seqGraph[fromState].getCharIndex(), 0,
                                                         graph->params.maxInsertion);
                        //0->insertion, 1-> match, 2,3...->deletions
                        ind_prec transitionIndex = (toState - matchoff)/(graph->params.maxInsertion+1);

                        if(toState <= endState){

                            prob_prec newViterbi = (*prevViterbi)[fromState] +
                            graph->transitionProbs[fromState][transitionIndex] + graph->maxEmissionProbs[toState].first;

                            itVit = curViterbi->find(toState);
                            if(itVit == curViterbi->end() || (*curViterbi)[toState] < newViterbi){
                                //regular viterbi here
                                (*curViterbi)[toState] = newViterbi;
                                backtrace[curTime][toState] = fromState;
                                (*nextTrSet).insert(toState);
                            }
                        }
                    }

                    curStateTransitions.clear();
                }
            }

            std::fill_n(allowedParentStates.begin(), allowedParentStates.size(), false);
            curTrSet->clear();

            //find the most possible states that will contribute to the viterbi of the next time
            std::priority_queue<std::pair<prob_prec, ind_prec>, std::vector<std::pair<prob_prec, ind_prec> >,
            std::greater<std::pair<prob_prec, ind_prec> > > maxValues;
            for(std::set<ind_prec>::iterator itSet = (*nextTrSet).begin(); itSet != (*nextTrSet).end(); ++itSet) {
                if(maxValues.size() < (size_t)graph->params.viterbiFilterSize)
                    maxValues.push(std::pair<prob_prec, ind_prec>((*curViterbi)[*itSet], *itSet));
                else if(maxValues.top().first < (*curViterbi)[*itSet]){
                    maxValues.pop();
                    maxValues.push(std::pair<prob_prec, ind_prec>((*curViterbi)[*itSet], *itSet));
                }
            }

            while(!maxValues.empty()){
                allowedParentStates.operator[](maxValues.top().second - startState) = true;
                maxValues.pop();
            }

            if(allowedParentStates[totalStates-1] &&
               (maxFinalStateTime == 0 || (*curViterbi)[endState] > maxFinalStateValue)){
                maxFinalStateTime = curTime;
                maxFinalStateValue = (*curViterbi)[endState];
            }

            //stop condition: we've reached final state and max valued final state is "way" left behind the current time
            if(maxFinalStateTime > 0 && curTime > maxFinalStateTime + 20) finished = true;
            else curTime++;

            //swaps
            tmpTrSet = curTrSet; curTrSet = nextTrSet; nextTrSet = tmpTrSet;
            tmpViterbi = curViterbi; curViterbi = prevViterbi; prevViterbi = tmpViterbi;
        }

        curTrSet->clear();
        nextTrSet->clear();
        allowedParentStates.clear();
        delete curTrSet; delete nextTrSet;
        delete curViterbi; delete prevViterbi;

        if(maxFinalStateTime > 0){
            resize(decodedOut[curSequence], maxFinalStateTime);
            ind_prec lastIndexUsed = endState;
            for(ind_prec curChar = maxFinalStateTime; curChar > 0; --curChar){
                ind_prec curIndex = backtrace[curChar][lastIndexUsed];
                decodedOut[curSequence][curChar-1] = graph->maxEmissionProbs[curIndex].second;
                lastIndexUsed = curIndex;
            }
        }

        {//block for obtaining new index for the next contig chunk
            std::lock_guard<std::mutex> lk(indexMutex);
            curSequence = seqIndex++;
        }
    }
}
