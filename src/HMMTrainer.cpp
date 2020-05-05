/** @file HMMTrainer.cpp
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/
#include "HMMTrainer.h"
#include <limits>

HMMTrainer::HMMTrainer(){graph = NULL;}
HMMTrainer::HMMTrainer(HMMGraph* graph):graph(graph){}
HMMTrainer::~HMMTrainer(){}

void HMMTrainer::calculateFB(std::vector<seqan::BamFileIn>& alignmentSetsIn,
                             std::vector<seqan::BamAlignmentRecord>& curRecords,
                             const std::vector<seqan::FaiIndex>& readFAIs, unsigned thread){
    if(graph == NULL) return;
    
    //Starting F/B calculations in separate threads (i.e., for each read alignment)
    for(size_t curSet = 0; curSet < alignmentSetsIn.size(); ++curSet){
        if(atEnd(alignmentSetsIn[curSet]) || curRecords[curSet].rID != graph->contigId) continue;
        
        bool shouldPolish = true;
        while(shouldPolish){
            ind_prec readIndex = 0; //for threads
            std::vector<std::thread> threads;
            std::vector<Read> reads;
            
            shouldPolish = fillBuffer(alignmentSetsIn[curSet], readFAIs[curSet], curRecords[curSet], reads,
                                      graph->contigId, graph->params.mapQ, 50000);
            
            //alignedReads.size will be 0 if there is no more alignment record left to read in the current file
            for(unsigned i = 0; i < thread && i < (unsigned)reads.size(); ++i)
                threads.push_back(std::thread(&HMMTrainer::fbThreadPool, this, std::ref(reads), std::ref(readIndex)));

            //buffer is cleared here. every thread needs to wait before the buffer gets reloaded again
            for(size_t i = 0; i < threads.size(); ++i) threads[i].join();
        }
    }//F/B calculation is now done
}

void HMMTrainer::maximizeEM(){
    if(graph == NULL) return;
    
    for(ind_prec curState = 0; curState < graph->numberOfStates; ++curState){
        if(graph->seqGraph[curState].isLastInsertionState())
            //for the last insertion state, the insertion probs change so that it wont have insertion transition
            graph->preCalculatedTransitionProbs[1] = graph->params.matchTransition + graph->params.insertionTransition;
        
        //if this state ever processed then its probs may need to be updated
        if(graph->stateProcessedCount[curState] > 0){
            for(ind_prec curTransition = (graph->seqGraph[curState].isLastInsertionState())?1:0;
                curTransition < graph->numOfTransitionsPerState; ++curTransition){

                graph->transitionProbs[curState][curTransition] = (graph->transitionProcessedCount[curState][curTransition] > 0)?
                graph->transitionProbs[curState][curTransition]/graph->transitionProcessedCount[curState][curTransition]:
                graph->preCalculatedTransitionProbs[curTransition];
            }

            graph->emissionProbs[curState][A] /= graph->stateProcessedCount[curState];
            graph->emissionProbs[curState][T] /= graph->stateProcessedCount[curState];
            graph->emissionProbs[curState][G] /= graph->stateProcessedCount[curState];
            graph->emissionProbs[curState][C] /= graph->stateProcessedCount[curState];
        }else{ //initial probs to be set unless this state has been processed
            for(ind_prec curTr = (graph->seqGraph[curState].isLastInsertionState())?1:0;
                curTr < graph->numOfTransitionsPerState; ++curTr)
                graph->transitionProbs[curState][curTr] = graph->preCalculatedTransitionProbs[curTr];

            graph->emissionProbs[curState][A] = graph->seqGraph[curState].getEmissionProb('A', graph->params);
            graph->emissionProbs[curState][T] = graph->seqGraph[curState].getEmissionProb('T', graph->params);
            graph->emissionProbs[curState][G] = graph->seqGraph[curState].getEmissionProb('G', graph->params);
            graph->emissionProbs[curState][C] = graph->seqGraph[curState].getEmissionProb('C', graph->params);
        }

        if(graph->seqGraph[curState].isLastInsertionState())
            graph->preCalculatedTransitionProbs[1] = graph->params.matchTransition;
    }

    //@dynamic init
    graph->maxEmissionProbs = new std::pair<prob_prec, char>[graph->numberOfStates];
    for(unsigned curState = 0; curState < graph->numberOfStates; ++curState){
        if(graph->emissionProbs[curState][A] >= std::max(graph->emissionProbs[curState][T],
                                                         std::max(graph->emissionProbs[curState][G],
                                                                  graph->emissionProbs[curState][C]))){
            graph->maxEmissionProbs[curState] = std::make_pair(graph->emissionProbs[curState][A], 'A');
        }else if(graph->emissionProbs[curState][T] >= std::max(graph->emissionProbs[curState][A],
                                                        std::max(graph->emissionProbs[curState][G],
                                                                 graph->emissionProbs[curState][C]))){
            graph->maxEmissionProbs[curState] = std::make_pair(graph->emissionProbs[curState][T], 'T');
        }else if(graph->emissionProbs[curState][G] >= std::max(graph->emissionProbs[curState][A],
                                                        std::max(graph->emissionProbs[curState][T],
                                                                 graph->emissionProbs[curState][C]))){
            graph->maxEmissionProbs[curState] = std::make_pair(graph->emissionProbs[curState][G], 'G');
        }else{
            graph->maxEmissionProbs[curState] = std::make_pair(graph->emissionProbs[curState][C], 'C');
        }
    }
}

void HMMTrainer::fbThreadPool(const std::vector<Read>& reads, ind_prec& readIndex){
    if(graph == NULL) return;
    
    std::vector<TransitionInfoNode> curStateTransitions;
    //@dynamic init
    prob_prec* curStateTransitionLikelihood = new prob_prec[graph->numOfTransitionsPerState];
    //@dynamic init
    prob_prec* curStateEmissionProbs = new prob_prec[totalNuc];
    
    ind_prec curRead = 0; //0-based. refers to the current read id
    {//block for obtaining the next aligned read
        std::lock_guard<std::mutex> lk(indexMutex);
        curRead = readIndex++;
    }
    
    while(curRead < (ind_prec)reads.size()){
        
        //to which character should correction extend at maximum
        //decide whether this is a short read (fragment) or long
        ind_prec readLength = (ind_prec) seqan::length(reads[curRead].read);
        ind_prec maxTransition = readLength + readLength/20 + 1;
        if(readLength < 500) maxTransition = readLength + readLength/3 + 1;

        ind_prec maxDistanceOnAssembly = std::min(reads[curRead].pos+maxTransition, graph->contigLength);
        //states prior to this wont be processed. offset value is to ignore these states
        ind_prec offset = MATCH_OFFSET(reads[curRead].pos, -1, graph->params.maxInsertion);
        //maximum number of states to be processed
        ind_prec fbMatrixSize = MATCH_OFFSET(maxDistanceOnAssembly, 1, graph->params.maxInsertion) - offset;

        if(fbMatrixSize > 0){
            ind_prec j; //j value in i,j transitions
            //@dynamic init
            prob_prec** forwardMatrix = new prob_prec*[readLength]; //forwardMatrix[time][state]
            for(ind_prec i = 0; i < readLength; i++){
                //@dynamic init
                forwardMatrix[i] = new prob_prec[fbMatrixSize];
                std::fill_n(forwardMatrix[i], fbMatrixSize, 0);
            }
            //@dynamic init
            prob_prec** backwardMatrix = new prob_prec*[readLength]; //backwardMatrix[time][state]
            for(ind_prec i = 0; i < readLength; i++){
                //@dynamic init
                backwardMatrix[i] = new prob_prec[fbMatrixSize];
                std::fill_n(backwardMatrix[i], fbMatrixSize, 0);
            }
            
            ind_prec startForBackward = fillForwardMatrix(forwardMatrix, toCString(reads[curRead].read), offset,
                                                          MATCH_OFFSET(maxDistanceOnAssembly, 1,
                                                                       graph->params.maxInsertion), readLength);

            if(startForBackward != 0 && graph->seqGraph[startForBackward].getCharIndex() > reads[curRead].pos){
                startForBackward = MATCH_OFFSET(graph->seqGraph[startForBackward].getCharIndex(), 1,
                                                graph->params.maxInsertion);

                if(fillBackwardMatrix(backwardMatrix, toCString(reads[curRead].read), startForBackward, offset,
                                      readLength)){
                    //updating probabilities wrt the f/b matrices computed just now
                    for(ind_prec curState = INSERTION_OFFSET(reads[curRead].pos, -1, 1, graph->params.maxInsertion);
                        curState < startForBackward; ++curState){

                        if(graph->seqGraph[curState].isLastInsertionState())
                            //for the last insertion state, the insertion probs change
                            graph->preCalculatedTransitionProbs[1] = graph->params.matchTransition +
                                                                     graph->params.insertionTransition;

                        if(curState-offset < fbMatrixSize){
                            ind_prec matchoff = MATCH_OFFSET(graph->seqGraph[curState].getCharIndex(), 0,
                                                        graph->params.maxInsertion);
                            std::fill_n(curStateTransitionLikelihood, graph->numOfTransitionsPerState, 0);
                            std::fill_n(curStateEmissionProbs, totalNuc, 0);
                            insertNewForwardTransitions(&curStateTransitions, graph->seqGraph[curState],
                                                        graph->params.maxDeletion, graph->params.maxInsertion);

                            for(ind_prec t = 0; t < readLength; ++t){
                                //transition probabilities
                                if(t < readLength-1){
                                    for(size_t curTr = 0; curTr < curStateTransitions.size(); ++curTr){
                                        if(curStateTransitions.at(curTr).toState - offset < fbMatrixSize){
                                            j = curStateTransitions[curTr].toState;
                                            //0->insertion, 1-> match, 2,3...->deletions
                                            ind_prec transitionIndex = (j - matchoff)/(graph->params.maxInsertion+1);
                                            curStateTransitionLikelihood[transitionIndex] +=
                                            forwardMatrix[t][curState-offset]*
                                            graph->preCalculatedTransitionProbs[transitionIndex]*
                                            graph->seqGraph[j].getEmissionProb(reads[curRead].read[t+1],
                                                                               graph->params)*backwardMatrix[t+1][j-offset];

                                        }
                                    }
                                }

                                //emission probabilities
                                char emitChar = (reads[curRead].read[t] != 'N')?
                                reads[curRead].read[t]:
                                (graph->seqGraph[curState].isMatchState())?
                                graph->seqGraph[curState].getNucleotide():'\0';
                                Nucleotide chosenNuc = (emitChar == 'A' || emitChar == 'a')?A:
                                (emitChar == 'T' || emitChar == 't')?T:
                                (emitChar == 'G' || emitChar == 'g')?G:
                                (emitChar == 'C' || emitChar == 'c')?C:totalNuc;
                                if(chosenNuc < totalNuc)
                                    curStateEmissionProbs[chosenNuc] +=
                                    forwardMatrix[t][curState-offset]*backwardMatrix[t][curState-offset];
                            }
                            curStateTransitions.clear();
                            prob_prec totalEmissionProbs = curStateEmissionProbs[A] + curStateEmissionProbs[T] +
                            curStateEmissionProbs[G] + curStateEmissionProbs[C];
                            prob_prec totalTransitionLikelihoods = 0;
                            prob_prec processedTransitionProb = 0;
                            for(ind_prec i = (graph->seqGraph[curState].isLastInsertionState())?1:0;
                                i < graph->numOfTransitionsPerState; ++i){
                                if(curStateTransitionLikelihood[i] > 0 || i == 0){
                                    totalTransitionLikelihoods += curStateTransitionLikelihood[i];
                                    processedTransitionProb += graph->preCalculatedTransitionProbs[i];
                                }
                            }

                            if(totalEmissionProbs != 0 && curState < startForBackward){
                                if(totalTransitionLikelihoods != 0){
                                    {//block for updating the transition probs and the transition proccessed count
                                        std::lock_guard<std::mutex> lk(transitionProbMutex);

                                        for(ind_prec i = (graph->seqGraph[curState].isLastInsertionState())?1:0;
                                            i < graph->numOfTransitionsPerState; ++i){
                                            if(curStateTransitionLikelihood[i] > 0 || i == 0){
                                                graph->transitionProbs[curState][i] +=
                                                (curStateTransitionLikelihood[i]/totalTransitionLikelihoods)*
                                                processedTransitionProb;
                                                graph->transitionProcessedCount[curState][i]++;
                                            }
                                        }
                                    }
                                }

                                {//block for updating the emission probs and the state processed count
                                    std::lock_guard<std::mutex> lk(emissionProbMutex);

                                    graph->emissionProbs[curState][A] += curStateEmissionProbs[A]/totalEmissionProbs;
                                    graph->emissionProbs[curState][T] += curStateEmissionProbs[T]/totalEmissionProbs;
                                    graph->emissionProbs[curState][G] += curStateEmissionProbs[G]/totalEmissionProbs;
                                    graph->emissionProbs[curState][C] += curStateEmissionProbs[C]/totalEmissionProbs;

                                    graph->stateProcessedCount[curState]++;
                                }
                            }
                        }

                        //for the last insertion state, the insertion probs change so that it wont have insertion
                        //transition. putting it back to normal now
                        if(graph->seqGraph[curState].isLastInsertionState())
                            graph->preCalculatedTransitionProbs[1] = graph->params.matchTransition;
                    }
                }
            }

            for(ind_prec i = 0; i < readLength; ++i) { delete[] backwardMatrix[i]; delete[] forwardMatrix[i];}
            delete[] backwardMatrix;
            delete[] forwardMatrix;
        }
        
        {//block for obtaining the next aligned read
            std::lock_guard<std::mutex> lk(indexMutex);
            curRead = readIndex++;
        }
    }

    delete[] curStateTransitionLikelihood;
    delete[] curStateEmissionProbs;
}

ind_prec HMMTrainer::fillForwardMatrix(prob_prec** forwardMatrix, const char* read, const ind_prec startPosition,
                                       ind_prec maxDistanceOnContig, const ind_prec readLength){
    
    if(startPosition == 0 || startPosition >= maxDistanceOnContig) return 0;
    if(graph->numberOfStates < maxDistanceOnContig)
        maxDistanceOnContig = graph->numberOfStates;

    prob_prec maxPrec = std::numeric_limits<prob_prec>::max() / 100;
    
    //which transitions requested from the previous time
    std::vector<TransitionInfoNode>* curTrSet = new std::vector<TransitionInfoNode>;
    //which transitions to be made for the next time
    std::vector<TransitionInfoNode>* nextTrSet = new std::vector<TransitionInfoNode>;
    std::vector<TransitionInfoNode>* tmpTrSet;

    bool* allowedParentStates = new bool[maxDistanceOnContig-startPosition+1];
    std::vector<bool> hasStateBeenProcessedBefore;
    hasStateBeenProcessedBefore.resize(maxDistanceOnContig-startPosition+1);
    std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);
    std::fill_n(allowedParentStates, maxDistanceOnContig-startPosition+1, false);

    //1-initialization (t = 1)
    ind_prec curTime = 0; //represents the current time (1...T)
    insertNewForwardTransitions(curTrSet, graph->seqGraph[startPosition], graph->params.maxDeletion,
                                graph->params.maxInsertion);
    for(size_t curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
        ind_prec matchoff = MATCH_OFFSET(graph->seqGraph[curTrSet->at(curTransition).from].getCharIndex(), 0,
                                         graph->params.maxInsertion);
        //0->insertion, 1-> match, 2,3...->deletions
        ind_prec transitionIndex = (curTrSet->at(curTransition).toState - matchoff)/(graph->params.maxInsertion+1);
        //@IMPORTANT: check maxPrec here.
        forwardMatrix[0][curTrSet->at(curTransition).toState - startPosition] += maxPrec *
        graph->preCalculatedTransitionProbs[transitionIndex] *
        graph->seqGraph[curTrSet->at(curTransition).toState].getEmissionProb(read[curTime], graph->params);

        insertNewForwardTransitions(nextTrSet, graph->seqGraph[curTrSet->at(curTransition).toState],
                                    graph->params.maxDeletion, graph->params.maxInsertion);
    }
    curTrSet->clear();

    //find the most likely states that should be allowed to make the next transitions
    findMaxValues(forwardMatrix[0], allowedParentStates, startPosition, maxDistanceOnContig, graph->params.filterSize);

    tmpTrSet = curTrSet;
    curTrSet = nextTrSet;
    nextTrSet = tmpTrSet;

    //2-recursion (1 < t <= T)
    while (curTime < readLength-1 && !curTrSet->empty()){
        curTime++;
        for(size_t curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
            const TransitionInfoNode& frontTr = curTrSet->at(curTransition);
            if(allowedParentStates[frontTr.from - startPosition] && frontTr.toState < maxDistanceOnContig){
                ind_prec matchoff = MATCH_OFFSET(graph->seqGraph[frontTr.from].getCharIndex(), 0, graph->params.maxInsertion);
                //0->insertion, 1-> match, 2,3...->deletions
                ind_prec transitionIndex = (frontTr.toState - matchoff)/(graph->params.maxInsertion+1);
                forwardMatrix[curTime][frontTr.toState-startPosition] +=
                    forwardMatrix[curTime-1][frontTr.from-startPosition]*
                    graph->preCalculatedTransitionProbs[transitionIndex]*
                    graph->seqGraph[frontTr.toState].getEmissionProb(read[curTime], graph->params);

                if(!hasStateBeenProcessedBefore[frontTr.toState - startPosition]){
                    insertNewForwardTransitions(nextTrSet, graph->seqGraph[frontTr.toState], graph->params.maxDeletion,
                                                graph->params.maxInsertion);
                    hasStateBeenProcessedBefore[frontTr.toState - startPosition] = true;
                }
            }
        }
        curTrSet->clear();

        std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);
        std::fill_n(allowedParentStates, maxDistanceOnContig-startPosition+1, false);

        //find the most likely states that should be allowed to make the next transitions
        findMaxValues(forwardMatrix[curTime], allowedParentStates, startPosition, maxDistanceOnContig,
                      graph->params.filterSize);
        tmpTrSet = curTrSet;
        curTrSet = nextTrSet;
        nextTrSet = tmpTrSet;
    }

    ind_prec max = 0;
    for(ind_prec curStateForMax = 0; curStateForMax < maxDistanceOnContig-startPosition+1; ++curStateForMax){
        if(allowedParentStates[curStateForMax] &&
           (max == 0 || forwardMatrix[curTime][curStateForMax] > forwardMatrix[curTime][max])){
            max = curStateForMax;
        }
    }

    delete curTrSet; delete nextTrSet; delete[] allowedParentStates;

    return max + startPosition;
}

bool HMMTrainer::fillBackwardMatrix(prob_prec** backwardMatrix, const char* read, ind_prec startPosition,
                                    ind_prec maxDistanceOnContig, const ind_prec readLength){
    
    if(maxDistanceOnContig == 0 || maxDistanceOnContig >= startPosition) return false;
    if(graph->numberOfStates < startPosition)
        startPosition = graph->numberOfStates;
    
    prob_prec maxPrec = std::numeric_limits<prob_prec>::max() / 100;
    
    //which transitions requested from the previous time
    std::vector<TransitionInfoNode>* curTrSet = new std::vector<TransitionInfoNode>;
    //which transitions to be made for the next time
    std::vector<TransitionInfoNode>* nextTrSet = new std::vector<TransitionInfoNode>;
    std::vector<TransitionInfoNode>* tmpTrSet;

    bool* allowedParentStates = new bool[startPosition - maxDistanceOnContig + 1];
    std::vector<bool> hasStateBeenProcessedBefore;
    hasStateBeenProcessedBefore.resize(startPosition - maxDistanceOnContig + 1);
    std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);
    std::fill_n(allowedParentStates, startPosition - maxDistanceOnContig + 1, false);

    //1-initialization
    ind_prec curTime = readLength-1; //@curTime value is 0-based. So 0th index is the first character [T....1]
    insertNewBackwardTransitions(curTrSet, graph->seqGraph[startPosition], graph->params.maxDeletion,
                                 graph->params.maxInsertion);
    for(size_t curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
        if(curTrSet->at(curTransition).toState > maxDistanceOnContig){
            ind_prec matchoff = MATCH_OFFSET(graph->seqGraph[curTrSet->at(curTransition).toState].getCharIndex(), 0,
                                        graph->params.maxInsertion);
            //0->insertion, 1-> match, 2,3...->deletions
            ind_prec transitionIndex = (curTrSet->at(curTransition).from - matchoff)/(graph->params.maxInsertion+1);
            //@IMPORTANT: check maxPrec here.
            backwardMatrix[curTime][curTrSet->at(curTransition).toState - maxDistanceOnContig] += maxPrec *
            graph->preCalculatedTransitionProbs[transitionIndex];

            insertNewBackwardTransitions(nextTrSet, graph->seqGraph[curTrSet->at(curTransition).toState],
                                         graph->params.maxDeletion, graph->params.maxInsertion);
        }
    }
    curTrSet->clear();

    findMaxValues(backwardMatrix[curTime], allowedParentStates, maxDistanceOnContig, startPosition,
                  graph->params.filterSize);
    tmpTrSet = curTrSet;
    curTrSet = nextTrSet;
    nextTrSet = tmpTrSet;

    //2-recursion
    while (curTime > 0 && !curTrSet->empty()){
        curTime--;
        for(size_t curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
            const TransitionInfoNode& frontTr = curTrSet->at(curTransition);
            if(allowedParentStates[frontTr.from - maxDistanceOnContig] && frontTr.toState > maxDistanceOnContig) {
                ind_prec matchoff = MATCH_OFFSET(graph->seqGraph[frontTr.toState].getCharIndex(), 0,
                                                 graph->params.maxInsertion);
                //0->insertion, 1-> match, 2,3...->deletions
                ind_prec transitionIndex = (frontTr.from - matchoff)/(graph->params.maxInsertion+1);
                backwardMatrix[curTime][frontTr.toState - maxDistanceOnContig] +=
                    backwardMatrix[curTime+1][frontTr.from - maxDistanceOnContig]*
                    graph->preCalculatedTransitionProbs[transitionIndex]*
                    graph->seqGraph[frontTr.from].getEmissionProb(read[curTime+1], graph->params);

                if(!hasStateBeenProcessedBefore[frontTr.toState - maxDistanceOnContig]){
                    insertNewBackwardTransitions(nextTrSet, graph->seqGraph[frontTr.toState], graph->params.maxDeletion,
                                                 graph->params.maxInsertion);
                    hasStateBeenProcessedBefore[frontTr.toState - maxDistanceOnContig] = true;
                }
            }
        }
        curTrSet->clear();

        std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);
        std::fill_n(allowedParentStates, startPosition - maxDistanceOnContig + 1, false);
        findMaxValues(backwardMatrix[curTime], allowedParentStates, maxDistanceOnContig, startPosition,
                      graph->params.filterSize);
        tmpTrSet = curTrSet;
        curTrSet = nextTrSet;
        nextTrSet = tmpTrSet;
    }

    delete curTrSet; delete nextTrSet; delete[] allowedParentStates;

    return true;
}
