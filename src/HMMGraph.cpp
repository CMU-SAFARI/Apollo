/** @file HMMGraph.cpp
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/

#include "HMMGraph.h"

HMMGraph::HMMGraph( const HMMParameters& parameters, const seqan::Dna5String& contig,
                   int32_t contigId):params(parameters), contig(contig), contigId(contigId){
    isParametersSet = true;
    isDataFreed = true;
}

HMMGraph::HMMGraph(){
    isParametersSet = false;
    isDataFreed = true;
}

HMMGraph::~HMMGraph(){
    this->freeGraph();
}

void HMMGraph::setParameters(const HMMParameters& parameters){
    params = parameters;
    isParametersSet = true;
}

bool HMMGraph::getIsParametersSet(){
    return isParametersSet;
}

void HMMGraph::setContig(const seqan::Dna5String& contig, int32_t contigId){
    this->contig = contig;
    this->contigId = contigId;
}

void HMMGraph::buildGraph(){
    
    if(!isParametersSet || length(contig) <= 0) return;
    this->freeGraph();
    
    contigLength = (ind_prec)length(contig);
    numberOfStates = GRAPH_SIZE(contigLength, params.maxInsertion);
    seqGraph = new SeqNode[numberOfStates];

    //constructing the pHMM graph
    seqGraph[0] = SeqNode(0, 0, params.maxInsertion, '\0', contig[0]);
    for(cnt_prec curIn = 1; curIn <= params.maxInsertion; ++curIn)
        seqGraph[curIn] = SeqNode(curIn, 0, params.maxInsertion, '\0', contig[0]);

    ind_prec matchOffset;
    //states in between the start and end states
    for(ind_prec curCharacter = 1; curCharacter < contigLength; ++curCharacter){
        matchOffset = MATCH_OFFSET(curCharacter, 0, params.maxInsertion);
        seqGraph[matchOffset] = SeqNode(matchOffset, curCharacter, params.maxInsertion, contig[curCharacter-1],
                                        contig[curCharacter]);
        for(cnt_prec curIn = 1; curIn <= params.maxInsertion; ++curIn)
            seqGraph[matchOffset+curIn] = SeqNode(matchOffset+curIn, curCharacter,params.maxInsertion,
                                                  contig[curCharacter-1], contig[curCharacter]);
    }
    
    //last character before the end state
    matchOffset = MATCH_OFFSET(contigLength, 0, params.maxInsertion);
    seqGraph[matchOffset] = SeqNode(matchOffset, contigLength, params.maxInsertion, contig[contigLength-1], '\0');
    for(cnt_prec curIn = 1; curIn <= params.maxInsertion; ++curIn)
        seqGraph[matchOffset+curIn] = SeqNode(matchOffset+curIn, contigLength, params.maxInsertion,
                                              contig[contigLength-1], '\0');

    //end state
    matchOffset = END_STATE(contigLength, params.maxInsertion);
    seqGraph[matchOffset] = SeqNode(matchOffset, contigLength+1, params.maxInsertion, '\0', '\0');
    std::vector<TransitionInfoNode> stateTransitions; //prefetching all possible transitions from a state to another
    insertNewForwardTransitions(&stateTransitions, seqGraph[0], params.maxDeletion, params.maxInsertion);
    numOfTransitionsPerState = (ind_prec)stateTransitions.size();

    //pre calculated initial transition values for any state in hmm
    preCalculatedTransitionProbs = new prob_prec[numOfTransitionsPerState];
    for(size_t curTr = 0; curTr < numOfTransitionsPerState; ++curTr)
        preCalculatedTransitionProbs[curTr]=seqGraph[0].transitionProbFromThisNode(stateTransitions[curTr].toState,
                                                                                   params);
    
    //array of transitionProbs[from state][to state]
    transitionProbs = new prob_prec*[numberOfStates];
    //array of emissionProbs[state][character index]
    emissionProbs = new prob_prec*[numberOfStates];
    //how many times a state has been processed [state]
    stateProcessedCount = new cnt_prec[numberOfStates];
    std::fill_n(stateProcessedCount, numberOfStates, 0);
    //how many times a transition has been processed [state][transition]
    transitionProcessedCount = new cnt_prec*[numberOfStates];
    for(ind_prec curState = 0; curState < numberOfStates; ++curState){
        transitionProcessedCount[curState] = new cnt_prec[numOfTransitionsPerState];
        std::fill_n(transitionProcessedCount[curState], numOfTransitionsPerState, 0);
        transitionProbs[curState] = new prob_prec[numOfTransitionsPerState];
        emissionProbs[curState] = new prob_prec[totalNuc];
        std::fill_n(transitionProbs[curState], numOfTransitionsPerState, 0);
        std::fill_n(emissionProbs[curState], totalNuc, 0);
    }
    
    isDataFreed = false;
}

void HMMGraph::freeGraph(){
    if(isDataFreed) return;
    
    for(ind_prec i = 0; i < numberOfStates; ++i){
        delete[] transitionProbs[i]; delete[] emissionProbs[i]; delete[] transitionProcessedCount[i];
    }
    
    delete[] seqGraph;
    delete[] preCalculatedTransitionProbs;
    delete[] transitionProbs;
    delete[] emissionProbs;
    delete[] stateProcessedCount;
    delete[] transitionProcessedCount;
    delete[] maxEmissionProbs;
    
    isDataFreed = true;
}

