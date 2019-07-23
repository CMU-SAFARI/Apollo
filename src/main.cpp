/** @file main.cpp
 * @brief One sentence brief
 *
 * More details
 * In multiple lines
 * Copyright © 2019 Can Fırtına. All rights reserved.
 *
 * @author Can Firtina
 * @bug No known bug
 */

#include <string>
#include <math.h>
#include <set>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include <map>
#include <thread>
#include "CommandLineParser.h"
#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/sequence.h>
#include <seqan/sequence/string_base.h>
#include <seqan/file.h>
#include <seqan/seq_io.h>
#include <seqan/bam_io.h>
#include <seqan/bam_io/bam_index_bai.h>
#include <mutex>
#include <condition_variable>

#define MATCH_OFFSET(char, offset, insize) ((char+offset)*(insize+1))
#define INSERTION_OFFSET(char, offset, innumber, insize) ((char+offset)*(insize+1) + innumber) //innumber 1 based
#define GRAPH_SIZE(size, insize) ((size+1)*(insize+1) + 1)
#define END_STATE(size, insize) (GRAPH_SIZE(size, insize)-1)

std::mutex indexMutex;
std::mutex emissionProbMutex;
std::mutex transitionProbMutex;

enum Nucleotide {A, T, G, C, totalNuc};

/** @brief Parameters to be used in profile hidden Markov model graph
 *
 */
struct HMMParameters{

public:

    HMMParameters(unsigned int filterSize, unsigned int viterbiFilterSize, unsigned int maxDeletion, unsigned int
                  maxInsertion, unsigned int batchSize, double matchTransition, double insertionTransition,
                  double deletionTransitionFactor, double matchEmission):
    filterSize(filterSize), viterbiFilterSize(viterbiFilterSize), maxDeletion(maxDeletion),
    maxInsertion(maxInsertion), batchSize(batchSize), matchTransition(matchTransition), insertionTransition(insertionTransition),
    deletionTransitionFactor(deletionTransitionFactor), matchEmission(matchEmission){
        deletionTransition = 1.000 - (matchTransition + insertionTransition);
        mismatchEmission = (double)(1 - matchEmission)/3.00;
        insertionEmission = (double)1/3.00; //total nucleotide = 4; Emission prob for each except one
    }

    HMMParameters(const HMMParameters& cpy):
    filterSize(cpy.filterSize), viterbiFilterSize(cpy.viterbiFilterSize), maxDeletion(cpy.maxDeletion),
    maxInsertion(cpy.maxInsertion), batchSize(cpy.batchSize), matchTransition(cpy.matchTransition),
    insertionTransition(cpy.insertionTransition), deletionTransition(cpy.deletionTransition),
    deletionTransitionFactor(cpy.deletionTransitionFactor), matchEmission(cpy.matchEmission),
    mismatchEmission(cpy.mismatchEmission), insertionEmission(cpy.insertionEmission){}

    unsigned filterSize;
    unsigned viterbiFilterSize;
    unsigned maxDeletion;
    unsigned maxInsertion;
    unsigned batchSize;
    double matchTransition;
    double insertionTransition;
    double deletionTransition;
    double deletionTransitionFactor;
    double matchEmission;
    double mismatchEmission;
    double insertionEmission;
};

/** @brief Represents a state in the profile hidden Markov model graph. Total size is usually 16 bytes
 *
 */
struct SeqNode{
public:

    SeqNode(){}
    SeqNode(unsigned int index, unsigned int charIndex, unsigned int insize, char nuc, char nextNuc):
    index(index), charIndex(charIndex), nuc(nuc), isMatch((index%(insize+1) == 0)?true:false),
    isLastInsertion((insize > 0 && index%(insize+1) == insize)?true:false), nextNuc(nextNuc){}

    /** @brief Emission probability calculation
     *
     *  @param character Basepair to calculate the emission probability
     *  @param parameters Calculation is based on the specified parameters
     *  @return Emission probability of the given character (basepair) [0-1]
     */
    double getEmissionProb(const char& character, const HMMParameters& parameters){

        if(isMatch){//match state: either match or mismatch emission probability
            if(character == nuc || character == 'N')
                return parameters.matchEmission;
            return parameters.mismatchEmission;
        }else if(character == 'N' || character == nextNuc){
            //insertion state: no emission for 'N' or for the next character in contig
            return 0.00;
        }

        return parameters.insertionEmission; //insertion state
    }

    /** @brief Transition probability calculation from this state to the specified `toIndex` state
     *
     *  @param toIndex Specifies which state to transit from this state
     *  @param parameters Calculation is based on the specified parameters
     *  @return Transition probability from this state to the state indicated with `toIndex`
     */
    double transitionProbFromThisNode(const unsigned int& toIndex, HMMParameters& parameters){

        if(MATCH_OFFSET(charIndex, 1, parameters.maxInsertion) == toIndex){
            if(isLastInsertion){
                //match transition prob for last insertion state
                return parameters.matchTransition + parameters.insertionTransition;
            }
            return parameters.matchTransition; //match transition prob
        }

        if(index+1 == toIndex) return parameters.insertionTransition; //insertion transition

        //deletion transition calculations: normalized polynomial distribution
        int count = 0;
        int start = 1;
        double transitionProb = 0;
        for(int curDel = parameters.maxDeletion+1; curDel > 1; --curDel){
            if(MATCH_OFFSET(charIndex, curDel, parameters.maxInsertion) == toIndex)
                transitionProb = parameters.deletionTransition*start;
            count+=start;
            start*=parameters.deletionTransitionFactor;
        }

        return (count==0)?0:transitionProb/(double)count;
    }

    bool isMatchState() const { return isMatch; }
    bool isLastInsertionState() const { return isLastInsertion; }
    unsigned int getIndex() const { return index;}
    char getNucleotide() const { return nuc;}
    unsigned int getCharIndex() const { return charIndex; }

private:
    unsigned int index; //what is the index of this state in the graph
    unsigned int charIndex; //which character it corresponds to in the sequencing read
    char nuc; //basepair in the charIndex
    bool isMatch; //is a match state?
    bool isLastInsertion;
    char nextNuc; //basepair in the next position
};

struct Read{
public:
    Read(const seqan::Dna5String& readSeq, unsigned int pos, const seqan::String<seqan::CigarElement<> >& cigar):
    pos(pos){
        this->read = readSeq;
        unsigned int curIndex = 0;
        for(unsigned int i = 0; i < length(cigar); ++i){
            char type = cigar[i].operation;
            if(type == 'S') seqan::erase(read, curIndex, curIndex+cigar[i].count);
            else if(type == 'M' || type == 'I') curIndex += cigar[i].count;
        }
    }
    
    seqan::String<char> read;
    unsigned int pos;
};

/*
 * Can be used as directed edge between from (outgoing edge source), and curState (incoming edge source).
 * operator< implemented so that std::set can evaulate this data structure
 */
struct TransitionInfoNode{

    int from, toState;
    TransitionInfoNode(int from, int toState):from(from), toState(toState){}

    bool operator==(const TransitionInfoNode& rhs) const{
        return from == rhs.from && toState == rhs.toState;
    }

    bool operator<(const TransitionInfoNode& rhs) const{
        return ((from < rhs.from) || (from == rhs.from && (toState < rhs.toState)));
    }
};

/*
 * Starting from startValues, until endValues, finds the indices that has the greatest values in values array and puts
 * in maxValues array. If there cannot be more than maxValuesSize different indices, rest of maxValues left as -1
 * Return parameter: maxValues array
 * Return value:
 */
template <typename T>
void findMaxValues(const T* values, bool* selectedIndices, const int startValues, const int endValues,
                     const size_t maxValuesSize){

    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int> >,
                        std::greater<std::pair<double, int> > > maxValues;

    for(int curState = 0; curState < endValues - startValues; ++curState) {
        if(maxValues.size() < maxValuesSize)
            maxValues.push(std::pair<double, int>(values[curState], curState));
        else if(maxValues.top().first < values[curState]){
            maxValues.pop();
            maxValues.push(std::pair<double, int>(values[curState], curState));
        }
    }

    while(!maxValues.empty()){
        selectedIndices[maxValues.top().second] = true;
        maxValues.pop();
    }
}

/*
 * Using the information provided with node and numberOfDeletions, inserts the transitions that should be made afterward
 * If you are going to change the transition structure, change it from there. These are imaginary edges in the graph.
 */
void insertNewForwardTransitions(std::vector<TransitionInfoNode>* transitionSet, const SeqNode& node,
                                 const int numberOfDeletions, const int maxInsertion){
    //next insertion
    if(!node.isLastInsertionState())
        transitionSet->push_back(TransitionInfoNode(node.getIndex(), node.getIndex()+1));

    //match and deletions
    for(int curOffset = 1; curOffset <= numberOfDeletions+1; ++curOffset)
        transitionSet->push_back(TransitionInfoNode(node.getIndex(), MATCH_OFFSET(node.getCharIndex(), curOffset,
                                                                                  maxInsertion)));
}

/*
 */
int fillForwardMatrix(SeqNode* graph, HMMParameters parameters, double* calculatedTransitionProbs,
                      double** forwardMatrix, const char* read, int const startPosition, int maxDistanceOnContig,
                      const int contigSize, int const readLength){

    if(startPosition < 0 || startPosition >= maxDistanceOnContig) return -1;
    if((int)GRAPH_SIZE(contigSize, parameters.maxInsertion) < maxDistanceOnContig)
        maxDistanceOnContig = GRAPH_SIZE(contigSize, parameters.maxInsertion);

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
    int curTime = 0; //represents the current time (1...T)
    insertNewForwardTransitions(curTrSet, graph[startPosition], parameters.maxDeletion, parameters.maxInsertion);
    for(size_t curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
        int matchoff = MATCH_OFFSET(graph[curTrSet->at(curTransition).from].getCharIndex(), 0, parameters.maxInsertion);
        //0->insertion, 1-> match, 2,3...->deletions
        int transitionIndex = (curTrSet->at(curTransition).toState - matchoff)/(parameters.maxInsertion+1);
        forwardMatrix[0][curTrSet->at(curTransition).toState - startPosition] += calculatedTransitionProbs[transitionIndex]*
        graph[curTrSet->at(curTransition).toState].getEmissionProb(read[curTime], parameters);

        insertNewForwardTransitions(nextTrSet, graph[curTrSet->at(curTransition).toState], parameters.maxDeletion,
                                    parameters.maxInsertion);
    }
    curTrSet->clear();

    //find the most likely states that should be allowed to make the next transitions
    findMaxValues(forwardMatrix[0], allowedParentStates, startPosition, maxDistanceOnContig, parameters.filterSize);

    tmpTrSet = curTrSet;
    curTrSet = nextTrSet;
    nextTrSet = tmpTrSet;

    //2-recursion (1 < t <= T)
    while (curTime < readLength-1 && !curTrSet->empty()){
        curTime++;
        for(size_t curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
            const TransitionInfoNode& frontTr = curTrSet->at(curTransition);
            if(allowedParentStates[frontTr.from - startPosition] && frontTr.toState < maxDistanceOnContig) {
                int matchoff = MATCH_OFFSET(graph[frontTr.from].getCharIndex(), 0, parameters.maxInsertion);
                //0->insertion, 1-> match, 2,3...->deletions
                int transitionIndex = (frontTr.toState - matchoff)/(parameters.maxInsertion+1);
                forwardMatrix[curTime][frontTr.toState-startPosition] +=
                    forwardMatrix[curTime-1][frontTr.from-startPosition]*
                    calculatedTransitionProbs[transitionIndex]*
                    graph[frontTr.toState].getEmissionProb(read[curTime],parameters);

                if(!hasStateBeenProcessedBefore[frontTr.toState - startPosition]){
                    insertNewForwardTransitions(nextTrSet, graph[frontTr.toState], parameters.maxDeletion,
                                                parameters.maxInsertion);
                    hasStateBeenProcessedBefore[frontTr.toState - startPosition] = true;
                }
            }
        }
        curTrSet->clear();

        std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);
        std::fill_n(allowedParentStates, maxDistanceOnContig-startPosition+1, false);

        //find the most likely states that should be allowed to make the next transitions
        findMaxValues(forwardMatrix[curTime], allowedParentStates, startPosition, maxDistanceOnContig,
                      parameters.filterSize);
        tmpTrSet = curTrSet;
        curTrSet = nextTrSet;
        nextTrSet = tmpTrSet;
    }

    int max = -1;
    for(int curStateForMax = 0; curStateForMax < maxDistanceOnContig-startPosition+1; ++curStateForMax){
        if(allowedParentStates[curStateForMax] &&
           (max == -1 || forwardMatrix[curTime][curStateForMax] > forwardMatrix[curTime][max])){
            max = curStateForMax;
        }
    }

    delete curTrSet; delete nextTrSet; delete[] allowedParentStates;

    return max + startPosition;
}

/*
 * Using the information provided with node and numberOfDeletions, inserts the transitions that should be made after
 * this node.
 * If you are going to change the transition structure, change it from there. These are imaginary edges in the graph.
 */
void insertNewBackwardTransitions(std::vector<TransitionInfoNode>* transitionSet, const SeqNode& node,
                                  const int numberOfDeletions, const int maxInsertion){

    if(node.getCharIndex() == 0) return;

    if(node.isMatchState()){
        for(int curBeforeOffset = 0; curBeforeOffset <= numberOfDeletions; ++curBeforeOffset){
            int offset = -1*curBeforeOffset - 1;
            transitionSet->push_back(TransitionInfoNode(node.getIndex(), MATCH_OFFSET(node.getCharIndex(), offset,
                                                                                      maxInsertion))); //deletion
            for(int curInsertion = 1; curInsertion <= maxInsertion; ++curInsertion)
                transitionSet->push_back(TransitionInfoNode(node.getIndex(),
                                                            INSERTION_OFFSET(node.getCharIndex(),offset, curInsertion,
                                                                             maxInsertion)));//match
        }
    }else{
        transitionSet->push_back(TransitionInfoNode(node.getIndex(), node.getIndex()-1));//match
    }
}

/*
 */
bool fillBackwardMatrix(SeqNode* graph, HMMParameters parameters, double* calculatedTransitionProbs,
                        double** backwardMatrix, const char* read, int startPosition, int maxDistanceOnContig,
                        const int contigSize, const int readLength){

    if(maxDistanceOnContig < 0 || maxDistanceOnContig >= startPosition) return false;
    if((int) GRAPH_SIZE(contigSize, parameters.maxInsertion) < startPosition)
        startPosition = GRAPH_SIZE(contigSize, parameters.maxInsertion);

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
    int curTime = readLength-1; //@curTime value is 0-based. So 0th index is the first character [T....1]
    insertNewBackwardTransitions(curTrSet, graph[startPosition], parameters.maxDeletion, parameters.maxInsertion);
    for(size_t curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
        if(curTrSet->at(curTransition).toState > maxDistanceOnContig){
            int matchoff = MATCH_OFFSET(graph[curTrSet->at(curTransition).toState].getCharIndex(), 0,
                                        parameters.maxInsertion);
            //0->insertion, 1-> match, 2,3...->deletions
            int transitionIndex = (curTrSet->at(curTransition).from - matchoff)/(parameters.maxInsertion+1);
            backwardMatrix[curTime][curTrSet->at(curTransition).toState - maxDistanceOnContig] +=
            calculatedTransitionProbs[transitionIndex];

            insertNewBackwardTransitions(nextTrSet, graph[curTrSet->at(curTransition).toState], parameters.maxDeletion,
                                         parameters.maxInsertion);
        }
    }
    curTrSet->clear();

    findMaxValues(backwardMatrix[curTime], allowedParentStates, maxDistanceOnContig, startPosition,
                  parameters.filterSize);
    tmpTrSet = curTrSet;
    curTrSet = nextTrSet;
    nextTrSet = tmpTrSet;

    //2-recursion
    while (curTime > 0 && !curTrSet->empty()){
        curTime--;
        for(size_t curTransition = 0; curTransition < curTrSet->size(); ++curTransition){
            const TransitionInfoNode& frontTr = curTrSet->at(curTransition);
            if(allowedParentStates[frontTr.from - maxDistanceOnContig] && frontTr.toState > maxDistanceOnContig) {
                int matchoff = MATCH_OFFSET(graph[frontTr.toState].getCharIndex(), 0, parameters.maxInsertion);
                //0->insertion, 1-> match, 2,3...->deletions
                int transitionIndex = (frontTr.from - matchoff)/(parameters.maxInsertion+1);
                backwardMatrix[curTime][frontTr.toState - maxDistanceOnContig] +=
                    backwardMatrix[curTime+1][frontTr.from - maxDistanceOnContig]*
                    calculatedTransitionProbs[transitionIndex]*
                    graph[frontTr.from].getEmissionProb(read[curTime+1], parameters);

                if(!hasStateBeenProcessedBefore[frontTr.toState - maxDistanceOnContig]){
                    insertNewBackwardTransitions(nextTrSet, graph[frontTr.toState], parameters.maxDeletion,
                                                 parameters.maxInsertion);
                    hasStateBeenProcessedBefore[frontTr.toState - maxDistanceOnContig] = true;
                }
            }
        }
        curTrSet->clear();

        std::fill_n(hasStateBeenProcessedBefore.begin(), hasStateBeenProcessedBefore.size(), false);

        std::fill_n(allowedParentStates, startPosition - maxDistanceOnContig + 1, false);

        findMaxValues(backwardMatrix[curTime], allowedParentStates, maxDistanceOnContig, startPosition,
                      parameters.filterSize);
        tmpTrSet = curTrSet;
        curTrSet = nextTrSet;
        nextTrSet = tmpTrSet;
    }

    delete curTrSet; delete nextTrSet; delete[] allowedParentStates;

    return true;
}

//offset is an index for state, not a basepair index
void calculateViterbiPool(SeqNode* graph, HMMParameters parameters, double** transitionProbs,
                          std::pair<double, char>* emissionProbs, unsigned int seqLength, unsigned int numberOfStates,
                          unsigned int& seqIndex, seqan::String<seqan::Dna5String>& sequence){

    unsigned int curSequence = 0; //0-based. refers to the current piece of sequence
    {//block for obtaining new index for the next contig chunk
        std::lock_guard<std::mutex> lk(indexMutex);
        curSequence = seqIndex++;
    }

    while(curSequence < length(sequence)){

        std::vector<std::map<int, int> > backtrace; //index
        unsigned int startState = MATCH_OFFSET(curSequence*seqLength, 0, parameters.maxInsertion);
        unsigned int endState = MATCH_OFFSET(((curSequence+1)* seqLength) + 1, 0, parameters.maxInsertion);

        if(endState > numberOfStates-1) endState = numberOfStates-1;
        unsigned int totalStates = endState - startState + 1;

        std::vector<bool> allowedParentStates; //index - offset
        allowedParentStates.resize(totalStates);
        std::fill_n(allowedParentStates.begin(), allowedParentStates.size(), false);

        std::set<int>* curTrSet = new std::set<int>(); //index
        std::set<int>* nextTrSet = new std::set<int>(); //index
        std::set<int>* tmpTrSet; //index

        std::map<int, double>* prevViterbi = new std::map<int, double>(); //index
        std::map<int, double>* curViterbi = new std::map<int, double>(); //index
        std::map<int, double>* tmpViterbi; //index
        std::map<int, double>::iterator itVit;

        //initialization step t = 1 (curTime = 0)
        std::vector<TransitionInfoNode> curStateTransitions; //pushing real index
        insertNewForwardTransitions(&curStateTransitions, graph[startState], parameters.maxDeletion, parameters.maxInsertion);
        backtrace.push_back(std::map<int, int >());

        unsigned int curTime = 0;
        prevViterbi->clear();
        while(!curStateTransitions.empty()){

            int next = curStateTransitions.begin()->toState;
            int transitionIndex = (next - startState)/(parameters.maxInsertion+1); //0->insertion, 1-> match, 2,3...->deletions

            if(next <= (int) endState && transitionProbs[startState][transitionIndex] > 0.001){
                (*prevViterbi)[next] = log10(transitionProbs[startState][transitionIndex]) + log10(emissionProbs[next].first);

                curTrSet->insert(next);
                allowedParentStates[next-startState] = true;
                backtrace[curTime][next] = curStateTransitions.begin()->from; //offset
            }
            curStateTransitions.erase(curStateTransitions.begin());
        }
        curStateTransitions.clear();
        curTime++;

        bool finished = false;
        int maxFinalStateTime = 0;
        double maxFinalStateValue = INT_MIN;
        while(!finished && curTime < (parameters.maxInsertion+1)*seqLength){

            //new storage for the current time
            curViterbi->clear();
            backtrace.push_back(std::map<int, int>());

            for(std::set<int>::iterator itSet = (*curTrSet).begin(); itSet != curTrSet->end(); ++itSet){
                //calculate the viterbi values from this state to the states that it has transitions
                int fromState = (*itSet);
                //if this state is one of the major ones from the previous time
                if(allowedParentStates[fromState-startState]){

                    insertNewForwardTransitions(&curStateTransitions, graph[fromState], parameters.maxDeletion,
                                                parameters.maxInsertion);
                    for(size_t toStateIndex = 0; toStateIndex < curStateTransitions.size(); ++toStateIndex){

                        int toState = curStateTransitions[toStateIndex].toState;
                        int matchoff = MATCH_OFFSET(graph[fromState].getCharIndex(), 0, parameters.maxInsertion);
                        //0->insertion, 1-> match, 2,3...->deletions
                        int transitionIndex = (toState - matchoff)/(parameters.maxInsertion+1);

                        if(toState <= (int) endState && transitionProbs[fromState][transitionIndex] > 0.001){

                            double newViterbi = (*prevViterbi)[fromState] +
                            log10(transitionProbs[fromState][transitionIndex]) + log10(emissionProbs[toState].first);

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
            std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int> >,
            std::greater<std::pair<double, int> > > maxValues;
            for(std::set<int>::iterator itSet = (*nextTrSet).begin(); itSet != (*nextTrSet).end(); ++itSet) {
                //!!! changed it was maxValues.size() < parameters.filterSize
                if(maxValues.size() < parameters.viterbiFilterSize)
                    maxValues.push(std::pair<double, int>((*curViterbi)[*itSet], *itSet));
                else if(maxValues.top().first < (*curViterbi)[*itSet]){
                    maxValues.pop();
                    maxValues.push(std::pair<double, int>((*curViterbi)[*itSet], *itSet));
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
            if(maxFinalStateTime > 0 && (int) curTime > maxFinalStateTime + 20) finished = true;
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
            resize(sequence[curSequence], maxFinalStateTime);
            int lastIndexUsed = endState;
            for(int curChar = maxFinalStateTime-1; curChar >= 0; --curChar){
                int curIndex = backtrace[curChar+1][lastIndexUsed];
                sequence[curSequence][curChar] = emissionProbs[curIndex].second;
                lastIndexUsed = curIndex;
            }
        }

        {//block for obtaining new index for the next contig chunk
            std::lock_guard<std::mutex> lk(indexMutex);
            curSequence = seqIndex++;
        }
    }
}

void backtraceWithViterbi(SeqNode* graph, HMMParameters parameters, double** transitionProbs,
                         std::pair<double, char>* emissionProbs, unsigned int numberOfStates, unsigned int thread,
                         unsigned int contigSize, seqan::String<seqan::Dna5String>& sequence){

    if(contigSize > 0){
        std::vector<std::thread> threads;
        unsigned int seqIndex = 0;
        unsigned int batchSize = (parameters.batchSize == 0)?contigSize:parameters.batchSize;
        resize(sequence, ((contigSize-1)/batchSize) + 1);

        for(unsigned int i = 0; i < thread; ++i){
            threads.push_back(std::thread(calculateViterbiPool, graph, parameters, transitionProbs, emissionProbs,
                                          batchSize, numberOfStates, std::ref(seqIndex), std::ref(sequence)));
        }
        for(unsigned int i = 0; i < thread; ++i) threads[i].join();
    }
}

void calculateFBPool(HMMParameters parameters, const std::vector<Read>& reads, unsigned int assemblySize, SeqNode* sequencingGraph, double** transitionProbs, double** emissionProbs, int* stateProcessedCount, int** transitionProcessedCount, unsigned int& readIndex){
    
    std::vector<TransitionInfoNode> curStateTransitions;
    std::vector<TransitionInfoNode> stateTransitions; //all possible transitions from a state to another
    insertNewForwardTransitions(&stateTransitions, sequencingGraph[0], parameters.maxDeletion, parameters.maxInsertion);
    unsigned int numOfTransitions = (unsigned int) stateTransitions.size();
    //@dynamic init array of pre calculated initial transition values for any state in hmm
    double* calculatedTransitionProbs = new double[numOfTransitions];
    for(unsigned int curTr = 0; curTr < numOfTransitions; ++curTr)
        calculatedTransitionProbs[curTr]=sequencingGraph[0].transitionProbFromThisNode(stateTransitions[curTr].toState,
                                                                                       parameters);

    //@dynamic init
    double* curStateTransitionLikelihood = new double[numOfTransitions];
    //@dynamic init
    double* curStateEmissionProbs = new double[totalNuc];
    
    unsigned int curRead = 0; //0-based. refers to the current read id
    {//block for obtaining the next aligned read
        std::lock_guard<std::mutex> lk(indexMutex);
        curRead = readIndex++;
    }
    
    while(curRead < reads.size()){
        
        //to which character should correction extend at maximum
        unsigned int readLength = (unsigned int) seqan::length(reads[curRead].read);
        unsigned int maxTransition = readLength + readLength/20 + 1;
        if(readLength < 500) maxTransition = readLength + readLength/3 + 1;

        unsigned int maxDistanceOnAssembly = std::min(reads[curRead].pos+maxTransition, assemblySize);
        //states prior to this wont be processed. offset value is to ignore these states
        unsigned int offset = MATCH_OFFSET(reads[curRead].pos, -1, parameters.maxInsertion);
        //maximum number of states to be processed
        unsigned int fbMatrixSize = MATCH_OFFSET(maxDistanceOnAssembly, 1, parameters.maxInsertion) - offset;

        if(fbMatrixSize > 0){
            int j; //j value in i,j transitions
            //@dynamic init
            double** forwardMatrix = new double*[readLength]; //forwardMatrix[time][state]
            for(unsigned int i = 0; i < readLength; i++){
                //@dynamic init
                forwardMatrix[i] = new double[fbMatrixSize];
                std::fill_n(forwardMatrix[i], fbMatrixSize, 0);
            }
            //@dynamic init
            double** backwardMatrix = new double*[readLength]; //backwardMatrix[time][state]
            for(unsigned int i = 0; i < readLength; i++){
                //@dynamic init
                backwardMatrix[i] = new double[fbMatrixSize];
                std::fill_n(backwardMatrix[i], fbMatrixSize, 0);
            }
            int startForBackward = fillForwardMatrix(sequencingGraph, parameters, calculatedTransitionProbs,
                                                     forwardMatrix, toCString(reads[curRead].read), offset,
                                                     MATCH_OFFSET(maxDistanceOnAssembly, 1, parameters.maxInsertion),
                                                     assemblySize, readLength);

            if(startForBackward != -1 && sequencingGraph[startForBackward].getCharIndex() > reads[curRead].pos){
                startForBackward = MATCH_OFFSET(sequencingGraph[startForBackward].getCharIndex(), 1,
                                                parameters.maxInsertion);

                if(fillBackwardMatrix(sequencingGraph, parameters, calculatedTransitionProbs, backwardMatrix,
                                      toCString(reads[curRead].read), startForBackward, offset,
                                      assemblySize, readLength)){
                    //updating probabilities wrt the f/b matrices computed just now
                    for(int curState = INSERTION_OFFSET(reads[curRead].pos, -1, 1, parameters.maxInsertion);
                        curState < startForBackward; ++curState){

                        if(sequencingGraph[curState].isLastInsertionState())
                            //for the last insertion state, the insertion probs change
                            calculatedTransitionProbs[1] = parameters.matchTransition + parameters.insertionTransition;

                        if(curState-offset < fbMatrixSize){
                            int matchoff = MATCH_OFFSET(sequencingGraph[curState].getCharIndex(), 0,
                                                        parameters.maxInsertion);
                            std::fill_n(curStateTransitionLikelihood, numOfTransitions, 0.0);
                            std::fill_n(curStateEmissionProbs, totalNuc, 0.0);
                            insertNewForwardTransitions(&curStateTransitions, sequencingGraph[curState],
                                                        parameters.maxDeletion, parameters.maxInsertion);

                            for(unsigned int t = 0; t < readLength; ++t){
                                //transition probabilities
                                if(t < readLength-1){
                                    for(size_t curTr = 0; curTr < curStateTransitions.size(); ++curTr){
                                        if(curStateTransitions.at(curTr).toState - offset < fbMatrixSize){
                                            j = curStateTransitions[curTr].toState;
                                            //0->insertion, 1-> match, 2,3...->deletions
                                            int transitionIndex = (j - matchoff)/(parameters.maxInsertion+1);
                                            curStateTransitionLikelihood[transitionIndex] +=
                                            forwardMatrix[t][curState-offset]*
                                            calculatedTransitionProbs[transitionIndex]*
                                            sequencingGraph[j].getEmissionProb(reads[curRead].read[t+1],
                                                                               parameters)*backwardMatrix[t+1][j-offset];

                                        }
                                    }
                                }

                                //emission probabilities
                                char emitChar = (reads[curRead].read[t] != 'N')?
                                reads[curRead].read[t]:
                                (sequencingGraph[curState].isMatchState())?
                                sequencingGraph[curState].getNucleotide():'\0';
                                Nucleotide chosenNuc = (emitChar == 'A' || emitChar == 'a')?A:
                                (emitChar == 'T' || emitChar == 't')?T:
                                (emitChar == 'G' || emitChar == 'g')?G:
                                (emitChar == 'C' || emitChar == 'c')?C:totalNuc;
                                if(chosenNuc < totalNuc)
                                    curStateEmissionProbs[chosenNuc] +=
                                    forwardMatrix[t][curState-offset]*backwardMatrix[t][curState-offset];
                            }
                            curStateTransitions.clear();
                            double totalEmissionProbs = curStateEmissionProbs[A] + curStateEmissionProbs[T] +
                            curStateEmissionProbs[G] + curStateEmissionProbs[C];
                            double totalTransitionLikelihoods = 0;
                            double processedTransitionProb = 0;
                            for(unsigned int i = (sequencingGraph[curState].isLastInsertionState())?1:0;
                                i < numOfTransitions; ++i){
                                if(curStateTransitionLikelihood[i] > 0 || i == 0){
                                    totalTransitionLikelihoods += curStateTransitionLikelihood[i];
                                    processedTransitionProb += calculatedTransitionProbs[i];
                                }
                            }

                            if(totalEmissionProbs != 0 && curState < startForBackward){
                                if(totalTransitionLikelihoods != 0){
                                    {//block for updating the transition probs and the transition proccessed count
                                        std::lock_guard<std::mutex> lk(transitionProbMutex);

                                        for(unsigned int i = (sequencingGraph[curState].isLastInsertionState())?1:0;
                                            i < numOfTransitions; ++i){
                                            if(curStateTransitionLikelihood[i] > 0 || i == 0){
                                                transitionProbs[curState][i] +=
                                                (curStateTransitionLikelihood[i]/totalTransitionLikelihoods)*
                                                processedTransitionProb;
                                                transitionProcessedCount[curState][i]++;
                                            }
                                        }
                                    }
                                }

                                {//block for updating the emission probs and the state processed count
                                    std::lock_guard<std::mutex> lk(emissionProbMutex);

                                    emissionProbs[curState][A] += curStateEmissionProbs[A]/totalEmissionProbs;
                                    emissionProbs[curState][T] += curStateEmissionProbs[T]/totalEmissionProbs;
                                    emissionProbs[curState][G] += curStateEmissionProbs[G]/totalEmissionProbs;
                                    emissionProbs[curState][C] += curStateEmissionProbs[C]/totalEmissionProbs;

                                    stateProcessedCount[curState]++;
                                }
                            }
                        }

                        //for the last insertion state, the insertion probs change so that it wont have insertion
                        //transition. putting it back to normal now
                        if(sequencingGraph[curState].isLastInsertionState())
                            calculatedTransitionProbs[1] = parameters.matchTransition;
                    }
                }
            }

            for(unsigned int i = 0; i < readLength; ++i) { delete[] backwardMatrix[i]; delete[] forwardMatrix[i];}
            delete[] backwardMatrix;
            delete[] forwardMatrix;
        }
        
        {//block for obtaining the next aligned read
            std::lock_guard<std::mutex> lk(indexMutex);
            curRead = readIndex++;
        }
    }

    delete[] calculatedTransitionProbs;
    delete[] curStateTransitionLikelihood;
    delete[] curStateEmissionProbs;
}

bool fillBuffer(seqan::BamFileIn& alignmentFileIn, const seqan::FaiIndex& readsIndex, seqan::BamAlignmentRecord& record,
                std::vector<Read>& reads, const uint64_t& contigId, unsigned int mapQ, unsigned int size){
    
    reads.clear();
    while((int)contigId == record.rID && reads.size() < size){
        
        if(!hasFlagUnmapped(record) && record.mapQ >= mapQ){
            seqan::Dna5String alignedString;
            int qId;
            if(getIdByName(qId, readsIndex, record.qName)){
                readSequence(alignedString, readsIndex, qId);
                if(hasFlagRC(record)) reverseComplement(alignedString);
                if(length(alignedString) > 0 && record.beginPos >= 0)
                    reads.push_back(Read(alignedString, record.beginPos+1, record.cigar));
            }
        }
        
        if(!atEnd(alignmentFileIn)){readRecord(record, alignmentFileIn);}
        else{break;}
    }
    
    return (contigId == record.rID)&(!atEnd(alignmentFileIn));
}

/** @brief Polishes given contig by applying first Forward/Backward calculationg and than Viterbi calculation
 *
 *  @param parameters Parameters needed for HMM graphs. User provides most of the parameters for it.
 *  @param alignmentSetsIn List of bam/sam files
 *  @param curRecords List of alignment record pointers that points to the current unprocessed record in @alignmentSetsIn
 *  @param readFAIs List of FAI indices
 *  @param contig contig sequence to polish
 *  @param contigId Id of the contig. This must be the same in each of the alignment - read files pairs
 *  @param correctedContig An output of this function. If it can polish a contig, the polished contig is stored here
 *  @param mapQ Threshold for the minimum mapping quality to use in polishing
 *  @param thread Number of threads that can be used
 *  @return correctedContig See above.
 */
void polishContig(HMMParameters parameters, std::vector<seqan::BamFileIn>& alignmentSetsIn, std::vector<seqan::BamAlignmentRecord>& curRecords, const std::vector<seqan::FaiIndex>& readFAIs, const seqan::Dna5String& contig, const uint64_t& contigId, seqan::Dna5String& correctedContig, unsigned int mapQ, unsigned int thread){
    
    unsigned int contigLength = (unsigned int)length(contig);

    //@dynamic init. array
    unsigned int hmmGraphSize = GRAPH_SIZE(contigLength, parameters.maxInsertion);
    SeqNode* sequencingGraph = new SeqNode[hmmGraphSize];

    //constructing the pHMM graph
    sequencingGraph[0] = SeqNode(0, 0, parameters.maxInsertion, '\0', contig[0]);
    for(int curIn = 1; curIn <= (int)parameters.maxInsertion; ++curIn)
        sequencingGraph[curIn] = SeqNode(curIn, 0, parameters.maxInsertion, '\0', contig[0]);

    unsigned int matchOffset;
    //states in between the start and end states
    for(unsigned int curCharacter = 1; curCharacter < contigLength; ++curCharacter){
        matchOffset = MATCH_OFFSET(curCharacter, 0, parameters.maxInsertion);
        sequencingGraph[matchOffset] = SeqNode(matchOffset, curCharacter, parameters.maxInsertion,
                                                      contig[curCharacter-1], contig[curCharacter]);
        for(unsigned int curIn = 1; curIn <= parameters.maxInsertion; ++curIn)
            sequencingGraph[matchOffset+curIn] = SeqNode(matchOffset+curIn, curCharacter,parameters.maxInsertion,
                                                                contig[curCharacter-1], contig[curCharacter]);
    }
    //last character before the end state
    matchOffset = MATCH_OFFSET(contigLength, 0, parameters.maxInsertion);
    sequencingGraph[matchOffset] = SeqNode(matchOffset, contigLength, parameters.maxInsertion,
                                                  contig[contigLength-1], '\0');
    for(unsigned int curIn = 1; curIn <= parameters.maxInsertion; ++curIn)
        sequencingGraph[matchOffset+curIn] = SeqNode(matchOffset+curIn, contigLength, parameters.maxInsertion,
                                                              contig[contigLength-1], '\0');

    //end state
    matchOffset = END_STATE(contigLength, parameters.maxInsertion);
    sequencingGraph[matchOffset] = SeqNode(matchOffset, contigLength+1, parameters.maxInsertion, '\0', '\0');
    std::vector<TransitionInfoNode> stateTransitions; //all possible transitions from a state to another
    insertNewForwardTransitions(&stateTransitions, sequencingGraph[0], parameters.maxDeletion, parameters.maxInsertion);
    unsigned int numOfTransitions = (unsigned int) stateTransitions.size();

    //pre calculated initial transition values for any state in hmm
    //@dynamic init. array
    double* calculatedTransitionProbs = new double[numOfTransitions];
    for(unsigned int curTr = 0; curTr < numOfTransitions; ++curTr)
        calculatedTransitionProbs[curTr]=sequencingGraph[0].transitionProbFromThisNode(stateTransitions[curTr].toState,
                                                                                        parameters);
    //@dynamic init. array of transitionProbs[from state][to state]
    double** transitionProbs = new double*[hmmGraphSize];
    //@dynamic init. array of emissionProbs[state][character index]
    double** emissionProbs = new double*[hmmGraphSize];
    //@dynamic init: how many times a state has been processed [state]
    int* stateProcessedCount = new int[hmmGraphSize];
    std::fill_n(stateProcessedCount, hmmGraphSize, 0);
    //@dynamic init: how many times a transition has been processed [state][transition]
    int** transitionProcessedCount = new int*[hmmGraphSize];
    for(unsigned int curState = 0; curState < hmmGraphSize; ++curState){
        //@dynamic init
        transitionProcessedCount[curState] = new int[numOfTransitions];
        std::fill_n(transitionProcessedCount[curState], numOfTransitions, 0);
        //@dynamic init
        transitionProbs[curState] = new double[numOfTransitions];
        //@dynamic init
        emissionProbs[curState] = new double[totalNuc];
        std::fill_n(transitionProbs[curState], numOfTransitions, 0.0);
        std::fill_n(emissionProbs[curState], totalNuc, 0.0);
    }

    for(unsigned int curSet = 0; curSet < alignmentSetsIn.size(); ++curSet){
        if(atEnd(alignmentSetsIn[curSet]) || curRecords[curSet].rID != contigId) continue;
        
        bool shouldPolish = true;
        while(shouldPolish){
            unsigned int readIndex = 0; //for threads
            std::vector<std::thread> threads;
            std::vector<Read> reads;
            
            shouldPolish = fillBuffer(alignmentSetsIn[curSet], readFAIs[curSet], curRecords[curSet], reads, contigId, mapQ, 50000);
            
            //alignedReads.size will be 0 if there is no more alignment record left to read in the current file
            for(unsigned int i = 0; i < thread && i < reads.size(); ++i){
                threads.push_back(std::thread(calculateFBPool, parameters, std::ref(reads), contigLength, sequencingGraph, transitionProbs, emissionProbs, stateProcessedCount, transitionProcessedCount, std::ref(readIndex)));
            }

            //buffer is cleared here. every thread needs to wait before the buffer gets reloaded again
            for(unsigned int i = 0; i < threads.size(); ++i) threads[i].join();
        }
    }
    
    for(unsigned int curState = 0; curState < hmmGraphSize; ++curState){
        if(sequencingGraph[curState].isLastInsertionState())
            //for the last insertion state, the insertion probs change so that it wont have insertion transition
            calculatedTransitionProbs[1] = parameters.matchTransition + parameters.insertionTransition;

        if(stateProcessedCount[curState] > 0){ //if this state ever processed then its probs may need to be updated
            for(unsigned int curTransition = (sequencingGraph[curState].isLastInsertionState())?1:0;
                curTransition < numOfTransitions; ++curTransition){

                transitionProbs[curState][curTransition] = (transitionProcessedCount[curState][curTransition] > 0)?
                transitionProbs[curState][curTransition]/transitionProcessedCount[curState][curTransition]:
                calculatedTransitionProbs[curTransition];
            }

            emissionProbs[curState][A] /= stateProcessedCount[curState];
            emissionProbs[curState][T] /= stateProcessedCount[curState];
            emissionProbs[curState][G] /= stateProcessedCount[curState];
            emissionProbs[curState][C] /= stateProcessedCount[curState];
        }else{ //initial probs to be set unless this state has been processed
            for(unsigned int curTr = (sequencingGraph[curState].isLastInsertionState())?1:0; curTr < numOfTransitions;
                ++curTr)
                transitionProbs[curState][curTr] = calculatedTransitionProbs[curTr];

            emissionProbs[curState][A] = sequencingGraph[curState].getEmissionProb('A', parameters);
            emissionProbs[curState][T] = sequencingGraph[curState].getEmissionProb('T', parameters);
            emissionProbs[curState][G] = sequencingGraph[curState].getEmissionProb('G', parameters);
            emissionProbs[curState][C] = sequencingGraph[curState].getEmissionProb('C', parameters);
        }

        if(sequencingGraph[curState].isLastInsertionState())
            calculatedTransitionProbs[1] = parameters.matchTransition;
    }

    delete[] stateProcessedCount;

    //@dynamic init
    std::pair<double, char>* maxEmissionProbs = new std::pair<double, char>[hmmGraphSize];
    for(unsigned int curState = 0; curState < hmmGraphSize; ++curState){
        if(emissionProbs[curState][A] >= std::max(emissionProbs[curState][T],
                                                  std::max(emissionProbs[curState][G], emissionProbs[curState][C]))){
            maxEmissionProbs[curState] = std::make_pair(emissionProbs[curState][A], 'A');
        }else if(emissionProbs[curState][T] >= std::max(emissionProbs[curState][A],
                                                        std::max(emissionProbs[curState][G], emissionProbs[curState][C]))){
            maxEmissionProbs[curState] = std::make_pair(emissionProbs[curState][T], 'T');
        }else if(emissionProbs[curState][G] >= std::max(emissionProbs[curState][A],
                                                        std::max(emissionProbs[curState][T], emissionProbs[curState][C]))){
            maxEmissionProbs[curState] = std::make_pair(emissionProbs[curState][G], 'G');
        }else{
            maxEmissionProbs[curState] = std::make_pair(emissionProbs[curState][C], 'C');
        }
    }

    seqan::String<seqan::Dna5String> sequence;
    backtraceWithViterbi(sequencingGraph, parameters, transitionProbs, maxEmissionProbs, hmmGraphSize, thread,
                         contigLength, sequence);

    for(unsigned int i = 0; i < hmmGraphSize; ++i){
        delete[] transitionProbs[i]; delete[] emissionProbs[i]; delete[] transitionProcessedCount[i];
    }
    delete[] transitionProbs; delete[] emissionProbs; delete[] transitionProcessedCount; delete[] sequencingGraph;
    delete[] maxEmissionProbs; delete[] calculatedTransitionProbs;

    for(unsigned i = 0; i < length(sequence); ++i)
        append(correctedContig, sequence[i]);
}

/** @brief Polishes each contig in @assemblyFile reading the alignments from @alignmentFile
 *
 *  @param parameters Parameters needed for HMM graphs. User provides most of the parameters for it.
 *  @param assemblyFile Fasta/Fastq file that includes sequencing reads for each contig in assembly
 *  @param readSets List of Fasta/Fastq files that include reads that are used to align the reads in @assemblyFile
 *  @param alignmentSets List of SAM/BAM files that include alignments of the reads in @readsFile to reads in @assemblyFile
 *  @param outputFile Output file to write the both corrected/uncorrected reads (contigs).
 *  @param mapQ Threshold for mapping quality of an alignment.
 *  @param thread Number of threads that can be used
 *  @param shouldQuite If true, there is no standard output
 *  @return True, if there was no problem at any step. Otherwise, false with an appropriate output of the error.
 */
bool polish(HMMParameters parameters, seqan::String<char> assemblyFile, std::vector<seqan::String<char> > readSets,
            std::vector<seqan::String<char> > alignmentSets, seqan::String<char> outputFile, unsigned int mapQ,
            unsigned int thread, bool shouldQuite){

    //Index file for assembly fasta
    seqan::FaiIndex assemblyFAI;
    if(!build(assemblyFAI, toCString(assemblyFile))){ //read the assembly file and build the index file
        std::cerr << "ERROR: Could not build FAI index for file " << assemblyFile << ". Make sure the fasta file "
        << "is structured correctly for creating its index (i.e. all lines have to have same number of basepairs)."
        << std::endl;
        return false;
    }
    //how many contigs in the assembly
    uint64_t nContigs = seqan::numSeqs(assemblyFAI);
    
    //reading/creating the indices for each of the alignment and the read set pairs
    std::vector<seqan::BamAlignmentRecord> curRecords(alignmentSets.size());
    std::vector<seqan::BamFileIn> alignmentSetsIn(alignmentSets.size());
    std::vector<seqan::BamIndex<seqan::Bai> > baiIndices(alignmentSets.size()); //indexed alignment files
    std::vector<seqan::FaiIndex> readFAIs(alignmentSets.size()); //indexed read sets
    for(unsigned int i = 0; i < alignmentSets.size(); ++i){

        if (!open(alignmentSetsIn[i], toCString(alignmentSets[i]))){
            std::cerr << "ERROR: Could not open " << alignmentSets[i] << std::endl;
            return false;
        }
        try{
            seqan::BamHeader header;
            readHeader(header, alignmentSetsIn[i]);
            if(!atEnd(alignmentSetsIn[i])){
                readRecord(curRecords[i], alignmentSetsIn[i]);
            }
        }catch(seqan::Exception const & e){
            std::cerr << "ERROR: " << e.what() << std::endl;
            return false;
        }

        //fetching/building&saving FAI files for each of the read sets. UPDATE: Just building currently, not fetching
//        if(!seqan::open(readFAIs[i], toCString(readSets[i]))){
            if(!build(readFAIs[i], toCString(readSets[i]))){
                std::cerr << "ERROR: Could not build FAI index for file " << readSets[i] << ". Make sure the fasta "
                << "file is structured correctly for creating its index (i.e., all lines have the same number of "
                << "basepairs)." << std::endl; return false;
            }
//            if (!save(readFAIs[i])){
//                std::cerr << "ERROR: Could not write FAI index for file " << readSets[i] << " to the disk. Apollo "
//                <<"will still continue running.\n";
//            }
//        }
    }

    //output file for the polished/unpolished contigs (all of them in a single data set with the order preserved)
    std::fstream correctedReadStream;
    try{
        correctedReadStream.open(toCString(outputFile), std::fstream::out | std::fstream::trunc);
    }catch(std::ios_base::failure e){
        std::cerr << "Could not open " << outputFile << std::endl;
        return false;
    }

    if(!shouldQuite) std::cout << "Polishing has begun..." << std::endl;
    seqan::SeqFileOut correctedReadsOut(correctedReadStream, seqan::Fasta());
    
    //for each contig, read all the aligned reads, save them, and correct the contig if there is at least one
    //alignment record
    uint64_t contigId = 0;
    //polish each of the contig
    while(contigId < nContigs){
        
        //name of the current contig to polish
        seqan::CharString curContigName = (seqan::CharString)seqan::sequenceName(assemblyFAI, (unsigned int) contigId);
        seqan::Dna5String curSeq; //assembly contig
        seqan::readSequence(curSeq, assemblyFAI, (unsigned int)contigId);
        seqan::Dna5String corr;
        
        //polish if there is at least one alignment in any of the read-to-assembly sets
        bool shouldPolish = false;
        for(int curSet = 0; curSet < curRecords.size() && !shouldPolish; ++curSet){
            shouldPolish |= (curRecords[curSet].rID == contigId)&(!atEnd(alignmentSetsIn[curSet]));
        }
        
        //start polishing
        if(shouldPolish){
            polishContig(parameters, alignmentSetsIn, curRecords, readFAIs, curSeq, contigId, corr, mapQ, thread);
            if(length(corr) > 0) curSeq = corr;
        }else{
          std::cerr << "The contig with id " << toCString(curContigName) << " could not be polished because there is no read aligning to it. Original (i.e., unpolished) sequence will be reported." << std::endl;
        }
        
        writeRecord(correctedReadsOut, curContigName, curSeq);
        ++contigId;
    }

    correctedReadStream.close();
    if(!shouldQuite) std::cout << std::endl << "Results have been written under " << outputFile << std::endl;

    return true;
}

int main(int argc, const char** argv){

    CommandLineParser options;
    seqan::ArgumentParser::ParseResult parseRes = parseCommandOptions(options, argc, argv);

    if(parseRes != seqan::ArgumentParser::PARSE_OK){
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }

    HMMParameters parameters(options.filterSize, options.viterbiFilterSize, options.maxDeletion, options.maxInsertion,
                             options.batchSize, options.matchTransition, options.insertionTransition,
                             options.deletionTransitionFactor, options.matchEmission);

    if(!options.shouldQuite){
        std::cout << "Assembly: " << toCString(options.assembly) << std::endl <<
        "Pair of a set of reads and their alignments:" << std::endl;

        for(unsigned int i = 0; i < options.readSets.size(); ++i)
            std::cout << toCString(options.readSets.at(i)) << ", " << toCString(options.alignmentSets.at(i)) << std::endl;

        std::cout << "Output file: " << options.output << std::endl << "Min mapping quality: " << options.mapQ <<
        std::endl << "Filter size: " << parameters.filterSize << std::endl << "Viterbi filter size: " <<
        parameters.viterbiFilterSize << std::endl << "Viterbi batch size: " << parameters.batchSize << std::endl <<
        "Maximum insertion: " << parameters.maxInsertion << std::endl << "Maximum deletion: " << parameters.maxDeletion
        << std::endl << "Match transition probability: " << parameters.matchTransition << std::endl <<
        "Insertion transition probability: " << parameters.insertionTransition << std::endl <<
        "Deletion transition probability: " << parameters.deletionTransition << std::endl <<
        "Deletion transition factor: " << parameters.deletionTransitionFactor << std::endl <<
        "Match emission probability: " << parameters.matchEmission << std::endl << "Mismatch emission probability: " <<
        parameters.mismatchEmission << std::endl << "Insertion emission probability: " << parameters.insertionEmission
        << std::endl << "Max thread: " << options.maxThread << std::endl;
    }

    polish(parameters, options.assembly, options.readSets, options.alignmentSets, options.output, options.mapQ,
           options.maxThread, options.shouldQuite);

    return 0;
}
