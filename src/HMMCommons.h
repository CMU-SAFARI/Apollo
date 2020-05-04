/** @file HMMCommons.h
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/

#ifndef HMMCommons_h
#define HMMCommons_h

#include <stdio.h>
#include <seqan/bam_io.h>
#include <seqan/seq_io.h>

#define MATCH_OFFSET(char, offset, insize) ((char+offset)*(insize+1))
#define INSERTION_OFFSET(char, offset, innumber, insize) ((char+offset)*(insize+1) + innumber) //innumber 1 based
#define MATCH_DOFFSET(char, offset, insize) ((char-offset)*(insize+1))
#define INSERTION_DOFFSET(char, offset, innumber, insize) ((char-offset)*(insize+1) + innumber) //innumber 1 based
#define GRAPH_SIZE(size, insize) ((size+1)*(insize+1) + 1)
#define END_STATE(size, insize) (GRAPH_SIZE(size, insize)-1)

enum Nucleotide {A, T, G, C, totalNuc};

typedef double prob_prec; //precision of the probabilities
typedef int cnt_prec; //precision of the count arrays -- relative to the avg. depth of coverage
typedef unsigned ind_prec; //precision of the size of contigs -- relative to genome size

struct HMMParameters{
public:
    
    HMMParameters(){}
    
    HMMParameters(cnt_prec filterSize, cnt_prec viterbiFilterSize, cnt_prec maxDeletion, cnt_prec maxInsertion,
                  cnt_prec batchSize, cnt_prec mapQ, prob_prec matchTransition, prob_prec insertionTransition,
                  prob_prec deletionTransitionFactor, prob_prec matchEmission):
    filterSize(filterSize), viterbiFilterSize(viterbiFilterSize), maxDeletion(maxDeletion), maxInsertion(maxInsertion),
    batchSize(batchSize), mapQ(mapQ), matchTransition(matchTransition), insertionTransition(insertionTransition),
    deletionTransitionFactor(deletionTransitionFactor), matchEmission(matchEmission){
        deletionTransition = 1.000 - (matchTransition + insertionTransition);
        mismatchEmission = (prob_prec)(1 - matchEmission)/3.00;
        insertionEmission = (prob_prec)1/3.00; //total nucleotide = 4; Emission prob for each except one
    }

    HMMParameters(const HMMParameters& cpy):
    filterSize(cpy.filterSize), viterbiFilterSize(cpy.viterbiFilterSize), maxDeletion(cpy.maxDeletion),
    maxInsertion(cpy.maxInsertion), batchSize(cpy.batchSize), mapQ(cpy.mapQ), matchTransition(cpy.matchTransition),
    insertionTransition(cpy.insertionTransition), deletionTransition(cpy.deletionTransition),
    deletionTransitionFactor(cpy.deletionTransitionFactor), matchEmission(cpy.matchEmission),
    mismatchEmission(cpy.mismatchEmission), insertionEmission(cpy.insertionEmission){}
    
    HMMParameters& operator=(const HMMParameters& rhs){
        filterSize = rhs.filterSize; viterbiFilterSize = rhs.viterbiFilterSize; maxDeletion = rhs.maxDeletion;
        maxInsertion = rhs.maxInsertion; batchSize = rhs.batchSize; mapQ = rhs.mapQ;
        matchTransition = rhs.matchTransition; insertionTransition = rhs.insertionTransition;
        deletionTransition = rhs.deletionTransition; deletionTransitionFactor = rhs.deletionTransitionFactor;
        matchEmission = rhs.matchEmission; mismatchEmission = rhs.mismatchEmission;
        insertionEmission = rhs.insertionEmission;
        
        return *this;
    }
    
    friend std::ostream& operator<<(std::ostream& out, const HMMParameters& rhs){
        out << "Min mapping quality: " << rhs.mapQ << std::endl << "Filter size: " << rhs.filterSize << std::endl
        << "Viterbi filter size: " << rhs.viterbiFilterSize << std::endl << "Viterbi batch size: " << rhs.batchSize <<
        std::endl << "Maximum insertion: " << rhs.maxInsertion << std::endl << "Maximum deletion: " << rhs.maxDeletion
        << std::endl << "Match transition probability: " << rhs.matchTransition << std::endl <<
        "Insertion transition probability: " << rhs.insertionTransition << std::endl <<
        "Deletion transition probability: " << rhs.deletionTransition << std::endl <<
        "Deletion transition factor: " << rhs.deletionTransitionFactor << std::endl <<
        "Match emission probability: " << rhs.matchEmission << std::endl << "Mismatch emission probability: " <<
        rhs.mismatchEmission << std::endl << "Insertion emission probability: " << rhs.insertionEmission;
        
        return out;
    }
    
    cnt_prec filterSize;
    cnt_prec viterbiFilterSize;
    cnt_prec maxDeletion;
    cnt_prec maxInsertion;
    cnt_prec batchSize;
    cnt_prec mapQ;
    prob_prec matchTransition;
    prob_prec insertionTransition;
    prob_prec deletionTransition;
    prob_prec deletionTransitionFactor;
    prob_prec matchEmission;
    prob_prec mismatchEmission;
    prob_prec insertionEmission;
};

/** @brief Represents a state in the profile hidden Markov model graph. Total size is usually 16 bytes
 *
 */
struct SeqNode{
public:

    SeqNode(){}
    SeqNode(ind_prec index, ind_prec charIndex, cnt_prec insize, char nuc, char nextNuc):
    index(index), charIndex(charIndex), nuc(nuc), isMatch((index%(insize+1) == 0)?true:false),
    isLastInsertion((insize > 0 && index%(insize+1) == (ind_prec)insize)?true:false), nextNuc(nextNuc){}

    /** @brief Emission probability calculation
     *
     *  @param character Basepair to calculate the emission probability
     *  @param parameters Calculation is based on the specified parameters
     *  @return Emission probability of the given character (basepair) [0-1]
     */
    prob_prec getEmissionProb(const char& character, const HMMParameters& parameters){

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
    prob_prec transitionProbFromThisNode(const ind_prec& toIndex, HMMParameters& parameters){

        if(MATCH_OFFSET(charIndex, 1, parameters.maxInsertion) == toIndex){
            if(isLastInsertion){
                //match transition prob for last insertion state
                return parameters.matchTransition + parameters.insertionTransition;
            }
            return parameters.matchTransition; //match transition prob
        }

        if(index+1 == toIndex) return parameters.insertionTransition; //insertion transition

        //deletion transition calculations: normalized polynomial distribution
        prob_prec count = 0;
        prob_prec start = 1;
        prob_prec transitionProb = 0;
        for(cnt_prec curDel = parameters.maxDeletion+1; curDel > 1; --curDel){
            if(MATCH_OFFSET(charIndex, curDel, parameters.maxInsertion) == toIndex)
                transitionProb = parameters.deletionTransition*start;
            count+=start;
            start*=parameters.deletionTransitionFactor;
        }

        return (count==0)?0:transitionProb/(prob_prec)count;
    }

    bool isMatchState() const { return isMatch; }
    bool isLastInsertionState() const { return isLastInsertion; }
    char getNucleotide() const { return nuc;}
    ind_prec getIndex() const { return index;}
    ind_prec getCharIndex() const { return charIndex; }

private:
    ind_prec index; 
    ind_prec charIndex; //which character it corresponds to in the sequencing read: 1-based
    char nuc; //basepair in the charIndex
    bool isMatch; //is a match state?
    bool isLastInsertion;
    char nextNuc; //basepair in the next position
};

struct Read{
public:
    Read(const seqan::Dna5String& readSeq, ind_prec pos, const seqan::String<seqan::CigarElement<> >& cigar):
    pos(pos){
        this->read = readSeq;
        ind_prec curIndex = 0;
        for(ind_prec i = 0; i < length(cigar); ++i){
            char type = cigar[i].operation;
            if(type == 'S') seqan::erase(read, curIndex, curIndex+cigar[i].count);
            else if(type == 'M' || type == 'I') curIndex += cigar[i].count;
        }
    }
    
    seqan::String<char> read;
    ind_prec pos;
};

/*
 * Can be used as directed edge between from (outgoing edge source), and curState (incoming edge source).
 * operator< implemented so that std::set can evaulate this data structure
 */
struct TransitionInfoNode{

    ind_prec from, toState;
    TransitionInfoNode(ind_prec from, ind_prec toState):from(from), toState(toState){}

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
inline void findMaxValues(const T* values, bool* selectedIndices, const ind_prec startValues, const ind_prec endValues,
                          const size_t maxValuesSize){

    std::priority_queue<std::pair<prob_prec, ind_prec>, std::vector<std::pair<prob_prec, ind_prec> >,
                        std::greater<std::pair<prob_prec, ind_prec> > > maxValues;

    for(ind_prec curState = 0; curState < endValues - startValues; ++curState) {
        if(maxValues.size() < maxValuesSize)
            maxValues.push(std::pair<prob_prec, ind_prec>(values[curState], curState));
        else if(maxValues.top().first < values[curState]){
            maxValues.pop();
            maxValues.push(std::pair<prob_prec, ind_prec>(values[curState], curState));
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
inline void insertNewForwardTransitions(std::vector<TransitionInfoNode>* transitionSet, const SeqNode& node,
                                        const cnt_prec numberOfDeletions, const cnt_prec maxInsertion){
    //next insertion
    if(!node.isLastInsertionState())
        transitionSet->push_back(TransitionInfoNode(node.getIndex(), node.getIndex()+1));

    //match and deletions
    for(cnt_prec curOffset = 1; curOffset <= numberOfDeletions+1; ++curOffset)
        transitionSet->push_back(TransitionInfoNode(node.getIndex(), MATCH_OFFSET(node.getCharIndex(), curOffset,
                                                                                  maxInsertion)));
}

/*
* Using the information provided with node and numberOfDeletions, inserts the transitions that should be made after
* this node.
* If you are going to change the transition structure, change it from there. These are imaginary edges in the graph.
*/
inline void insertNewBackwardTransitions(std::vector<TransitionInfoNode>* transitionSet, const SeqNode& node,
                                         const cnt_prec numberOfDeletions, const cnt_prec maxInsertion){

    if(node.getCharIndex() == 0) return;

    if(node.isMatchState()){
        for(ind_prec offset = 1; offset <= (ind_prec)numberOfDeletions+1 && offset <= node.getCharIndex(); ++offset){
            
            transitionSet->push_back(TransitionInfoNode(node.getIndex(), MATCH_DOFFSET(node.getCharIndex(), offset,
                                                                                       maxInsertion))); //deletion
            
            for(cnt_prec curInsertion = 1; curInsertion <= maxInsertion; ++curInsertion){
                transitionSet->push_back(TransitionInfoNode(node.getIndex(),INSERTION_DOFFSET(node.getCharIndex(), offset,
                                                                                             curInsertion,maxInsertion)));
            }
        }
    }else{
        transitionSet->push_back(TransitionInfoNode(node.getIndex(), node.getIndex()-1));//match //can never be <0
    }
}

inline bool fillBuffer(seqan::BamFileIn& alignmentFileIn, const seqan::FaiIndex& readsIndex,
                       seqan::BamAlignmentRecord& record, std::vector<Read>& reads, const int32_t& contigId,
                       cnt_prec mapQ, ind_prec size){
    
    reads.clear();
    while(contigId == record.rID && reads.size() < size){
        
        if(!hasFlagUnmapped(record) && !hasFlagLast(record) && record.mapQ >= mapQ){
            seqan::Dna5String alignedString;
            unsigned qId = 0;
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

#endif /* HMMCommons_h */
