/** @file HMMGraph.h
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/

#ifndef HMMGraph_h
#define HMMGraph_h

#include <stdio.h>
#include "HMMCommons.h"

class HMMGraph{

    bool isDataFreed;
    bool isParametersSet;
    
public:
    
    HMMParameters params;
    SeqNode* seqGraph;
    seqan::Dna5String contig;
    seqan::Dna5String correctedContig;
    int32_t contigId;
    ind_prec numberOfStates;
    ind_prec contigLength;
    ind_prec numOfTransitionsPerState;
    prob_prec* preCalculatedTransitionProbs;
    prob_prec** transitionProbs;
    prob_prec** emissionProbs;
    cnt_prec* stateProcessedCount;
    cnt_prec** transitionProcessedCount;
    std::pair<prob_prec, char>* maxEmissionProbs;
    
    HMMGraph(const HMMParameters& parameters, const seqan::Dna5String& contig, int32_t contigId = 0);
    HMMGraph();
    ~HMMGraph();
    
    void setParameters(const HMMParameters& parameters);
    bool getIsParametersSet();
    
    void setContig(const seqan::Dna5String& contig, int32_t contigId = 0);
    void buildGraph();
    void freeGraph();
};

#endif /* HMMGraph_h */
