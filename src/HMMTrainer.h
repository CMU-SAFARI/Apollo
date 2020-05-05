/** @file HMMTrainer.h
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/
#ifndef HMMTrainer_h
#define HMMTrainer_h

#include <stdio.h>
#include "HMMGraph.h"

class HMMTrainer{
    
    void fbThreadPool(const std::vector<Read>& reads, ind_prec& readIndex);
    
    ind_prec fillForwardMatrix(prob_prec** forwardMatrix, const char* read, const ind_prec startPosition,
                               ind_prec maxDistanceOnContig, const ind_prec readLength);
    
    bool fillBackwardMatrix(prob_prec** backwardMatrix, const char* read, ind_prec startPosition,
                            ind_prec maxDistanceOnContig, const ind_prec readLength);
        
    HMMGraph* graph;
    std::mutex indexMutex;
    std::mutex emissionProbMutex;
    std::mutex transitionProbMutex;
    
public:
    HMMTrainer();
    HMMTrainer(HMMGraph* graph);
    ~HMMTrainer();
    
    void calculateFB(std::vector<seqan::BamFileIn>& alignmentSetsIn, std::vector<seqan::BamAlignmentRecord>& curRecords,
                     const std::vector<seqan::FaiIndex>& readFAIs, unsigned thread);
    void maximizeEM();
};
#endif /* HMMTrainer_h */
