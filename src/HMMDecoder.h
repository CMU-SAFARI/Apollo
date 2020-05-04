/** @file HMMDecoder.h
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/

#ifndef HMMDecoder_h
#define HMMDecoder_h

#include <stdio.h>
#include "HMMGraph.h"

class HMMDecoder{
    
    //offset is an index for state, not a basepair index
    void viterbiThreadPool(ind_prec seqLength, ind_prec& seqIndex, seqan::String<seqan::Dna5String>& decodedOut);
    
    HMMGraph* graph;
    std::mutex indexMutex;
    
public:
    
    HMMDecoder();
    HMMDecoder(HMMGraph* graph);
    ~HMMDecoder();
    
    void backtrace(unsigned threadCnt, seqan::String<seqan::Dna5String>& decodedOut);
};

#endif /* HMMDecoder_h */
