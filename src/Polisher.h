/** @file Polisher.h
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/

#ifndef Polisher_h
#define Polisher_h

#include <stdio.h>
#include "HMMGraph.h"

class Polisher{
    
    HMMParameters hmmParams;
    bool isParametersSet;
    
public:
    Polisher();
    ~Polisher();
    
    void setParameters(cnt_prec filterSize, cnt_prec viterbiFilterSize, cnt_prec maxDeletion, cnt_prec maxInsertion,
                       cnt_prec batchSize, cnt_prec chunkSize, cnt_prec mapQ, prob_prec matchTransition,
                       prob_prec insertionTransition, prob_prec deletionTransitionFactor, prob_prec matchEmission);
    
    void setParameters(const HMMParameters& parameters);
    
    /** @brief Polishes each contig in @assemblyFile reading the alignments from @alignmentFile
    *
    *  @param assemblyFile Fasta/Fastq file that includes sequencing reads for each contig in assembly
    *  @param readSets List of Fasta/Fastq files that include reads that are used to align the reads in @assemblyFile
    *  @param alignmentSets List of SAM/BAM files that include alignments of the reads in @readsFile to reads in @assemblyFile
    *  @param outputFile Output file to write the both corrected/uncorrected reads (contigs)
    *  @param thread Number of threads that can be used
    *  @param shouldQuite If true, there is no standard output
    *  @return True, if there was no problem at any step. Otherwise, false with an appropriate output of the error.
    */
    bool polish(seqan::String<char> assemblyFile, std::vector<seqan::String<char> > readSets,
                std::vector<seqan::String<char> > alignmentSets, seqan::String<char> outputFile, unsigned thread,
                bool shouldQuite);
    
};
#endif
