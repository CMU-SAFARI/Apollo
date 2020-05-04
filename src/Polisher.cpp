/** @file Polisher.cpp
* @brief One sentence brief
*
* More details
* In multiple lines
* Copyright © 2020 SAFARI
*
* @author Can Firtina
* @bug No known bug
*/

#include "Polisher.h"
#include "HMMTrainer.h"
#include "HMMDecoder.h"

Polisher::Polisher():isParametersSet(false){}
Polisher::~Polisher(){}

void Polisher::setParameters(cnt_prec filterSize, cnt_prec viterbiFilterSize, cnt_prec maxDeletion,cnt_prec maxInsertion,
                             cnt_prec batchSize, cnt_prec mapQ, prob_prec matchTransition, prob_prec insertionTransition,
                             prob_prec deletionTransitionFactor, prob_prec matchEmission){
    
    hmmParams = HMMParameters(filterSize, viterbiFilterSize, maxDeletion, maxInsertion, batchSize, mapQ, matchTransition,
                              insertionTransition, deletionTransitionFactor, matchEmission);
    isParametersSet = true;
}

void Polisher::setParameters(const HMMParameters& parameters){
    
    hmmParams = parameters;
    isParametersSet = true;
}

bool Polisher::polish(seqan::String<char> assemblyFile, std::vector<seqan::String<char> > readSets,
                      std::vector<seqan::String<char> > alignmentSets, seqan::String<char> outputFile, unsigned thread,
                      bool shouldQuite){
    
    if(!isParametersSet){
        std::cout << "Error: Parameters are not set properly. Will not polish." << std::endl;
        return false;
    }
    
    if(!shouldQuite){
        std::cout << "Output file: " << outputFile << std::endl << hmmParams << std::endl << "Max thread: " << thread
        << std::endl;
    }
    
    //Index file for assembly fasta
    seqan::FaiIndex assemblyFAI;
    if(!build(assemblyFAI, toCString(assemblyFile))){ //read the assembly file and build the index file
        
        std::cerr << "ERROR: Could not build FAI index for file " << assemblyFile << ". Make sure the FASTA file "
        << "is structured correctly (e.g., all lines must include same number of bases)." << std::endl;
        
        return false;
    }
    
    //Number of contigs in the assembly file
    uint64_t nContigs = seqan::numSeqs(assemblyFAI);
    
    //reading/creating the indices for each of the alignment and the read set pairs
    std::vector<seqan::BamAlignmentRecord> curRecords(alignmentSets.size());
    std::vector<seqan::BamFileIn> alignmentSetsIn(alignmentSets.size());
    std::vector<seqan::BamIndex<seqan::Bai> > baiIndices(alignmentSets.size()); //indexed alignment files
    std::vector<seqan::FaiIndex> readFAIs(alignmentSets.size()); //indexed read sets
    for(size_t i = 0; i < alignmentSets.size(); ++i){

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
        }catch(seqan::Exception const& e){
            std::cerr << "ERROR: " << e.what() << std::endl;
            return false;
        }

        if(!build(readFAIs[i], toCString(readSets[i]))){
            std::cerr << "ERROR: Could not build FAI index for file " << readSets[i] << ". Make sure the FASTA file "
            << "is structured correctly (e.g., all lines must include same number of bases)." << std::endl;
            
            return false;
        }
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
    
    int32_t contigId = 0;
    //Polishing each contig
    while(contigId < (int32_t)nContigs){
        
        //name of the current contig to polish
        seqan::CharString curContigName = (seqan::CharString)seqan::sequenceName(assemblyFAI, (unsigned)contigId);
        seqan::CharString polishedPrefix = "polished_";
        seqan::Dna5String curSeq; //assembly contig
        seqan::readSequence(curSeq, assemblyFAI, (unsigned)contigId);
        seqan::Dna5String correctedContig;
        
        //Contig can only be polished if there is an alignment to that contig.
        //shouldPolish is false if there is no alignment
        bool shouldPolish = false;
        for(size_t curSet = 0; curSet < curRecords.size() && !shouldPolish; ++curSet){
            shouldPolish |= (curRecords[curSet].rID == contigId)&(!atEnd(alignmentSetsIn[curSet]));
        }
        
        //We confirmed there is an alignment for the current contig to be polished. Starting polishing.
        if(shouldPolish){
            
            HMMGraph *graph = new HMMGraph(hmmParams, curSeq, contigId);
            graph->buildGraph();
            
            HMMTrainer trainer(graph);
            trainer.calculateFB(alignmentSetsIn, curRecords, readFAIs, thread);
            trainer.maximizeEM();
            
            seqan::String<seqan::Dna5String> decodedOut;
            HMMDecoder decoder(graph);
            decoder.backtrace(thread, decodedOut);
            
            for(ind_prec i = 0; i < (ind_prec)length(decodedOut); ++i)
                append(correctedContig, decodedOut[i]);
            
            delete graph;
            
            if(length(correctedContig) > 0){
                curSeq = correctedContig;
                append(polishedPrefix, curContigName);
                curContigName = polishedPrefix;
            }
        }else{ //If there is no alignment to a contig, report that the original (unpolished) contig will be produced.
          std::cerr << "The contig with id " << toCString(curContigName) << " could not be polished because there is no read aligning to it. Original (i.e., unpolished) sequence will be reported." << std::endl;
        }
        
        //write polished/unpolished contig to the output file
        writeRecord(correctedReadsOut, curContigName, curSeq);
        ++contigId; //iterate to next contig to be polished
    }

    correctedReadStream.close();
    if(!shouldQuite) std::cout << std::endl << "Results have been written under " << outputFile << std::endl;

    return true;
}
