/** @file main.cpp
 * @brief One sentence brief
 *
 * More details
 * In multiple lines
 * Copyright © 2020 SAFARI
 *
 * @author Can Firtina
 * @bug No known bug
 */

#include <stdio.h>
#include "CommandLineParser.h"
#include "Polisher.h"

int main(int argc, const char** argv){

    CommandLineParser options;
    seqan::ArgumentParser::ParseResult parseRes = parseCommandOptions(options, argc, argv);

    if(parseRes != seqan::ArgumentParser::PARSE_OK){
        return parseRes == seqan::ArgumentParser::PARSE_ERROR;
    }
    
    if(!options.shouldQuite){
        std::cout << "Assembly: " << toCString(options.assembly) << std::endl <<
        "Pair of a set of reads and their alignments:" << std::endl;
        for(unsigned int i = 0; i < options.readSets.size(); ++i)
            std::cout << toCString(options.readSets.at(i)) << ", " << toCString(options.alignmentSets.at(i)) << std::endl;
    }
    
    Polisher* polisher = new Polisher();
    polisher->setParameters(options.filterSize, options.viterbiFilterSize, options.maxDeletion,
                            options.maxInsertion, options.batchSize, options.chunkSize, options.mapQ,
                            options.matchTransition, options.insertionTransition, options.deletionTransitionFactor,
                            options.matchEmission);
    
    
    polisher->polish(options.assembly, options.readSets, options.alignmentSets, options.output, options.maxThread,
                     options.shouldQuite);

    delete polisher;
    
    return 0;
}
