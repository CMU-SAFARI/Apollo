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
    
    Polisher* polisher = new Polisher();
    polisher->setParameters(options.filterSize, options.viterbiFilterSize, options.maxDeletion, options.maxInsertion,
                            options.batchSize, options.mapQ, options.matchTransition, options.insertionTransition,
                            options.deletionTransitionFactor, options.matchEmission);
    
    
    polisher->polish(options.assembly, options.readSets, options.alignmentSets, options.output, options.maxThread,
                     options.shouldQuite);

    delete polisher;
    
    return 0;
}
