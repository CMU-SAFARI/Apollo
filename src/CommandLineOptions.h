/** @file CommandLineOptions.h
 * @brief Parses command line options.
 *
 * It is used to parse command line options for both preprocessing and correction step.
 * Copyright © 2018 Can Firtina. All rights reserved.
 *
 * @author Can Firtina
 * @bug No bug currently
 */

#ifndef COMMAND_LINE_OPTIONS_H_
#define COMMAND_LINE_OPTIONS_H_

#include <iostream>
#include <vector>
#include <seqan/arg_parse.h>

#define CORRECT_FLAG 2

/** @brief Struct for holding command line options and their values.
 *
 *
 *  @param CommandLineOptions All variables are defined below
 *  @see parseCommandOptions()
 *  @see parseInitialCommand()
 */
struct CommandLineOptions
{
    CommandLineOptions():
    filterSize(100), viterbiFilterSize(5), maxDeletion(10), maxInsertion(3), batchSize(5000), matchTransition(0.85),
    insertionTransition(0.1), shouldQuite(false), matchEmission(0.97), maxThread(1), mapQ(0),
    deletionTransitionFactor(2.5){}
    
    unsigned int filterSize;
    unsigned int viterbiFilterSize;
    unsigned int maxDeletion;
    unsigned int maxInsertion;
    unsigned int batchSize;
    double matchTransition;
    
    double insertionTransition;
    bool shouldQuite;
    double matchEmission;
    unsigned int maxThread;
    unsigned int mapQ;
    double deletionTransitionFactor;
    
    seqan::String<char> assembly;
    seqan::String<char> alignment;
    std::vector<seqan::String<char> > reads;
    seqan::String<char> output;
};

/** @brief Parse values in order to run either preprocessing step or correction step.
 *
 *  @param options Stores parsed values
 *  @param argc Number of arguments specified while running Apollo
 *  @param argv Argument values array
 *  @param phase Preprocessing or Correction phase?
 *  @see parseInitialCommand()
 *  @return seqan::ArgumentParser::PARSE_OK if everything went well
 */
seqan::ArgumentParser::ParseResult
parseCommandOptions(CommandLineOptions& options, int argc, char const **argv){
    
    using namespace std;
    seqan::ArgumentParser parser("Apollo: A Profile HMM-based assembly polishing algorithm for long reads");
    
    setVersion(parser, "0.1");
    setDate(parser, "April 2018");
    
    
    addOption(parser, seqan::ArgParseOption("a", "assembly", "fast{a,q} file which contains the assembly constructed "
                                            "using provided long reads", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
    setRequired(parser, "a");
    
    addOption(parser, seqan::ArgParseOption("l", "long", "Long sequencing reads", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
    setRequired(parser, "l");
    
    addOption(parser, seqan::ArgParseOption("m", "alignment", "{s,b}am file which contains alignments of long reads"
                                            " for the assembly.", seqan::ArgParseArgument::INPUT_FILE, "FILE"));
    setRequired(parser, "m");
    
    addOption(parser, seqan::ArgParseOption("o", "output", "Output file to write the resulting reads",
                                            seqan::ArgParseArgument::OUTPUT_FILE, "FILE"));
    setRequired(parser, "o");
    
    addOption(parser, seqan::ArgParseOption("q", "mapq", "Minimum mapping quality for a long-assembly alignment "
                                            "to use in correction. Note that if multiple alignment specified, "
                                            "aligners report these mapping quality as 0",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(getOption(parser, "mapq"), options.mapQ);
    seqan::setMinValue(parser, "mapq", "0");
    seqan::setMaxValue(parser, "mapq", "255");
    
    addOption(parser, seqan::ArgParseOption("f", "filter", "Filter size that allows calculation of at most mf "
                                            "many most probable transitions in each time step. This parameter is "
                                            "directly proportional to running time.",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(getOption(parser, "filter"), options.filterSize);
    seqan::setMinValue(parser, "filter", "1");
    
    addOption(parser, seqan::ArgParseOption("v", "vf", "Filter size for the Viterbi algorithm that allows calculation"
                                            " of at most vf many most probable states in each time step. This "
                                            "parameter is directly proportional to running time.",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(getOption(parser, "vf"), options.viterbiFilterSize);
    seqan::setMinValue(parser, "vf", "1");
    
    addOption(parser, seqan::ArgParseOption("i", "maxi", "Maximum number of insertions in a row. This "
                                            "parameter is directly proportional to the running time.",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(getOption(parser, "maxi"), options.maxInsertion);
    seqan::setMinValue(parser, "maxi", "0");
    
    addOption(parser, seqan::ArgParseOption("d", "maxd", "Maximum number of deletions in a row. This "
                                            "parameter is directly proportional to the running time.",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(getOption(parser, "maxd"), options.maxDeletion);
    seqan::setMinValue(parser, "maxd", "0");
    
    addOption(parser, seqan::ArgParseOption("tm", "mtransition", "Initial transition probability to a match "
                                            "state. See --itransition as well.",
                                            seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    setDefaultValue(getOption(parser, "tm"), options.matchTransition);
    seqan::setMinValue(parser, "tm", "0");
    seqan::setMaxValue(parser, "tm", "1");
    
    addOption(parser, seqan::ArgParseOption("ti", "itransition", "Initial transition probability to a "
                                            "insertion state. Note that: deletion transition probability = "
                                            "1 - (matchTransition + insertionTransition)",
                                            seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    setDefaultValue(getOption(parser, "ti"), options.insertionTransition);
    seqan::setMinValue(parser, "ti", "0");
    seqan::setMaxValue(parser, "ti", "1");
    
    addOption(parser, seqan::ArgParseOption("df", "dfactor", "Factor for the polynomial distribution to calculate each "
                                            "deletion transition. Higher value favors less deletions.",
                                            seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    setDefaultValue(getOption(parser, "df"), options.deletionTransitionFactor);
    seqan::setMinValue(parser, "df", "0");
    
    addOption(parser, seqan::ArgParseOption("em", "memission", "Initial emission probability of a match to a "
                                            "reference. Note that: mismatch emission probability = "
                                            "(1-matchEmission)/3", seqan::ArgParseArgument::DOUBLE, "FLOAT"));
    setDefaultValue(getOption(parser, "em"), options.matchEmission);
    seqan::setMinValue(parser, "em", "0");
    seqan::setMaxValue(parser, "em", "1");
    
    addOption(parser, seqan::ArgParseOption("b", "batch", "Number of basepairs that Viterbi decodes per thread. "
                                            "Setting it to zero will decode the entire contig in a single thread",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(getOption(parser, "b"), options.batchSize);
    seqan::setMinValue(parser, "b", "0");
    
    addOption(parser, seqan::ArgParseOption("t", "thread", "Number of threads to use",
                                            seqan::ArgParseArgument::INTEGER, "INT"));
    setDefaultValue(getOption(parser, "t"), options.maxThread);
    seqan::setMinValue(parser, "t", "1");
    
    addOption(parser, seqan::ArgParseOption("nv", "noVerbose", "Apollo runs quitely with no informative output"));
    
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);
    if (res == seqan::ArgumentParser::PARSE_OK){
        
        getOptionValue(options.assembly, parser, "a");
        getOptionValue(options.alignment, parser, "m");
        seqan::CharString shortFileStr;
        getOptionValue(shortFileStr, parser, "l");
        options.reads.push_back(shortFileStr);
        getOptionValue(options.output, parser, "o");
        getOptionValue(options.mapQ, parser, "q");
        getOptionValue(options.filterSize, parser, "f");
        getOptionValue(options.viterbiFilterSize, parser, "vf");
        getOptionValue(options.maxInsertion, parser, "i");
        getOptionValue(options.maxDeletion, parser, "d");
        getOptionValue(options.matchTransition, parser, "tm");
        getOptionValue(options.insertionTransition, parser, "ti");
        getOptionValue(options.deletionTransitionFactor, parser, "df");
        getOptionValue(options.matchEmission, parser, "em");
        getOptionValue(options.batchSize, parser, "b");
        getOptionValue(options.maxThread, parser, "t");
        options.shouldQuite = isSet(parser, "nv");
        
        if(options.matchTransition + options.insertionTransition > 1){
            std::cerr << "ERROR: (matchTransition + insertionTransition) cannot be larger than 1 whereas the sum "
            << "is now: " << options.matchTransition + options.insertionTransition << std::endl;
            return seqan::ArgumentParser::PARSE_ERROR;
        }
    }
    
    return res;
}

#endif
