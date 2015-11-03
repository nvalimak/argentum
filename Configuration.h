#ifndef _Configuration_h_
#define _Configuration_h_
#include "default.h"
#include "InputReader.h"
#include "TreeDistance.h"
#include <string>
#include <iostream>

/**
 * Command line option parser
 */
class Configuration
{
public:
    Configuration(int, char **);
    void print_short_usage();
    void print_usage();

    std::string prog;
    std::string inputfile;
    std::string dotfile;
    InputReader::input_format_t inputformat;
    TreeDistance::tree_distance_t treedistance;
    TreeDistance::distance_scaling_t distancescaling;
    unsigned nrows;
    unsigned skip;
    bool newick;
    bool verbose;
    bool debug;
    unsigned debug_interval;
    bool good;

private:
    bool parse(int, char **);
    Configuration();
};
#endif
