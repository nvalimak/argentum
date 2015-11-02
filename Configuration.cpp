#include "Configuration.h"
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <getopt.h>
using namespace std;

int atoi_min(char const *value, int min, char const *parameter, char const *name)
{
    std::istringstream iss(value);
    int i;
    char c;
    if (!(iss >> i) || iss.get(c))
    {
        cerr << "error: argument of " << parameter << " must be of type <int>, and greater than or equal to " << min << endl
             << "Check " << name << " --help' for more information." << endl;
        std::exit(1);
    }

    if (i < min)
    {
        cerr << "error: argument of " << parameter << " must be greater than or equal to " << min << endl
             << "Check `" << name << " --help' for more information." << endl;
        std::exit(1);
    }
    return i;
}

void Configuration::print_short_usage()
{
    cerr << "See `" << prog << " --help` for command-line options." << endl;
    exit(1);
}

void Configuration::print_usage()
{
    cerr << "usage: " << prog << " [options]" << endl
         << " -i,--input <file>    Input filename" << endl
         << " -V,--VCF             Input format is VCF" << endl
         << " -s,--scrm            Input format is SCRM text (without trees)" << endl
         << " -S,--plaintext       Input format is simple text" << endl
         << " -r,--rewind          Rewind test" << endl
//         << " -z,--skip <n>        Skip first n sites of input" << endl
         << " -n,--nrows <n>       Process n sites of input" << endl
         << " -d,--dot <file>      Graphwiz DOT output filename" << endl
         << " -v,--verbose         Verbose output" << endl;
}

Configuration::Configuration(int argc, char ** argv)
    : prog(argv[0]), inputfile(""), dotfile(""), inputformat(InputReader::input_unset), nrows(~0u),
      skip(~0u), newick(false), verbose(false), debug(false), debug_interval(2000), good(true)
{
    good = parse(argc, argv);
}

bool Configuration::parse(int argc, char ** argv)
{
    const int debug_skip_opt = 256;
    
    static struct option long_options[] =
        {
            {"input",     required_argument, 0, 'i'},
            {"vcf",       no_argument,       0, 'V'},
            {"scrm",      no_argument,       0, 's'},
            {"plaintext", no_argument,       0, 'S'},
            {"newick",    no_argument,       0, 'N'},
            {"nrows",     required_argument, 0, 'n'},
//            {"skip",      required_argument, 0, 'z'}, // TODO
            {"dot",       required_argument, 0, 'd'},
            {"verbose",   no_argument,       0, 'v'},
            {"debug",     no_argument,       0, 'D'},
            {"debug-skip",required_argument, 0, debug_skip_opt},
            {"help",      no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "i:VsSNn:d:vDh",
                            long_options, &option_index)) != -1) 
    {
        switch(c) 
        {
        case 'i':
            inputfile = string(optarg); break;
        case 'V': 
            inputformat = InputReader::input_vcf; break;
        case 's': 
            inputformat = InputReader::input_scrm; break;
        case 'S': 
            inputformat = InputReader::input_plaintext; break;
        case 'N':
            newick = true; break;
        case 'n':
            nrows = atoi_min(optarg, 1, "-n,--nrows", argv[0]); break;
//        case 'z':
//            skip = atoi_min(optarg, 1, "-S,--skip", argv[0]); break;
        case 'd':
            dotfile = string(optarg); break;
        case 'v':
            verbose = true; break;
        case 'D':
            debug = true; break;
        case debug_skip_opt:
            debug_interval = atoi_min(optarg, 1, "--debug-skip", argv[0]); break;
        case 'h':
            print_usage(); exit(1);
        default:
            cerr << "error: invalid parameter; see -h for usage!" << endl;
            return false;
        }
    }

    if (inputformat == InputReader::input_unset)
    {
        cerr << "error: input format needs to be specified (-s,--plaintext or -V,--vcf)" << endl;
        return false;
    }

    if (inputfile == "")
    {
        cerr << "error: input file needs to be specified (-i, --input)" << endl;
        return false;
    }
    return true;
}
