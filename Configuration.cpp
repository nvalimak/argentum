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
         << " -s,--scrm            Input format is SCRM text (with trees)" << endl
         << " -S,--plaintext       Input format is simple text" << endl
         << " --forward            Run just the forward run (without prediction)" << endl
         << " --backward-forward   Run with backward prediction" << endl
         << " --forward-backward-forward   Run with forward and backward prediction" << endl
         << " --distance <type>    Use distance type <leaves> or <depth> (default: none)" << endl
         << " --scaling <type>     Use distance scaling of <log> or <expm> (default: none)" << endl
         << " -n,--nrows <n>       Process n sites of input" << endl
         << " --newick             Output Newick trees (to standard output)" << endl
         << " --enumerate          Output parent-child ranges (to standard output)" << endl
         << " -d,--dot <file>      Graphwiz DOT output filename" << endl
         << " -v,--verbose         Verbose output" << endl;
}

Configuration::Configuration(int argc, char ** argv)
    : prog(argv[0]), inputfile(""), dotfile(""), inputformat(InputReader::input_unset),
      treedistance(TreeDistance::distance_unset), distancescaling(TreeDistance::scaling_none), nrows(~0u),
      newick(false), enumerate(false), no_prediction(false), scrm_prediction(false), verbose(false), debug(false), debug_interval(2000), good(true)
{
    good = parse(argc, argv);
}

bool Configuration::parse(int argc, char ** argv)
{
    enum long_opt_names { debug_skip_opt = 256, no_prediction_opt, scrm_prediction_opt, forward_opt, backward_forward_opt, forward_backward_forward_opt, distance_opt, scaling_opt };
    
    static struct option long_options[] =
        {
            {"input",     required_argument, 0, 'i'},
            {"vcf",       no_argument,       0, 'V'},
            {"scrm",      no_argument,       0, 's'},
            {"plaintext", no_argument,       0, 'S'},
            {"newick",    no_argument,       0, 'N'},
            {"enumerate", no_argument,       0, 'E'},
            {"no-prediction", no_argument,   0, no_prediction_opt},
            {"scrm-prediction", no_argument, 0, scrm_prediction_opt},
            {"forward",   no_argument,       0, forward_opt},
            {"backward-forward", no_argument,0, backward_forward_opt},
            {"forward-backward-forward", no_argument, 0, forward_backward_forward_opt},
            {"distance",  required_argument, 0, distance_opt},
            {"scaling",   required_argument, 0, scaling_opt},
            {"nrows",     required_argument, 0, 'n'},
            {"dot",       required_argument, 0, 'd'},
            {"verbose",   no_argument,       0, 'v'},
            {"debug",     no_argument,       0, 'D'},
            {"debug-skip",required_argument, 0, debug_skip_opt},
            {"help",      no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "i:VsSNEn:d:vDh",
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
        case 'E':
            enumerate = true; break;
        case no_prediction_opt:
            no_prediction = true; break;
        case scrm_prediction_opt:
            scrm_prediction = true; break;
        case 'n':
            nrows = atoi_min(optarg, 1, "-n,--nrows", argv[0]); break;
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

    if ((!no_prediction && !scrm_prediction) || (no_prediction == scrm_prediction))
    {
        cerr << "error: specify either --no-prediction or --scrm-prediction, but not both" << endl;
        return false;
    }

    if (enumerate && newick)
    {
        cerr << "error: specify either --newick or --enumerate, but not both" << endl;
        return false;
    }
    return true;
}
