#include "InputReader.h"
#include "Tree.h"
#include "TreeControllerSimple.h"
#include <sstream>
#include <cstdlib>
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

void print_usage(char const * prgrm)
{
    cerr << "usage: " << prgrm << " [options]" << endl
         << " -i,--input <file>    Input filename" << endl
         << " -V,--VCF             Input format is VCF" << endl
         << " -s,--scrm            Input format is SCRM text (without trees)" << endl
         << " -S,--plaintext       Input format is simple text" << endl
         << " -r,--rewind          Rewind test" << endl
         << " -z,--skip <n>        Skip first n sites of input" << endl
         << " -n,--nrows <n>       Process n sites of input" << endl
         << " -d,--dot <file>      Graphwiz DOT output filename" << endl
         << " -v,--verbose         Verbose output" << endl;
}

int main(int argc, char ** argv)
{
    string inputfile = "", dotfile = "";
    InputReader::input_format_t inputformat = InputReader::input_unset;
    unsigned nrows = ~0u;
    unsigned skip = ~0u;
    bool rewind = false;
    bool verbose = false;
    bool debug = false;
    unsigned debug_interval = 1000;
    
    static struct option long_options[] =
        {
            {"input",     required_argument, 0, 'i'},
            {"vcf",       no_argument,       0, 'V'},
            {"scrm",      no_argument,       0, 's'},
            {"plaintext", no_argument,       0, 'S'},
            {"rewind",    no_argument,       0, 'r'},
            {"nrows",     required_argument, 0, 'n'},
            {"skip",      required_argument, 0, 'z'},
            {"dot",       required_argument, 0, 'd'},
            {"verbose",   no_argument,       0, 'v'},
            {"debug",     required_argument, 0, 'D'},
            {"help",      no_argument,       0, 'h'},
            {0, 0, 0, 0}
        };
    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "i:VsSrn:z:d:vD:h",
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
        case 'r':
            rewind = true; break;
        case 'n':
            nrows = atoi_min(optarg, 1, "-n,--nrows", argv[0]); break;
        case 'z':
            skip = atoi_min(optarg, 1, "-S,--skip", argv[0]); break;
        case 'd':
            dotfile = string(optarg); break;
        case 'v':
            verbose = true; break;
        case 'D':
            debug = true; debug_interval = atoi_min(optarg, 1, "-D,--debug", argv[0]); break;
        case 'h':
            print_usage(argv[0]); return 1;
        default:
            cerr << "error: invalid parameter; see -h for usage!" << endl;
            return 1;
        }
    }

    if (inputformat == InputReader::input_unset)
    {
        cerr << "error: input format needs to be specified (-s,--plaintext or -V,--vcf)" << endl;
        return 1;
    } 
    
    unsigned step = 0;
    InputColumn ic;
    InputReader *inputr = InputReader::build(inputformat, inputfile);
    if (!inputr->good())
    {
        cerr << "error: could not read input file '" << inputfile << "'" << endl;
        return 1;
    }

    /**
     * Main loop over input data
     */
    if (!inputr->next(ic))
    {
        cerr << "error: empty input file?!" << endl;
        return 1;
    }
    PointerTree tree(ic);
    TreeControllerSimple tc(tree, debug, dotfile);
    vector<InputColumn> columns;
    unsigned prev_hsize = 0;
    unsigned prev_reused = 0;
    unsigned prev_stashed = 0;
    unsigned prev_relocates = 0;
    if (skip != ~0u)
        do
        {
            step ++;
            if (rewind)
                columns.push_back(ic);
        } while (inputr->next(ic) && step < skip);

    do
    {
        tc.process(ic, step);
        if (rewind)
            columns.push_back(ic);
        
        if (!dotfile.empty())
            tree.outputDOT(dotfile, step);

        if (verbose && step % debug_interval == 0)
        {
            cerr << "at step " << step << ", history = " << tree.historySize() << " (" << tree.historySize()-prev_hsize
                 << ", reused " << prev_reused << "+" << tree.reusedHistoryEvents()-prev_reused << ", stashed " << prev_stashed << "+" << tree.numberOfStashed()-prev_stashed
                 << ", relocates " << prev_relocates << "+" << tree.numberOfRelocates()-prev_relocates
                 << "), tree size = " << tree.nnodes() << " (" << tree.size() << " leaves, "
                 << tc.countGhostBranches(tree.root()) << " ghostbranch, " << tc.countUnaryGhosts(tree.root()) << " unaryghosts, "
                 << tc.countBranchingGhosts(tree.root()) << " branchingghost, " << tc.countActive(tree.root()) << " active)" << endl;
            prev_hsize = tree.historySize();
            prev_reused = tree.reusedHistoryEvents();
            prev_stashed = tree.numberOfStashed();
            prev_relocates = tree.numberOfRelocates();
        }
        /*        if (debug && step >178 && step < 200)
        {
            cerr << "row " << step << ": ";
            for (InputColumn::iterator it = ic.begin(); it != ic.end(); ++it)
                cerr << (int)*it;
            cerr << endl;

            }*/
        ++step;
    } while (inputr->next(ic) && step < nrows);

    if (verbose && rewind)
        cerr << "Finished after " << step << " steps of input. Rewinding " << tree.historySize() << " history events..." << endl;

    assert (!rewind || columns.size() == step || columns.size() == step-skip);
    tc.deFloatAll(tree.root());
    unsigned h = columns.size();
    while (h-- > 0 || (skip != ~0u && h > skip))
    {
        if (verbose && h % 1000 == 0)
            cerr << "at " << h << "/" << step << ", history = " << tree.historySize() << ", tree size = " << tree.nnodes() << " (" << tree.size() << " leaves, "
                 << tc.countGhostBranches(tree.root()) << " ghostbranch, " << tc.countUnaryGhosts(tree.root()) << " unaryghosts, "
                 << tc.countBranchingGhosts(tree.root()) << " branchingghost, " << tc.countActive(tree.root()) << " active)" << endl;
        
        tc.rewind(columns[h], h);                      
    }
    
    if (verbose)
        cerr << "main-simple done. Finished after " << step << " steps of input." << endl;
    
    // Clean up
    delete inputr; inputr = 0;
    return 0;
}
