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
    InputReader *inputr = InputReader::build(inputformat, inputfile, nrows);
    if (!inputr->good())
    {
        cerr << "error: could not read input file '" << inputfile << "'" << endl;
        return 1;
    }

    /**
     * Forward loop over input data
     */
    if (!inputr->next(ic))
    {
        cerr << "error: empty input file?!" << endl;
        return 1;
    }
    PointerTree forward_tree(ic);
    TreeControllerSimple forward_tc(forward_tree, debug, dotfile);
    unsigned prev_hsize = 0;
    unsigned prev_reused = 0;
    unsigned prev_stashed = 0;
    unsigned prev_relocates = 0;
    if (skip != ~0u)
        do
            step ++;
        while (inputr->next(ic) && step < skip);

    do
    {
        forward_tc.process(ic, step);
       
        if (!dotfile.empty())
            forward_tree.outputDOT("forward_" + dotfile, step);

        if (verbose && step % debug_interval == 0)
        {
            cerr << "at forward step " << step << ", history = " << forward_tree.historySize() << " (" << forward_tree.historySize()-prev_hsize
                 << ", reused " << prev_reused << "+" << forward_tree.reusedHistoryEvents()-prev_reused << ", stashed " << prev_stashed << "+" << forward_tree.numberOfStashed()-prev_stashed
                 << ", relocates " << prev_relocates << "+" << forward_tree.numberOfRelocates()-prev_relocates
                 << "), forward_tree size = " << forward_tree.nnodes() << " (" << forward_tree.size() << " leaves, "
                 << forward_tc.countGhostBranches(forward_tree.root()) << " ghostbranch, " << forward_tc.countUnaryGhosts(forward_tree.root()) << " unaryghosts, "
                 << forward_tc.countBranchingGhosts(forward_tree.root()) << " branchingghost, " << forward_tc.countActive(forward_tree.root()) << " active)" << endl;
            prev_hsize = forward_tree.historySize();
            prev_reused = forward_tree.reusedHistoryEvents();
            prev_stashed = forward_tree.numberOfStashed();
            prev_relocates = forward_tree.numberOfRelocates();
        }
        /*        if (debug && step >178 && step < 200) TODO Cleanup
        {
            cerr << "row " << step << ": ";
            for (InputColumn::iterator it = ic.begin(); it != ic.end(); ++it)
                cerr << (int)*it;
            cerr << endl;

            }*/
        ++step;
    } while (step < nrows && inputr->next(ic));

    if (verbose && rewind)
        cerr << "Forward scan done after " << step << " steps of input. Rewinding " << forward_tree.historySize() << " history events..." << endl;

    /**
     * Backward loop over the input data
     */
    assert (!rewind || inputr->size() == step || inputr->size() == step-skip);
    forward_tc.deFloatAll(forward_tree.root());
    unsigned h = inputr->size();
    PointerTree backward_tree(inputr->col(h-1));
    TreeControllerSimple backward_tc(backward_tree, debug, dotfile.empty() ? "" : "backward");
    prev_hsize = 0;
    prev_reused = 0;
    prev_stashed = 0;
    prev_relocates = 0;
    while (h > 0 && (skip == ~0u || h > skip))
    {
        h--;
        LeafDistance dist(ic.size(), 0);
        forward_tc.distance(dist, inputr->col(h));
        forward_tc.rewind(inputr->col(h), h);
        
        backward_tc.process(inputr->col(h), inputr->size()-h-1, dist);
        if (!dotfile.empty())
            backward_tree.outputDOT("backward_" + dotfile, inputr->size()-h-1);

        if (verbose && h % debug_interval == 0)
        {
            //cerr << "rewind forward at " << h << "/" << step << ", history = " << forward_tree.historySize() << ", forward_tree size = " << forward_tree.nnodes() << " (" << forward_tree.size() << " leaves, "
            //     << forward_tc.countGhostBranches(forward_tree.root()) << " ghostbranch, " << forward_tc.countUnaryGhosts(forward_tree.root()) << " unaryghosts, "
            //     << forward_tc.countBranchingGhosts(forward_tree.root()) << " branchingghost, " << forward_tc.countActive(forward_tree.root()) << " active)" << endl;
            cerr << "at backward step " << step << ", history = " << backward_tree.historySize() << " (" << backward_tree.historySize()-prev_hsize
                 << ", reused " << prev_reused << "+" << backward_tree.reusedHistoryEvents()-prev_reused << ", stashed " << prev_stashed << "+" << backward_tree.numberOfStashed()-prev_stashed
                 << ", relocates " << prev_relocates << "+" << backward_tree.numberOfRelocates()-prev_relocates
                 << "), backward_tree size = " << backward_tree.nnodes() << " (" << backward_tree.size() << " leaves, "
                 << backward_tc.countGhostBranches(backward_tree.root()) << " ghostbranch, " << backward_tc.countUnaryGhosts(backward_tree.root()) << " unaryghosts, "
                 << backward_tc.countBranchingGhosts(backward_tree.root()) << " branchingghost, " << backward_tc.countActive(backward_tree.root()) << " active)" << endl;
            prev_hsize = backward_tree.historySize();
            prev_reused = backward_tree.reusedHistoryEvents();
            prev_stashed = backward_tree.numberOfStashed();
            prev_relocates = backward_tree.numberOfRelocates();
        }
    }

    if (verbose && rewind)
        cerr << "Backward scan done after " << backward_tree.historySize() << " backward history events..." << endl;

    /**
     * Forward loop over the input data
     */
    backward_tc.deFloatAll(backward_tree.root());
    h = 0;
    PointerTree final_tree(inputr->col(0));
    TreeControllerSimple final_tc(final_tree, debug, dotfile.empty() ? "" : "final");
    prev_hsize = 0;
    prev_reused = 0;
    prev_stashed = 0;
    prev_relocates = 0;
    while (h < inputr->size()) // TODO FIXME --skip support
    {
        LeafDistance dist(ic.size(), 0);
        backward_tc.distance(dist, inputr->col(h));
        backward_tc.rewind(inputr->col(h), inputr->size()-h-1);
        
        final_tc.process(inputr->col(h), h, dist);
        if (!dotfile.empty())
            final_tree.outputDOT("final_" + dotfile, h);

        if (verbose && h % debug_interval == 0)
        {
            //cerr << "rewind backward at " << h << "/" << step << ", history = " << backward_tree.historySize() << ", backward_tree size = " << backward_tree.nnodes() << " (" << backward_tree.size() << " leaves, "
            //     << backward_tc.countGhostBranches(backward_tree.root()) << " ghostbranch, " << backward_tc.countUnaryGhosts(backward_tree.root()) << " unaryghosts, "
            //     << backward_tc.countBranchingGhosts(backward_tree.root()) << " branchingghost, " << backward_tc.countActive(backward_tree.root()) << " active)" << endl;
            cerr << "at final step " << h << ", history = " << final_tree.historySize() << " (" << final_tree.historySize()-prev_hsize
                 << ", reused " << prev_reused << "+" << final_tree.reusedHistoryEvents()-prev_reused << ", stashed " << prev_stashed << "+" << final_tree.numberOfStashed()-prev_stashed
                 << ", relocates " << prev_relocates << "+" << final_tree.numberOfRelocates()-prev_relocates
                 << "), final_tree size = " << final_tree.nnodes() << " (" << final_tree.size() << " leaves, "
                 << final_tc.countGhostBranches(final_tree.root()) << " ghostbranch, " << final_tc.countUnaryGhosts(final_tree.root()) << " unaryghosts, "
                 << final_tc.countBranchingGhosts(final_tree.root()) << " branchingghost, " << final_tc.countActive(final_tree.root()) << " active)" << endl;
            prev_hsize = final_tree.historySize();
            prev_reused = final_tree.reusedHistoryEvents();
            prev_stashed = final_tree.numberOfStashed();
            prev_relocates = final_tree.numberOfRelocates();
        }
        h++;
    }

    
    if (verbose)
        cerr << "All done. Finished after " << step << " steps of input. In total " << final_tree.historySize() << " final history events." << endl;
    
    // Clean up
    delete inputr; inputr = 0;
    return 0;
}
