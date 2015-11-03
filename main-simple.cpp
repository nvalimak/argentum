#include "Configuration.h"
#include "InputReader.h"
#include "Tree.h"
#include "TreeControllerSimple.h"
#include <sstream>
#include <cstdlib>
#include <getopt.h>
using namespace std;


/**
 * Extract trees (without extra information)
 */
void extract(InputReader &ir, direction_t direction, PointerTree &tree, TreeControllerSimple &tc, Configuration const &config, string const &dotprefix, bool newick)
{
    unsigned prev_hsize = 0;
    unsigned prev_reused = 0;
    unsigned prev_stashed = 0;
    unsigned prev_relocates = 0;
    unsigned step = 0;
    if (direction == direction_backward)
        step = ir.size();
    
    do
    {
        if (direction == direction_backward)
            --step;

        unsigned increasingStep = direction == direction_forward ? step : ir.size()-step-1;
        tc.process(ir.col(step), increasingStep);

        if (!config.dotfile.empty())
            tree.outputDOT(dotprefix + "_" + config.dotfile, increasingStep);

        if (newick)
            tree.outputNewick("-", ir.position(step));

        if (config.verbose && increasingStep % config.debug_interval == 0)
        {
            cerr << "at step " << increasingStep << "/" << ir.size()
                 << " (history " << prev_hsize << "+" << tree.historySize()-prev_hsize
                 << ", reused " << prev_reused << "+" << tree.reusedHistoryEvents()-prev_reused
                 << ", stashed " << prev_stashed << "+" << tree.numberOfStashed()-prev_stashed
                 << ", relocates " << prev_relocates << "+" << tree.numberOfRelocates()-prev_relocates
                 << "), tree size = " << tree.nnodes() << " (" << tree.size() << " leaves, "
                 << tc.countGhostBranches(tree.root()) << " ghostbranch, " << tc.countUnaryGhosts(tree.root()) << " unaryghosts, "
                 << tc.countBranchingGhosts(tree.root()) << " branchingghost, " << tc.countActive(tree.root()) << " active)" << endl;
            prev_hsize = tree.historySize();
            prev_reused = tree.reusedHistoryEvents();
            prev_stashed = tree.numberOfStashed();
            prev_relocates = tree.numberOfRelocates();
        }

        if (direction == direction_forward)
            ++step;
    } while (0 < step && step < ir.size());
}

/**
 * Extract trees (with given extra information)
 */
void extract(InputReader &ir, direction_t direction, PointerTree &tree, TreeControllerSimple &tc, Configuration const &config, string const &dotprefix, TreeControllerSimple &predictor_tc, bool newick)
{
    unsigned prev_hsize = 0;
    unsigned prev_reused = 0;
    unsigned prev_stashed = 0;
    unsigned prev_relocates = 0;
    LeafDistance dist(ir.front().size(), 0);
    unsigned step = 0;
    if (direction == direction_backward)
        step = ir.size();
    
    do
    {
        if (direction == direction_backward)
            --step;

        unsigned increasingStep = direction == direction_forward ? step : ir.size()-step-1;
        unsigned decreasingStep = ir.size()-increasingStep-1;
        predictor_tc.distance(dist, ir.col(step));
        predictor_tc.rewind(ir.col(step), decreasingStep);
        
        tc.process(ir.col(step), increasingStep, dist);

        if (!config.dotfile.empty())
            tree.outputDOT(dotprefix + "_" + config.dotfile, step);

        if (newick)
            tree.outputNewick("-", ir.position(step));

        if (config.verbose && increasingStep % config.debug_interval == 0)
        {
            cerr << "at step " << increasingStep << "/" << ir.size()
                 << " (history " << prev_hsize << "+" << tree.historySize()-prev_hsize
                 << ", reused " << prev_reused << "+" << tree.reusedHistoryEvents()-prev_reused
                 << ", stashed " << prev_stashed << "+" << tree.numberOfStashed()-prev_stashed
                 << ", relocates " << prev_relocates << "+" << tree.numberOfRelocates()-prev_relocates
                 << "), tree size = " << tree.nnodes() << " (" << tree.size() << " leaves, "
                 << tc.countGhostBranches(tree.root()) << " ghostbranch, " << tc.countUnaryGhosts(tree.root()) << " unaryghosts, "
                 << tc.countBranchingGhosts(tree.root()) << " branchingghost, " << tc.countActive(tree.root()) << " active)" << endl;
            prev_hsize = tree.historySize();
            prev_reused = tree.reusedHistoryEvents();
            prev_stashed = tree.numberOfStashed();
            prev_relocates = tree.numberOfRelocates();
        }

        if (direction == direction_forward)
            ++step;
    } while (0 < step && step < ir.size());
}


int main(int argc, char ** argv)
{
    Configuration config(argc, argv);
    if (!config.good)
        config.print_short_usage();
    if (config.verbose)
        cerr << "Reading input file " << config.inputfile << endl;
    InputReader *inputr = InputReader::build(config.inputformat, config.inputfile, config.nrows);
    if (!inputr->good())
    {
        cerr << "error: could not read input file " << config.inputfile << endl;
        return 1;
    }
    if (config.verbose)
        cerr << "Input size is " << inputr->col(0).size() << " sequences and " << inputr->size() << " sites." << endl;

    /**
     * Forward loop over the input
     */
    PointerTree forward_tree(inputr->front());
    TreeControllerSimple forward_tc(forward_tree, config.debug, config.dotfile.empty() ? "" : "forward");
    extract(*inputr, direction_forward, forward_tree, forward_tc, config, "forward", config.newick);
    cerr << "Forward scan done after " << inputr->size() << " steps of input. In total " << forward_tree.historySize() << " history events." << endl;

    /**
     * Backward loop over the input data
     */
    forward_tc.deFloatAll(forward_tree.root());
    PointerTree backward_tree(inputr->back());
    TreeControllerSimple backward_tc(backward_tree, config.debug, config.dotfile.empty() ? "" : "backward");
    extract(*inputr, direction_backward, backward_tree, backward_tc, config, "backward", forward_tc, false);
    cerr << "Backward scan done after " << inputr->size() << " steps of input. In total " << backward_tree.historySize() << " history events." << endl;

    /**
     * Final loop over the input data
     */
    backward_tc.deFloatAll(backward_tree.root());
    PointerTree final_tree(inputr->front());
    TreeControllerSimple final_tc(final_tree, config.debug, config.dotfile.empty() ? "" : "final");
    extract(*inputr, direction_forward, final_tree, final_tc, config, "final", backward_tc, false);
    cerr << "Final scan done after " << inputr->size() << " steps of input. In total " << final_tree.historySize() << " history events." << endl;
    
    // Clean up
    delete inputr; inputr = 0;
    return 0;
}
