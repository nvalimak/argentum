#include "Configuration.h"
#include "InputReader.h"
#include "Tree.h"
#include "TreeController.h"
#include "TreeDistance.h"
#include "TreeEnumerator.h"
#include <sstream>
#include <cstdlib>
#include <getopt.h>
using namespace std;

/**
 * Extract trees (without extra information)
 */
unsigned extract(InputReader &ir, direction_t direction, PointerTree &tree, TreeController &tc, Configuration const &config, string const &dotprefix)
{
    unsigned total_intnodes = 0;
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
        tc.assignLabels(ir.col(step));
        tc.process(ir.col(step), increasingStep);

        if (!config.dotfile.empty())
            tree.outputDOT(dotprefix + "_" + config.dotfile, increasingStep);
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
        total_intnodes += tc.countInternalNodes();
        if (direction == direction_forward)
            ++step;
    } while (0 < step && step < ir.size());
    return total_intnodes;
}

/**
 * Extract trees (with given rewind tree)
 */
unsigned extract(InputReader &ir, direction_t direction, PointerTree &tree, TreeController &tc, Configuration const &config, string const &dotprefix, TreeController &predictor_tc, PointerTree &predictor_tree, TreeDistance &predictor_dist, bool newick)
{
    unsigned total_intnodes = 0;
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
        unsigned decreasingStep = ir.size()-increasingStep-1;
        tc.assignLabels(ir.col(step));
        predictor_tc.assignLabels(ir.col(step));
        if (newick)
            predictor_tree.outputNewick("-", ir.position(step));
        predictor_tc.rewind(ir.col(step), decreasingStep, 0);
        tc.process(ir.col(step), increasingStep, predictor_dist);

        if (!config.dotfile.empty())
            tree.outputDOT(dotprefix + "_" + config.dotfile, increasingStep);
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
        total_intnodes += tc.countInternalNodes();
        if (direction == direction_forward)
            ++step;        
    } while (0 < step && step < ir.size());
    return total_intnodes;
}

/**
 * Extract trees (with given Newick tree)
 */
unsigned extract(InputReader &ir, direction_t direction, PointerTree &tree, TreeController &tc, Configuration const &config, string const &dotprefix, NewickTree &predictor_tree, NewickDistance &predictor_dist)
{
    unsigned total_intnodes = 0;
    unsigned prev_hsize = 0;
    unsigned prev_reused = 0;
    unsigned prev_stashed = 0;
    unsigned prev_relocates = 0;
    unsigned step = 0;
    unsigned pred_base = 0;
    unsigned cur_base = 0;    
    if (direction == direction_backward)
        step = ir.size();
    do
    {
        if (direction == direction_backward)
            --step;

        unsigned increasingStep = direction == direction_forward ? step : ir.size()-step-1;

        /**
         * debug print 
         *
        cerr << " col at step " << step << ", predt_base " << pred_base << ", cur_base " << cur_base << ": ";
        for (unsigned i = 0; i < ir.col(step).size(); ++i)
            cerr << (int)ir.col(step)[i];
        cerr << endl;
        */
        cur_base += ir.position(step);
        while (cur_base > pred_base)
        {
            predictor_tree.next();
            pred_base += predictor_tree.validForBases();
        }
        tc.assignLabels(ir.col(step));
        predictor_tree.assignLabels(ir.col(step));
        tc.process(ir.col(step), increasingStep, predictor_dist);

        if (!config.dotfile.empty())
            tree.outputDOT(dotprefix + "_" + config.dotfile, increasingStep);
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
        total_intnodes += tc.countInternalNodes();
        if (direction == direction_forward)
            ++step;
    } while (0 < step && step < ir.size());
    return total_intnodes;
}

/**
 * Rewind and extract Newick trees
 */
unsigned output_newick(InputReader &ir, direction_t direction, PointerTree &tree, TreeController &tc)
{
    unsigned total_intnodes = 0;
    unsigned step = 0;
    if (direction == direction_backward)
        step = ir.size();
    do
    {
        if (direction == direction_backward)
            --step;

        unsigned increasingStep = direction == direction_forward ? step : ir.size()-step-1;
        unsigned decreasingStep = ir.size()-increasingStep-1;
        tc.assignLabels(ir.col(step));
        tree.outputNewick("-", ir.position(step));
        tc.rewind(ir.col(step), decreasingStep, 0);

        total_intnodes += tc.countInternalNodes();
        if (direction == direction_forward)
            ++step;
    } while (0 < step && step < ir.size());
    return total_intnodes;
}

/**
 * Rewind and enumerate parent-child ranges (--enumerate)
 */
unsigned output_parent_child(InputReader &ir, direction_t direction, PointerTree &tree, TreeController &tc, TreeEnumerator &te)
{
    unsigned total_intnodes = 0;
    unsigned step = 0;
    if (direction == direction_backward)
        step = ir.size();
    do
    {
        if (direction == direction_backward)
            --step;
        
        unsigned increasingStep = direction == direction_forward ? step : ir.size()-step-1;
        unsigned decreasingStep = ir.size()-increasingStep-1;
        tc.assignLabels(ir.col(step));
        tc.rewind(ir.col(step), decreasingStep, &te);

        total_intnodes += tc.countInternalNodes();
        if (direction == direction_forward)
            ++step;
    } while (0 < step && step < ir.size());
    return total_intnodes;
}

void backward_forward(Configuration &config, InputReader &inputr)
{
    /**
     * Backward loop over the input data
     */
    PointerTree backward_tree(inputr.back());
    TreeController backward_tc(backward_tree, config.debug, config.dotfile.empty() ? "" : "backward");
    extract(inputr, direction_backward, backward_tree, backward_tc, config, "backward");
    cerr << "Backward scan done after " << inputr.size() << " steps of input. In total " << backward_tree.historySize() << " history events." << endl;

    /**
     * Forward loop over the input
     */
    backward_tc.deTagAll(backward_tree.root());
    PointerTree forward_tree(inputr.front());
    PointerTreeDistance predictor_dist(backward_tree, forward_tree);
    TreeController forward_tc(forward_tree, config.debug, config.dotfile.empty() ? "" : "forward");
    extract(inputr, direction_forward, forward_tree, forward_tc, config, "forward", backward_tc, backward_tree, predictor_dist, config.newick);
    cerr << "Forward scan done after " << inputr.size() << " steps of input. In total " << forward_tree.historySize() << " history events." << endl;
}

void scrm_forward(Configuration &config, InputReader &inputr)
{
    /**
     * Init SCRM tree
     */
    NewickTree newick_tree(config.inputfile);
    assert(newick_tree.nleaves() == inputr.col(0).size());
    
    /**
     * Forward loop over the input
     */
    PointerTree forward_tree(inputr.front());
    NewickDistance predictor_dist(newick_tree, forward_tree);
    TreeController forward_tc(forward_tree, config.debug, config.dotfile.empty() ? "" : "forward");
    unsigned intnodes = extract(inputr, direction_forward, forward_tree, forward_tc, config, "forward", newick_tree, predictor_dist);
    cerr << "Forward scan done after " << inputr.size() << " steps of input. In total " << forward_tree.historySize() << " history events." << endl;
    cerr << "Avg. internal nodes: " << (double)intnodes/inputr.size()
         << ", reused " << forward_tree.reusedHistoryEvents()
         << ", reused from root " << forward_tree.reusedRootHistoryEvents()
         << ", stashed " << forward_tree.numberOfStashed()
         << ", relocates " << forward_tree.numberOfRelocates()
         << ", one cuts " << forward_tc.numberOfOneCuts() << endl;
        
    /**
     * Output newick with rewind
     */
    if (config.newick)        
    {
        cerr << "warning: output in backward order, make sure to pipe it through `tac` or `tail -r` to reverse the line order!" << endl;
        unsigned intnodes = output_newick(inputr, direction_backward, forward_tree, forward_tc);
        cerr << "Avg. internal nodes: " << (double)intnodes/inputr.size() << endl;
    }
}

void nopred_forward(Configuration &config, InputReader &inputr)
{
    /**
     * Forward loop over the input
     */
    PointerTree forward_tree(inputr.front());
    TreeController forward_tc(forward_tree, config.debug, config.dotfile.empty() ? "" : "forward");
    unsigned intnodes = extract(inputr, direction_forward, forward_tree, forward_tc, config, "forward");
    cerr << "Forward scan done after " << inputr.size() << " steps of input. In total " << forward_tree.historySize() << " history events." << endl;
    cerr << "Avg. internal nodes: " << (double)intnodes/inputr.size()
         << ", reused " << forward_tree.reusedHistoryEvents()
         << ", reused from root " << forward_tree.reusedRootHistoryEvents()
         << ", stashed " << forward_tree.numberOfStashed()
         << ", relocates " << forward_tree.numberOfRelocates()
         << ", one cuts " << forward_tc.numberOfOneCuts() << endl;
    
    /**
     * Output newick with rewind
     */
    if (config.newick)
    {
        cerr << "warning: output in backward order, make sure to pipe it through `tac` or `tail -r` to reverse the line order!" << endl;
        unsigned intnodes = output_newick(inputr, direction_backward, forward_tree, forward_tc);
        cerr << "Avg. internal nodes: " << (double)intnodes/inputr.size() << endl;
    }
    /**
     * Enumerate parent-child ranges
     */
    else if (config.enumerate)
    {
        cerr << "Gathering parent-child ranges (--enumerate)..." << endl;
        TreeEnumerator te;
        unsigned intnodes = output_parent_child(inputr, direction_backward, forward_tree, forward_tc, te);
        cerr << "Outputting parent-child ranges (--enumerate)..." << endl;
        te.output();
        cerr << "Avg. internal nodes: " << (double)intnodes/inputr.size() << endl;
    }
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

    //backward_forward(config, *inputr);
    if (config.scrm_prediction)
        scrm_forward(config, *inputr);
    else if (config.no_prediction)
        nopred_forward(config, *inputr);

    // Clean up
    delete inputr; inputr = 0;
    return 0;
}
