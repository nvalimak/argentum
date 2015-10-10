#include "TreeController.h"
#include <iostream>
#include <utility>
#include <cassert>
using namespace std;

void TreeController::process(InputColumn const &ic, unsigned step)
{
    PointerTree::PointerNode &root = t.root();
    reduce(&root, ic);
    if (debug && !dotfile.empty())
        t.outputDOT("reduced", step);
    resolveNonBinary(&root);
    if (debug && !dotfile.empty())
        t.outputDOT("resolvednonbinaryr", step);

    if (debug)
        t.validate();

    recombine.clear();
    recombine.reserve(2 * t.size());
    collectRecombine(&root);
    recombineSubtrees(true, false); // Keep history; do not keep parent counts
    // Note: Recombine can mess up 0-1-count values in the upward path FIXME
    //if (debug)
    //    t.validate(); // Does not pass assert FIXME
}

// Reduce the tree (updates all the 0-1-count values in the given tree)
pair<unsigned, unsigned> TreeController::reduce(PointerTree::PointerNode *pn, InputColumn const &ic)
{
    assert(!pn->leaf() || pn->leafId() != ~0u);
    if (pn->ghostbranch())
        return make_pair(0,0);

    if (pn->hasShortcut())
        return reduce(pn->descendantShortcut(), ic);
    
    if (pn->leaf())
    {
        // Sets the new column leaf label
        pn->setLabel(ic[pn->leafId()]);
        return make_pair(pn->nZeros(), pn->nOnes());
    }
    if (pn->size() == 0)
    {
        // New ghost branch
        assert(pn->numberOfRefs() > 0);
        assert(pn->nZeros() > 0 || pn->nOnes() > 0);
        pn->nZeros(0); // Becomes a ghost
        pn->nOnes(0);
        return make_pair(0,0);
    }

    unsigned nzeros = 0;
    unsigned nones = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        pair<unsigned, unsigned> tmp = reduce(*it, ic);
        nzeros += tmp.first;
        nones += tmp.second;
    }
    pn->nZeros(nzeros);
    pn->nOnes(nones);
    return make_pair(nzeros, nones);
}

// Resolve non-binary nodes locally
// These modifications are not tracked in history events
void TreeController::resolveNonBinary(PointerTree::PointerNode *pn)
{
    if (pn->leaf() || pn->reduced() || pn->ghostbranch())
        return;

    if (pn->hasShortcut())
    {
        resolveNonBinary(pn->descendantShortcut());
        return;
    }

    // Depth first
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        resolveNonBinary(*it);

    // Modify local structure
    recombine.clear();
    unsigned nones = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if ((*it)->reduced() && (*it)->reducedLabel() == 1)
        {
            nones += (*it)->nOnes();
            recombine.push_back(*it);
        }
    if (recombine.size() > 1)
    {
        recombineSubtrees(false, true); // Not included in history events; does keep parent counter 
        (*(recombine.begin()))->parent()->nOnes(nones); // The new parent must be a reduced node
        (*(recombine.begin()))->parent()->nZeros(0);    // with 'nones' number of ones.
    }
}

// Collects 'recombine' to contain all reduced branches whose symbol == '1'
void TreeController::collectRecombine(PointerTree::PointerNode *pn)
{
    if (pn->ghostbranch())
        return;
    if (pn->leaf())
    {
        if (pn->reducedLabel() == 1)
            recombine.push_back(pn);
        return;
    }

    if (pn->hasShortcut())
    {
        collectRecombine(pn->descendantShortcut());
        return;
    }

    bool unarypath = false;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        if ((*it)->nZeros() == pn->nZeros() && (*it)->nOnes() == pn->nOnes())
            unarypath = true;
    }
    if (!unarypath && pn->reduced())
    {
        if (pn->reducedLabel() == 1)
            recombine.push_back(pn);
        return;
    }

    // Assert: 'pn' is an unary ghost node, or a nonreducible node
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        collectRecombine(*it);
}

// Relocates all selected subtrees to the first selected subtree
void TreeController::recombineSubtrees(bool keephistory, bool keepparentcounts)
{
    if (recombine.size() < 2)
        return; // One subtree cannot be relocated to itself
    
    vector<PointerTree::PointerNode *>::iterator it = recombine.begin();
    PointerTree::PointerNode * dest = *it; // Choose first subree as the destination
    dest = t.createDest(dest);
    while (++it != recombine.end())
        t.relocate(*it, dest, keephistory, keepparentcounts); // Relocate all selected 'recombine' subtrees to destination
}

// Counts the number of active nodes
unsigned TreeController::countActive(PointerTree::PointerNode *pn)
{
    unsigned g = 0;
    if (pn->leaf())
        return 0;
    if (pn->ghostbranch())
        return 0;

    if (pn->hasShortcut())
        return countActive(pn->descendantShortcut());
    g++;    
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        g += countActive(*it);
    return g;
}


// Counts the number of ghostbranch nodes (1 or 0 children)
unsigned TreeController::countGhostBranches(PointerTree::PointerNode *pn)
{
    unsigned g = 0;
    if (pn->leaf())
        return 0;
    if (pn->ghostbranch())
        g += 1;

    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        g += countGhostBranches(*it);
    return g;
}

// Counts the number of non-ghostbranch unary-nodes (1 or 0 children)
unsigned TreeController::countUnaryGhosts(PointerTree::PointerNode *pn)
{
    unsigned g = 0;
    if (pn->leaf())
        return 0;
    if (pn->ghostbranch())
        return 0;

    if (pn->size() == 1)
        g++;

    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        g += countUnaryGhosts(*it);
    return g;
}

// Counts the number of non-ghostbranch unary-nodes (1 or 0 children)
unsigned TreeController::countBranchingGhosts(PointerTree::PointerNode *pn)
{
    unsigned g = 0;
    if (pn->leaf())
        return 0;
    if (pn->ghostbranch())
        return 0;

    if (pn->size() == 1)
        for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
            g += countBranchingGhosts(*it);
    else
    {
        unsigned nonghosts = 0;
        for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        {
            if (!(*it)->ghostbranch())
                nonghosts++;
            g += countBranchingGhosts(*it);
        }
        if (nonghosts==1)
            g++;
    }
    return g;
}

/**
 * REMOVE
InputLabel TreeController::reduceGhosts(PointerTree::PointerNode *pn)
{
    assert(!pn->leaf() || pn->leafId() != ~0u);
    if (pn->ghostbranch())
        return PointerTree::ghostbranch;
    
    if (pn->leaf())
        return 0;

    if (pn->size() == 0)
    {
        // Ghost leaf node
        //cerr << "found ghost branch" << endl;
        assert(pn->numberOfRefs() > 0);
        pn->setReduced(PointerTree::ghostbranch);
        return PointerTree::ghostbranch;
    }
    
    unsigned nnonghostbranch = 0;
    PointerTree::PointerNode *nonghost_chld = 0;

    PointerTree::PointerNode::iterator it = pn->begin();
    InputLabel il = reduceGhosts(*it); // Label of first child
    if (il != PointerTree::ghostbranch)
    {
        nonghost_chld = *it;
        ++nnonghostbranch;
    }
    
    while (++it != pn->end())
    {
        // Process rest of the children
        InputLabel tmp = reduceGhosts(*it);
        if (tmp != PointerTree::ghostbranch)
        {
            nonghost_chld = *it;
            ++nnonghostbranch;
        }
        if (il != tmp)
            il = PointerTree::nonreducible;
    }
    // All children labels are ghostbranch
    if (il == PointerTree::ghostbranch)
        pn->setReduced(il);
    // Check if there's only one "real" descendant.
    else if (nnonghostbranch == 1)
    {
        // Update the descendant shortcut
        if (nonghost_chld->hasShortcut())
            pn->descendantShortcut(nonghost_chld->descendantShortcut());
        else
            if (nonghost_chld->size() < 2 && !nonghost_chld->leaf())
                pn->descendantShortcut(nonghost_chld);
    }
        
    return il;
    } REMOVE */

