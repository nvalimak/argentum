#include "TreeController.h"
#include <iostream>
#include <utility>
#include <cassert>
using namespace std;

void TreeController::process(InputColumn const &ic, unsigned step)
{
    PointerTree::PointerNode *root = t.root();
    reduce(root, ic);
    if (debug && !dotfile.empty())
        t.outputDOT("reduced", step);
    resolveNonBinary(root);
    if (debug && !dotfile.empty())
        t.outputDOT("resolvednonbinaryr", step);

    if (debug)
        t.validate();

    recombine.clear();
    recombine.reserve(2 * t.size());
    recombineStrategy(root);

    t.unstash(); // Recover stashed subtrees (inserted to the root)
    if (debug)
        t.validate();
}

// Reduce the tree (updates all the 0-1-count values in the given tree)
pair<int, int> TreeController::reduce(PointerTree::PointerNode *pn, InputColumn const &ic)
{
    assert(!pn->leaf() || pn->leafId() != ~0u);
    if (pn->ghostbranch())
        return make_pair(0,0);

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

    int nzeros = 0;
    int nones = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        pair<int, int> tmp = reduce(*it, ic);
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

    // Depth first
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        resolveNonBinary(*it);

    // Modify local structure
    recombine.clear();
    int nones = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if ((*it)->reduced() && (*it)->reducedLabel() == 1)
        {
            nones += (*it)->nOnes(); // FIXME If *it is a unary ghost, gets merged with other reduced subtrees
            recombine.push_back(*it);
        }
    if (recombine.size() > 1)
    {
        recombineSubtrees(false, true); // Not included in history events; does keep parent counter 
        (*(recombine.begin()))->parentPtr()->nOnes(nones); // The new parent must be a reduced node
        (*(recombine.begin()))->parentPtr()->nZeros(0);    // with 'nones' number of ones.
    }
}

void TreeController::findReduced(PointerTree::PointerNode *pn, InputLabel il)
{
    if (pn->ghostbranch())
        return;
    if (pn->leaf())
    {
        if (pn->reducedLabel() == il)
            recombine.push_back(pn);
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
        if (pn->reducedLabel() == il)
            recombine.push_back(pn);
        return;
    }

    // Assert: 'pn' is an unary ghost node, or a nonreducible node
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        findReduced(*it, il);
}

/**
 * Recombine strategy
 *
 * For each subtree:
 *   1) Relocates all non-consistent 0-reduced branches to the root
 *   2) Relocates all non-consistent 1-reduced branches below the same node
 */
pair<int,int> TreeController::recombineStrategy(PointerTree::PointerNode *pn)
{
    if (pn->ghostbranch())
        return make_pair(0,0);       // Assert: Ghost branch, ignored here
    if (pn->leaf())
    {
        if (pn->reducedLabel() == 1) // Assert: Leaf node
            return make_pair(0,1);
        else
            return make_pair(1,0);
    }
    
    bool unarypath = false;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if ((*it)->nZeros() == pn->nZeros() && (*it)->nOnes() == pn->nOnes())
            unarypath = true;

    if (!unarypath && pn->reduced())
    {
        if (pn->reducedLabel() == 1) // Assert: Reduced internal node that is not a unary ghost node
            return make_pair(0,1);
        else
            return make_pair(1,0);
    }
    
    // Assert: 'pn' is an unary ghost node, or a nonreducible node
    int nzeroreduced = 0;
    int nonereduced = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        pair<int,int> tmp = recombineStrategy(*it);
        nzeroreduced += tmp.first;
        nonereduced += tmp.second;
    }       
    if (unarypath)
        return make_pair(nzeroreduced, nonereduced);

    // Assert: pn is a nonreducible node
    if (nzeroreduced < nonereduced)
    {
        // Cut the zero-reduced subtrees and relocate under root
        recombine.clear();
        findReduced(pn, 0);
        for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
            t.stash(*it, true, false);

        pn->nZeros(0);
        return make_pair(0,1); // This subtree becomes reduced
    }
    if (nzeroreduced == nonereduced)
    {
        // Balanced tree; no operation at this step
        return make_pair(nzeroreduced, nonereduced);
    }
    // Assert: nonereduced < nzeroreduced
    if (nonereduced == 1)
        return make_pair(nzeroreduced, nonereduced);
    
    // Assert: nonereduced >= 2 && nonereduced < nzeroreduced
    recombine.clear();
    findReduced(pn, 1);
    recombineSubtrees(true, false);
    return make_pair(nzeroreduced, 1);
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

