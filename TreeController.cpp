#include "TreeController.h"
#include <iostream>
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

    recombine.clear();
    recombine.reserve(2 * t.size());
    collectRecombine(&root);
    recombineSubtrees(true);

    reduceGhosts(&root);
}

// Reduce the tree (updates all the 'reduce' values in the given tree)
// Note: Returns PointerTree::nonreducible if the subtree cannot be reduced.
// Ghost nodes inherit the child node's status
InputLabel TreeController::reduce(PointerTree::PointerNode *pn, InputColumn const &ic)
{
    assert(!pn->leaf() || pn->leafId() != ~0u);
    if (pn->ghostbranch())
        return PointerTree::ghostbranch;

    if (pn->hasShortcut())
        return reduce(pn->descendantShortcut(), ic);
    
    pn->setReduced(PointerTree::nonreducible);
    if (pn->leaf())
    {
        // Leaves are always reducable
        pn->setReduced(ic[pn->leafId()]);
        return ic[pn->leafId()];
    }
    if (pn->size() == 0)
    {
        // Ghost branch
        assert(pn->numberOfRefs() > 0);
        pn->setReduced(PointerTree::ghostbranch);
        return PointerTree::ghostbranch;
    }

    unsigned nnonghostbranch = 0;
    PointerTree::PointerNode::iterator it = pn->begin();
    InputLabel il = reduce(*it, ic); // Label of first child
    if (il != PointerTree::ghostbranch)
        nnonghostbranch ++;
    while (++it != pn->end())
    {
        // Process rest of the children
        InputLabel tmp = reduce(*it, ic);
        if (tmp != PointerTree::ghostbranch)
            nnonghostbranch ++;
        if (il == PointerTree::ghostbranch)
            il = tmp;
        if (il != tmp && tmp != PointerTree::ghostbranch)
            il = PointerTree::nonreducible;
    }
    // All nonghost children labels are equal to il (or il == nonreducible)
    pn->setReduced(il);
    if (nnonghostbranch == 1)
    {
        // Ghost node
        pn->setReduced(PointerTree::ghostnode);
        return PointerTree::ghostnode;
    }
    return il;
}


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
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if ((*it)->reducedLabel() == 1)
            recombine.push_back(*it);
    if (recombine.size() > 1)
    {
        recombineSubtrees(false); // Not included in history events
        (*(recombine.begin()))->parent()->setReduced(1); // The new parent must be a reduced node
    }
}

// Collects 'recombine' to contain all reduced branches whose symbol == '1'
void TreeController::collectRecombine(PointerTree::PointerNode *pn)
{
    if (pn->ghostbranch() || pn->reduced())
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
    
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        collectRecombine(*it);
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
            if ((*it)->reducedLabel() != PointerTree::ghostnode && (*it)->reducedLabel() != PointerTree::ghostbranch)
                nonghosts++;
            g += countBranchingGhosts(*it);
        }
        if (nonghosts==1)
            g++;
    }
    return g;
}


// Relocates all selected subtrees to the first selected subtree
void TreeController::recombineSubtrees(bool keephistory)
{
    if (recombine.size() < 2)
        return; // One subtree cannot be relocated to itself
    
    vector<PointerTree::PointerNode *>::iterator it = recombine.begin();
    PointerTree::PointerNode * dest = *it; // Choose first subree as the destination
    dest = t.createDest(dest);
    while (++it != recombine.end())
        t.relocate(*it, dest, keephistory); // Relocate all selected 'recombine' subtrees to destination
}
