#include "TreeController.h"
#include <iostream>
#include <utility>
#include <set>
#include <map>
#include <cassert>
#include <cstdlib> // std::abort()
using namespace std;

void TreeController::process(InputColumn const &ic, unsigned step_)
{
    step = step_;
    PointerTree::PointerNode *root = t.root();
    reduce(root, ic);
    if (debug && !dotfile.empty())
        t.outputDOT("reduced", step);

    resolveNonBinary(root);
    PointerTree::clearNonBranchingInternalNodes();
    if (debug && !dotfile.empty())
        t.outputDOT("resolvednonbinaryr", step);
    if (debug)
        t.validate();

    recombine.clear();
    findReduced(root, 1);
    updatedThisStep = recombine;
    
    recombine.clear();
    recombine.reserve(2 * t.size());
    pair<int,int> checksum = recombineStrategy(root);
    assert (checksum.second == 1);
        
    
    t.unstash(); // Recover stashed subtrees (inserted to the root)

    for (vector<PointerTree::PointerNode *>::iterator it = updatedThisStep.begin(); it != updatedThisStep.end(); ++it)
        if (!(*it)->leaf())
            (*it)->previousUpdate(step);
    PointerTree::clearNonBranchingInternalNodes();
}


void TreeController::rewind(InputColumn const &ic, unsigned step_)
{
    step = step_;
    PointerTree::PointerNode *root = t.root();
    reduce(root, ic);

    recombine.clear();
    findReduced(root, 1);
    unsigned nones = recombine.size();
    recombine.clear();
    findReduced(root, 0);
    unsigned nzeros = recombine.size();
    if (!dotfile.empty() && debug)
        t.outputDOT("prerewind", step);

    if (nones > 1 && nzeros > 1)
    {
        unsigned rootzeros = 0;
        for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
            if ((*it)->parentPtr() == t.root())
                ++rootzeros;
        recombine.clear();
        findReduced(root, 1);
        unsigned rootones = 0;
        PointerTree::PointerNode *nonrootparent = 0;
        for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
            if ((*it)->parentPtr() == t.root())
                ++rootones;
            else
                nonrootparent = (*it)->parentPtr();

        unsigned nonrparentones = 0;
        for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
            if ((*it)->parentPtr() == nonrootparent)
                ++nonrparentones;
        
        cerr << "assert failed at TreeController::rewind(): step = " << step
             << "    nzeros = " << nzeros << ", nones = " << nones << endl
             << " root zeros = " << rootzeros << ", root ones = " << rootones << endl
             << " nonroot ones = " << nonrparentones << endl;
        
        abort();
    }

    t.rewind(step);
    PointerTree::clearNonBranchingInternalNodes();
    if (!dotfile.empty() && debug)
        t.outputDOT("rewind", step);

    if (debug)
        t.validate();
}

/**
 * Reduce the tree
 *
 * Updates all the 
 *     0-1-count values in the given tree, and
 *     removes tags from binary nodes' children
 */
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
    unsigned floating = 0;
    unsigned nonghosts = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        pair<int, int> tmp = reduce(*it, ic);
        nzeros += tmp.first;
        nones += tmp.second;
        if (tmp.first + tmp.second > 0)
        {
            ++ nonghosts;
            if ((*it)->floating())
                ++ floating;
        }
    }
    pn->nZeros(nzeros);
    pn->nOnes(nones);
    if (!pn->root() && nonghosts <= 2 && floating > 0)
        for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
            (*it)->floating(false);

    return make_pair(nzeros, nones);
}

/**
 * Resolve non-binary nodes locally
 *
 * These modifications are not tracked in history events.
 * Only non-tagged subtrees are considered.
 */
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
        if (!(*it)->floating() && (*it)->reduced() && (*it)->reducedLabel() == 1)
        {
            nones += (*it)->nOnes(); // FIXME If *it is a unary ghost, gets merged with other reduced subtrees
            recombine.push_back(*it);
        }
    if (recombine.size() > 1)
        recombineNonBinarySubtrees(nones, false, true); // Not included in history events; does keep parent counter 
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
    if (!unarypath && pn->reduced() && !pn->root())
    {
        if (pn->reducedLabel() == il)
            recombine.push_back(pn);
        return;
    }

    // Assert: 'pn' is an unary ghost node; or a nonreducible node; or the root node
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        findReduced(*it, il);
}

/**
 * Recombine strategy
 *
 * For each subtree, bottom-up:
 *   1) Relocates all non-consistent 0-reduced branches to the root
 *   2) Relocates all non-consistent 1-reduced branches below the same node
 *
 * Non-floating nodes create a recombination if prune-and-regraphed.
 * A floating node does not create new recombinations as long as it is regraphed
 * under a subtree that is older than the previous prune-operation of the floater.
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

    if (!unarypath && pn->reduced() && !pn->root())
    {
        if (pn->reducedLabel() == 1) // Assert: Reduced internal node that is not a unary ghost node
            return make_pair(0,1);
        else
            return make_pair(1,0);
    }
    
    // Assert: 'pn' is an unary ghost node, or a nonreducible node; or the root node
    int nzeroreduced = 0;
    int nonereduced = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); pn->nodeId() != PointerTree::nonreserved && it != pn->end(); ++it)
    {
        pair<int,int> tmp = recombineStrategy(*it);
        nzeroreduced += tmp.first;
        nonereduced += tmp.second;
    }
    if (pn->nodeId() == PointerTree::nonreserved)
        return make_pair(0,0);
    if (unarypath && !pn->root())
        return make_pair(nzeroreduced, nonereduced);

    // Assert: pn is a nonreducible node
    if (!pn->root() && nzeroreduced < nonereduced)
    {
        // Cut the zero-reduced subtrees and relocate under root
        recombine.clear();
        findReduced(pn, 0);
        for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
            t.stash(*it, step, true, false);

        pn->nZeros(0);
        pn->previousUpdate(step);
        return make_pair(0,1); // This subtree becomes reduced
    }
    if (!pn->root() && nzeroreduced == nonereduced)
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
    recombineSubtrees(pn, true, false);
    return make_pair(nzeroreduced, 1);
}


// Relocates all selected subtrees under new internal node
void TreeController::recombineNonBinarySubtrees(unsigned nones, bool keephistory, bool keepparentcounts)
{
    if (recombine.size() < 2)
        return; // One subtree cannot be relocated to itself

    // Choose the largest, non-floating subtree as destination
    unsigned msize = 0;
    PointerTree::PointerNode *mpn = 0;
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if (!(*it)->floating() && (*it)->nZeros() + (*it)->nOnes() > msize)
        {
            msize = (*it)->nZeros()+(*it)->nOnes();
            mpn = *it;
        }
    if (mpn == 0) // All floating; choose the youngest floater as destination
    {
        msize = 0;
        for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
            if (t.getPreviousEventStep(*it) >= msize)
            {
                msize = t.getPreviousEventStep(*it);
                mpn = *it;
            }
    }
    
    // Count the number of floaters younger than the target subtree
    unsigned dest_updated = mpn->previousUpdate();
    unsigned nfloaters = 0;
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if ((*it)->floating())
        {
            unsigned flt_created = t.getPreviousEventStep(*it);
            if (flt_created > dest_updated)
                nfloaters++;
        }
       
    PointerTree::PointerNode * dest = mpn;
    if (mpn->leaf() || (nfloaters < recombine.size()-1 && !mpn->floating()))
        dest = t.createDest(mpn, step); // Create new internal node if leaf, or not all siblings are floaters and destination is not floater
    if (dest->floating())
        keephistory = true;
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if (*it != mpn)
            t.relocate(*it, dest, step, step, keephistory, keepparentcounts); // Relocate all selected 'recombine' subtrees to destination
    
    dest->nOnes(nones);
    dest->nZeros(0);
    if (!dest->parentPtr()->root() && dest->parentPtr()->size() <= 2)
        for (PointerTree::PointerNode::iterator it = dest->parentPtr()->begin(); it != dest->parentPtr()->end(); ++it)
            (*it)->floating(false);
    if (dest->size() <= 2)
        for (PointerTree::PointerNode::iterator it = dest->begin(); it != dest->end(); ++it)
            (*it)->floating(false);
    if (dest!=mpn && mpn->size() <= 2)
        for (PointerTree::PointerNode::iterator it = mpn->begin(); it != mpn->end(); ++it)
            (*it)->floating(false);
}

// Relocates all selected subtrees to the largest selected subtree
void TreeController::recombineSubtrees(PointerTree::PointerNode *subtree_root, bool keephistory, bool keepparentcounts)
{
    if (recombine.size() < 2)
        return; // One subtree cannot be relocated to itself

    // Choose the largest, non-floating subtree as destination
    unsigned msize = 0;
    PointerTree::PointerNode *mpn = 0;
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if (!(*it)->floating() && (*it)->nZeros()+(*it)->nOnes() > msize)
        {
            msize = (*it)->nZeros()+(*it)->nOnes();
            mpn = *it;
        }

    if (mpn == 0)
        mpn = recombine.front(); // All floating, take first
    
    map<PointerTree::PointerNode *, unsigned> upward_path; // Lowest common ancestor candidates between mpn and rest of the subtrees
    PointerTree::PointerNode *tmp = mpn;
    unsigned largest_age = mpn->previousUpdate();
    while (tmp != subtree_root)
    {
        tmp = tmp->parentPtr();
        if (tmp->previousUpdate() > largest_age)
            largest_age = tmp->previousUpdate();
        upward_path[tmp] = largest_age;
    }
    
    PointerTree::PointerNode * dest = mpn;
    if (mpn->leaf())
        dest = t.createDest(mpn, PointerTree::nohistory);
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if (*it != mpn)
        {
            unsigned lca_age = 0;
            if (upward_path.count((*it)->parentPtr()) == 0)
                (*it)->floating(false); // Not floating for this dest.
            else
                lca_age = upward_path[(*it)->parentPtr()];

            t.relocate(*it, dest, lca_age, step, keephistory, keepparentcounts); // Relocate all selected 'recombine' subtrees to destination
            (*it)->floating(true);
        }
    if (mpn->leaf())
        dest->previousUpdate(step);
}

// Update the float-flags in the whole tree
void TreeController::deFloatAll(PointerTree::PointerNode *pn)
{
    pn->floating(false);
    if (pn->leaf())
        return;

    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        deFloatAll(*it);
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
        return 1;

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

