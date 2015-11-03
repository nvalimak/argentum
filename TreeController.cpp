#include "TreeController.h"
#include <iostream>
#include <utility>
#include <set>
#include <map>
#include <cassert>
#include <cstdlib> // std::abort()
#include <cmath>
using namespace std;

/**
 * Process the next input column (without any additional information)
 */
void TreeController::process(InputColumn const &ic, unsigned step_)
{
    step = step_;
    PointerTree::PointerNode *root = t.root();
    reduce(root, ic);
    if (debug && !dotfile.empty())
        t.outputDOT(dotfile + "_reduced", step);

    resolveNonBinary(root);
    t.clearNonBranchingInternalNodes();
    if (debug && !dotfile.empty())
        t.outputDOT(dotfile + "_resolvednonbinaryr", step);
    recombine.clear();
    findReduced(root, 1);
    updatedThisStep = recombine;
    
    recombine.clear();
    recombine.reserve(2 * t.size());
    pair<int,int> checksum = recombineStrategy(root);
    assert (checksum.second <= 1);
        
    t.unstash(); // Recover stashed subtrees (inserted to the root)

    recombine.clear();
    findReduced(t.root(), 1);
    recombineSubtrees(t.root(), true, false);

    for (vector<PointerTree::PointerNode *>::iterator it = updatedThisStep.begin(); it != updatedThisStep.end(); ++it)
        if (!(*it)->leaf())
            (*it)->previousUpdate(step);
    t.clearNonBranchingInternalNodes();
    if (debug)
        t.validate();
}

/**
 * Process the next input column (with given distance information)
 */
void TreeController::process(InputColumn const &ic, unsigned step_, LeafDistance const &dist)
{
    step = step_;
    PointerTree::PointerNode *root = t.root();
    reduce(root, ic);
    // FIXME
    TreeController::computeMaxDists(t.root());
    
    if (debug && !dotfile.empty())
        t.outputDOT(dotfile + "_reduced", step);

    resolveNonBinary(root);
    t.clearNonBranchingInternalNodes();
    if (debug && !dotfile.empty())
        t.outputDOT(dotfile + "_resolvednonbinaryr", step);
    recombine.clear();
    findReduced(root, 1);
    updatedThisStep = recombine;
    
    recombine.clear();
    recombine.reserve(2 * t.size());
    pair<int,int> checksum = recombineStrategy(root, dist);
    assert (checksum.second <= 1);

    t.unstash(); // Recover stashed subtrees (inserted to the root)
    
    recombine.clear();
    findReduced(t.root(), 1);
    recombineSubtrees(t.root(), true, false, dist);

    for (vector<PointerTree::PointerNode *>::iterator it = updatedThisStep.begin(); it != updatedThisStep.end(); ++it)
        if (!(*it)->leaf())
            (*it)->previousUpdate(step);
    t.clearNonBranchingInternalNodes();
    if (debug)
        t.validate();
}

/*void distByTraversal(PointerTree::PointerNode * pn, PointerTree::PointerNode * src, unsigned sup, unsigned sdown, LeafDistance &d)
{
    if (pn->ghostbranch())
        return;
    if (pn->leaf())
    {
        if (pn->reducedLabel() == 1)
            d[pn->leafId()] = max(sup,sdown);
        return;
    }

    bool unarypath = false;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if ((*it)->nZeros() == pn->nZeros() && (*it)->nOnes() == pn->nOnes())
            unarypath = true;

    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if (*it != src)
            distByTraversal(*it, pn, sup, unarypath ? sdown : sdown+1, d);
            
    if (!pn->root() && pn->parentPtr() != src)
        distByTraversal(pn->parentPtr(), pn, unarypath ? sup : sup+1, sdown, d);
        }*/


// Assert: maxdepth values are set
void TreeController::distanceByTraversal(PointerTree::PointerNode * pn, PointerTree::PointerNode * src, unsigned mdepth, LeafDistance &d)
{
    if (pn->ghostbranch())
        return;
    if (pn->leaf())
    {
        if (pn->reducedLabel() == 1)
            d[pn->leafId()] = mdepth;
        return;
    }

    //bool unarypath = false;
    //for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    //    if ((*it)->nZeros() == pn->nZeros() && (*it)->nOnes() == pn->nOnes())
    //        unarypath = true;

    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if (*it != src)
            TreeController::distanceByTraversal(*it, pn, max(pn->maxDepth(), mdepth), d);
            
    if (!pn->root() && pn->parentPtr() != src)
        TreeController::distanceByTraversal(pn->parentPtr(), pn, mdepth, d);
}

/* depth version */
unsigned TreeController::computeMaxDists(PointerTree::PointerNode * pn)
{
    if (pn->ghostbranch())
        return 0;
    if (pn->leaf())
    {
        if (pn->reducedLabel() == 1)
            return 1;
        return 1;
    }

    bool unarypath = false;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if ((*it)->nZeros() == pn->nZeros() && (*it)->nOnes() == pn->nOnes())
            unarypath = true;

    unsigned maxd = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        unsigned d = computeMaxDists(*it);
        if (maxd < d)
            maxd = d;
    }

    pn->maxDepth(maxd);
    if (unarypath)
        return maxd;
    return (maxd>0 ? maxd+1 : maxd);
}

void TreeController::distance(LeafDistance &dist, InputColumn const & ic)
{
    reduce(t.root(), ic);
    TreeController::distance(t, dist, ic);
}

void TreeController::distance(PointerTree &t, LeafDistance &dist, InputColumn const & ic)
{    
    TreeController::computeMaxDists(t.root());
    
    for (size_t l = 0; l < t.size(); ++l) // TODO make linear time
    {
        LeafDistance d(dist.size(), 0);
        PointerTree::PointerNode * pn = t.leaf(l);
        if (pn->reducedLabel() != 0)
            continue;
        TreeController::distanceByTraversal(pn->parentPtr(), pn, 0, d);
        dist[l] = 0;
        for (size_t j = 0; j < d.size(); ++j)
            if (t.leaf(j)->reducedLabel() == 1)
                dist[l] += d[j]; //exp(-d[j]);
    }
}

void TreeController::rewind(InputColumn const &ic, unsigned step_)
{
    step = step_;
    PointerTree::PointerNode *root = t.root();    
    reduce(root, ic);

    
    if (!dotfile.empty() && debug)
        t.outputDOT(dotfile + "_prerewind", step);
    t.prerewind(step);
    
    recombine.clear();
    findReduced(root, 1);
    unsigned nones = recombine.size();
    recombine.clear();
    findReduced(root, 0);
    unsigned nzeros = recombine.size();
    if (nones > 1 && nzeros > 1)
    {
        // TODO cleanup debug code
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
    t.clearNonBranchingInternalNodes();
    if (!dotfile.empty() && debug)
        t.outputDOT(dotfile + "_rewind", step);

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
        if ((*it)->reduced() && (*it)->reducedLabel() == 1)
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
        updatedThisStep.push_back(pn);
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
    //  recombine.clear();
    //  findReduced(pn, 1);
    //  recombineSubtrees(pn, true, false); // FIXME cleanup done at root
    return make_pair(nzeroreduced, 0);
}


/**
 * Recombine strategy (with given distance information)
 */
pair<int,int> TreeController::recombineStrategy(PointerTree::PointerNode *pn, LeafDistance const &dist)
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
        updatedThisStep.push_back(pn);
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
    //    recombine.clear();
    //    findReduced(pn, 1); // FIXME Cleanup will be recombined at root
    //    recombineSubtrees(pn, true, false, dist);
    return make_pair(nzeroreduced, 0);
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
        if (!(*it)->floating() && (*it)->nZeros() + (*it)->nOnes() > (int)msize)
        {
            msize = (*it)->nZeros()+(*it)->nOnes();
            mpn = *it;
        }
    if (mpn == 0) // All floating; choose the oldest floater as destination
    {
        msize = ~0u;
        for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
            if (t.getPreviousEventStep(*it) <= msize)
            {
                msize = t.getPreviousEventStep(*it);
                mpn = *it;
            }
    }
    
    // Count the number of floaters younger than the target subtree
    assert(mpn != 0);
    unsigned dest_updated = mpn->previousUpdate();
    if (mpn->previousUpdate() == PointerTree::nohistory)
        dest_updated = 0;
    unsigned nfloaters = 0;
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if ((*it)->floating() && *it != mpn)
        {
            unsigned flt_created = t.getPreviousEventStep(*it);
            if (flt_created >= dest_updated)
                nfloaters++;
        }
       
    PointerTree::PointerNode * dest = mpn;
    //if (mpn->leaf() || nfloaters < recombine.size()-1)
        dest = t.createDest(mpn, PointerTree::nohistory); // Create new internal node if leaf, or not all siblings are floaters
        //else
        //cerr << "skipped dest creation for " << mpn->nodeId() << ", floating = " << mpn->floating() << ", leaf = " << mpn->leaf() << endl
        //     << "  nfloaters = " << nfloaters << ", recomb size = " << recombine.size() << endl;
    
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if (*it != mpn)
            t.relocate(*it, dest, step, step, keephistory, keepparentcounts); // Relocate all selected 'recombine' subtrees to destination
    
    dest->nOnes(nones);
    dest->nZeros(0);
    updatedThisStep.push_back(dest);
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
    int msize = 0;
    PointerTree::PointerNode *mpn = 0;
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if (!(*it)->floating() && (*it)->nZeros()+(*it)->nOnes() > msize)
        {
            msize = (*it)->nZeros()+(*it)->nOnes();
            mpn = *it;
        }

    if (mpn == 0)
        mpn = recombine.front(); // All floating, take first
    
    map<PointerTree::PointerNode *, pair<unsigned,bool> > upward_path; // Lowest common ancestor candidates between mpn and rest of the subtrees
    PointerTree::PointerNode *tmp = mpn;
    unsigned largest_age = mpn->previousUpdate();
    bool floating_lca = mpn->floating();
    if (largest_age == PointerTree::nohistory)
        largest_age = 0;
    while (tmp != subtree_root)
    {
        tmp = tmp->parentPtr();
        upward_path[tmp] = make_pair(largest_age, floating_lca); // Largest age not including the LCA
        if (tmp->previousUpdate() > largest_age && tmp->previousUpdate() != PointerTree::nohistory) 
            largest_age = tmp->previousUpdate();
        if (tmp->floating())
            floating_lca = true;
    }

    PointerTree::PointerNode * dest = mpn;
    if (mpn->leaf())
       dest = t.createDest(mpn, PointerTree::nohistory);
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if (*it != mpn)
        {
            pair<unsigned,bool> lca_age(0, false);
            if (upward_path.count((*it)->parentPtr()) == 0)
                (*it)->floating(false); //lca_age = make_pair(largest_age, false);
            else
                lca_age = upward_path[(*it)->parentPtr()];
            if (lca_age.second)
                (*it)->floating(false); // Not floating for this dest.
            
            //,cerr << "relocating node " << (*it)->nodeId() << ", to " << dest->nodeId() << ", lca_age = " << lca_age << endl;
            t.relocate(*it, dest, lca_age.first, step, keephistory, keepparentcounts); // Relocate all selected 'recombine' subtrees to destination
            (*it)->floating(true);
        }
    if (mpn->leaf())
        updatedThisStep.push_back(dest);
}

// Relocates all selected subtrees to the largest selected subtree (minimizing the difference to the given distance values)
void TreeController::recombineSubtrees(PointerTree::PointerNode *subtree_root, bool keephistory, bool keepparentcounts, LeafDistance const &dist)
{
    if (recombine.size() < 2)
        return; // One subtree cannot be relocated to itself

    // Choose the destination subtree
    double eval = -1;
    int msize = 0;
    PointerTree::PointerNode *mpn = 0;
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if (!(*it)->floating())
        {
            double e = evaluateDistance(t, recombine, subtree_root, *it, dist);
            if (printEval)
                cerr << "eval = " << e << " for dest = " << (*it)->nodeId() << ", isleaf = " << (*it)->leaf() << endl;
            if (eval == -1 || e < eval || (e == eval && msize < (*it)->nZeros()+(*it)->nOnes()))
            {
                eval = e;
                msize = (*it)->nZeros()+(*it)->nOnes();
                mpn = *it;
            }
        }
    
    if (mpn == 0)
        mpn = recombine.front(); // All floating, take first
    
    map<PointerTree::PointerNode *, pair<unsigned,bool> > upward_path; // Lowest common ancestor candidates between mpn and rest of the subtrees
    PointerTree::PointerNode *tmp = mpn;
    unsigned largest_age = mpn->previousUpdate();
    bool floating_lca = mpn->floating();
    if (largest_age == PointerTree::nohistory)
        largest_age = 0;
    while (tmp != subtree_root)
    {
        tmp = tmp->parentPtr();
        upward_path[tmp] = make_pair(largest_age, floating_lca); // Largest age not including the LCA
        if (tmp->previousUpdate() > largest_age && tmp->previousUpdate() != PointerTree::nohistory) 
            largest_age = tmp->previousUpdate();
        if (tmp->floating())
            floating_lca = true;
    }

    PointerTree::PointerNode * dest = mpn;
    if (mpn->leaf())
        dest = t.createDest(mpn, PointerTree::nohistory);
    for (vector<PointerTree::PointerNode *>::iterator it = recombine.begin(); it != recombine.end(); ++it)
        if (*it != mpn)
        {
            pair<unsigned,bool> lca_age(0, false);
            if (upward_path.count((*it)->parentPtr()) == 0)
                (*it)->floating(false); // Not floating for this dest.
            else
                lca_age = upward_path[(*it)->parentPtr()];
            if (lca_age.second)
                (*it)->floating(false); // Not floating for this dest.
            
            //,cerr << "relocating node " << (*it)->nodeId() << ", to " << dest->nodeId() << ", lca_age = " << lca_age << endl;
            t.relocate(*it, dest, lca_age.first, step, keephistory, keepparentcounts); // Relocate all selected 'recombine' subtrees to destination
            (*it)->floating(true);
        }
    if (mpn->leaf())
        updatedThisStep.push_back(dest);
}

void TreeController::distanceByTraversal(PointerTree::PointerNode * pn, PointerTree::PointerNode const *src, PointerTree::PointerNode const *dest, set<PointerTree::PointerNode *> &oner, unsigned mdepth, LeafDistance &d)
{
    if (pn == dest)
    {
        unsigned mm = mdepth;
        for (set<PointerTree::PointerNode *>::iterator it = oner.begin(); it != oner.end(); ++it)
            if (mm < (*it)->maxDepth())
                mm = (*it)->maxDepth();
        
        for (set<PointerTree::PointerNode *>::iterator it = oner.begin(); it != oner.end(); ++it)
        {
            if (*it != dest || dest->leaf())
                TreeController::distanceByTraversal(*it, (*it)->parentPtr(), mm+1, d);
            else
                TreeController::distanceByTraversal(*it, (*it)->parentPtr(), mm, d);
        }
        return;
    }
    if (oner.count(pn))
        return; // Subtree not here...
    
    if (pn->ghostbranch())
        return;
    if (pn->leaf())
    {
        if (pn->reducedLabel() == 1)
        {
            assert(false);
            d[pn->leafId()] = mdepth;
        }
        return;
    }

    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if (*it != src)
            TreeController::distanceByTraversal(*it, pn, dest, oner, max(pn->maxDepth(), mdepth), d);

    if (!pn->root() && pn->parentPtr() != src)
        TreeController::distanceByTraversal(pn->parentPtr(), pn, dest, oner, max(pn->maxDepth(), mdepth), d);
}


double TreeController::evaluateDistance(PointerTree &t, std::vector<PointerTree::PointerNode *> const &recombine, PointerTree::PointerNode *subtree_root, PointerTree::PointerNode *dest, LeafDistance const & dist)
{
    set<PointerTree::PointerNode *> oner(recombine.begin(), recombine.end());
        
    LeafDistance distForDest(dist.size(), 0); 
    
    for (size_t l = 0; l < t.size(); ++l) // TODO make linear time
    {
        PointerTree::PointerNode * pn = t.leaf(l);
        if (pn->reducedLabel() != 0)
            continue;
        LeafDistance d(dist.size(), 0);
        TreeController::distanceByTraversal(pn->parentPtr(), pn, dest, oner, 0, d);
        
        distForDest[l] = 0;
        for (size_t j = 0; j < d.size(); ++j)
            if (t.leaf(j)->reducedLabel() == 1)
                distForDest[l] += d[j]; //exp(-d[j]);
    }
    
    double diff = 0;
    for (size_t l = 0; l < t.size(); ++l) // TODO make linear time
    {
        PointerTree::PointerNode * pn = t.leaf(l);
        if (pn->reducedLabel() != 0)
            continue;
        double d = 0;
        if (false)
            cerr << "dist[" << pn->leafId() << "] = " <<  dist[pn->leafId()] << ", distForDest[" << pn->leafId() << ", " << dest->nodeId() << "] = " << distForDest[pn->leafId()] << endl;
        d = distForDest[pn->leafId()] - dist[pn->leafId()];
        diff += abs(d)/max( distForDest[pn->leafId()], dist[pn->leafId()]);
    }
    return diff;
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

