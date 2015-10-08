#include "TreeController.h"
#include <iostream>
#include <cassert>
using namespace std;

void TreeController::process(InputColumn const &ic, unsigned step)
{
    PointerTree::PointerNode &root = t.root();
    reduce(&root, ic);
    if (debug)
        t.outputDOT("reduced", step); // Debug
    cerr << "to resolve nonbinary at step = " << step << endl;
    resolveNonBinary(&root);
    if (debug)
        t.outputDOT("resolvednonbinaryr", step); // Debug

    // FIXME Faster to reduce already during resolveNonBinary()
    reduce(&root, ic);

    recombine.clear();
    recombine.reserve(2 * t.size());
    collectRecombine(&root);
    cerr << "recombine at step = " << step << endl;
    recombineSubtrees(true);
}

// Reduce the tree (updates all the 'reduce' values in the given tree)
// Note: Returns PointerTree::nonreducible if the subtree cannot be reduced.
InputLabel TreeController::reduce(PointerTree::PointerNode *pn, InputColumn const &ic)
{
    assert(!pn->leaf() || pn->leafId() != ~0u);
    pn->setReduced(PointerTree::nonreducible);
    if (pn->leaf())
    {
        // Leaves are always reducable
        pn->setReduced(ic[pn->leafId()]);
        return ic[pn->leafId()];
    }
    if (pn->size() == 0)
    {
        // Ghost leaf node
        assert(pn->numberOfRefs() > 0);
        pn->setReduced(PointerTree::ghostnode);
        return PointerTree::ghostnode;
    }
    
    PointerTree::PointerNode::iterator it = pn->begin();
    InputLabel il = reduce(*it, ic); // Label of first child
    while (++it != pn->end())
    {
        // Process rest of the children
        InputLabel tmp = reduce(*it, ic);
        if (il == PointerTree::ghostnode) // skip ghost nodes
            il = tmp;
        if (il != tmp && tmp != PointerTree::ghostnode)
            il = PointerTree::nonreducible;
    }
    // All children labels are equal to il
    if (il != PointerTree::nonreducible)
        pn->setReduced(il);
    return il;
}

// Resolve non-binary nodes locally
// These modifications are not tracked in history events
void TreeController::resolveNonBinary(PointerTree::PointerNode *pn)
{
    if (pn->leaf())
        return;
    if (pn->reduced())
        return;

    // Depth first
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        resolveNonBinary(*it);

    // Modify local structure
    recombine.clear();
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if ((*it)->reducedLabel() == 1)
            recombine.push_back(*it);
    if (recombine.size() > 1)
        recombineSubtrees(false); // Not included in history events
}

// Collects 'recombine' to contain all reduced branches whose symbol == '1'
void TreeController::collectRecombine(PointerTree::PointerNode *pn)
{
    if (pn->reduced())
    {
        if (pn->reducedLabel() == 1)
            recombine.push_back(pn);
        return;
    }
    
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        collectRecombine(*it);
}

// Relocates all selected subtrees to the first selected subtree
void TreeController::recombineSubtrees(bool keephistory)
{
    cerr << "recombine size = " << recombine.size() << endl;
    if (recombine.size() < 2)
        return; // One subtree cannot be relocated to itself
    
    vector<PointerTree::PointerNode *>::iterator it = recombine.begin();
    PointerTree::PointerNode * dest = *it; // Choose first subree as the destination
    dest = t.createDest(dest);
    while (++it != recombine.end())
        t.relocate(*it, dest, keephistory); // Relocate all selected 'recombine' subtrees to destination
}
