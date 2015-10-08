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
    recombine.clear();
    recombine.reserve(2 * t.size());
    collectRecombine(&root);
    recombineSubtrees();
}

// Reduce the tree (updates all the 'reduce' values in the given tree)
// Note: Returns PointerTree::non_reducable if the subtree cannot be reduced.
InputLabel TreeController::reduce(PointerTree::PointerNode *pn, InputColumn const &ic)
{
    assert(!pn->leaf() || pn->leafId() != ~0u);
    pn->setReduced(false, PointerTree::nonreducible);
    if (pn->leaf())
    {
        // Leaves are always reducable
        pn->setReduced(true, ic[pn->leafId()]);
        return ic[pn->leafId()];
    }

    assert(pn->size() >= 2); // Assume at least two child nodes
    
    PointerTree::PointerNode::iterator it = pn->begin();
    InputLabel il = reduce(*it, ic); // Label of first child
    while (++it != pn->end())
    {
        // Process rest of the children
        InputLabel tmp = reduce(*it, ic);
        if (il != tmp)
            il = PointerTree::nonreducible;
    }
    // All children labels are equal to il
    if (il != PointerTree::nonreducible)
        pn->setReduced(true, il);
    return il;
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
void TreeController::recombineSubtrees()
{
    cerr << "recombine size = " << recombine.size() << endl;
    if (recombine.size() < 2)
        return; // One subtree cannot be relocated to itself
    
    vector<PointerTree::PointerNode *>::iterator it = recombine.begin();
    PointerTree::PointerNode * dest = *it; // Choose first subree as the destination
    dest = t.createDest(dest);
    while (++it != recombine.end())
        t.relocate(*it, dest); // Relocate all selected 'recombine' subtrees to destination
}
