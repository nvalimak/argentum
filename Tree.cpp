#include "Tree.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
using namespace std;

// Default flag values
const LeafId PointerTree::nonleaf = (LeafId)~0u;
const unsigned PointerTree::nohistory = ~0u;
const unsigned PointerTree::unknown = ~0u;
size_t PointerTree::N = 0; // Total number of nodes

PointerTree::PointerTree(InputColumn const &ic)
    : r(), n(ic.size()), history(), validationReachable()
{
    history.reserve(HISTORY_INIT_SIZE);
    LeafId id = 0;
    for (InputColumn::const_iterator it = ic.begin(); it != ic.end(); ++it)
        r.insert(new PointerNode(id++, 0.0, &r));

    PointerTree::N = ic.size()+1; // Initial size of the tree
    validationReachable.resize(n);
}

/**
 * Create a new internal node above the given node
 */
PointerTree::PointerNode * PointerTree::createDest(PointerNode *pn)
{
    ++N;
    PointerNode *parent = pn->parent();
    assert (parent != 0);
    parent->erase(pn);
    TreeDepth d = (parent->depth() - pn->depth())/2;
    PointerNode * dest = new PointerNode(d, parent);
    parent->insert(dest);
    dest->insert(pn);
    pn->setParent(dest);
    dest->nZeros(pn->nZeros());
    dest->nOnes(pn->nOnes());
    return dest;
}

/**
 * Truncate the given internal node if there are no active references to it.
 */
void PointerTree::clearNonBranchingInternalNode(PointerNode *pn)
{
    assert(pn->size() < 2);
    assert(!pn->root());
    if (pn->numberOfRefs() > 0)
    {
        // Cannot be removed since one or more references appear in history
        if (pn->size() == 1)
        {
            // Update the descendant shortcut
            PointerNode *only_chld = *(pn->begin());
            if (only_chld->hasShortcut())
                pn->descendantShortcut(only_chld->descendantShortcut());
            else
                if (only_chld->size() < 2 && !only_chld->leaf())
                    pn->descendantShortcut(only_chld);
        }
        return;
    }
    
    // Safe to truncate 'pn' out from the tree
    PointerNode *parent = pn->parent();
    parent->erase(pn);
    if (pn->size() == 1)
    {
        // Rewire the only child of pn:
        PointerNode *only_chld = *(pn->begin());
        pn->erase(only_chld);      // Release the only child from pn
        parent->insert(only_chld); // Rewire pn's parent to point to the only child
        only_chld->setParent(parent);
    }
    // Now safe to delete 'pn' out from the tree
    PointerTree::N--;
    assert(pn->size() == 0); // pn does not own any objects
    delete pn;
}

/**
 * Relocate subtree 'pn' as a new child of 'dest'.
 * Adds an event into the history queue.
 * Cleans up the source subtree.
 */
void PointerTree::relocate(PointerNode *pn, PointerNode *dest, bool keephistory, bool keepparentcounts)
{
    PointerNode *src = pn->parent();
    src->erase(pn);
    if (!keepparentcounts)
    {
        src->addZeros(-(pn->nZeros()));
        src->addOnes(-(pn->nOnes()));
    }

    // Update history event queue
    if (keephistory)
    {
//              pn->previousEvent(history.size());
// FIXME        history.push_back(Event(src, pn));
    }
    
    // Clean source subtree if needed
    if (src->size() == 1 && !src->root())
        PointerTree::clearNonBranchingInternalNode(src);

    dest->insert(pn);
    pn->setParent(dest);
    dest->addZeros(pn->nZeros());
    dest->addOnes(pn->nOnes());
}

/**
 * Rewind the given event
 */
void PointerTree::rewind(Event &e)
{
    PointerNode *dest = e.getSource();
    PointerNode *pn = e.getNode();
    PointerNode *src = pn->parent();
    dest->insert(pn);
    dest->addZeros(pn->nZeros());
    dest->addOnes(pn->nOnes());
    pn->setParent(dest);

    // Clean source subtree if needed
    src->erase(pn);
    if (src->size() == 1 && !src->root())
        PointerTree::clearNonBranchingInternalNode(src);

    e.rewind(); // May delete src (truncated if pn is the only child)
}

/**
 * Validate the integrity of the tree structure
 */
pair<unsigned,unsigned> validate(PointerTree::PointerNode *pn, PointerTree::PointerNode *p, vector<bool> &validationReachable)
{
    assert (pn->nZeros() != PointerTree::unknown);
    assert (pn->nOnes() != PointerTree::unknown);
    assert (p == pn->parent());
    if (pn->leaf())
    {
        assert (pn->leafId() != PointerTree::nonleaf);
        assert (pn->reduced());
        assert (pn->numberOfRefs() == 0);
	validationReachable[pn->leafId()] = true;
        if (pn->reducedLabel() == 1)
            return make_pair(0,1);
        else
            return make_pair(1,0);
    }
    if (pn->hasShortcut())
        return validate(pn->descendantShortcut(), pn->descendantShortcut()->parent(), validationReachable);
    
    if (!pn->root())
    {
        assert(pn->size() > 1 || pn->numberOfRefs() > 0);
    }
    unsigned nzeros = 0;
    unsigned nones = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        pair<unsigned,unsigned> tmp = validate(*it, pn, validationReachable);
        nzeros += tmp.first;
        nones += tmp.second;
    }
    assert(nzeros == pn->nZeros());
    assert(nones == pn->nOnes());
    return make_pair(nzeros, nones);
}
void PointerTree::validate()
{
    for (unsigned i = 0; i<n; ++i)
	validationReachable[i]=false;
    pair<unsigned,unsigned> count = ::validate(&r, 0, validationReachable);
    unsigned reachable=0;
    for (unsigned i = 0; i<n; ++i)
	if(validationReachable[i])
            ++reachable;
    assert(reachable == n);
    assert(count.first + count.second == n); // Number of ones and zeros at root
}

/**
 * Graphical output. Requires Graphwiz DOT.
 */
unsigned outputDOT(PointerTree::PointerNode *pn, unsigned id, ostream &of)
{    
    unsigned cur_id = id;
    if (pn->leaf())
        return id;
    
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        of << " n" << cur_id << " -> n" << ++id  << endl;
        of << " n" << id << " [label=\"";
        if ((*it)->leaf())
            of << (*it)->nZeros() << "v" << (*it)->nOnes() << "\",shape=box";
        else
        {
            of << (*it)->depth() << " " << (*it)->nZeros() << "v" << (*it)->nOnes();
            if ((*it)->reduced())
                of << "R" << (int) (*it)->reducedLabel();
            of << "\"";
        }
        of << "]" << endl;
        id = outputDOT(*it, id, of);
    }
    return id;
}
void PointerTree::outputDOT(string const &filename, unsigned step)
{
    /**
     * Example graph:
     *  digraph G {
     *     main -> parse -> execute;
     *     main -> init;
     *  }
     */
    ostream *of = 0;
    if (filename == "-")
        of = &std::cout;
    else
    {
        char fn[256];
        snprintf(fn, 256, "%s.%u.dot", filename.c_str(), step);
        of = new ofstream (fn);
    }
    (*of) << "digraph G {" << endl;
    ::outputDOT(&r, 0, *of);
    (*of) << "}" << endl;
    if (of != &std::cout)
        delete of;
}
