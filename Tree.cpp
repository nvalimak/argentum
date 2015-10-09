#include "Tree.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
using namespace std;

// Default flag values
const InputLabel PointerTree::nonreducible = (InputLabel)~0u;
const InputLabel PointerTree::ghostnode = ((InputLabel)~0u)-1;
const InputLabel PointerTree::ghostbranch = ((InputLabel)~0u)-2;
const LeafId PointerTree::nonleaf = (LeafId)~0u;
const unsigned PointerTree::nohistory = ~0u;
size_t PointerTree::N = 0;

PointerTree::PointerTree(InputColumn const &ic)
    : r(), n(ic.size()), history()
{
    history.reserve(HISTORY_INIT_SIZE);
    LeafId id = 0;
    for (InputColumn::const_iterator it = ic.begin(); it != ic.end(); ++it)
        r.insert(new PointerNode(id++, 0.0, &r));
    PointerTree::N = ic.size()+1;
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
        pn->erase(only_chld);      // Release the only child from src
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
void PointerTree::relocate(PointerNode *pn, PointerNode *dest, bool keephistory)
{
    PointerNode *src = pn->parent();
    src->erase(pn);

    // Update history event queue
    if (keephistory)
    {
        pn->previousEvent(history.size());
        history.push_back(Event(src, pn));
    }
    
    // Clean source subtree if needed
    if (src->size() == 1 && !src->root())
        PointerTree::clearNonBranchingInternalNode(src);

    dest->insert(pn);
    pn->setParent(dest);
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

    // Clean source subtree if needed
    if (src->size() == 1 && !src->root())
        PointerTree::clearNonBranchingInternalNode(src);

    dest->insert(pn);
    pn->setParent(dest);

    e.rewind(); // May delete dest (truncated if pn is the only child)
}

/**
 * Validate the integrity of the tree structure
 */
void validate(PointerTree::PointerNode *pn, PointerTree::PointerNode *p)
{
    assert (p == pn->parent());
    if (pn->leaf())
    {
        assert (pn->leafId() != PointerTree::nonleaf);
        assert (pn->numberOfRefs() == 0);
        return;
    }
    
    if (!pn->root())
    {
        assert(pn->size() > 1 || pn->numberOfRefs() > 0);
    }    
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        validate(*it, pn);
}
void PointerTree::validate()
{
    ::validate(&r, 0);
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
            of << (int) (*it)->reducedLabel() << "\",shape=box";
        else
        {
            of << (*it)->depth() << " " << (int)(*it)->reducedLabel();
            if ((*it)->reduced())
                of << "R";
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
