#include "Tree.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
using namespace std;

// Default flag values
const NodeId PointerTree::nonreserved = (NodeId)~0u;
const unsigned PointerTree::nohistory = ~0u;
const int PointerTree::unknown = (int)-1000;

// Static variables (FIXME clean up)
size_t PointerTree::N = 0; // Total number of nodes
vector<PointerTree::PointerNode *> PointerTree::nodes = std::vector<PointerTree::PointerNode *>();
set<PointerTree::PointerNode *> PointerTree::nonbranching = set<PointerTree::PointerNode *>();
vector<NodeId> PointerTree::vacant = vector<NodeId>();
NodeId PointerTree::nextVacant = 0;

PointerTree::PointerTree(InputColumn const &ic)
    : r(0), n(ic.size()), stashv(), history(), validationReachable()
{
    assert (nodes.size() == 0); // Only one instance of this class is allowed
    history.reserve(HISTORY_INIT_SIZE);
    if (BUFFER_INIT_SIZE > ic.size()*2)
        nodes.resize(BUFFER_INIT_SIZE);
    else
        nodes.resize(ic.size()*2);
        
    size_t id = 0;
    nodes[id++] = new PointerNode(1.0, PointerTree::nonreserved);
    for (InputColumn::const_iterator it = ic.begin(); it != ic.end(); ++it)
    {
        nodes[id] = new PointerNode(id, 0.0, r);
        nodes[r]->insert(id++);
    }
    nextVacant = id;
    for (; id < nodes.size(); ++id)
        nodes[id] = new PointerNode();

    PointerTree::N = ic.size()+1; // Initial size of the tree
    validationReachable.resize(n);
}

PointerTree::~PointerTree()
{
    for (size_t i = 0; i < nodes.size(); ++i)
        delete nodes[i];
}

/**
 * Create a new internal node above the given node
 */
PointerTree::PointerNode * PointerTree::createDest(PointerNode *pn, unsigned step)
{
    ++N;
    NodeId parent = pn->parent();
    assert (parent != PointerTree::nonreserved);
    nodes[parent]->erase(pn);
    TreeDepth d = (nodes[parent]->depth() - pn->depth())/2;
    NodeId dest = createNode(d, parent, step);
    nodes[parent]->insert(dest);
    nodes[dest]->insert(pn);
    pn->setParent(dest);
    nodes[dest]->nZeros(pn->nZeros());
    nodes[dest]->nOnes(pn->nOnes());
    return nodes[dest];
}

/**
 * Truncate the given internal node if there are no active references to it.
 */
void PointerTree::clearNonBranchingInternalNode(PointerNode *pn)
{
    assert(pn->size() < 2);
    assert(!pn->root());
    if (pn->numberOfRefs() > 0)        
        return; // Cannot be removed since one or more references appear in history

    nonbranching.insert(pn);
}

void PointerTree::clearNonBranchingInternalNodes()
{
    for (set<PointerNode *>::iterator it = nonbranching.begin(); it != nonbranching.end(); ++it)
    {
        if ((*it)->nodeId() == PointerTree::nonreserved)
            continue;
        PointerNode *pn = *it;
        // Safe to truncate 'pn' out from the tree
        while (pn->size() < 2 && pn->numberOfRefs() == 0)
        {
            NodeId parent = pn->parent();        
            assert(parent != PointerTree::nonreserved);
            nodes[parent]->erase(pn);
            if (pn->size() == 1)
            {
                // Rewire the only child of pn:
                PointerNode * only_chld = *(pn->begin());
                pn->erase(only_chld);      // Release the only child from pn
                nodes[parent]->insert(only_chld); // Rewire pn's parent to point to the only child
                only_chld->setParent(parent);
            }
            // Now safe to delete 'pn' out from the tree
            PointerTree::N--;
            assert(pn->size() == 0); // pn does not own any objects
            discardNode(pn->nodeId());
            pn = nodes[parent];
        }
    }
    nonbranching.clear();
}
    
/**
 * Relocate subtree 'pn' as a new child of 'dest'.
 * Adds an event into the history queue.
 * Cleans up the source subtree.
 */
void PointerTree::relocate(PointerNode *pn, PointerNode *dest, unsigned destPreviousUpdate, unsigned step, bool keephistory, bool keepparentcounts)
{
    NodeId src = pn->parent();
    nodes[src]->erase(pn);

    if (!keepparentcounts)
        propagateUpwardCounts(nodes[src], -(pn->nZeros()), -(pn->nOnes()));

    // Update history event queue
    if (keephistory)
        if (!pn->floating() || getPreviousEventStep(pn) < destPreviousUpdate)
            history.push_back(Event(nodes[src], pn, step, history.size()));
    
    // Clean source subtree if needed
    if (nodes[src]->size() == 1 && !nodes[src]->root())
        PointerTree::clearNonBranchingInternalNode(nodes[src]);

    dest->insert(pn);
    pn->setParent(dest);
    if (!keepparentcounts)
        propagateUpwardCounts(dest, pn->nZeros(), pn->nOnes());
}

void PointerTree::propagateUpwardCounts(PointerNode *pn, int nzeros, int nones)
{
    while (true) {
        pn->addZeros(nzeros);
        pn->addOnes(nones);
        if (pn->parent() == PointerTree::nonreserved)
            break;
        pn = nodes[pn->parent()];
    }
}

/**
 * Stash subtree 'pn' (to be later inserted at root)
 * Adds an event into the history queue.
 * Cleans up the source subtree.
 */
void PointerTree::stash(PointerNode *pn, unsigned step, bool keephistory, bool keepparentcounts)
{
    stashv.push_back(pn);
    NodeId src = pn->parent();
    nodes[src]->erase(pn);
    if (!keepparentcounts)
        propagateUpwardCounts(nodes[src], -(pn->nZeros()), -(pn->nOnes()));

    // Update history event queue
    if (keephistory)
        history.push_back(Event(nodes[src], pn, step, history.size()));

    // Clean source subtree if needed
    if (nodes[src]->size() == 1 && !nodes[src]->root())
        PointerTree::clearNonBranchingInternalNode(nodes[src]);
}

/**
 * Unstash all subtrees (insert them at root)
 */
void PointerTree::unstash()
{
    PointerNode *dest = nodes[r];
    for (vector<PointerNode *>::iterator it = stashv.begin(); it != stashv.end(); ++it)
    {
        PointerNode *pn = *it;
        dest->insert(pn);
        pn->setParent(dest);
        pn->floating(true);
        dest->addZeros(pn->nZeros());
        dest->addOnes(pn->nOnes());        
    }
    stashv.clear();
}

/**
 * Rewind the given event
 */
void PointerTree::rewind(Event &e)
{
    PointerNode *dest = e.getSource();
    PointerNode *pn = e.getNode();
    if (pn->nodeId() == PointerTree::nonreserved)
    {
        e.rewind();
        return;
    }
    PointerNode *src = pn->parentPtr();
    if (src == dest)
    {
        e.rewind();
        return;
    } 
    dest->insert(pn);
    propagateUpwardCounts(dest, pn->nZeros(), pn->nOnes());
    propagateUpwardCounts(src, -(pn->nZeros()), -(pn->nOnes()));

    pn->setParent(dest);

    // Clean source subtree if needed
    src->erase(pn);
    if (src->size() <= 1 && !src->root())
        PointerTree::clearNonBranchingInternalNode(src);

    e.rewind(); // May delete src (truncated if pn is the only child)
}

/**
 * Rewind the given site, if it is next in the history
 */
void PointerTree::rewind(unsigned h)
{
    if (history.empty())
        return;
    Event &e = history.back();
    while (e.getStep() == h)
    {
        rewind(e);
        history.pop_back();
        if (history.empty())
            return;
        e = history.back();
    }
}

/**
 * Validate the integrity of the tree structure
 */
pair<int, int> validate(PointerTree::PointerNode *pn, PointerTree::PointerNode *p, vector<bool> &validationReachable)
{
    assert (pn->root() || pn->leaf() || pn->previousUpdate() != PointerTree::nohistory);
    if (! (pn->root() || pn->ghostbranch() || !pn->floating() || p->size() > 2 || p->root()))
        cerr << "at node id = " << pn->nodeId() << " nzeros = " << pn->nZeros() << " nones = " << pn->nOnes() << ", size = " << pn->size() << ", floating = " << pn->floating() << endl
             << "with parent id = " << p->nodeId() << " nzeros = " << p->nZeros() << " nones = " << p->nOnes() << ", size = " << p->size() << endl;
    assert (pn->root() || pn->ghostbranch() || !pn->floating() || p->size() > 2 || p->root()); // floating node cannot attach to a bi-/unary node (except root)
    assert (pn->nZeros() != PointerTree::unknown);
    assert (pn->nOnes() != PointerTree::unknown);
    assert (pn->root() || p == pn->parentPtr());
    if (pn->leaf())
    {
        assert (pn->leafId() != PointerTree::nonreserved);
        assert (pn->reduced());
        assert (pn->numberOfRefs() == 0);
	validationReachable[pn->leafId()] = true;
        if (pn->reducedLabel() == 1)
            return make_pair(0,1);
        else
            return make_pair(1,0);
    }
    
    if (!pn->root())
    {
        if (!(pn->size() > 1 || pn->numberOfRefs() > 0))
            cerr << "at node = " << pn->nodeId() << ", size = " << pn->size() << endl;
        assert(pn->size() > 1 || pn->numberOfRefs() > 0);
    }
    int nzeros = 0;
    int nones = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        pair<int,int> tmp = validate(*it, pn, validationReachable);
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
    pair<int,int> count = ::validate(nodes[r], 0, validationReachable);
    int reachable=0;
    for (unsigned i = 0; i<n; ++i)
	if(validationReachable[i])
            ++reachable;
    assert((unsigned)reachable == n);
    assert(count.first + count.second == (int)n); // Number of ones and zeros at root
}

/**
 * Graphical output. Requires Graphwiz DOT.
 */
unsigned PointerTree::outputDOT(PointerNode *pn, unsigned id, ostream &of)
{    
    unsigned cur_id = id;
    if (pn->leaf())
        return id;
    
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        of << " n" << cur_id << " -> n" << ++id  << endl;
        of << " n" << id << " [label=\"";
        of << (*it)->nodeId();
        if ((*it)->previousUpdate() != PointerTree::nohistory)
            of << "/" << (*it)->previousUpdate();
        of << " ";
        if ((*it)->floating())
        {
            assert(PointerTree::history.size() > (*it)->previousEvent());
            of << "F" << PointerTree::history[(*it)->previousEvent()].getStep() << " ";
        }
        if ((*it)->leaf())
            of << (*it)->nZeros() << "v" << (*it)->nOnes() << "\",shape=box";
        else
        {
            of << " " << (*it)->nZeros() << "v" << (*it)->nOnes();
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
    outputDOT(nodes[r], 0, *of);
    (*of) << "}" << endl;
    if (of != &std::cout)
        delete of;
}
