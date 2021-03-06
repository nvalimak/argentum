#include "Tree.h"
#include "TreeDistance.h"
#include "TreeEnumerator.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
using namespace std;

// Default flag values
const NodeId PointerTree::nonreserved = (NodeId)~0u;
const unsigned PointerTree::nohistory = ~0u;
const int PointerTree::unknown = (int)-1000;

PointerTree::PointerTree(InputColumn const &ic)
    : nodes(), N(0), r(0), n(ic.size()), leaves(), stashv(), nstashed(0), nrelocate(0), history(), validationReachable(),
      reusedRootHistoryEvent(0), reusedHistoryEvent(0),
      nonbranching(), vacant(), nextVacant(0), nextUid(0)
{
    history.reserve(HISTORY_INIT_SIZE);
    vacant.reserve(HISTORY_INIT_SIZE);
    if (nodes.size() == 0)
    {
        // Initialize the static node buffer
        if (BUFFER_INIT_SIZE > ic.size()*2)
            nodes.resize(BUFFER_INIT_SIZE);
        else
            nodes.resize(ic.size()*2);

        for (size_t id = 0; id < nodes.size(); ++id)
            nodes[id] = new PointerNode(this);
    }

    NodeId id = 0;
    leaves.reserve(ic.size());
    r = createNode(1.0, PointerTree::nonreserved, PointerTree::nohistory, true); // New root
    for (InputColumn::const_iterator it = ic.begin(); it != ic.end(); ++it)
    {
        NodeId c = createNode(0.0, r, PointerTree::nohistory, false, id++);
        nodes[r]->insert(c);
        leaves.push_back(nodes[c]);
    }
    validationReachable.resize(n);
}

PointerTree::~PointerTree()
{
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        delete nodes[i]; nodes[i] = 0;
    }
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
        while (!pn->root() && pn->size() < 2 && pn->numberOfRefs() == 0)
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

void PointerTree::clearNonBranchingInternalNodes(unsigned step, TreeEnumerator *te)
{
    for (set<PointerNode *>::iterator it = nonbranching.begin(); it != nonbranching.end(); ++it)
    {
        if ((*it)->nodeId() == PointerTree::nonreserved)
            continue;
        PointerNode *pn = *it;
        // Safe to truncate 'pn' out from the tree
        while (!pn->root() && pn->size() < 2 && pn->numberOfRefs() == 0)
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
    ++nrelocate;
    NodeId src = pn->parent();
    nodes[src]->erase(pn);

    if (!keepparentcounts)
        propagateUpwardCounts(nodes[src], -(pn->nZeros()), -(pn->nOnes()));

    // Update history event queue
    if (keephistory)
    {
        if (!pn->tagged())
            history.push_back(Event(nodes[src], pn, step, history.size(), false));
        else if (getPreviousEventStep(pn) < destPreviousUpdate) // Assert: tagged
            history.push_back(Event(nodes[src], pn, step, history.size(), false));
        else
        {
            // Assert: tagged, reused history event
            if (src == r)
                reusedRootHistoryEvent ++;
            else
                reusedHistoryEvent ++;
        }
    }
        
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

void propagateStashed(PointerTree::PointerNode *pn, bool value)
{
    pn->stashed(value);
    if (pn->leaf())
        return;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        propagateStashed(*it, value);
} 

/**
 * Stash subtree 'pn' (to be later inserted at root)
 * Adds an event into the history queue.
 * Cleans up the source subtree.
 */
void PointerTree::stash(PointerNode *pn, unsigned step, bool keephistory, bool keepparentcounts)
{
    ++nstashed;
    ::propagateStashed(pn, true);
    pn->stashedFrom(pn->parent());
    stashv.push_back(pn);
    NodeId src = pn->parent();
    nodes[src]->erase(pn);
    if (!keepparentcounts)
        propagateUpwardCounts(nodes[src], -(pn->nZeros()), -(pn->nOnes()));

    // Update history event queue
    if (keephistory)
        history.push_back(Event(nodes[src], pn, step, history.size(), true));
    
    // Clean source subtree if needed
    if (nodes[src]->size() == 1 && !nodes[src]->root())
        PointerTree::clearNonBranchingInternalNode(nodes[src]);
}

/**
 * Unstash all subtrees (insert them on upward path to root)
 */
void PointerTree::unstash()
{
    for (vector<PointerNode *>::iterator it = stashv.begin(); it != stashv.end(); ++it)
    {
        PointerNode *pn = *it;
        PointerNode *dest = findUnstashUpwardDestination(nodes[pn->stashedFrom()]);
        dest->insert(pn);
        pn->setParent(dest);
        pn->tagged(true);
        propagateUpwardCounts(dest, pn->nZeros(), pn->nOnes());
        ::propagateStashed(pn, false);
    }
    stashv.clear();
}

// chooses the first non-reduced parent
PointerTree::PointerNode * PointerTree::findUnstashUpwardDestination(PointerNode *pn)
{
    if (pn->root())
        return pn;
    if (!pn->reduced())
        return pn;    
    return findUnstashUpwardDestination(pn->parentPtr());
}


// Unstash all subtrees (according to the given tree distance)
void PointerTree::unstash(TreeDistance &dist)
{
    vector<PointerNode *> destv(stashv.size(), 0);
    for (size_t i = 0; i < stashv.size(); ++i)
        destv[i] = findUnstashDestination(nodes[(stashv[i])->stashedFrom()], dist);

    for (size_t i = 0; i < stashv.size(); ++i)
    {
        PointerNode *pn = stashv[i];
        PointerNode *dest = destv[i];
        dest->insert(pn);
        pn->setParent(dest);
        pn->tagged(true);
        propagateUpwardCounts(dest, pn->nZeros(), pn->nOnes());
        //dest->addZeros(pn->nZeros());
        //dest->addOnes(pn->nOnes());        
        ::propagateStashed(pn, false);
    }
    stashv.clear();
}


pair<PointerTree::PointerNode *, double> evalZeroSkeleton(PointerTree::PointerNode *dest, PointerTree::PointerNode *best_dest, double best_eval, PointerTree::PointerNode *eval_pn, TreeDistance &dist)
{
    if (dest->reduced() && dest->reducedLabel() == 1)
        return make_pair(best_dest, best_eval);

    bool unarypath = false;
    for (PointerTree::PointerNode::iterator it = dest->begin(); it != dest->end(); ++it)
        if ((*it)->nZeros() == dest->nZeros() && (*it)->nOnes() == dest->nOnes())
            unarypath = true;

    double e = -1.0;
    if (!unarypath)
        e = dist.evalZeroSkeleton(dest, eval_pn);
    if (best_eval < e)
    {
        best_eval = e;
        best_dest = dest;
    }        
    for (PointerTree::PointerNode::iterator it = dest->begin(); it != dest->end(); ++it)
    {
        pair<PointerTree::PointerNode *, double> tmp = evalZeroSkeleton(*it, best_dest, best_eval, eval_pn, dist);
        if (best_eval < tmp.second)
        {
            best_eval = tmp.second;
            best_dest = tmp.first;
        }
    }
    return make_pair(best_dest, best_eval);
}

PointerTree::PointerNode * PointerTree::findUnstashDestination(PointerNode *pn, TreeDistance &dist)
{
    pair<PointerNode *, double> eval = ::evalZeroSkeleton(nodes[r], 0, -1.0, pn, dist);

    if (eval.first == 0) // All zero-branches were stashed
        eval.first = nodes[r];
     
    return eval.first;        
}

/**
 * Pre-rewind the given event
 *
 * Forces the tree become "tree consistent" at step h by relocating
 * all 0-branch events under the root (as a temporary measure).
 * All these 0-branches are then properly rewinded later in rewind(h).
 */
void PointerTree::prerewind(unsigned h, TreeEnumerator *te)
{
    if (history.empty())
        return;

    //cerr << "--- prerewind AT STEP " << h << " ";
    //te->debugPrint();
        
    for (vector<Event>::reverse_iterator it = history.rbegin(); it != history.rend() && it->getStep() == h; ++it)
        if (it->allZeros())
        {
            PointerNode *oldparent = it->getNode()->parentPtr();
            if (te && !oldparent->root())
            {
                //std::cerr << "!!!!!!!!!! prerewind closeChild called oldparent = " << oldparent->uniqueId() << " (ghost=" << oldparent->ghostbranch() << ",unary=" << oldparent->isUnary() << "), it = " << it->getNode()->uniqueId() << std::endl;
                te->closeChild(oldparent, it->getNode(), h+1, true);
                te->insertChild(nodes[r], it->getNode(), ~0u);
            }
            relocate(it->getNode(), nodes[r], 0, 0, false, false);
            if (te && !oldparent->root() && oldparent->isUnary())
                te->truncate(oldparent, h);
            else if (te && oldparent->ghostbranch())
                te->truncateGhost(oldparent, h);
        }
}

/**
 * Rewind the given event
 */
void PointerTree::rewind(Event &e, TreeEnumerator *te)
{
    PointerNode *dest = e.getSource();
    PointerNode *pn = e.getNode();
    //std::cerr << "!!!!!!!!!! rewind pn " << pn->uniqueId() << " under parent " << pn->parentPtr()->uniqueId() << " to dest " << dest->uniqueId() << std::endl;
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

    if (te && dest->isUnary())
        te->splitUnary(dest, e.getStep()); // Was unary
    else if (te && dest->ghostbranch())
        te->insertGhost(dest, e.getStep()); // Was ghost
    
    dest->insert(pn);
    propagateUpwardCounts(dest, pn->nZeros(), pn->nOnes());
    propagateUpwardCounts(src, -(pn->nZeros()), -(pn->nOnes()));
    pn->setParent(dest);
    if (te)
    {
        //if (!e.allZeros()) FIXME With old allzero
        //{
            //std::cerr << "rewind closeChild called src = " << src->uniqueId() << std::endl;
        te->closeChild(src, pn, e.getStep()+1, true);
            //}
        te->insertChild(dest, pn, e.getStep());
    }
    
    // Clean source subtree if needed
    src->erase(pn);
    if (te && src->isUnary()) //if (te && !e.allZeros() && src->isUnary())
        te->truncate(src, e.getStep()); // Becomes unary
    else if (te && src->ghostbranch()) // FIXME else if (te && !e.allZeros() && src->ghostbranch())
        te->truncateGhost(src, e.getStep()); // Becomes ghost

    if (src->size() <= 1 && !src->root())
        PointerTree::clearNonBranchingInternalNode(src);

    e.rewind(); // May delete src (truncated if pn is the only child)
}

/**
 * Rewind the given site, if it is next in the history
 */
void PointerTree::rewind(unsigned h, TreeEnumerator *te)
{
    //cerr << "--- AT STEP " << h << " ";
    //te->debugPrint();
    if (history.empty())
        return;
    // Process all rooted subtrees first
    for (vector<Event>::reverse_iterator it = history.rbegin(); it != history.rend() && it->getStep() == h; ++it)
        if (it->allZeros())
            rewind(*it, te);

    // Process rest of the subtress
    Event &e = history.back();
    while (e.getStep() == h)
    {
        if (!e.allZeros())
            rewind(e, te);
        history.pop_back();
        if (history.empty())
            break;
        e = history.back();
    }
}

/**
 * Update the max depth of each subtree
 */
unsigned PointerTree::updateMaxDists()
{
    return updateMaxDists(nodes[r])-1;
}
unsigned PointerTree::updateMaxDists(PointerTree::PointerNode * pn)
{
    if (pn->ghostbranch())
        return 0;
    if (pn->leaf())
    {
        if (pn->reducedLabel() == 1)
            return 0;
        return 1; // 0-leaves
    }

    bool unarypath = false;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if ((*it)->nZeros() == pn->nZeros() && (*it)->nOnes() == pn->nOnes())
            unarypath = true;

    unsigned maxd = 0;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        unsigned d = updateMaxDists(*it);
        if (maxd < d)
            maxd = d;
    }

    pn->maxDepth(maxd);
    if (unarypath || maxd == 0)
        return maxd;
    return maxd + 1;
}

/**
 * Validate the integrity of the tree structure
 */
pair<int, int> validate(PointerTree::PointerNode *pn, PointerTree::PointerNode *p, vector<bool> &validationReachable)
{
    assert (pn->root() || pn->leaf() || pn->previousUpdate() != PointerTree::nohistory);
    /*if (! (pn->root() || pn->ghostbranch() || !pn->tagged() || p->size() > 2 || p->root()))
        cerr << "at node id = " << pn->nodeId() << " nzeros = " << pn->nZeros() << " nones = " << pn->nOnes() << ", size = " << pn->size() << ", tagged = " << pn->tagged() << endl
             << "with parent id = " << p->nodeId() << " nzeros = " << p->nZeros() << " nones = " << p->nOnes() << ", size = " << p->size() << endl;
             assert (pn->root() || pn->ghostbranch() || !pn->tagged() || p->size() > 2 || p->root()); // tagged node cannot attach to a bi-/unary node (except root)*/
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
void PointerTree::outputDOT(PointerNode *pn, unsigned id, ostream &of)
{
    if (pn->leaf())
        return;
    if (pn->ghostbranch())
        return;

    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        if ((*it)->ghostbranch())
            continue;
        bool unarypath = false;
        for (PointerTree::PointerNode::iterator itt = (*it)->begin(); itt != (*it)->end(); ++itt)
            if ((*itt)->nZeros() == (*it)->nZeros() && (*itt)->nOnes() == (*it)->nOnes())
                unarypath = true;
    
        if (unarypath)
        {
            outputDOT(*it, id, of);
            continue;
        }
        of << " n" << id << " -> n" << (*it)->uniqueId() << endl;
        of << " n" << (*it)->uniqueId() << " [label=\"";
        if ((*it)->leaf())
            of << (*it)->leafId() + 1;
        else
            of << (*it)->uniqueId();
        if ((*it)->previousUpdate() != PointerTree::nohistory)
            of << "/" << (*it)->previousUpdate();
        of << " ";
        if ((*it)->tagged())
        {
            assert((*it)->previousEvent() == PointerTree::nohistory || PointerTree::history.size() >= (*it)->previousEvent());
            of << "F";
            if ((*it)->previousEvent() != PointerTree::nohistory)
                of << PointerTree::history[(*it)->previousEvent()].getStep();
            of << " ";
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
        outputDOT(*it, (*it)->uniqueId(), of);
    }
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

/**
 * Newick output
 */
void PointerTree::outputNewick(PointerNode *pn, unsigned depth, ostream &of)
{
    if (pn->leaf())
    {
        of << pn->leafId()+1 << ":" << depth << ".0";
        return;
    }
    if (pn->ghostbranch())
        return;

    bool unarypath = false;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if ((*it)->nZeros() == pn->nZeros() && (*it)->nOnes() == pn->nOnes())
            unarypath = true;

    if (!unarypath)
        of << "(";
    bool first = true;
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
    {
        if ((*it)->ghostbranch())
            continue;
    
        if (unarypath)
        {
            outputNewick(*it, depth, of);
            return;
        }
        
        if (!first)
            of << ",";
        first = false;
        outputNewick(*it, depth+1, of);
    }
    of << ")";
    if (!pn->root())
        of << ":" << depth << ".0";
}
void PointerTree::outputNewick(string const &filename, unsigned position)
{
    /**
     * Example tree:
     *  [66]((3:0.0368455,1:0.0368455):1.40615,(4:0.234725,2:0.234725):1.20828);
     */
    ostream *of = 0;
    if (filename == "-")
        of = &std::cout;
    else
    {
        char fn[256];
        snprintf(fn, 256, "%s.%u.newick", filename.c_str(), position);
        of = new ofstream (fn);
    }
    (*of) << "[" << position << "]";
    outputNewick(nodes[r], 0, *of);
    (*of) << ";" << endl;
    if (of != &std::cout)
        delete of;
}

