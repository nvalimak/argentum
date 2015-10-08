#include "Tree.h"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cassert>
using namespace std;

const InputLabel PointerTree::nonreducible = (InputLabel)~0u;

PointerTree::PointerTree(InputColumn const &ic)
    : r(), n(ic.size())
{
    LeafId id = 0;    
    for (InputColumn::const_iterator it = ic.begin(); it != ic.end(); ++it)
        r.insert(new PointerNode(id++, 0.0, &r));
}


/**
 * Create a new internal node above the given node
 */
PointerTree::PointerNode * PointerTree::createDest(PointerNode *pn)
{
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
 * Relocate subtree 'pn' as a new child of 'dest'.
 * Cleans up the source subtree.
 */
void PointerTree::relocate(PointerNode *pn, PointerNode *dest)
{
    PointerNode *src = pn->parent();
    src->erase(pn);

    // Clean source subtree if needed
    if (src->size() == 1 && !src->root())
    {
        // Discard 'src' as non-branching internal node
        PointerNode *parent = src->parent();
        PointerNode *only_chld = *(src->begin());
        parent->erase(src);
        src->erase(only_chld);     // Release the only child from src
        parent->insert(only_chld); // Rewire source's parent to point to the only child
        only_chld->setParent(parent);
        assert(src->size() == 0);
        delete src;
    }

    dest->insert(pn);
    pn->setParent(dest);
}

void validate(PointerTree::PointerNode *pn, PointerTree::PointerNode *p)
{
    assert (p == pn->parent());
    if (pn->leaf())
    {
        assert (pn->leafId() != (LeafId)~0u);
        return;
    }
    
    assert (pn->root() || pn->size() > 1);
    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        validate(*it, pn);
}
        

void PointerTree::validate()
{
    ::validate(&r, 0);
}

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
            of << (*it)->depth();
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
    
}
