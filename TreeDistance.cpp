#include "TreeDistance.h"
#include <iostream>
#include <utility>
#include <set>
#include <map>
#include <cassert>
#include <cstdlib> // std::abort()
#include <cmath>
using namespace std;

/*TreeDistance* TreeDistance::build(PointerTree &source, PointerTree &target, tree_distance_t tdist, distance_scaling_t dscaling)
{
    return new TreeDistance(source, target);
    }*/


inline double distanceScaling(double d)
{
    return exp(-d); // FIXME
}

inline double distanceDifference(double t, double s)
{
    double d = t - s;
    return d * d;
}

void distanceNormalization(LeafDistance &ld)
{
    double mean = 0;
    for (LeafDistance::iterator it = ld.begin(); it != ld.end(); ++it)
        mean += *it;
    mean /= ld.size();

    double sd = 0;
    for (LeafDistance::iterator it = ld.begin(); it != ld.end(); ++it)
        sd += (*it - mean) * (*it - mean);
    sd /= ld.size();

    for (LeafDistance::iterator it = ld.begin(); it != ld.end(); ++it)
        *it = (*it - mean) / sd;    
}

PointerTreeDistance::PointerTreeDistance(PointerTree &source_, PointerTree &target_)
    : TreeDistance(target_), source(source_)
{
    // NOP
}

void PointerTreeDistance::recomputeDistances(InputColumn const & ic)
{
    source.updateMaxDists();
    if (sourcedist.empty())
        sourcedist.resize(ic.size(), -1.0);
    
    for (size_t l = 0; l < source.size(); ++l) // TODO make linear time
    {
        LeafDistance d(sourcedist.size(), 0);
        PointerTree::PointerNode * pn = source.leaf(l);
        if (pn->reducedLabel() != 0)
            continue;
        distanceByTraversal(pn->parentPtr(), pn, 0, d);
        sourcedist[l] = 0;
        for (size_t j = 0; j < d.size(); ++j)
            if (source.leaf(j)->reducedLabel() == 1)
                sourcedist[l] += ::distanceScaling(d[j]);
    }
}

// Assert: maxdepth values are set
void TreeDistance::distanceByTraversal(PointerTree::PointerNode * pn, PointerTree::PointerNode * src, unsigned mdepth, LeafDistance &d)
{
    mdepth = max(pn->maxDepth(), mdepth);
    if (pn->ghostbranch())
        return;
    if (pn->leaf())
    {
        if (pn->reducedLabel() == 1)
            d[pn->leafId()] = mdepth;
        return;
    }

    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if (*it != src)
            distanceByTraversal(*it, pn, mdepth, d);
            
    if (!pn->root() && pn->parentPtr() != src)
        distanceByTraversal(pn->parentPtr(), pn, mdepth, d);
}

void TreeDistance::distanceByTraversal(PointerTree::PointerNode * pn, PointerTree::PointerNode const *src, PointerTree::PointerNode const *dest, set<PointerTree::PointerNode *> &oner, unsigned mdepth, LeafDistance &d)
{
    mdepth = max(pn->maxDepth(), mdepth);
    if (pn == dest)
    {
        unsigned mm = mdepth;
        for (set<PointerTree::PointerNode *>::iterator it = oner.begin(); it != oner.end(); ++it)
            if (mm < (*it)->maxDepth())
                mm = (*it)->maxDepth();
        
        for (set<PointerTree::PointerNode *>::iterator it = oner.begin(); it != oner.end(); ++it)
        {
            if (*it != dest || dest->leaf())
                distanceByTraversal(*it, (*it)->parentPtr(), mm+1, d);
            else
                distanceByTraversal(*it, (*it)->parentPtr(), mm, d);
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
            assert (0); // This is obsolete code
            d[pn->leafId()] = mdepth;
        }
        return;
    }

    for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        if (*it != src)
            distanceByTraversal(*it, pn, dest, oner, mdepth, d);

    if (!pn->root() && pn->parentPtr() != src)
        distanceByTraversal(pn->parentPtr(), pn, dest, oner, mdepth, d);
}


double TreeDistance::evaluateDistance(std::vector<PointerTree::PointerNode *> const &recombine, PointerTree::PointerNode *subtree_root, PointerTree::PointerNode *dest)
{
    set<PointerTree::PointerNode *> oner(recombine.begin(), recombine.end());
    LeafDistance distForDest(sourcedist.size(), 0); 
    
    for (size_t l = 0; l < target.size(); ++l) // TODO make linear time
    {
        PointerTree::PointerNode * pn = target.leaf(l);
        if (pn->reducedLabel() != 0)
            continue;
        LeafDistance d(sourcedist.size(), 0);
        distanceByTraversal(pn->parentPtr(), pn, dest, oner, 0, d);
        
        distForDest[l] = 0;
        for (size_t j = 0; j < d.size(); ++j)
            if (target.leaf(j)->reducedLabel() == 1)
                distForDest[l] += ::distanceScaling(d[j]);
    }

    ::distanceNormalization(distForDest);
    
    double diff = 0;
    for (size_t l = 0; l < target.size(); ++l) // TODO make linear time
    {
        PointerTree::PointerNode * pn = target.leaf(l);
        if (pn->reducedLabel() != 0)
            continue;
        if (false)
            cerr << "dist[" << pn->leafId() << "] = " <<  sourcedist[pn->leafId()] << ", distForDest[" << pn->leafId() << ", " << dest->nodeId() << "] = " << distForDest[pn->leafId()] << endl;
        diff += ::distanceDifference(distForDest[pn->leafId()], sourcedist[pn->leafId()]);
    }
    return diff;
}

NewickDistance::NewickDistance(NewickTree &source_, PointerTree &target_)
    : TreeDistance(target_), source(source_)
{
    // NOP
}

// Assert: maxdepth values are set
void NewickDistance::distanceByTraversal(NewickTree::Node * pn, NewickTree::Node * src, unsigned mdepth, LeafDistance &d)
{
    mdepth = max(pn->mdepth, mdepth);
    if (pn->leaf)
    {
        if (pn->llabel == 1)
            d[pn->lid] = mdepth;
        return;
    }

    for (NewickTree::children_set_t::iterator it = pn->ch.begin(); it != pn->ch.end(); ++it)
        if (*it != src)
            distanceByTraversal(*it, pn, mdepth, d);
            
    if (pn->parent != 0 && pn->parent != src)
        distanceByTraversal(pn->parent, pn, mdepth, d);
}

void NewickDistance::recomputeDistances(InputColumn const & ic)
{
    source.updateMaxDists();
    if (sourcedist.empty())
        sourcedist.resize(ic.size(), 0);
    
    for (size_t l = 0; l < source.nleaves(); ++l) // TODO make linear time
    {
        LeafDistance d(sourcedist.size(), 0);
        NewickTree::Node * pn = source.leaf(l);
        if (pn->llabel != 0)
            continue;
        distanceByTraversal(pn->parent, pn, 0, d);
        sourcedist[l] = 0;
        for (size_t j = 0; j < d.size(); ++j)
            if (source.leaf(j)->llabel == 1)
                sourcedist[l] += ::distanceScaling(d[j]);
    }

    ::distanceNormalization(sourcedist);
}
