#ifndef _TreeDistance_h_
#define _TreeDistance_h_
#include "default.h"
#include "Tree.h"
#include "NewickTree.h"
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

/**
 * Base class for distance computation
 *
 * Use TreeDistance::build() to construct an instance.
 */
class TreeDistance
{
public:
    enum tree_distance_t { distance_unset, distance_depth, distance_leaves };
    enum distance_scaling_t { scaling_none, scaling_exp_minus, scaling_log };
   
    // Constructor
    //static TreeDistance* build(PointerTree &, PointerTree &, tree_distance_t, distance_scaling_t);
    
    // Compute distances to be used for prediction at the given input column
    virtual void recomputeDistances(InputColumn const &) = 0;
    // Evaluate the distance for the given destination node on the target tree
    double evaluateDistance(std::vector<PointerTree::PointerNode *> const &, PointerTree::PointerNode *, PointerTree::PointerNode *);

    virtual ~TreeDistance()
    { }
protected:
    TreeDistance(PointerTree &target_)
        : target(target_), sourcedist()
    { }
            
    void distanceByTraversal(PointerTree::PointerNode *, PointerTree::PointerNode *, unsigned, LeafDistance &);
    void distanceByTraversal(PointerTree::PointerNode *, PointerTree::PointerNode const *, PointerTree::PointerNode const *, std::set<PointerTree::PointerNode *> &, unsigned, LeafDistance &);
    unsigned rootDistanceByTraversal(PointerTree::PointerNode *, PointerTree::PointerNode const *, std::set<PointerTree::PointerNode *> &);

    PointerTree &target;
    LeafDistance sourcedist;
};

/**
 * Distance computation for PointerTrees
 */
class PointerTreeDistance : public TreeDistance
{
public:
    PointerTreeDistance(PointerTree &, PointerTree &);

    // Compute distances to be used for prediction at the given input column
    virtual void recomputeDistances(InputColumn const &);
private:
    PointerTree &source;
};

/**
 * Distance computation for Newick trees
 */
class NewickDistance : public TreeDistance
{
public:
    // Constructor
    //static TreeDistance* build(NewickTree &, PointerTree &, tree_distance_t, distance_scaling_t);
    NewickDistance(NewickTree &, PointerTree &);

    // Compute distances to be used for prediction at the given input column
    virtual void recomputeDistances(InputColumn const &);

    virtual ~NewickDistance()
    { }
private:

    void distanceByTraversal(NewickTree::Node *, NewickTree::Node *, unsigned, LeafDistance &);
    
    NewickTree &source;
};
#endif
