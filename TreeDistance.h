#ifndef _TreeDistance_h_
#define _TreeDistance_h_
#include "default.h"
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
    //    static InputReader* build(tree_distance_t, distance_scaling_t);
};
#endif
