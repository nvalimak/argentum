#ifndef _TreeController_h_
#define _TreeController_h_
#include "default.h"
#include "Tree.h"
#include <vector>
#include <string>
#include <utility>

/**
 * Controller to handle tree updates
 */
class TreeController
{
public:
    TreeController (PointerTree &t_, bool debug_, std::string dotfile_)
        : t(t_), recombine(), updatedThisStep(), step(0), debug(debug_), dotfile(dotfile_), printEval(false)
    { }
    void process(InputColumn const &, unsigned);
    void process(InputColumn const &, unsigned, LeafDistance const &);
    void rewind(InputColumn const &, unsigned);
    
    // Debug functionality
    unsigned countGhostBranches(PointerTree::PointerNode *);
    unsigned countUnaryGhosts(PointerTree::PointerNode *);
    unsigned countBranchingGhosts(PointerTree::PointerNode *);
    unsigned countActive(PointerTree::PointerNode *);
    void deFloatAll(PointerTree::PointerNode *);

    void distance(LeafDistance &, InputColumn const &);
    static void distance(PointerTree &t, LeafDistance &, InputColumn const &);
    static unsigned computeMaxDists(PointerTree::PointerNode * pn);
    static void distanceByTraversal(PointerTree::PointerNode *, PointerTree::PointerNode *, unsigned, LeafDistance &);
    static void distanceByTraversal(PointerTree::PointerNode *, PointerTree::PointerNode const *, PointerTree::PointerNode const *, std::set<PointerTree::PointerNode *> &, unsigned, LeafDistance &);
    static double evaluateDistance(PointerTree &t, std::vector<PointerTree::PointerNode *> const &, PointerTree::PointerNode *, PointerTree::PointerNode *, LeafDistance const &);

private:
    std::pair<int, int> reduce(PointerTree::PointerNode *, InputColumn const &);
    void resolveNonBinary(PointerTree::PointerNode *);
    std::pair<int,int> recombineStrategy(PointerTree::PointerNode *);
    std::pair<int,int> recombineStrategy(PointerTree::PointerNode *, LeafDistance const &);
    void recombineSubtrees(PointerTree::PointerNode *, bool, bool);
    void recombineSubtrees(PointerTree::PointerNode *, bool, bool, LeafDistance const &);
    void recombineNonBinarySubtrees(unsigned, bool, bool);
    void findReduced(PointerTree::PointerNode *, InputLabel);
        
    PointerTree &t;
    std::vector<PointerTree::PointerNode *> recombine;
    std::vector<PointerTree::PointerNode *> updatedThisStep;
    unsigned step;
    bool debug;
    std::string dotfile;
    bool printEval;
    
    TreeController();
    // No copy constructor or assignment
    TreeController(TreeController const&);
    TreeController& operator = (TreeController const&);
};
#endif
