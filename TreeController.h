#ifndef _TreeController_h_
#define _TreeController_h_
#include "default.h"
#include "Tree.h"
#include "TreeDistance.h"
#include "TreeEnumerator.h"
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
        : t(t_), recombine(), updatedThisStep(), step(0), debug(debug_), dotfile(dotfile_), printEval(false), nOnesCut(0)
    { }
    void process(InputColumn const &, unsigned);
    void process(InputColumn const &, unsigned, TreeDistance &);
    void rewind(InputColumn const &, unsigned, TreeEnumerator *);
    void assignLabels(InputColumn const &);
    void deTagAll(PointerTree::PointerNode *);

    // Debug functionality
    unsigned countGhostBranches(PointerTree::PointerNode *);
    unsigned countUnaryGhosts(PointerTree::PointerNode *);
    unsigned countBranchingGhosts(PointerTree::PointerNode *);
    unsigned countActive(PointerTree::PointerNode *);
    unsigned countInternalNodes();
    unsigned numberOfOneCuts()
    { return nOnesCut; }
private:
    std::pair<int, int> reduce(PointerTree::PointerNode *, InputColumn const &);
    void resolveNonBinary(PointerTree::PointerNode *);
    std::pair<int,int> recombineStrategy(PointerTree::PointerNode *);
    void recombineSubtrees(PointerTree::PointerNode *, bool, bool);
    void recombineSubtrees(PointerTree::PointerNode *, bool, bool, TreeDistance &);
    void recombineNonBinarySubtrees(unsigned, bool, bool);
    void findReduced(PointerTree::PointerNode *, InputLabel);
    
    PointerTree &t;
    std::vector<PointerTree::PointerNode *> recombine;
    std::vector<PointerTree::PointerNode *> updatedThisStep;
    unsigned step;
    bool debug;
    std::string dotfile;
    bool printEval;
    unsigned nOnesCut;
    
    TreeController();
    // No copy constructor or assignment
    TreeController(TreeController const&);
    TreeController& operator = (TreeController const&);
};
#endif
