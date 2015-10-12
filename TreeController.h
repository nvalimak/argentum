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
        : t(t_), recombine(), debug(debug_), dotfile(dotfile_)
    { }
    void process(InputColumn const &, unsigned);

    // Debug functionality
    unsigned countGhostBranches(PointerTree::PointerNode *);
    unsigned countUnaryGhosts(PointerTree::PointerNode *);
    unsigned countBranchingGhosts(PointerTree::PointerNode *);
    unsigned countActive(PointerTree::PointerNode *);
        
private:
    std::pair<int, int> reduce(PointerTree::PointerNode *, InputColumn const &);
    void resolveNonBinary(PointerTree::PointerNode *);
    std::pair<int,int> recombineStrategy(PointerTree::PointerNode *);
    void recombineSubtrees(bool, bool);
    void findReduced(PointerTree::PointerNode *, InputLabel);

    PointerTree &t;
    std::vector<PointerTree::PointerNode *> recombine;
    bool debug;
    std::string dotfile;
    
    TreeController();
    // No copy constructor or assignment
    TreeController(TreeController const&);
    TreeController& operator = (TreeController const&);
};
#endif
