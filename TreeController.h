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
        
protected:
    std::pair<unsigned, unsigned> reduce(PointerTree::PointerNode *, InputColumn const &);
//    InputLabel reduceGhosts(PointerTree::PointerNode *);
    void resolveNonBinary(PointerTree::PointerNode *);
    void collectRecombine(PointerTree::PointerNode *);
    void recombineSubtrees(bool, bool);
    PointerTree &t;
    std::vector<PointerTree::PointerNode *> recombine;
    bool debug;
    std::string dotfile;
private:
    TreeController();
    // No copy constructor or assignment
    TreeController(TreeController const&);
    TreeController& operator = (TreeController const&);
};
#endif
