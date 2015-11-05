#ifndef _NewickTree_h_
#define _NewickTree_h_
#include "default.h"
#include <string>
#include <vector>
#include <set>
#include <cassert>
#include <iostream>

/**
 * Data structure to read in SCRM trees
 */
class NewickTree
{
public:
    class Node;
    typedef std::set<Node *> children_set_t;
    class Node
    {
    public:
        Node()
            : ch(), parent(0), lid(-1), llabel(-1), leaf(false), size(0), mdepth(0)
        { }
        Node(Node *p)
            : ch(), parent(p), lid(-1), llabel(-1), leaf(false), size(0), mdepth(0)
        {
            p->ch.insert(this); // Insert as child to *p
        }
        
        ~Node()
            { for (children_set_t::iterator it = ch.begin(); it != ch.end(); ++it) delete *it; ch.clear(); }
        children_set_t ch; // Children
        Node *parent; 
        int lid; // Leaf ID
        int llabel; // Leaf label
        bool leaf; // Is leaf?
        unsigned size; // Total number of subtree leaves
        unsigned mdepth; // Max depth for this subtree
    };

    // Distance computation for newickcomp.cpp
    void distanceByTraversal(Node *, Node *, Node *, unsigned, std::vector<unsigned> &);
    void subtreeDistance(std::vector<unsigned> &);
    NewickTree(std::string const &i);
    ~NewickTree();

    Node * root()
    { return root_; }
    unsigned validForBases()
    { return valid; }
    unsigned nleaves()
    { return root_->size; }
    Node * leaf(size_t i)
    { return leaves[i]; }
    
    void assignLabels(InputColumn const &);
    unsigned updateMaxDists();
    
    bool next();
    bool good()
    { return good_; }
    
    // Debuging output (Graphwiz DOT format)
    void outputDOT(std::string const &, unsigned);
private:
    void parse(std::string const &);
    size_t parse(std::string const &, size_t, Node *);
    int updateSizes(Node *);
    unsigned outputDOT(Node *, std::ostream &, unsigned);
    void collectLeaves(Node *);
    unsigned updateMaxDists(Node *);

    Node *root_;
    std::vector<Node *> leaves;
    unsigned valid; // Tree is valid for this many sites
    std::istream *fp;
    bool good_;
};
#endif
