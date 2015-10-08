#ifndef _Tree_h_
#define _Tree_h_
#include "default.h"
#include <unordered_set>

/**
 * Base class for trees
 */
class Tree
{
public:
    class Node
    {
    public:
        virtual bool leaf() const = 0;

        virtual ~Node()
        { }
        // Default copy constructor and assignment
    protected:
        Node()
        { }
    };

    virtual ~Tree()
    { }
protected:
    Tree()
    { }
private:
    // No copy constructor or assignment
    Tree(Tree const&);
    Tree& operator = (Tree const&);
};

/**
 * Pointer-based implementation of the tree
 */
class PointerTree : public Tree
{
public: // To be private? 
    class PointerNode : public Tree::Node
    {
    public:
        PointerNode() // Root node constructor
            : ch(), id(~0u), d(1.0), p(0), reduce(false)
        { }
        // Internal node constructor
        PointerNode(TreeDepth d_, PointerNode *p_)
            : ch(), id(~0u), d(d_), p(p_), reduce(false)
        { }
        // Leaf node constructor
        PointerNode(LeafId id_, TreeDepth d_, PointerNode *p_)
            : ch(), id(id_), d(d_), p(p_), reduce(false)
        { }

        inline bool root() const
        { return p == 0; }
        inline bool leaf() const
        { return ch.size() == 0; }
        inline LeafId leafId() const
        { return id; }
        inline std::unordered_set<PointerNode *>::size_type size() const
        { return ch.size(); }
        inline TreeDepth depth() const
        { return d; }
        inline PointerNode * parent() const
        { return p; }
        inline bool reduced() const
        { return reduce; }
        inline bool reducedLabel() const
        { return reducel; }
        inline void setReduced(bool b)
        { reduce = b; }
        inline void setReduced(bool b, InputLabel il)
        { reduce = b; reducel = il; }
        
        // Iterator to child nodes
        typedef std::unordered_set<PointerNode *>::iterator iterator;
        typedef std::unordered_set<PointerNode *>::const_iterator const_iterator;
        inline iterator begin()
        { return ch.begin(); }
        inline const_iterator begin() const
        { return ch.begin(); }
        inline iterator end()
        { return ch.end(); }
        inline const_iterator end() const
        { return ch.end(); }

        inline void setParent(PointerNode * p_)
        { p = p_; }
        // Modify the set of children        
        inline void erase(PointerNode *pn)
        { ch.erase(pn); /* Note: PointerNode instance still in memory. */ }
        inline void insert(PointerNode *pn)
        { ch.insert(pn); }

        virtual ~PointerNode()
        {
            // Free child nodes; assumes that all PointerNodes are owned by this class
            for (iterator it = ch.begin(); it != ch.end(); ++it)
                delete *it;
        }
        
        // Default copy constructor and assignment
    private:
        std::unordered_set<PointerNode *> ch;
        LeafId id;
        TreeDepth d;
        PointerNode * p;
        bool reduce;
        InputLabel reducel;
    };

public:
    PointerTree(InputColumn const &);

    void validate();
    void outputDOT(std::string const &, unsigned);

    PointerNode & root()
    { return r; }
    std::size_t size()
    { return n; }

    PointerNode * createDest(PointerNode *);
    void relocate(PointerNode *, PointerNode *);
    
    // Marks nonreducible subtrees
    static const InputLabel nonreducible;
private:
    PointerNode r;
    std::size_t n;
};
#endif
