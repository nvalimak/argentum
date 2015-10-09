#ifndef _Tree_h_
#define _Tree_h_
#include "default.h"
#include <unordered_set>
#include <iostream>
/**
 * Base class for trees
 *
 * To be removed? - This is an uncessecary class in current implementation. 
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
public:
    /**
     * Reperesentation for a node (both internal and leaf).
     *
     * Leaf nodes are static, and the pointer to a leaf node never changes.
     *
     * Internal nodes get created and destroyed during the computation.
     * Once created, the pointer to an existing internal node does not change.
     * Root node never gets deleted and is static.
     * Internal nodes can get deleted (or truncated) if they have less than two children 
     *   *and* there are no active history references to the internal node.
     * Thus, an internal node is not removed as long as there is one or more history
     * events that refer to it.
     */
    class PointerNode : public Tree::Node
    {
    public:
        PointerNode() // Root node constructor
            : ch(), id(PointerTree::nonleaf), d(1.0), p(0), reducel(PointerTree::nonreducible), nrefs(0), preve(PointerTree::nohistory), descentNonGhost(0)
        { }
        // Internal node constructor
        PointerNode(TreeDepth d_, PointerNode *p_)
            : ch(), id(PointerTree::nonleaf), d(d_), p(p_), reducel(PointerTree::nonreducible), nrefs(0), preve(PointerTree::nohistory), descentNonGhost(0)
        { }
        // Leaf node constructor
        PointerNode(LeafId id_, TreeDepth d_, PointerNode *p_)
            : ch(), id(id_), d(d_), p(p_), reducel(PointerTree::nonreducible), nrefs(0), preve(PointerTree::nohistory), descentNonGhost(0)
        { }

        // General accessors
        inline bool root() const
        { return p == 0; }
        inline bool leaf() const
        { return id != PointerTree::nonleaf; }
        inline LeafId leafId() const
        { return id; }
        inline std::unordered_set<PointerNode *>::size_type size() const
        { return ch.size(); }
        inline TreeDepth depth() const
        { return d; }
        inline unsigned previousEvent() const
        { return preve; }
        inline void previousEvent(unsigned e)
        { preve = e; }
        
        // Accessors to recuded tree representation
        inline bool ghostbranch() const
        { return reducel == PointerTree::ghostbranch; }
        inline bool reduced() const
        { return reducel < PointerTree::ghostbranch; }
        inline InputLabel reducedLabel() const
        { return reducel; }
        inline void setReduced(InputLabel il)
        { reducel = il; }
        
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

        // Parent accessors
        inline PointerNode * parent() const
        { return p; }
        inline void setParent(PointerNode * p_)
        { p = p_; }

        // Descendant shortcut
        inline bool hasShortcut() const
            { return false ; }  //descentNonGhost != 0; }
        inline PointerNode * descendantShortcut() const
        { return descentNonGhost; }
        inline void descendantShortcut(PointerNode *dng)
        { descentNonGhost = dng; }
        
        // Accessors to keep track of history references
        inline void addRef()
        { ++nrefs; }
        /**
         * Decreases the number of outside references to *this.
         * Take care when calling this function for pointer P:
         * if number of references becomes 0, P may get deleted here.
         */
        inline void removeRef()
        {
            --nrefs;
            // If no other references remain, clean the tree
            if (nrefs == 0 && ch.size() < 2 && p != 0)
                PointerTree::clearNonBranchingInternalNode(this);
            // Note: this object can commit suicide ('delete this') above!
        }
        inline unsigned numberOfRefs() const
        { return nrefs; }
        
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
        std::unordered_set<PointerNode *> ch; // Children
        LeafId id;          // Leaf identifier (unique)
        TreeDepth d;        // Given depth
        PointerNode * p;    // Parent node
        InputLabel reducel; // Label of the reduced subtree
        unsigned nrefs;     // Number of active history references.
        unsigned preve;     // Previous event number (refers to the history vector)
        PointerNode *descentNonGhost; // Short-cut descent to next non-ghost node below.
    };

    /**
     * Representation for relocation event that relocates 'pn' from under 'src'.
     */
    class Event
    {
    public:
        Event(PointerNode *src_, PointerNode *pn_)
            : src(src_), pn(pn_)
        {
            src->addRef(); // Increase the number of references to src.
        }

        /**
         * Set new source for the event
         * Reference counters need to be updated.
         */
        void rewire(PointerNode *new_src)
        {
            src->removeRef();  // Decrease the number of references; may delete src!
            src = 0;
            new_src->addRef(); // Increase the number of references to new_src.
            src = new_src;
        }

        PointerNode * getSource() const
        { return src; }
        PointerNode * getNode() const
        { return pn; }

        
        void rewind()
        {
            src->removeRef(); // Decrease the number of references; may delete src!
            src = 0; // Cannot be referenced anymore
        }
        
        ~Event ()
        { /* Event does not own any of the object instances */ }

        // Default copy constructor and assignment
    private:
        Event();
        PointerNode *src; // Source node
        PointerNode *pn;  // Subtree that was relocated
    };
    
public:
    PointerTree(InputColumn const &);

    void validate();
    void outputDOT(std::string const &, unsigned);

    // Helpful accessors
    PointerNode & root()
    { return r; }
    std::size_t size()
    { return n; }
    std::size_t nnodes()
    { return N; }

    // Tree modification (used in the class TreeController)
    PointerNode * createDest(PointerNode *);
    void relocate(PointerNode *, PointerNode *, bool);
    void rewind(Event &);
    
    // Clean nonbranching internal node, if possible
    static void clearNonBranchingInternalNode(PointerNode *);
    
    // Flags for nonreducible subtrees etc.
    static const InputLabel nonreducible;
    static const InputLabel ghostnode;
    static const InputLabel ghostbranch;
    static const LeafId nonleaf;
    static const unsigned nohistory;
private:
    PointerNode r;
    std::size_t n; // Number of leaves
    static std::size_t N; // Number of nodes
    std::vector<Event> history;
};
#endif
