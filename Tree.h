#ifndef _Tree_h_
#define _Tree_h_
#include "default.h"
#include <set>
#include <vector>
#include <list>
#include <iostream>
#include <cassert>

/**
 * Pointer-based implementation of the tree
 */
class PointerTree
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
    class PointerNode
    {
        typedef std::set<NodeId> children_set_t;
    public:        
        PointerNode() // Empty node constructor
            : ch(), id(PointerTree::nonreserved), lf(false), d(1.0), p(PointerTree::nonreserved),
              nzeros(PointerTree::unknown), nones(PointerTree::unknown), nrefs(0),
              prevupdate(PointerTree::nohistory), flt(false), preve(PointerTree::nohistory)
        { }
        // Root node constructor
        PointerNode(TreeDepth d_, NodeId p_)
            : ch(), id(0), lf(false), d(d_), p(p_),
              nzeros(PointerTree::unknown), nones(PointerTree::unknown), nrefs(0),
              prevupdate(PointerTree::nohistory), flt(false), preve(PointerTree::nohistory)
        { }
        // Leaf node constructor
        PointerNode(LeafId id_, TreeDepth d_, NodeId p_)
            : ch(), id(id_), lf(true), d(d_), p(p_),
              nzeros(PointerTree::unknown), nones(PointerTree::unknown), nrefs(0),
              prevupdate(PointerTree::nohistory), flt(false), preve(PointerTree::nohistory)
        { }

        // General accessors
        inline bool root() const
        { return id == 0; }
        inline bool leaf() const
        { return lf; }
        inline NodeId nodeId() const
        { return id; }
        inline LeafId leafId() const
        { return id-1; /* Root is id==0 */ } 
        inline children_set_t::size_type size() const
        { return ch.size(); }
        inline TreeDepth depth() const
        { return d; }

        // History tracking
        inline bool hasPreviousEvent() const
        { return preve != PointerTree::nohistory; }
        inline unsigned previousEvent() const
        { return preve; }
        inline void previousEvent(unsigned e)
        { preve = e; }
        inline unsigned previousUpdate() const
        { return prevupdate; }
        inline void previousUpdate(unsigned e)
        { prevupdate = e; }
        inline bool floating() const
        { return flt; }
        inline void floating(bool flt_)
        { flt = flt_; }
        
        // Accessors to recuded tree representation
        inline bool ghostbranch() const
        { return (nones == 0 && nzeros == 0); }
        inline bool reduced() const
        { return (nones == 0 || nzeros == 0) && nones+nzeros>0; }
        inline InputLabel reducedLabel() const
        { return (nones > 0 && nzeros == 0); }
        inline void nZeros(int zeros_)
        { nzeros = zeros_; }
        inline void nOnes(int nones_)
        { nones = nones_; }
        inline void addZeros(int zeros_)
        { nzeros += zeros_; }
        inline void addOnes(int nones_)
        { nones += nones_; }
        inline int nZeros() const
        { return nzeros; }
        inline int nOnes() const
        { return nones; }
        void setLabel(InputLabel il)
        {
            nzeros = 0; nones = 0;
            if (il == 1)
                nones ++;
            else
                nzeros ++;
        } 
        
        // Iterator to child nodes
        class iterator 
        {            
            children_set_t::iterator p;
        public:
            iterator(children_set_t::iterator x) :p(x) {}
            iterator(const iterator& mit) : p(mit.p) {}
            iterator& operator++() {++p;return *this;}
            iterator operator++(int) {iterator tmp(*this); operator++(); return tmp;}
            bool operator==(const iterator& rhs) {return p==rhs.p;}
            bool operator!=(const iterator& rhs) {return p!=rhs.p;}
            PointerNode *& operator*() {return PointerTree::nodes[*p];}
        };

        inline iterator begin()
        { return iterator(ch.begin()); }
        inline iterator end()
        { return iterator(ch.end()); }

        // Parent accessors
        inline NodeId parent() const
        { return p; }
        inline PointerNode * parentPtr() const
        { return PointerTree::nodes[p]; }
        inline void setParent(NodeId p_)
        { p = p_; }
        inline void setParent(PointerNode *pn)
        { p = pn->nodeId(); }
    
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
            assert(nrefs > 0);
            --nrefs;
            // If no other references remain, clean the tree
            if (nrefs == 0 && ch.size() < 2 && p != 0)
                PointerTree::clearNonBranchingInternalNode(PointerTree::nodes[id]);
            // Note: this object can commit suicide ('delete this') above!
        }
        inline unsigned numberOfRefs() const
        { return nrefs; }
        
        // Modify the set of children        
        inline void erase(NodeId pn)
        { ch.erase(pn); /* Note: PointerNode instance still in memory. */ }
        inline void insert(NodeId pn)
        { ch.insert(pn); }
        inline void erase(PointerNode *pn)
        { ch.erase(pn->nodeId()); }
        inline void insert(PointerNode *pn)
        { ch.insert(pn->nodeId()); }

        inline void reset()
        {
            id = PointerTree::nonreserved;
            assert (lf == false);
            d = 1.0;
            p = PointerTree::nonreserved;
            nzeros  = PointerTree::unknown;
            nones = PointerTree::unknown;
            assert (nrefs  == 0);
            prevupdate = PointerTree::nohistory;
            flt = false;
            preve = PointerTree::nohistory;
            assert (ch.size() == 0);
            assert(nrefs == 0);
        }

        // Internal node reset
        inline void reset(NodeId id_, TreeDepth d_, NodeId p_, unsigned step)
        {
            assert(id == PointerTree::nonreserved);
            assert(p == PointerTree::nonreserved);
            id = id_;
            d = d_;
            p = p_;
            assert (lf == false);
            assert (nzeros == PointerTree::unknown);
            assert (nones == PointerTree::unknown);
            assert (nrefs == 0);
            prevupdate = step;
            flt = false;
            assert (preve = PointerTree::nohistory);
            assert (ch.size() == 0);
        }
        
        virtual ~PointerNode()
        { }
        
        // Default copy constructor and assignment
    private:
        children_set_t ch;  // Children
        NodeId id;          // Node identifier (unique)
        bool lf;            // Leaf node?
        TreeDepth d;        // Given depth
        NodeId p;           // Parent node
        int nzeros;         // Number of zero labels in the subtree
        int nones;          // Number of one labels in the subtree
        unsigned nrefs;     // Number of active history references.
        unsigned prevupdate;// Previous update to the subtree
        bool flt;           // Floating node (flagged)
        unsigned preve;     // Previous event number (refers to the history vector)
    };

    /**
     * Representation for relocation event that relocates 'pn' from under 'src' at column 'step'.
     */
    class Event
    {
    public:
        Event(PointerNode *src_, PointerNode *pn_, unsigned step_, unsigned histid)
            : src(src_), pn(pn_), step(step_), preve(PointerTree::nohistory)
        {
            src->addRef(); // Increase the number of references to src.
            if (pn->hasPreviousEvent())
                preve = pn->previousEvent();
            pn->previousEvent(histid);
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
        unsigned getStep() const
        { return step; }
       
        void rewind()
        {
            src->removeRef(); // Decrease the number of references; may delete src!
            src = 0;          // src cannot be referenced anymore
            pn->previousEvent(preve); // Restore previous event (or PointerTree::nohistory)
        }
        
        ~Event ()
        { /* Event does not own any of the object instances */ }

        // Default copy constructor and assignment
    private:
        Event();
        PointerNode * src; // Source node
        PointerNode * pn;  // Subtree that was relocated
        unsigned step;     // Column that initiated this event
        unsigned preve;    // Previous event number for 'src' (refers to the history vector)
    };

public:
    explicit PointerTree(InputColumn const &);
    virtual ~PointerTree();
    
    // Helpful accessors
    PointerNode * root()
    { return nodes[r]; }
    std::size_t size()
    { return n; }
    std::size_t nnodes()
    { return N; }

    inline unsigned getPreviousEventStep(PointerNode *pn)
    { return history[pn->previousEvent()].getStep(); }
    
    // Tree modification (used in the class TreeController)
    PointerNode * createDest(PointerNode *, unsigned);
    void relocate(PointerNode *, PointerNode *, unsigned, bool, bool);
    void stash(PointerNode *, unsigned, bool, bool);
    void unstash();
    void rewind(Event &);
    
    // Clean nonbranching internal node, if possible
    static void clearNonBranchingInternalNode(PointerNode *);

    // Debugging
    void validate();
    void outputDOT(std::string const &, unsigned);
    
    // Flags
    static const NodeId nonreserved;
    static const unsigned nohistory;
    static const int unknown;

    // FIXME Public?
    static std::vector<PointerNode *> nodes; // Buffer for node allocation    
    static std::size_t N; // Number of nodes
private:
    PointerTree();
    // No copy constructor or assignment
    PointerTree(PointerTree const&);
    PointerTree& operator = (PointerTree const&);
   
    static NodeId createNode(TreeDepth d_, NodeId p_, unsigned step)
    {
        if (vacant.empty())
        {
            if (nextVacant >= nodes.size())
            {                
                nodes.resize(nodes.size()*2);
                for (size_t i = nodes.size()/2; i < nodes.size(); ++i)
                    nodes[i] = new PointerNode();                
            }
            nodes[nextVacant]->reset(nextVacant, d_, p_, step);
            return nextVacant++;
        }
        NodeId id = vacant.back();
        vacant.pop_back();
        nodes[id]->reset(id, d_, p_, step);
        return id;
    }
    
    static void discardNode(NodeId id)
    {
        assert(nodes[id]->numberOfRefs() == 0);
        assert(nodes[id]->size() == 0);
        nodes[id]->reset();
        vacant.push_back(id);
    }

    void propagateUpwardCounts(PointerNode *, int, int);

    unsigned outputDOT(PointerNode *, unsigned, std::ostream &);
    
    NodeId r;
    std::size_t n; // Number of leaves
    std::vector<PointerTree::PointerNode *> stashv;
    std::vector<Event> history;
    std::vector<bool> validationReachable;
    static std::vector<NodeId> vacant;
    static NodeId nextVacant;
};
#endif
