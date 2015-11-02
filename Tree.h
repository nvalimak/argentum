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
        explicit PointerNode(PointerTree *t_) // Empty node constructor
            : t(t_), ch(), id(PointerTree::nonreserved), lfId(PointerTree::nonreserved), d(1.0), p(PointerTree::nonreserved),
              nzeros(PointerTree::unknown), nones(PointerTree::unknown), nrefs(0),
              prevupdate(PointerTree::nohistory), flt(false), preve(PointerTree::nohistory), root_(false), stashed_(false),
              mdepth(0)
        { }

        // General accessors
        inline bool root() const
        { return root_; }
        inline bool leaf() const
        { return (lfId != PointerTree::nonreserved); }
        inline NodeId nodeId() const
        { return id; }
        inline LeafId leafId() const
        { return lfId; } 
        inline children_set_t::size_type size() const
        { return ch.size(); }
        inline TreeDepth depth() const
        { return d; }
        inline bool stashed() const
        { return stashed_; }
        inline void stashed(bool b)
        { stashed_ = b; }
        inline unsigned maxDepth() const
        { return mdepth; }
        inline void maxDepth(unsigned u)
        { mdepth = u; }

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
            PointerTree *t;
            children_set_t::iterator p;
        public:
            iterator(PointerTree *t_, children_set_t::iterator x) :t(t_), p(x) {}
            iterator(const iterator& mit) : p(mit.p) {}
            iterator& operator++() {++p;return *this;}
            iterator operator++(int) {iterator tmp(*this); operator++(); return tmp;}
            bool operator==(const iterator& rhs) {return p==rhs.p;}
            bool operator!=(const iterator& rhs) {return p!=rhs.p;}
            PointerNode *& operator*() {return t->nodes[*p];}
        };

        inline iterator begin()
        { return iterator(t, ch.begin()); }
        inline iterator end()
        { return iterator(t, ch.end()); }

        // Parent accessors
        inline NodeId parent() const
        { return p; }
        inline PointerNode * parentPtr() const
        { return t->nodes[p]; }
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
            if (nrefs == 0 && ch.size() < 2) // && p != 0)
                t->clearNonBranchingInternalNode(t->nodes[id]);
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
            assert (lfId == PointerTree::nonreserved);
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
            assert(root_ == false);
            assert (stashed_ == false);
        }

        // Internal node reset
        inline void reset(NodeId id_, NodeId lfId_, TreeDepth d_, NodeId p_, unsigned step, bool r)
        {
            assert(id == PointerTree::nonreserved);
            assert(p == PointerTree::nonreserved);
            id = id_;
            d = d_;
            p = p_;
            assert (lfId == PointerTree::nonreserved);
            lfId = lfId_;
            assert (nzeros == PointerTree::unknown);
            assert (nones == PointerTree::unknown);
            assert (nrefs == 0);
            prevupdate = step;
            flt = false;
            assert (preve = PointerTree::nohistory);
            assert (ch.size() == 0);
            root_ = r;
            assert (stashed_ == false);
            mdepth = 0;
        }
        
        virtual ~PointerNode()
        { }
        
        // Default copy constructor and assignment
    private:
        PointerNode(); // No default constructor
        
        PointerTree *t;
        children_set_t ch;  // Children
        NodeId id;          // Node identifier (unique)
        NodeId lfId;        // Leaf node?
        TreeDepth d;        // Given depth
        NodeId p;           // Parent node
        int nzeros;         // Number of zero labels in the subtree (reduced)
        int nones;          // Number of one labels in the subtree (reduced)
        unsigned nrefs;     // Number of active history references.
        unsigned prevupdate;// Previous update to the subtree
        bool flt;           // Floating node (flagged)
        unsigned preve;     // Previous event number (refers to the history vector)
        bool root_;
        bool stashed_; 
        unsigned mdepth;    // Subtree max depth to leaf (in non-reduced tree)
    };

    /**
     * Representation for relocation event that relocates 'pn' from under 'src' at column 'step'.
     */
    class Event
    {
    public:
        Event(PointerNode *src_, PointerNode *pn_, unsigned step_, unsigned histid, bool allzeros_)
            : src(src_), pn(pn_), step(step_), preve(PointerTree::nohistory), allzeros(allzeros_)
        {
            src->addRef(); // Increase the number of references to src.
            if (!pn->leaf())
                pn->addRef(); // Increase the number of references to pn.
            if (pn->hasPreviousEvent())
                preve = pn->previousEvent();
            pn->previousEvent(histid);
            //std::cerr << "adding event to node = " << pn->nodeId() << ", from src = " << src->nodeId() << ", step = " << step << ", preve = " << preve << std::endl;
        }

        /**
         * Rewire new pointer for the event
         * Reference counters need to be updated.
         */
        void rewire(PointerNode *new_pn)
        {
            if (!pn->leaf())
                pn->removeRef();  // Decrease the number of references; may delete src!
            pn = 0;
            new_pn->addRef(); // Increase the number of references to new_src.
            pn = new_pn;            
        }
        
        PointerNode * getSource() const
        { return src; }
        PointerNode * getNode() const
        { return pn; }
        unsigned getStep() const
        { return step; }

        bool allZeros() const
        { return allzeros; }

        void rewind()
        {
            src->removeRef(); // Decrease the number of references; may delete src!
            src = 0;          // src cannot be referenced anymore
            pn->previousEvent(preve); // Restore previous event (or PointerTree::nohistory)
            if (!pn->leaf())
                pn->removeRef(); // Decrease the number of references; may delete pn!
            pn = 0;
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
        bool allzeros;     // Subtree was all zeros at column 'step' 
    };

public:
    explicit PointerTree(InputColumn const &);
    virtual ~PointerTree();
    
    // Helpful accessors
    PointerNode * root()
    { return nodes[r]; }
    std::size_t size() const
    { return n; }
    std::size_t nnodes() const
    { return N; }
    std::size_t historySize() const
    { return history.size(); }
    PointerNode * leaf(std::size_t i)
    { return leaves[i]; }
    
    inline unsigned getPreviousEventStep(PointerNode *pn)
    { return history[pn->previousEvent()].getStep(); }
    
    // Tree modification (used in the class TreeController)
    PointerNode * createDest(PointerNode *, unsigned);
    void relocate(PointerNode *, PointerNode *, unsigned, unsigned, bool, bool);
    void stash(PointerNode *, unsigned, bool, bool);
    void unstash();
    void rewind(Event &);
    void rewind(unsigned);
    void prerewind(unsigned);

    void rewire(PointerNode *pn, PointerNode *new_pn)
    {
        assert (pn->hasPreviousEvent());
        history[pn->previousEvent()].rewire(new_pn);
    }
    
    // Clean nonbranching internal node, if possible
    void clearNonBranchingInternalNode(PointerNode *);
    void clearNonBranchingInternalNodes();

    // Debugging
    unsigned reusedHistoryEvents() const
    { return reusedHistoryEvent; }
    unsigned numberOfStashed() const
    { return nstashed; }
    unsigned numberOfRelocates() const
    { return nrelocate; }
    void validate();
    void outputDOT(std::string const &, unsigned);
    void outputNewick(std::string const &, unsigned);
    
    // Flags
    static const NodeId nonreserved;
    static const unsigned nohistory;
    static const int unknown;

    // FIXME Public?
    std::vector<PointerNode *> nodes; // Buffer for node allocation    
    std::size_t N; // Number of nodes
private:
    PointerTree();
    // No copy constructor or assignment
    PointerTree(PointerTree const&);
    PointerTree& operator = (PointerTree const&);
   
    NodeId createNode(TreeDepth d_, NodeId p_, unsigned step, bool root = false, NodeId lfId = PointerTree::nonreserved)
    {
        if (vacant.empty())
        {
            if (nextVacant >= nodes.size())
            {                
                nodes.resize(nodes.size()*2);
                for (size_t i = nodes.size()/2; i < nodes.size(); ++i)
                    nodes[i] = new PointerNode(this);                
            }
            nodes[nextVacant]->reset(nextVacant, lfId, d_, p_, step, root);
            return nextVacant++;
        }
        NodeId id = vacant.back();
        vacant.pop_back();
        nodes[id]->reset(id, lfId, d_, p_, step, root);
        return id;
    }
    
    void discardNode(NodeId id)
    {
        assert(nodes[id]->numberOfRefs() == 0);
        assert(nodes[id]->size() == 0);
        nodes[id]->reset();
        vacant.push_back(id);
    }

    void propagateUpwardCounts(PointerNode *, int, int);

    void outputDOT(PointerNode *, unsigned, std::ostream &);
    void outputNewick(PointerNode *, unsigned, std::ostream &);
    
    NodeId r;
    std::size_t n; // Number of leaves
    std::vector<PointerTree::PointerNode *> leaves;
    std::vector<PointerTree::PointerNode *> stashv;
    std::size_t nstashed;
    std::size_t nrelocate;
    std::vector<Event> history;
    std::vector<bool> validationReachable;
    unsigned reusedHistoryEvent;
    std::set<PointerTree::PointerNode *> nonbranching;
    std::vector<NodeId> vacant;
    NodeId nextVacant;
};
#endif
