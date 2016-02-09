#ifndef _TreeEnumerator_h_
#define _TreeEnumerator_h_
#include "default.h"
#include "InputReader.h"
#include "Tree.h"
#include <iostream>
#include <cassert>
#include <map>
#include <utility>

//#define TE_DEBUG_PRINT

class TreeEnumerator
{
private:
    struct ARGChild
    {
        NodeId id;      // Unique ID of the node; 0 == root
        std::vector<unsigned> mutations;
        unsigned lRange;  // Range corresponds to VCF row numbers
        unsigned rRange; 
        bool recomb; // Was cut at position rRange due to recomb
        ARGChild()
            : id(0), mutations(), lRange(0), rRange(0), recomb(false)
        { }
        ARGChild(NodeId id_, std::vector<unsigned> &mut_, unsigned lRange_, unsigned rRange_)
            : id(id_), mutations(mut_), lRange(lRange_), rRange(rRange_), recomb(false)
        { }
    };

public:
    TreeEnumerator()
        : insertBuffer(), rightPos(), mutationPos(), leftRightPos(), init(true)
    { }

    bool initialize() const
    { return init; }
    void initialize(PointerTree &tree, unsigned step)
    {
        init = false;
        initRightPos(tree.root(), tree.root()->uniqueId(), step);
    }

    ~TreeEnumerator()
    {
        if (rightPos.size() > 0)
            std::cerr << "TreeEnumerator()::~TreeEnumerator() warning: size of rightPos at destructor = " << rightPos.size() << std::endl;
    }
    
    void insertChild(PointerTree::PointerNode *pn, PointerTree::PointerNode *cpn, unsigned rstep)
    {
        if (cpn->ghostbranch())
            return;
        if (!cpn->leaf() && cpn->isUnary())
        {
            while (!cpn->leaf() && cpn->isUnary())
            {
                PointerTree::PointerNode *unarypn = 0;
                for (PointerTree::PointerNode::iterator it = cpn->begin(); it != cpn->end(); ++it)
                    if ((*it)->nZeros() == cpn->nZeros() && (*it)->nOnes() == cpn->nOnes())
                        unarypn = *it;
                cpn = unarypn;
            }                        
        }
        unsigned cid = cpn->uniqueId();

        PointerTree::PointerNode *r = pn;
        while (!r->root() && (r->isUnary() || r->ghostbranch()))
            r = r->parentPtr();
        
        rightPos[r->uniqueId()][cid] = rstep;        
    }

    void splitUnary(PointerTree::PointerNode *pn, PointerTree::PointerNode *dest, PointerTree::PointerNode *r, unsigned step)
    {
        //std::cerr << "recursive splitUnary called pn = " << pn->uniqueId() << ", dest = " << dest->uniqueId() << " (pn ghost=" << pn->ghostbranch() << ",pn unary=" << pn->isUnary() << ")" << std::endl;
        if (pn->leaf())
        {
            assert(rightPos[r->uniqueId()].count(pn->uniqueId()) > 0);
            closeChild(r, pn, step+1);
            rightPos[dest->uniqueId()][pn->uniqueId()] = step;
            return;
        }
        if (pn->ghostbranch())
            return;

        if (pn->isUnary())
        {
            for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
                splitUnary(*it, dest, r, step);
            return;
        }

        // Assert: pn is non-unary node that points to r
        if (rightPos[r->uniqueId()].count(pn->uniqueId())  == 0)
        {
            debugPrint();
            std::cerr << "r = " << r->uniqueId() << ", pn = " << pn->uniqueId() << std::endl;
        }
        assert(rightPos[r->uniqueId()].count(pn->uniqueId()) > 0);
        forceCloseChild(r, pn, step+1);
        rightPos[dest->uniqueId()][pn->uniqueId()] = step;
    }
    void splitUnary(PointerTree::PointerNode *dest, unsigned step)
    {
#ifdef TE_DEBUG_PRINT
        std::cerr << "splitUnary called dest = " << dest->uniqueId() << " (ghost=" << dest->ghostbranch() << ",unary=" << dest->isUnary() << ")" << std::endl;
#endif
        if (dest->root())
            return;
        assert (dest->isUnary());
        PointerTree::PointerNode *r = dest;
        while (!r->root() && r->isUnary())
            r = r->parentPtr();

        if (rightPos.count(r->uniqueId()) == 0)
        {
            debugPrint();
            std::cerr << "r == " << r->uniqueId() << std::endl;
        }
        assert (rightPos.count(r->uniqueId()) > 0);

        // Assert: r is the lowest common non-unary node
        //std::cerr << "using non-unary r = " << r->uniqueId() << " instead of " << dest->uniqueId() << std::endl;

        if (!dest->root())
            rightPos[r->uniqueId()][dest->uniqueId()] = step;
        
        for (PointerTree::PointerNode::iterator it = dest->begin(); it != dest->end(); ++it)
            splitUnary(*it, dest, r, step);
    }

    void insertGhost(PointerTree::PointerNode *dest, unsigned step)
    { 
#ifdef TE_DEBUG_PRINT
        std::cerr << "insertGhost called dest = " << dest->uniqueId() << " (ghost=" << dest->ghostbranch() << ",unary=" << dest->isUnary() << ")" << std::endl;
#endif
        assert (dest->ghostbranch());
        while (dest->ghostbranch() && !dest->root())
            dest = dest->parentPtr();

        if (dest->isUnary())
            splitUnary(dest, step); // Unary node "dest" becomes non-unary
    }
    
    void truncate(PointerTree::PointerNode *pn, PointerTree::PointerNode *dest, PointerTree::PointerNode *r, unsigned step)
    {
        if (pn->leaf())
        {
            assert(rightPos[dest->uniqueId()].count(pn->uniqueId()) > 0);
            closeChild(dest, pn, step+1);
            rightPos[r->uniqueId()][pn->uniqueId()] = step;
            return;
        }
        if (pn->ghostbranch())
            return;

        if (pn->isUnary())
        {
            for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
                truncate(*it, dest, r, step);
            return;
        }

        // Assert: pn is non-unary node
        assert (rightPos[dest->uniqueId()].count(pn->uniqueId()) > 0);
        closeChild(dest, pn, step+1);
        rightPos[r->uniqueId()][pn->uniqueId()] = step;
    }
    void truncate(PointerTree::PointerNode *dest, unsigned step)
    {
#ifdef TE_DEBUG_PRINT
        std::cerr << "truncate called dest = " << dest->uniqueId() << std::endl;
#endif
        if (dest->root())
            return;
        PointerTree::PointerNode *r = dest;
        while (!r->root() && r->isUnary())
            r = r->parentPtr();

        assert (rightPos.count(r->uniqueId()) > 0);
        assert (rightPos[r->uniqueId()].count(dest->uniqueId()) > 0);
        forceCloseChild(r, dest, step+1);
               
        // Assert: r is the lowest common non-unary node
        for (PointerTree::PointerNode::iterator it = dest->begin(); it != dest->end(); ++it)
            truncate(*it, dest, r, step);
    }

    void truncateGhost(PointerTree::PointerNode *dest, unsigned step)
    {
#ifdef TE_DEBUG_PRINT
        std::cerr << "truncateGhost called dest = " << dest->uniqueId() << std::endl;
#endif
        assert (dest->ghostbranch());
        while (dest->ghostbranch() && !dest->root())
            dest = dest->parentPtr();
        if (dest->root())
            return;

        if (dest->isUnary())
            truncate(dest, step); // node "dest" became unary
    }
    
    void insertMutation(PointerTree::PointerNode *dest, NodeId cid, unsigned step)
    {
#ifdef TE_DEBUG_PRINT
        NodeId uid = dest->uniqueId();
        std::cerr << "insert mutation for " << uid << " " << cid << " for time " << step << std::endl;
#endif
        mutationPos[cid].push_back(step);
    }

    PointerTree::PointerNode * findNonUnaryChild(PointerTree::PointerNode *pn)
    {
        while (pn->size() > 0 && pn->isUnary())
        {
            PointerTree::PointerNode *new_pn = 0;
            unsigned checksum = 0;
            for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
                if (!(*it)->ghostbranch())
                {
                    checksum ++;
                    new_pn = *it;
                }
            assert(checksum == 1);
            assert(new_pn != 0);
            pn = new_pn;
        }
        return pn;
    }
        
    void closeChild(PointerTree::PointerNode *pn, PointerTree::PointerNode *cpn, unsigned lstep, bool recomb = false)
    {
#ifdef TE_DEBUG_PRINT
        std::cerr << "closeChild called for parent " << pn->uniqueId() << " child " << cpn->uniqueId() << std::endl;
#endif
        // Find highest non-unary node below cpn
        if (cpn->isUnary())
            cpn = findNonUnaryChild(cpn);
        
        assert(!cpn->isUnary());
        NodeId cid = cpn->uniqueId();

        // Find lowest node that cid points to
        NodeId uid = pn->uniqueId();
        while (rightPos.count(uid) == 0 || rightPos[uid].count(cid) == 0)
        {
            if (pn->root())
                break;
            pn = pn->parentPtr();
            uid = pn->uniqueId();
        }

        assert (rightPos.count(uid) > 0);
        assert (rightPos[uid].count(cid) > 0);
        
        unsigned rstep = rightPos[uid][cid];
        rightPos[uid].erase(cid);
        if (rightPos[uid].size() == 0)
            rightPos.erase(uid);

        if (rstep != ~0u && lstep <= rstep)
        {
            leftRightPos[uid].push_back(ARGChild(cid, mutationPos[cid], lstep, rstep));
            if (recomb)
                leftRightPos[uid].back().recomb = true;
            mutationPos[cid].clear();
        }
    }

    void forceCloseChild(PointerTree::PointerNode *pn, PointerTree::PointerNode *cpn, unsigned lstep)
    {
#ifdef TE_DEBUG_PRINT
        std::cerr << "forceCloseChild called for parent " << pn->uniqueId() << " child " << cpn->uniqueId() << std::endl;
#endif
        NodeId cid = cpn->uniqueId();

        // Find lowest node that cid points to
        NodeId uid = pn->uniqueId();
        while (rightPos.count(uid) == 0 || rightPos[uid].count(cid) == 0)
        {
            if (pn->root())
                break;
            pn = pn->parentPtr();
            uid = pn->uniqueId();
        }

        //std::cerr << "closeChild determined parent " << uid << " child " << cid << std::endl;
        assert (rightPos.count(uid) > 0);
        assert (rightPos[uid].count(cid) > 0);
        
        unsigned rstep = rightPos[uid][cid];
        rightPos[uid].erase(cid);
        if (rightPos[uid].size() == 0)
            rightPos.erase(uid);

        if (rstep != ~0u && lstep <= rstep)
        {
            leftRightPos[uid].push_back(ARGChild(cid, mutationPos[cid], lstep, rstep));
            mutationPos[cid].clear();
        }
    }

    
    /**
     * Naive text format output (TODO: replace with binary output format)
     */
    void output(unsigned nleaves, InputReader const &inputr)
    {
        // Flush open ranges
        for (std::map<NodeId,std::map<NodeId,unsigned> >::iterator it = rightPos.begin(); it != rightPos.end(); ++it)
            for (std::map<NodeId,unsigned>::iterator itt = it->second.begin(); itt != it->second.end(); ++itt)
            {
                unsigned rstep = rightPos[it->first][itt->first];
                if (rstep != ~0u)
                    leftRightPos[it->first].push_back(ARGChild(itt->first, mutationPos[itt->first], 0, rstep));
            }
        rightPos.clear();
    
        
        // Print header:   number of leaves           largest key value                      number of keys
        std::cout << "ARGraph " << nleaves << ' '  << leftRightPos.rbegin()->first << ' ' << leftRightPos.size() << '\n';
        
        // Print data rows
        for (std::map<NodeId,std::vector<struct ARGChild> >::iterator it = leftRightPos.begin(); it != leftRightPos.end(); ++it)
        {
            NodeId uid = it->first;
            unsigned nmut = 0;
            for (std::vector<struct ARGChild>::iterator itt = it->second.begin(); itt != it->second.end(); ++itt)
                nmut += itt->mutations.size();
                
            std::cout << "parent " << uid << ' ' << it->second.size() << ' ' << nmut << '\n';

            // Output ranges
            for (std::vector<struct ARGChild>::iterator itt = it->second.begin(); itt != it->second.end(); ++itt)
                std::cout << "child " << itt->id << ' ' << itt->lRange << ' ' << itt->rRange << ' ' << inputr.position(itt->lRange) << ' ' << inputr.position(itt->rRange) << ' ' << itt->recomb  << '\n';
            // Output mutations
            for (std::vector<struct ARGChild>::iterator itt = it->second.begin(); itt != it->second.end(); ++itt)
            {
                std::vector<unsigned> &muts = itt->mutations;
                for (std::vector<unsigned>::iterator ittt = muts.begin(); ittt != muts.end(); ++ittt)
                    std::cout << "mutation " << itt->id << ' ' << *ittt << ' ' << inputr.position(*ittt) << '\n';
            }
        }
        
    }

    void debugPrint()
    {
        std::cerr << "debug print called:" << std::endl;
        for (std::map<NodeId,std::map<NodeId,unsigned> >::iterator it = rightPos.begin(); it != rightPos.end(); ++it)
            for (std::map<NodeId,unsigned>::iterator itt = it->second.begin(); itt != it->second.end(); ++itt)
            {
                unsigned rstep = rightPos[it->first][itt->first];
                std::cerr << "parent " << it->first << " child " << itt->first << " from rstep " << rstep << std::endl;
            }
    }
    
private:
    void initRightPos(PointerTree::PointerNode *pn, NodeId parent_uid, unsigned step)
    {
        if (pn->leaf())
        {
            rightPos[parent_uid][pn->uniqueId()] = step;
            return;
        }
        if (pn->ghostbranch())
            return;

        {
            bool unarypath = false;
            PointerTree::PointerNode *unarypn = 0;
            for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
                if ((*it)->nZeros() == pn->nZeros() && (*it)->nOnes() == pn->nOnes())
                {
                    unarypn = *it;
                    unarypath = true;
                }
            
            if (unarypath)
            {
                initRightPos(unarypn, parent_uid, step);
                return;
            }
        }
        
        if (!pn->root())
            rightPos[parent_uid][pn->uniqueId()] = step;
        for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
            initRightPos(*it, pn->uniqueId(), step);
    }

    std::vector<std::pair<PointerTree::PointerNode *,NodeId> > insertBuffer;
    std::map<NodeId,std::map<NodeId,unsigned> > rightPos; // keeps track of right positions (position of child insertion)
    std::map<NodeId,std::vector<unsigned> > mutationPos; // keeps track of mutation positions
    //std::map<NodeId,std::vector<std::pair<NodeId,std::pair<std::vector<unsigned>,std::pair<unsigned,unsigned> > > > > leftRightPos; // final list of mutations and [left,right] positions (positions where active child)
    std::map<NodeId,std::vector<struct ARGChild> > leftRightPos;
    bool init;
};
#endif
