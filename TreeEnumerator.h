#ifndef _TreeEnumerator_h_
#define _TreeEnumerator_h_
#include "default.h"
#include "Tree.h"
#include <iostream>
#include <cassert>
#include <map>
#include <utility>

class TreeEnumerator
{
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
        debugPrint();
    }

    ~TreeEnumerator()
    {
        if (rightPos.size() > 0)
            std::cerr << "TreeEnumerator()::~TreeEnumerator() warning: size of rightPos at destructor = " << rightPos.size() << std::endl;
    }
    
    void insertChild(PointerTree::PointerNode *pn, NodeId cid, unsigned rstep)
    {
        //insertBuffer.push_back(std::make_pair(pn, cid));

        PointerTree::PointerNode *r = pn;
        while (!r->root() && r->isUnary())
            r = r->parentPtr();
        
        rightPos[r->uniqueId()][cid] = rstep;        
    }

    void flushInserts(unsigned rstep)
    {
        /*for (std::vector<std::pair<PointerTree::PointerNode *,NodeId> >::iterator it = insertBuffer.begin(); it != insertBuffer.end(); ++it)
        {
            PointerTree::PointerNode *pn = it->first;
            bool unarypath;
            {
                unarypath = false;
                for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
                    if ((*it)->nZeros() == pn->nZeros() && (*it)->nOnes() == pn->nOnes())
                        unarypath = true;
                if (unarypath)
                    pn = pn->parentPtr();
            } while (unarypath && !pn->root());            
            rightPos[pn->uniqueId()][it->second] = rstep;
        }
        insertBuffer.clear();*/
    }

    void splitUnary(PointerTree::PointerNode *pn, PointerTree::PointerNode *dest, PointerTree::PointerNode *r, unsigned step)
    {
        if (pn->leaf())
        {
            assert(rightPos[r->uniqueId()].count(pn->uniqueId()) > 0);
            closeChild(r, pn->uniqueId(), step+1);
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
        assert(rightPos[r->uniqueId()].count(pn->uniqueId()) > 0);
        closeChild(r, pn->uniqueId(), step+1);
        rightPos[dest->uniqueId()][pn->uniqueId()] = step;
        
        /*for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        {
            if ((*it)->ghostbranch())
                continue;
            assert (rightPos[pn->uniqueId()].count((*it)->uniqueId()) > 0);
            closeChild(pn, (*it)->uniqueId(), step+1);
            rightPos[dest->uniqueId()][(*it)->uniqueId()] = step;
            }*/
    }
    void splitUnary(PointerTree::PointerNode *dest, unsigned step)
    {
        std::cerr << "splitUnary called dest = " << dest->uniqueId() << std::endl;
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
        
        for (PointerTree::PointerNode::iterator it = dest->begin(); it != dest->end(); ++it)
            splitUnary(*it, dest, r, step);
    }


    
    void truncate(PointerTree::PointerNode *pn, PointerTree::PointerNode *dest, PointerTree::PointerNode *r, unsigned step)
    {
        if (pn->leaf())
        {
            //debugPrint();
            //std::cerr << "truncate r " << r->uniqueId() << " for pn " << pn->uniqueId() << " and dest " << dest->uniqueId() << std::endl;
            assert(rightPos[dest->uniqueId()].count(pn->uniqueId()) > 0);
            closeChild(dest, pn->uniqueId(), step+1);
            rightPos[r->uniqueId()][pn->uniqueId()] = step;
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

        // Assert: pn is non-unary node
        assert (rightPos[dest->uniqueId()].count(pn->uniqueId()) > 0);
        closeChild(dest, pn->uniqueId(), step+1);
        rightPos[r->uniqueId()][pn->uniqueId()] = step;
        
        /*for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
        {
            assert (rightPos[dest->uniqueId()].count((*it)->uniqueId()) > 0);
            closeChild(pn, (*it)->uniqueId(), step+1);
            rightPos[r->uniqueId()][(*it)->uniqueId()] = step;
            }*/
    }
    void truncate(PointerTree::PointerNode *dest, unsigned step)
    {
        std::cerr << "truncate called dest = " << dest->uniqueId() << std::endl;
        if (dest->root())
            return;
        PointerTree::PointerNode *r = dest;
        while (!r->root() && r->isUnary())
            r = r->parentPtr();

        assert (rightPos.count(r->uniqueId()) > 0);
        // Assert: r is the lowest common non-unary node

        for (PointerTree::PointerNode::iterator it = dest->begin(); it != dest->end(); ++it)
            truncate(*it, dest, r, step);
    }
    
    void insertMutation(NodeId uid, NodeId cid, unsigned step)
    {
        mutationPos[uid][cid] = step;
    }

    /*void closeChild(NodeId uid, NodeId cid, unsigned lstep)
    {
        assert (rightPos.count(uid) > 0);
        assert (rightPos[uid].count(cid) > 0);

        unsigned rstep = rightPos[uid][cid];
        rightPos[uid].erase(cid);
        if (rightPos[uid].size() == 0)
            rightPos.erase(uid);
        
        leftRightPos[uid].push_back(std::make_pair(cid, std::make_pair(lstep, rstep)));
        }*/

    void closeChild(PointerTree::PointerNode *pn, NodeId cid, unsigned lstep)
    {
        std::cerr << "closeChild called for parent " << pn->uniqueId() << " child " << cid << std::endl;
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
        
        leftRightPos[uid].push_back(std::make_pair(cid, std::make_pair(lstep, rstep)));
    }

    
    /**
     * Naive text format output (TODO: replace with binary output format)
     */
    void output()
    {
        // Flush open ranges
        for (std::map<NodeId,std::map<NodeId,unsigned> >::iterator it = rightPos.begin(); it != rightPos.end(); ++it)
            for (std::map<NodeId,unsigned>::iterator itt = it->second.begin(); itt != it->second.end(); ++itt)
            {
                unsigned rstep = rightPos[it->first][itt->first];
                leftRightPos[it->first].push_back(std::make_pair(itt->first, std::make_pair(0, rstep)));
            }
        rightPos.clear();

        
        for (std::map<NodeId,std::vector<std::pair<NodeId,std::pair<unsigned,unsigned> > > >::iterator it = leftRightPos.begin(); it != leftRightPos.end(); ++it)
        {
            NodeId uid = it->first;
            unsigned nmut = 0;
            if (mutationPos.count(uid))
                nmut = mutationPos[uid].size();
            std::cout << "parent " << uid << " " << it->second.size() << " " << nmut << '\n';

            // Output ranges
            for (std::vector<std::pair<NodeId,std::pair<unsigned,unsigned> > >::iterator itt = it->second.begin(); itt != it->second.end(); ++itt)
            {
                NodeId cid = itt->first;
                std::pair<unsigned,unsigned> range = itt->second;
                std::cout << "child " << cid << " " << range.first << " " << range.second << '\n';
            }
            // Output mutations
            if (nmut)
                for (std::map<NodeId,unsigned>::iterator itt = mutationPos[uid].begin(); itt != mutationPos[uid].end(); ++itt)
                {
                    NodeId cid = itt->first;
                    std::cout << "mutation " << cid << " " << itt->second << '\n';
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
        
        rightPos[parent_uid][pn->uniqueId()] = step;
        for (PointerTree::PointerNode::iterator it = pn->begin(); it != pn->end(); ++it)
            initRightPos(*it, pn->uniqueId(), step);
    }

    std::vector<std::pair<PointerTree::PointerNode *,NodeId> > insertBuffer;
    std::map<NodeId,std::map<NodeId,unsigned> > rightPos; // keeps track of right positions (position of child insertion)
    std::map<NodeId,std::map<NodeId,unsigned> > mutationPos; // keeps track of mutation positions
    std::map<NodeId,std::vector<std::pair<NodeId,std::pair<unsigned,unsigned> > > > leftRightPos; // final [left,right] positions (positions where active child)
    bool init;
};
#endif
