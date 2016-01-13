/**
 * Standalone example of reading --enumerate output from `main` and `main-simple`.
 *
 * This example shows how to read the text format and construct the ARGraph class from it.
 *
 * Notes:
 *
 * 1) If data is enumerated from a VCF, the position and range values refer to VCF row numbers. 
 *    If enumerating plain-text, the position and range values refer to column numbers. 
 * 
 * 2) All id values refer to unique IDs of the nodes (32bit unsigned int). The root node has always id 0.
 *    Leaf nodes are represented by id value in range from 1 to the total number of leaves. 
 *
 * Small usage example:
 *       ./main --input test.input --plaintext --verbose --debug --no-prediction --enumerate | ./enumerate-example
 *
 * For large data, you need to recompile all the software with compiler optimizations (`-O2 -DNDEBUG` recommended)
 * and use commandline options e.g.
 *
 *       ./main --input <(zcat large.vcf.gz) --vcf --verbose --no-prediction --enumerate | ./enumerate-example
 *
 */
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <cassert>
using namespace std;

// Note: disable this for large inputs
#define VALIDATE_STRUCTURE

typedef unsigned NodeId;
typedef unsigned Position;

// Example class structure to be constructed
class ARGraph
{
public:
    struct ARGchild
    {
        NodeId id;      // Unique ID of the node; 0 == root
        Position lRange;  // Range corresponds to VCF row numbers
        Position rRange; 
    };
    
    class ARNode
    {
    public:
        vector<struct ARGchild> child;
        vector<pair<NodeId,Position> > mutation; // id is the unique ID of the node; position is the VCF row number
    };

    ARGraph()
        : nodes(), nleaves(0), rRangeMax(0), ok_(true)
    {
        ok_ = construct();
    }

    bool ok() const
    { return ok_; }

    bool validate()
    {
        if (!ok())
            return false;

        // Validate the data structure
        // here
        assert (nleaves > 0);
#ifdef VALIDATE_STRUCTURE
        for (Position i = 0; i <= rRangeMax; i++)
            if (!traverseCol(i))
                return false;
#else
        cerr << "warning: VALIDATE_STRUCTURE was not defined; no validation performed" << endl;
#endif
        return true;
    }
    
private:

    bool inputError(unsigned nentries, string msg)
    {
        if (nentries == ~0u)
            cerr << "enumerate-example error: input failed while reading header " << msg << "." << endl;
        else
            cerr << "enumerate-example error: input failed at step " << msg << " after " << nentries << " input entries remain." << endl;
        return false;
    }
    
    // Returns false if there are problems reading the input
    bool construct()
    {
        // Parse the header row
        string s;
        cin >> s;
        if (s != "ARGraph")
            return inputError(~0u, "tag");
        cin >> nleaves;
        if (nleaves == 0)
            return inputError(~0u, "nleaves");

        unsigned largestid;
        cin >> largestid;
        if (largestid == 0)
            return inputError(~0u, "largestid");
        nodes.resize(largestid+1);
        
        unsigned nentries = 0;
        cin >> nentries;
        if (nentries == 0)
            return inputError(~0u, "nentries");
        
        // Parse the data rows (in total nentries elements)
        cerr << "Reading " << nentries << " entries with largest id of " << largestid << endl;
        while (nentries--)
        {
            // Read parent header
            cin >> s;
            if (s != "parent")
                return inputError(nentries, "parent");
            NodeId pid = 0;
            cin >> pid; // Parent id
            unsigned nchild = 0;
            cin >> nchild;
            if (nchild == 0)
                return inputError(nentries, "nchild");
            unsigned nmut = 0;
            cin >> nmut;

            // Parse the child ranges
            ARNode arnode;
            while (nchild--)
            {
                cin >> s;
                if (s != "child")
                    return inputError(nentries, "child-tag");
                struct ARGchild argchild;
                cin >> argchild.id; // Child id
                if (argchild.id == 0)
                    return inputError(nentries, "child-id"); // Root cannot appear as a child node
                cin >> argchild.lRange;
                cin >> argchild.rRange;
                // Keeps track of maximum position value
                rRangeMax = rRangeMax > argchild.rRange ? rRangeMax : argchild.rRange;
                arnode.child.push_back(argchild);
            }
            while (nmut--)
            {
                cin >> s;
                if (s != "mutation")
                    return inputError(nentries, "mutation-tag");
                NodeId cid = 0;
                cin >> cid; // Child id
                if (cid == 0)
                    return inputError(nentries, "mutation-id"); // Root cannot appear as a child node
                Position p = 0;
                cin >> p;
                arnode.mutation.push_back(make_pair(cid,p));
            }
            assert (nodes[pid].child.empty());
            assert (nodes[pid].mutation.empty());
            nodes[pid] = arnode;
        }
        return true;
    }    

#ifdef VALIDATE_STRUCTURE
    static bool childStructCompare(struct ARGchild &i, struct ARGchild &j)
    {
        if (i.lRange != j.lRange)
            return (i.lRange < j.lRange); // lRange values in increasing order
        return (i.rRange > j.rRange);     // rRange values in decreasing order
    }

    // Traverse the active tree at column i and check that all the leaf nodes are reachable
    bool traverseCol(Position i)
    {
        // Init
        set<NodeId> reachableLeaf;
        for (NodeId j = 1; j <= nleaves; ++j)
            reachableLeaf.insert(j);
        
        // Recursive traversal
        traverseCol(0, i, reachableLeaf);
        
        // Assert: all the leaves must be reachable
        if (reachableLeaf.size() != 0)
        {
            cerr << "error: validation failed for column i = " << i << endl;
            return false;
        }
        return true;        
    }

    void traverseCol(NodeId id, Position i, set<NodeId> &reachableLeaf)
    {
        // Handle leaf nodes
        if (1 <= id && id <= nleaves)
        {
            // Leaf must still exists in the array
            assert (reachableLeaf.count(id) > 0);
            // Remove leaf from the array to mark it found
            reachableLeaf.erase(id);
            return;
        }
        
        ARNode &arnode = nodes[id];
        for (vector<struct ARGchild>::iterator it = arnode.child.begin(); it != arnode.child.end(); ++it)
            if (it->lRange <= i && i <= it->rRange)
                traverseCol(it->id, i, reachableLeaf);
    }            
#endif
    
    std::vector<ARNode> nodes;
    unsigned nleaves; // Number of leaves
    Position rRangeMax; // Largest position value encountered
    bool ok_;
};

int main()
{
    // Read data from standard input
    ARGraph arg;
    if (!arg.ok())
    {
        cerr << "enumerate-example error: unable to read standard input and construct the ARGraph class!" << endl;
        return 1;
    }
    
    // Validity check
    if (!arg.validate())
    {
        cout << "ARGraph class validation failed" << endl;
        return 1;
    }
    
    cout << "ARGraph class constructed OK" << endl;
    return 0;
}
