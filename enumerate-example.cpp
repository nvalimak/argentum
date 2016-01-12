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
#include <utility>
#include <cassert>
using namespace std;

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
        : nodes(), nleaves(0), ok_(true)
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
        return true;
    }
    
private:

    bool inputError(string msg)
    {
        cerr << "enumerate-example error: input failed at step " << msg << " after reading " << nodes.size() << " input entries." << endl;
        return false;
    }
    
    // Returns false if there are problems reading the input
    bool construct()
    {
        // Parse the header row
        string s;
        cin >> s;
        if (s != "ARGraph")
            return inputError("header");
        cin >> nleaves;
        if (nleaves == 0)
            return inputError("nleaves");

        unsigned nnodes = 0;
        cin >> nnodes;
        if (nnodes == 0)
            return inputError("nnodes");
        
        // Parse the data rows (in total nnodes elements)
        while (nnodes--)
        {
            // Read parent header
            cin >> s;
            if (s != "parent")
                return inputError("parent");
            NodeId pid = 0;
            cin >> pid; // Parent id
            // Assumes that the root node is the first element
            if (pid == 0 && nodes.size() > 0)
                return inputError("pid");
            unsigned nchild = 0;
            cin >> nchild;
            if (nchild == 0)
                return inputError("nchild");
            unsigned nmut = 0;
            cin >> nmut;

            // Parse the child ranges
            ARNode arnode;
            while (nchild--)
            {
                cin >> s;
                if (s != "child")
                    return inputError("child-tag");
                struct ARGchild argchild;
                cin >> argchild.id; // Child id
                if (argchild.id == 0)
                    return inputError("child-id"); // Root cannot appear as a child node
                cin >> argchild.lRange;
                cin >> argchild.rRange;
                arnode.child.push_back(argchild);
            }
            while (nmut--)
            {
                cin >> s;
                if (s != "mutation")
                    return inputError("mutation-tag");
                NodeId cid = 0;
                cin >> cid; // Child id
                if (cid == 0)
                    return inputError("mutation-id"); // Root cannot appear as a child node
                Position p = 0;
                cin >> p;
                arnode.mutation.push_back(make_pair(cid,p));
            }
            if (nodes.size() <= pid)
                nodes.resize(pid+1);
            nodes[pid] = arnode;
        }
        return true;
    }    

    std::vector<ARNode> nodes;
    unsigned nleaves; // Number of leaves
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
