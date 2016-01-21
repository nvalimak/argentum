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
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <utility>
#include <algorithm>
#include <cassert>
using namespace std;

// Note: disable this for large inputs
#define VALIDATE_STRUCTURE

// Enable this to output range distributions (per node and total sum)
//#define OUTPUT_RANGE_DISTRIBUTION

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
        bool include;
        ARGchild()
            : id(0), lRange(0), rRange(0), include(true)
        { }
    };
    
    class ARNode
    {
    public:
        enum event_t { mutation_event = 0, insert_child_event, delete_child_event }; 
        vector<struct ARGchild> child;
        map<Position,NodeId> parent;
        vector<pair<NodeId,Position> > mutation; // id points to the child array; position is the VCF row number
#ifdef OUTPUT_RANGE_DISTRIBUTION
        map<unsigned,unsigned> rangeDist;
#endif

        unsigned childToCheck;
        double timestamp;
        bool inStack;
        ARNode()
            : child(), mutation(), childToCheck(0), timestamp(-1.0), inStack(false)
        {   }

        // Return parent node at step i
        NodeId getParent(Position i)
        {
            map<Position,NodeId>::iterator z = parent.lower_bound(i);            
            assert (z != parent.end());
            return z->second;
        }
        
        //Return value: first of each pair is event type: mutation - 0, inserte child - 1, delete child - 2; second of each pair is a pointer to event: if mutation - index in ARNode.mutation, otherwise index in ARNode.child
        vector<pair<event_t,unsigned> > getEvents()
        {
            vector<pair<event_t,unsigned> > events;
            for( size_t i = 0; i < mutation.size(); i++ )
                if ( child[mutation[i].first].include )//TODO - prevents from compiling due to referring to "ARGraph.nodes" instead of "ARNode.child"
                    events.push_back(std::make_pair(mutation_event, i));

            for( size_t i = 0; i < child.size(); i++ )
                if (child[i].include)
                {
                    events.push_back(std::make_pair(insert_child_event, i));
                    events.push_back(std::make_pair(delete_child_event, i));
                }
            std::sort(events.begin(), events.end(), compareEvents(this) );
            return events;
        }
        
        unsigned getPosition(const pair<event_t,unsigned> &ev)
        {
            switch (ev.first)
            {
            case mutation_event:
                return mutation[ev.second].second;  
            case insert_child_event:
                return child[ev.second].lRange; 
            case delete_child_event:
                return child[ev.second].rRange;
            default:
                assert (0);
            }
            return 0;
        }

        struct compareEvents : std::binary_function<pair<event_t,unsigned>,pair<event_t,unsigned>,bool>
        {
            compareEvents(ARNode * p)
                : instance(p)
            { }
            bool operator() (const pair<event_t,unsigned>& o1, const pair<event_t,unsigned>& o2) {
                return instance->compareEvents_(o1, o2);
            }
            ARNode * instance;
        };
        bool compareEvents_(const pair<event_t,unsigned> &e1, const pair<event_t,unsigned> &e2)
        {
            unsigned pos1, pos2;
            pos1 = getPosition(e1);
            pos2 = getPosition(e1);
            return pos1 < pos2;
        }
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

    void assignTimes()
    {
        std::vector<NodeId> nodeStack;
        nodeStack.push_back(0); //root node
        while( nodeStack.size() )
        {
            NodeId nodeRef = nodeStack.back();
            if( nodes[ nodeRef ].childToCheck == nodes[ nodeRef ].child.size() )
            {
                assignTime(nodeRef);
                nodes[ nodeRef ].inStack = false;
                nodeStack.pop_back();
            }
            else
            {
                NodeId childRef = nodes[ nodeRef ].child[ nodes[nodeRef].childToCheck ].id;
                nodes[nodeRef].childToCheck++;
                if ( nodes[childRef].inStack )
                {
                    size_t i = nodeStack.size();
                    do
                    {
                        i--;
                        NodeId n = nodeStack[i];
                        nodes[n].child[ nodes[n].childToCheck - 1 ].include = false;
                    } while(i > 0 && nodeStack[i] != childRef);
                }
                else
                {
                    nodeStack.push_back(childRef);
                    nodes[childRef].inStack = true;
                }
            }
        }
    }

    void assignParentPtrs()
    {
        for (vector<ARNode>::iterator it = nodes.begin(); it != nodes.end(); ++it)
            for (vector<struct ARGchild>::iterator itt = it->child.begin(); itt != it->child.end(); ++itt)
                nodes[itt->id].parent[itt->rRange] = std::distance(nodes.begin(), it);
    }

    pair<unsigned,unsigned> numberOfExcludedNodes()
    {
        unsigned nexcluded = 0;
        unsigned nedges = 0;
        for (vector<ARNode>::iterator it = nodes.begin(); it != nodes.end(); ++it)
            for (vector<struct ARGchild>::iterator itt = it->child.begin(); itt != it->child.end(); ++itt)
            {
                nedges ++;
                if (!(itt->include))
                    nexcluded ++;
            }

        return make_pair(nedges,nexcluded);
    }
    
    
#ifdef OUTPUT_RANGE_DISTRIBUTION
    void outputRangeDistributions()
    {
        map<unsigned,unsigned> totals;
        for (vector<ARNode>::iterator it = nodes.begin(); it != nodes.end(); ++it)
            for (map<unsigned,unsigned>::iterator itt = it->rangeDist.begin(); itt != it->rangeDist.end(); ++itt)
            {
                cout << "RANGE\t" << std::distance(nodes.begin(), it) << '\t' << itt->first << '\t' << itt->second << '\n';
                totals[itt->first] += itt->second;
            }
        for (map<unsigned,unsigned>::iterator itt = totals.begin(); itt != totals.end(); ++itt)
            cout << "RANGE\tTOTAL\t" << itt->first << '\t' << itt->second << '\n';
    }
#endif

    void outputTimes(NodeId x, NodeId y)
    {
        assert (x > 0);
        assert (y > 0);
        assert (x <= nleaves);
        assert (y <= nleaves);
        for (Position i = 0; i <= rRangeMax; i++)
        {
            double ts = getLCATime(i, x, y);
            cout << "TIME\t" << x << '\t' << y << '\t' << i << '\t' << ts << '\n';
        }
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
            map<NodeId,map<Position,pair<size_t,Position> > > pos;
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
                pos[argchild.id][argchild.rRange] = make_pair(arnode.child.size(), argchild.lRange);
                arnode.child.push_back(argchild);
#ifdef OUTPUT_RANGE_DISTRIBUTION
                arnode.rangeDist[argchild.rRange-argchild.lRange+1] += 1;
#endif
            }
            while (nmut--)
            {
                cin >> s;
                if (s != "mutation")
                    return inputError(nentries, "mutation-tag");
                NodeId cid = 0;
                cin >> cid; // Child id
                if (pos.count(cid) == 0)
                    return inputError(nentries, "mutation-id"); // Root cannot appear as a child node
                Position p = 0;
                cin >> p;
                map<Position,pair<size_t,Position> >::iterator it = pos[cid].lower_bound(p);
                assert (pos.count(cid) > 0);
                if (it == pos[cid].end())
                    it = pos[cid].begin();
                pair<size_t,Position> tmp = it->second;
                if (tmp.second <= p && p <= it->first)
                    arnode.mutation.push_back(std::make_pair(tmp.first, p));
            }
            // Assert: there is no overlap in node id's
            assert (nodes[pid].child.empty());
            assert (nodes[pid].mutation.empty());
            nodes[pid] = arnode;
        }
        return true;
    }    

#ifdef VALIDATE_STRUCTURE
    // Traverse the active tree at column i and check that all the leaf nodes are reachable
    bool traverseCol(Position i)
    {
        // Init
        set<NodeId> reachableLeaf;
        for (NodeId j = 1; j <= nleaves; ++j)
            reachableLeaf.insert(j);
        
        // Recursive traversal
        traverseCol(0, i, reachableLeaf, 0);
        
        // Assert: all the leaves must be reachable
        if (reachableLeaf.size() != 0)
        {
            cerr << "error: validation failed for column i = " << i << endl;
            return false;
        }
        return true;        
    }

    void traverseCol(NodeId id, Position i, set<NodeId> &reachableLeaf, NodeId active_parent)
    {
        ARNode &arnode = nodes[id];
        if (id > 0)
        {
            // Validate active parent ptr            
            assert(active_parent == arnode.getParent(i));
        }
        
        // Handle leaf nodes
        if (1 <= id && id <= nleaves)
        {
            // Leaf must still exists in the array
            assert (reachableLeaf.count(id) > 0);
            // Remove leaf from the array to mark it found
            reachableLeaf.erase(id);
            return;
        }
        
        for (vector<struct ARGchild>::iterator it = arnode.child.begin(); it != arnode.child.end(); ++it)
            if (it->lRange <= i && i <= it->rRange)
                traverseCol(it->id, i, reachableLeaf, id);
    }            
#endif


    void assignTime(NodeId nodeRef)
    {
        if (nodeRef == 0)
        {
            nodes[ nodeRef ].timestamp = -1;
            return;
        }
        if (nodeRef <= nleaves)
        {
            nodes[ nodeRef ].timestamp = 0;
            return;
        }
        vector<pair<ARNode::event_t,unsigned> > events = nodes[nodeRef].getEvents();
        if (events.size() == 0)
        {
            // include == false for all events under the node nodeRef
            nodes[ nodeRef ].timestamp = 0;
            return;
        }   
        nodes[nodeRef].timestamp = 0;
        double mu = 0.000001, rho = 0.000001;
        double A = 1, B = 0, C = 0, d = 0;
        assert(events.size() > 0);
        unsigned lRange = nodes[nodeRef].getPosition(events[0]), rRange = 0, La = 0, Lb = 0;
        std::set<NodeId> activeNodes;
        
        for(size_t i = 0; i < events.size(); i++)
        {
            switch (events[i].first)
            {
            case ARNode::mutation_event:
                rRange = nodes[nodeRef].getPosition( events[i] );
                La += (rRange - lRange)*activeNodes.size();
                d = 0;
                for (set<NodeId>::iterator it = activeNodes.begin(); it != activeNodes.end(); ++it)
                    d += nodes[ *it ].timestamp;
                Lb += (rRange - lRange) * d;
                    
                A += mu*La;
                B += mu*Lb;
                C += 1;
                lRange = rRange;
                La = 0;
                Lb = 0;
                break;
            
            case ARNode::insert_child_event:
                rRange = nodes[nodeRef].getPosition( events[i] );
                La += (rRange - lRange) * activeNodes.size();
                d = 0;
                for (set<NodeId>::iterator it = activeNodes.begin(); it != activeNodes.end(); ++it)
                    d += nodes[ *it ].timestamp;
                Lb += (rRange - lRange) * d;
                lRange = rRange;
                activeNodes.insert( nodes[ nodeRef ].child[events[i].second].id );
                break;
                
            case ARNode::delete_child_event:
                rRange = nodes[nodeRef].getPosition( events[i] );
                La += (rRange - lRange) * activeNodes.size();
                d = 0;
                for (set<NodeId>::iterator it = activeNodes.begin(); it != activeNodes.end(); ++it)
                    d += nodes[ *it ].timestamp;
                Lb += (rRange - lRange) * d;
                
                A += rho*La;
                B += rho*Lb;
                C += 1;
                activeNodes.erase( nodes[ nodeRef ].child[events[i].second].id );
                lRange = rRange;
                La = 0;
                Lb = 0;
                break;
                
            default:
                assert(0);
            }
        }
        
        nodes[ nodeRef ].timestamp = (B+C)/A;
        double maxTimestamp = 0;
        for(size_t i = 0; i < nodes[ nodeRef ].child.size(); i++)
        {
            double ts = nodes[nodes[ nodeRef ].child[i].id ].timestamp;
            maxTimestamp = maxTimestamp > ts ? maxTimestamp : ts;
        }
    }

    double getLCATime(Position i, NodeId x, NodeId y)
    {
        // Collect all nodes from x to root at position i
        assert (x > 0);
        assert (y > 0);
        set<NodeId> xParents;
        NodeId j = x;
        while (j != 0)
        {
            xParents.insert(j);
            j = nodes[j].getParent(i); // Move upwards in the tree
        }

        // Similarly, traverse from y to root at position i, and check if LCA is found
        j = y;
        while (j != 0)
        {
            if (xParents.count(j))
                return nodes[j].timestamp;
            j = nodes[j].getParent(i); // Move upwards in the tree
        }
        return nodes[0].timestamp;
    }
    
    std::vector<ARNode> nodes;
    unsigned nleaves; // Number of leaves
    Position rRangeMax; // Largest position value encountered
    bool ok_;
};

int atoi_min(char const *value, int min)
{
    std::istringstream iss(value);
    int i;
    char c;
    if (!(iss >> i) || iss.get(c))
    {
        cerr << "error: value " << value << " must be of type <int>, and greater than or equal to " << min << endl;
        std::exit(1);
    }

    if (i < min)
    {
        cerr << "error: value " << value << " must be greater than or equal to " << min << endl;
        std::exit(1);
    }
    return i;
}

map<unsigned,unsigned> init_pop_map(const char *fn)
{
    // Init population map
    ifstream ifs(fn);
    map<unsigned,unsigned> popm;
    map<unsigned,unsigned> popcount;
    if (!ifs.good())
    {
        cerr << "error: unable to read file " << fn << endl;
        return popm;
    }
    while (ifs.good())
    {
        unsigned lid = 0, p = 0;
        ifs >> lid;
        if (!ifs.good())
            break;
        ifs >> p;
        assert (p > 0);
        assert (lid > 0);
        popm[lid] = p;
        popcount[p]+=1;
    }
    cerr << "Read " << popcount.size() << " populations:" << endl;
    for (map<unsigned,unsigned>::iterator it = popcount.begin(); it != popcount.end(); ++it)
        cerr << "  population " << it->first << " with " << it->second << " leaf ids" << endl;
    return popm;
}

int main(int argc, char ** argv)
{
    if (argc != 4)
    {
        cerr << "usage: " << argv[0] << " [pairs.txt] [pop1] [pop2] < input > output" << endl;
        cerr << "  where" <<endl;
        cerr << "     pops_map.txt  - text file that lists pairs of <node id, pop id>" << endl;
        cerr << "     pop1          - Population to compare against" << endl;
        cerr << "     pop2          - another population/same population." << endl;
        cerr << "     input         - standard input (pipe from `./main --enumerate`)" << endl;
        cerr << "     output        - standard output" << endl;
        return 1;
    }

    map<unsigned,unsigned> popmap = init_pop_map(argv[1]);
    if (popmap.empty())
    {
        cerr << "error: unable to read pop map input" << endl;
        return 1;
    }

    unsigned popl = atoi_min(argv[2], 1);
    unsigned popr = atoi_min(argv[3], 1);
    cerr << "comparing pairs from pop " << popl << " vs " << popr << endl;
    
    // Read data from standard input
    ARGraph arg;
    if (!arg.ok())
    {
        cerr << "enumerate-example error: unable to read standard input and construct the ARGraph class!" << endl;
        return 1;
    }
    
    // Bidirectional tree
    arg.assignParentPtrs();
    
    // Validity check
    if (!arg.validate())
    {
        cout << "ARGraph class validation failed" << endl;
        return 1;
    }
    
    cerr << "ARGraph class constructed OK" << endl;

#ifdef OUTPUT_RANGE_DISTRIBUTION
    arg.outputRangeDistributions();
#endif
    
    cerr << "Assigning times..." << endl;
    arg.assignTimes();

    {
        pair<unsigned,unsigned> tmp = arg.numberOfExcludedNodes();
        cerr << "Number of edges = " << tmp.first << ", number of excluded edges = " << tmp.second << endl;
    }
    
    cerr << "Extracting times..." << endl;
    for(map<unsigned,unsigned>::iterator it = popmap.begin(); it != popmap.end(); ++it)
    {
        if (it->second != popl)
            continue;
        map<unsigned,unsigned>::iterator itt = popmap.begin();
        while (itt->second != popr && itt != popmap.end())
            ++itt;
        while (itt != popmap.end())
        {
            assert(it->second == popl);  // Compare exactly these populations
            assert(itt->second == popr);

            if (it->first != itt->first  && (popl != popr || it->first < itt->first))
                arg.outputTimes(it->first, itt->first);
            
            do ++itt;
            while (itt->second != popr && itt != popmap.end());
        }
    }
    
    return 0;
}
