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
#include <cmath>
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
	void PrintARG(){
		for (int i = 0; i < nodes.size(); i++){
			cerr << "id = " << i;
			for (int j = 0; j < nodes[i].child.size(); j++){
				cerr << "\t" << nodes[i].child[j].id;
			}
			cerr << endl;
		}
	}
	
    struct ARGchild
    {
        NodeId id;      // Unique ID of the node; 0 == root
        Position lRange;  // Range corresponds to VCF row numbers
        Position rRange; 
        bool include;
        bool recomb; // Was cut at position rRange due to recomb 
        ARGchild()
            : id(0), lRange(0), rRange(0), include(true), recomb(false)
        { }
    };
	
	struct Poly
	{
		double t;
		unsigned k;
		double e;
	};
	
	struct Root
	{
		bool success;
		double root;	
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

		bool edgesKnown;
		map<NodeId, pair<int, double> > edgesCh;
		map<NodeId, pair<int, double> > edgesP;
		vector<Poly> polynom;
		
        unsigned childToCheck;
        double timestamp;
        bool inStack;
		double probability;
		double weight;
        ARNode()
            : child(), mutation(), edgesKnown(false), childToCheck(0), timestamp(-1.0), inStack(false), probability(1.0)
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
                return child[ev.second].rRange + 1;
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
        : nodes(), knuthShuffle(), nleaves(0), rRangeMax(0), ok_(true), mu(0.00000001*300), rho(0.00000001*300)
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

    void iterateTimes(bool output = false)
    {
        // Initialize Knuth shuffle on first call to this method
        if (knuthShuffle.empty())
        {
            knuthShuffle.resize(nodes.size());
            for (unsigned i = 0; i < nodes.size(); ++i)
                knuthShuffle[i] = i;
        }
        // Knuth shuffle
        assert (knuthShuffle.size() == nodes.size());
        for (unsigned i = 1*0; i < nodes.size(); ++i)
        {
            unsigned j = rand() % (nodes.size()-1*0) + 1*0; // Not to Skip root
            unsigned tmp = knuthShuffle[i];
            knuthShuffle[i] = knuthShuffle[j]; // Swap values
            knuthShuffle[j] = tmp;
        }
	// timeUpdateReset(); // FIXME Not required anymore!?
		it_norm_abs = 0.0;
		it_norm_rel = 0.0;
		mean_abs_change = 0.0;
		mean_rel_change = 0.0;
        for (unsigned i = 1*0; i < nodes.size(); ++i)
            UpdateTime(knuthShuffle[i]);
		if (output){
			cerr << "Iteration max  norm: absolute = " << it_norm_abs << "\trelative = " << it_norm_rel << endl;
			cerr << "Iteration mean norm: absolute = " << mean_abs_change/nodes.size() << "\trelative = " << mean_rel_change/nodes.size() << endl;
		}
    }

    
	void timeRefine(unsigned iterationNumber){
		initializeEdges();
		for (int i = 0; i < iterationNumber; i++){
			updateTimes();
			if (i > 0 && i%100 == 0)
				cerr << i << " iterations performed." << endl;
		}
	}
	
	void timeUpdateReset(){ //Can be done in timeAssign()? TODO
		for (std::vector<ARNode>::iterator it = nodes.begin(); it != nodes.end(); it++){
			assert(!it->inStack);
			it->childToCheck = 0;
		}
	}

	void initializeEdges(){
		for (vector< ARNode >::iterator it = nodes.begin(); it != nodes.end(); ++it){
			for (vector< ARGchild >::iterator itt = it->child.begin(); itt != it->child.end(); ++itt){
				if (itt->rRange - itt->lRange + 1 == 0)
					continue;
				else
					it->edgesCh[ itt->id ].second += mu*(itt->rRange - itt->lRange + 1);
				if (itt->recomb)
					it->edgesCh[ itt->id ].first ++;
			}
			for (std::map<NodeId, std::pair<int, double> >::iterator itt = it->edgesCh.begin(); itt != it->edgesCh.end(); ++itt)
			{
				if (itt->second.second == 0){
					it->edgesCh.erase(itt);
					assert(false);
				}
			}
	//		extract mutations
			for (vector<pair<NodeId,Position> >::iterator itt = it->mutation.begin(); itt != it->mutation.end(); ++itt){
				if (it->edgesCh.find(itt->first) != it->edgesCh.end())
					it->edgesCh[ itt->first ].first++;
			}
			for (std::map<NodeId, std::pair<int, double> >::iterator itt = it->edgesCh.begin(); itt != it->edgesCh.end(); ++itt)
			{
				if (itt->second.second == 0){
					it->edgesCh.erase(itt);
					assert(false);
				}
			}
	//		add	information to child	
			for (vector< ARGchild >::iterator itt = it->child.begin(); itt != it->child.end(); ++itt){
				nodes[ itt->id ].edgesP[ it - nodes.begin() ] = it->edgesCh[ itt->id ];
			}
		}
	}

    void updateTimes()
    {
        std::vector<NodeId> nodeStack;
        nodeStack.push_back(0); //root node
		timeUpdateReset();
        while( nodeStack.size() )
        {
            NodeId nodeRef = nodeStack.back();
            if( nodes[ nodeRef ].childToCheck == nodes[ nodeRef ].child.size() )
            {
                UpdateTime(nodeRef);
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

	void ComputeWeights(){
		unsigned i;
		for (i = 0; i < nodes.size(); i++){
			nodes[i].probability = ComputeProbability(nodes[ i ].timestamp, i);
		}
		for (i = 0; i < nodes.size(); i++){
			nodes[i].weight = ComputeWeight(i);
		}
	}

    void outputTimes(NodeId x, NodeId y)
    {
        assert (x > 0);
        assert (y > 0);
        assert (x <= nleaves);
        assert (y <= nleaves);
		ComputeWeights();
        for (Position i = 0; i <= rRangeMax; i++)
        {
//            double ts = getLCATime(i, x, y);
			int nodeRef = getLCATime(i, x, y);
            cout << "TIME\t" << x << '\t' << y << '\t' << i << '\t' << nodes[nodeRef].timestamp << '\t' << nodes[nodeRef].weight << '\t' << nodes[nodeRef].probability  <<'\n';
//            cout << "TIME\t" << x << '\t' << y << '\t' << i << '\t' << ts <<'\n';
        }
    }

	int FindNextPoint(int n, NodeId nodeRef){
		while(nodes[nodeRef].polynom[n].k == 0 && n < nodes[nodeRef].polynom.size())
			n++;
		return n;
	}

    void outputDebug_(NodeId nodeRef, const char *outputPrefix)
    {
        // Output timestamps etc.
        {
            char fn[256];
            snprintf(fn, 256, "%s.%u.nodeRef.type.nodeId.timestamp.recomb.range.tsv", outputPrefix, nodeRef);
            ofstream of(fn);
            for (std::map<NodeId, std::pair<int, double> >::iterator it = nodes[nodeRef].edgesCh.begin(); it != nodes[nodeRef].edgesCh.end(); ++it){
                of << nodeRef << '\t' << "child" << '\t' << it->first << '\t' << nodes[ it->first ].timestamp << '\t' <<  it->second.first << '\t' << it->second.second << '\n';
			}
            for (std::map<NodeId, std::pair<int, double> >::iterator it = nodes[nodeRef].edgesP.begin(); it != nodes[nodeRef].edgesP.end(); ++it)
                of << nodeRef << '\t' << "parent" << '\t' << it->first << '\t' << nodes[ it->first ].timestamp << '\t' <<  it->second.first << '\t' << it->second.second << '\n';
        }
        // Output estimated values (Note: following is copy-paste from UpdateTime()
        {
            char fn[256];
            snprintf(fn, 256, "%s.%u.nodeRef.a.b.x.p.tsv", outputPrefix, nodeRef);
            ofstream of(fn);
            
			
/*			nodes[nodeRef].edgesCh[1].first = 1;
			nodes[nodeRef].edgesCh[1].second = 0.5;
			nodes[1].timestamp = 1;
			nodes[nodeRef].edgesCh[2].first = 1;
			nodes[nodeRef].edgesCh[2].second = 1.3;
			nodes[2].timestamp = 2;
			nodes[nodeRef].edgesCh[5].first = 1;
			nodes[nodeRef].edgesCh[5].second = 0.79;
			nodes[5].timestamp = 3;

			for (std::map<NodeId, std::pair<int, double> >::iterator it = nodes[nodeRef].edgesCh.begin(); it != nodes[nodeRef].edgesCh.end(); ++it){
				if (it->first != 1 && it->first != 2 && it->first != 5)
					nodes[nodeRef].edgesCh.erase(it->first);
			}
			nodes[nodeRef].edgesCh.erase(14);
			
			for (std::map<NodeId, std::pair<int, double> >::iterator it = nodes[nodeRef].edgesCh.begin(); it != nodes[nodeRef].edgesCh.end(); ++it){
				cerr << "ch " << it->first << "\t" << nodes[it->first].timestamp << "\t" << it->second.first << "\t" << it->second.second << endl;
			}
			
			nodes[nodeRef].edgesP[15].first = 0;
			nodes[nodeRef].edgesP[15].second = 2.0;
			nodes[15].timestamp = 4;
			nodes[nodeRef].edgesP[17].first = 1;
			nodes[nodeRef].edgesP[17].second = 3.7;
			nodes[17].timestamp = 5;


			for (std::map<NodeId, std::pair<int, double> >::iterator it = nodes[nodeRef].edgesP.begin(); it != nodes[nodeRef].edgesP.end(); ++it){
				cerr << "p " << it->first << "\t" << nodes[it->first].timestamp << "\t" << it->second.first << "\t" << it->second.second << endl;
			}*/
			
			
			
            // Init range vector; second bool is true for parent pointer
            vector<pair<NodeId,bool> > range;
            for (map<NodeId, pair<int, double> >::iterator it = nodes[ nodeRef ].edgesCh.begin(); it != nodes[ nodeRef ].edgesCh.end(); ++it)
                range.push_back(make_pair(it->first, false));
            for (map<NodeId, pair<int, double> >::iterator it = nodes[ nodeRef ].edgesP.begin(); it != nodes[ nodeRef ].edgesP.end(); ++it)
                range.push_back(make_pair(it->first, true));
            
            std::sort(range.begin(), range.end(), compareEvents(this));		

			GetPolynom(nodeRef, range);
			
			for (unsigned i = 0; i < nodes[nodeRef].polynom.size(); i++){
				cerr << nodes[nodeRef].polynom[i].t << "\t" << nodes[nodeRef].polynom[i].k << "\t" << nodes[nodeRef].polynom[i].e << endl;
			}
			
            for (size_t i = 0; i < nodes[nodeRef].polynom.size() - 1; ++i)
            {
                double a = nodes[nodeRef].polynom[i].t;
                double b = nodes[nodeRef].polynom[i+1].t;;
                double new_p = 0;
                double x = 0;
                Root root;
                root = RootBisection(a, b, nodeRef);
                if (!root.success){
                    of << nodeRef << '\t' << a << '\t' << b << '\t' << x << '\t' << "NA" << '\n';
                    continue;
                }
                x = root.root;
                new_p = ComputeProbability(x, nodeRef);
                of << nodeRef << '\t' << a << '\t' << b << '\t' << x << '\t' << new_p << '\n';
            }
            
            for (size_t i = 0; i < nodes[nodeRef].polynom.size(); ++i){
                if (nodes[nodeRef].polynom[i].k != 0)
                    continue;
                double x = nodes[nodeRef].polynom[i].t;
                double new_p = ComputeProbability(x, nodeRef);
                of << nodeRef << '\t' << "NA" << '\t' << "NA" << '\t' << x << '\t' << new_p << '\n';
            }
            
            // FIRST
            Root root = RootNewton(nodeRef, nodes[nodeRef].polynom[0].t, true);//true for -inf, false for +inf
            if (root.success){
                double x = root.root;
                double new_p = ComputeProbability(x, nodeRef);
                of << nodeRef << '\t' << "-inf" << '\t' << nodes[nodeRef].polynom[0].t << '\t' << x << '\t' << new_p << '\n';
            }
            // LAST 
            root = RootNewton(nodeRef, nodes[nodeRef].polynom.back().t, false);//true for -inf, false for +inf
			if (root.success){
				double x = root.root;
				double new_p = ComputeProbability(x, nodeRef);
				of << nodeRef << '\t' << nodes[nodeRef].polynom.back().t << '\t' << "inf" << '\t' << x << '\t' << new_p << '\n';
			}
        }
    }
	
	void GetPolynom(NodeId nodeRef, vector<pair<NodeId, bool> > &range){
		Poly poly;
		nodes[nodeRef].polynom.clear();
		nodes[nodeRef].polynom.reserve(range.size());
		for (unsigned i = 0; i < range.size(); i++){
			NodeId id = range[i].first;
			if (nodes[nodeRef].polynom.empty()){
				poly.t = nodes[id].timestamp;
				if (range[i].second){
					poly.k = nodes[nodeRef].edgesP[id].first;
					poly.e = nodes[nodeRef].edgesP[id].second;
				}
				else{
					poly.k = nodes[nodeRef].edgesCh[id].first;
					poly.e = nodes[nodeRef].edgesCh[id].second;
				}
				nodes[nodeRef].polynom.push_back(poly);
			}
			else{
				if (nodes[id].timestamp == nodes[nodeRef].polynom.back().t){
					if (range[i].second){
						nodes[nodeRef].polynom.back().k += nodes[nodeRef].edgesP[id].first;
						nodes[nodeRef].polynom.back().e += nodes[nodeRef].edgesP[id].second;
					}
					else{
						nodes[nodeRef].polynom.back().k += nodes[nodeRef].edgesCh[id].first;
						nodes[nodeRef].polynom.back().e += nodes[nodeRef].edgesCh[id].second;
					}
				}
				else{
					poly.t = nodes[id].timestamp;
					if (range[i].second){
						poly.k = nodes[nodeRef].edgesP[id].first;
						poly.e = nodes[nodeRef].edgesP[id].second;
					}
					else{
						poly.k = nodes[nodeRef].edgesCh[id].first;
						poly.e = nodes[nodeRef].edgesCh[id].second;
					}
					nodes[nodeRef].polynom.push_back(poly);
				}
			}
		}
	}

    NodeId findNextDebugNode(unsigned nEdges, NodeId nodeRef)
    {
        unsigned npolyk = 0;
        do
        {
            ++nodeRef;
            if (nodeRef >= nodes.size())
                break;
            vector<pair<NodeId,bool> > range;
            for (map<NodeId, pair<int, double> >::iterator it = nodes[ nodeRef ].edgesCh.begin(); it != nodes[ nodeRef ].edgesCh.end(); ++it)
                range.push_back(make_pair(it->first, false));
            for (map<NodeId, pair<int, double> >::iterator it = nodes[ nodeRef ].edgesP.begin(); it != nodes[ nodeRef ].edgesP.end(); ++it)
                range.push_back(make_pair(it->first, true));
            std::sort(range.begin(), range.end(), compareEvents(this));
	    GetPolynom(nodeRef, range);
            // check non-zero poly.k values
            npolyk = 0;
            for (vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it)
                if (it->k != 0)
                    ++npolyk;
        } while (npolyk != nEdges);
        return nodeRef;
    }

    void outputDebug(unsigned nEdges, const char *outputPrefix)
    {
        NodeId curNode = 0;
        do
        {
            curNode = findNextDebugNode(nEdges, curNode);
            if (curNode < nodes.size())
                outputDebug_(curNode, outputPrefix);
        } while (0); // Change to 1 to output *all* nodes with nEdges
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
                cin >> argchild.recomb;
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
        if (nodeRef == 0 && false)
        {
            nodes[ nodeRef ].timestamp = -1;
            return;
        }
        if (nodeRef <= nleaves && nodeRef != 0)
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
//        double mu = 0.00000001, rho = 0.00000001;
        double A = 0, B = 0, C = 0, d = 0, La = 0.0, Lb = 0.0;
        assert(events.size() > 0);
        unsigned lRange = nodes[nodeRef].getPosition(events[0]), rRange = 0;
        std::set<NodeId> activeNodes;
		bool stop = false;
        for(size_t i = 0; i < events.size(); i++)
        {
			if (stop)
				break;
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
				if (rRange == rRangeMax + 1){
					stop = true;
				}
                La += (rRange - lRange) * activeNodes.size();
                d = 0;
                for (set<NodeId>::iterator it = activeNodes.begin(); it != activeNodes.end(); ++it)
                    d += nodes[ *it ].timestamp;
                Lb += (rRange - lRange) * d;
                if (rRange == nodes[nodeRef].getPosition( events[ events.size() - 1 ] ))
					stop = true;
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
		if (A <= 0){
			cerr << "assignTime:\tnan produced at " << nodeRef << ", A = " << A << ", number of active nodes " << activeNodes.size() << ", num of events " << events.size() << endl;
			assert (A > 0);
		}
    }
	
	double ComputeProbability(double x, NodeId nodeRef, bool fast = false)// fast = true: cactorial is dropped for maxima computing FIXME - fast=false
	{
		double prob = 0;
		double lambda;
//		assert(fast);
		for (std::map<NodeId, std::pair<int, double> >::iterator it = nodes[nodeRef].edgesCh.begin(); it != nodes[nodeRef].edgesCh.end(); ++it){
			if (x == nodes[ it->first ].timestamp)
				continue;
			lambda = abs(x - nodes[ it->first ].timestamp )*it->second.second;
//			assert(lambda > 0);
//			cerr << "\tid = " << it->first << "\tx = " << x << "\tt = " << nodes[ it->first ].timestamp << "\te = " << it->second.second << endl;
//			prob = prob* pow(lambda, it->second.first)*exp(-lambda);
			prob += it->second.first*log( lambda) - lambda;
			if (!fast)
				prob = prob/Factorial(it->second.first);
		}
		for (std::map<NodeId, std::pair<int, double> >::iterator it = nodes[nodeRef].edgesP.begin(); it != nodes[nodeRef].edgesP.end(); ++it){
			if (x == nodes[ it->first ].timestamp)
				continue;
			lambda = abs(x - nodes[ it->first ].timestamp )*it->second.second;
//			cerr << "\tlambda = " << lambda << endl;
//			prob = prob* pow(lambda, it->second.first)*exp(-lambda);
			prob += it->second.first*log( lambda) - lambda;
			if (!fast)
				prob = prob/Factorial(it->second.first);
		}
//		prob = pow(prob, 1/(nodes[nodeRef].edgesCh.size()+nodes[nodeRef].edgesP.size() ) );
		return prob;
	}
	
	double ComputeWeight(NodeId nodeRef){
		double weight = 1, edge;
		double lambda;
		double x = nodes[nodeRef].timestamp;
		for (std::map<NodeId, std::pair<int, double> >::iterator it = nodes[nodeRef].edgesCh.begin(); it != nodes[nodeRef].edgesCh.end(); ++it){
			lambda = abs(x - nodes[ it->first ].timestamp )*it->second.second;
			edge = pow(lambda, it->second.first)*exp(-lambda)/Factorial(it->second.first);
			if (edge != 0)
				weight = weight * pow(edge, nodes[ it->first ].probability/edge);
			else
				weight = 0;
		}
		for (std::map<NodeId, std::pair<int, double> >::iterator it = nodes[nodeRef].edgesP.begin(); it != nodes[nodeRef].edgesP.end(); ++it){
			lambda = abs(x - nodes[ it->first ].timestamp )*it->second.second;
			edge = pow(lambda, it->second.first)*exp(-lambda)/Factorial(it->second.first);
			if (edge != 0)
				weight = weight * pow(edge, nodes[ it->first ].probability/edge);
			else
				weight = 0;
		}
//		assert(!isnan(weight));
		return weight;
	}
	
	double Polynomial(double x, NodeId nodeRef, bool side = true){//side true for left, false for right
		double s = 0, result = 0, product = 1;
		bool f = false;
		double norm = 1.0;
		if (nodes[nodeRef].polynom.back().t > 1)
			norm = nodes[nodeRef].polynom.back().t;
		for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
			if (it->k == 0)
				continue;
			if (x != it->t)
				product = product * ( x - it->t )/norm;
			else{
				product = product * it->k/norm;
				f = true;
			}
		}
		if (f)
			return product;
		for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
			if (it->k > 0)
				result += it->k*product/(x - it->t);	
			if (x == it->t && side)
				s -= it->e;
			else if (x == it->t && !side)
				s += it->e;
			else
				s += signum(x - it->t) * it->e;
		}
		double tmp = result;
		result -= product*s;
		if (isnan(result)){
			cerr << "RootNewton: node " << nodeRef << " produces nan." << endl;
			cerr << "RootNewton: prod = " << product << ", s = " << s << ", first term = " << tmp << endl;
			cerr << "RootNewton: x = " << x << ", side = " << side << endl;
			for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
				cerr << "\tk = " << it->k << "\te = " << it->e << "\tt = " << it->t << endl;
			}
			assert(false);
		}
//		assert(!isnan(result) );
		return result;
	}
	
	double PolynomialDerivative(double x, NodeId nodeRef, bool side = true){//true for -inf, false for +inf
		double s = 0, result = 0, product = 1, partSum, secondTerm = 0.0, tmp;
		bool f = false;
		unsigned nonZero = 0;
		for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
			s += it->e;
			if (it->k > 0)
				nonZero++;
			if (nonZero > 1)
				break;
		}
		if (nonZero < 2){
			assert(nonZero==1);
			if (!side)
				return -s;
			else
				return s;
		}
		double norm = 1.0;
		if (nodes[nodeRef].polynom.back().t > 1)
			norm = nodes[nodeRef].polynom.back().t;
		for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it)
		{
			if ( it->k == 0)
				continue;
			if (x != it->t)
				product = product * (x - it->t)/norm;
			else{
				product = product/norm;
				f = true;
			}
		}
		s = 0;
		for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
			s += it->e;
			if (it->k == 0)
				continue;
			if ( it->t == x )
				secondTerm += product;
			else if (!f)
				secondTerm += signum(x - it->t)*product / (x - it->t);
			partSum = 0.0;
			for ( vector<Poly>::iterator itt = nodes[nodeRef].polynom.begin(); itt != nodes[nodeRef].polynom.end(); ++itt){
				if (it->t == itt->t || itt->k == 0)
					continue;
				if ( f ){
					if (it->t != x && itt->t != x )
						continue;
					else
						tmp = x - it->t + x - itt->t;
				}
				else
					tmp = (x - it->t) * (x - itt->t);
				partSum += product/tmp;
				assert(tmp != 0);
			}
			result += it->k * partSum;
		}
		secondTerm = s*secondTerm;
		result -= secondTerm;
		return result;
	}

	Root RootNewton(NodeId nodeRef, double end, bool side = true){//true for -inf, false for +inf
		double epsilon = pow(10, -5); //TODO precision?
		double y;
		Root root;
		double x = end;
		unsigned degree = 0;
		time_t ts, tc;
		time(&ts);
		for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
			if (it->k > 0)
				degree++;
		}
		double coef = pow(-1, degree);
		if (side && coef*PolynomialDerivative(x, nodeRef, side ) >= 0){
			while ( coef*PolynomialDerivative(x, nodeRef, side ) >= 0 && Polynomial(x, nodeRef, side ) > epsilon){ //TODO step with the distance proportional to timestamps delta
				x -= 1000.0;
				time(&tc);
				if (abs(difftime(tc, ts) ) > 2){
					cerr << "node " << nodeRef << " takes too much time. Checkpoint 1." << endl;
					cerr << "RootNewton: x = " << x << ", pol(x) = " << Polynomial(x, nodeRef, side ) << ", der(x) = " << PolynomialDerivative(x, nodeRef, side ) << ", side = " << side << endl;
					for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
						cerr << "\tk = " << it->k << "\te = " << it->e << "\tt = " << it->t << endl;
					}
					assert(false);
				}
			}
		}
		else if (!side && PolynomialDerivative(x, nodeRef, side ) >= 0)
			while ( PolynomialDerivative(x, nodeRef, side ) >= 0 && Polynomial(x, nodeRef, side ) > epsilon){
				x += 1000.0;
				time(&tc);
				if (abs(difftime(tc, ts) ) > 2){
					cerr << "node " << nodeRef << " takes too much time. Checkpoint 2." << endl;
					cerr << "RootNewton: x = " << x << ", pol(x) = " << Polynomial(x, nodeRef, side ) << ", der(x) = " << PolynomialDerivative(x, nodeRef, side ) << ", side = " << side << endl;
					for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
						cerr << "\tk = " << it->k << "\te = " << it->e << "\tt = " << it->t << endl;
					}
					assert(false);
				}
			}
		y = Polynomial(x, nodeRef, side );
		if (y == 0){
			root.success = true;
			root.root = x;
			return root;
		}
		time(&ts);
		//cerr << "RootNewton: new x=" << x << "\ty=" << y << "\ty'=" << PolynomialDerivative(x, nodeRef, side ) << "\tside=" << side << endl;
		int counter = 0;
		int counterDebug = 0;
		while ( abs( y ) > epsilon && counter < 100)
		{
			counter++;
	        double orig_x = x;
	        double tmp = y / PolynomialDerivative(x, nodeRef, side );
			time(&tc);
	        x = x - tmp;
			if (abs(difftime(tc, ts) ) > 2){
//			if (nodeRef == 31587){
//				cerr << "RootNewton: node " << nodeRef << " takes too much time. Checkpoint 3." << endl;
				
				cerr << "RootNewton: x = " << x << ", pol(x) = " << y << ", der(x) = " << PolynomialDerivative(x, nodeRef, side ) << ", side = " << side << endl;
				cerr << x - orig_x << endl;
//				for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
//					cerr << "\tk = " << it->k << "\te = " << it->e << "\tt = " << it->t << endl;
//				}
				counterDebug ++;
				if (counterDebug == 20)
					assert(false);
			}
	        if (x == orig_x)
	        {
	            if (signum(tmp) > 0)  // minimal increment
	                x = nextafter(x, -numeric_limits<double>::infinity());
	            else
	                x = nextafter(x, numeric_limits<double>::infinity());
	            double orig_y = y;
	            y = Polynomial(x, nodeRef, side );
	            if (signum(orig_y) != signum(y))
	            {
	                root.success = true;
	                root.root = orig_x;
	                if ( (side && x <= end) || (!side && x >= end) )
	                    if ( abs(orig_y) > abs(y) )
	                        root.root = x;
	                return root;
	            }
	        }
	        else if (abs(x - orig_x) <= 0.00000001*abs(x) ){
				double orig_y = y;
				y = Polynomial(x, nodeRef, side );
				if ( signum(y) != signum(orig_y) )
					break;
			}
			else {
				y = Polynomial(x, nodeRef, side );
			}
        
	        if ( (side && x > end) || (!side && x < end) ){
	            root.success = false;
	            return root;
	        }
	        assert( orig_x != x );
		}

//		cerr << "\troot = " << x << ", precision = " << y << endl;
		root.success = true;
		root.root = x;
		return root;
	}
	

	Root RootBisection(double a, double b, NodeId nodeRef){
		Root root;
		root.success = false;
		if (a == b){
			return root;
		}
		if (a > b){
			a = a + b;
			b = a - b;
			a = a - b;
		}
		double x = (a + b)/2.0;
		double Fa, Fb, Fx;
		Fa = Polynomial(a, nodeRef, false);
		Fb = Polynomial(b, nodeRef);
		if ( signum(Fa) == signum(Fb) )
			return root;
		Fx = Polynomial(x, nodeRef);		
		double epsilon = min(abs( Fa - Fb )/1000000.0, 0.000000001);
		time_t ts, tc;
		time(&ts);
		while (abs( Fx ) > epsilon) {
			time(&tc);
			if (abs(difftime(tc, ts) ) > 2){
				cerr << "RootNewton: node " << nodeRef << " takes too much time." << endl;
				cerr << "RootNewton: x = " << x << ", Fa = " << Fa << ", Fb = " << Fb << endl;
				for ( vector<Poly>::iterator it = nodes[nodeRef].polynom.begin(); it != nodes[nodeRef].polynom.end(); ++it){
					cerr << "\tk = " << it->k << "\te = " << it->e << "\tt = " << it->t << endl;
				}
				assert(false);
			}
			if (signum(Fa) == signum (Fx) ){
				a = x;
				Fa = Fx;
			}
			else {
				b = x;
				Fb = Fx;
			}
			double orig_x = x;
			x = (a + b)/2.0;
			if (orig_x == x)
				break;
			Fx = Polynomial(x, nodeRef);
		}
//		cerr << "RootBisection: root = " << x << ", precision = " << Fx << endl;
		root.success = true;
		root.root = x;
		return root;
	}
	
	int Factorial(int k){
		int f = 1, i;
		for (i = 2; i < k + 1; i++)
			f = f * i;
		return f;
	}

	int signum(double x){//attention 0!
		if (x < 0)
			return -1;
		else
			return 1;
	}

    void UpdateTime(NodeId nodeRef)
    {
        if (nodeRef == 0 && false)
        {
            nodes[ nodeRef ].timestamp = -1;
            return;
        }
        if (nodeRef <= nleaves && nodeRef != 0)
        {
            nodes[ nodeRef ].timestamp = 0;
            return;
        }

		double A = 0, B = 0;
        int C = 0;
		
        for (map<NodeId, pair<int, double> >::iterator it = nodes[ nodeRef ].edgesCh.begin(); it != nodes[ nodeRef ].edgesCh.end(); ++it){
			A += it->second.second;
			B += it->second.second*nodes[it->first].timestamp;
			C += it->second.first;
		}
        for (map<NodeId, pair<int, double> >::iterator it = nodes[ nodeRef ].edgesP.begin(); it != nodes[ nodeRef ].edgesP.end(); ++it){
			if (it->first == 0)
				continue;
			A += it->second.second;
			B -= it->second.second*nodes[it->first].timestamp;
			C -= it->second.first;
		}
		double new_time = nodes[nodeRef].timestamp;
		assert(A != 0);
		if (A != 0)
			new_time = (B+C)/A;
		if (new_time < 0)
			new_time = 0;

		double abs_change = abs(nodes[nodeRef].timestamp - new_time);
		double rel_change = 0;
		if (nodes[nodeRef].timestamp != 0)
			rel_change = abs(nodes[nodeRef].timestamp - new_time)/nodes[nodeRef].timestamp;
		mean_abs_change += abs_change;
		mean_rel_change += rel_change;
		if ( abs_change > it_norm_abs )
			it_norm_abs = abs_change;
		if ( rel_change > it_norm_rel )
			it_norm_rel = rel_change;
        nodes[nodeRef].timestamp = new_time;
    }

    void UpdateTime1(NodeId nodeRef)
    {
        if (nodeRef == 0 && false)
        {
            nodes[ nodeRef ].timestamp = -1;
            return;
        }
        if (nodeRef <= nleaves && nodeRef != 0)
        {
            nodes[ nodeRef ].timestamp = 0;
            return;
        }

        vector<pair<NodeId,bool> > range;
        for (map<NodeId, pair<int, double> >::iterator it = nodes[ nodeRef ].edgesCh.begin(); it != nodes[ nodeRef ].edgesCh.end(); ++it)
            range.push_back(make_pair(it->first, false));
        for (map<NodeId, pair<int, double> >::iterator it = nodes[ nodeRef ].edgesP.begin(); it != nodes[ nodeRef ].edgesP.end(); ++it)
            range.push_back(make_pair(it->first, true));
        
        std::sort(range.begin(), range.end(), compareEvents(this));
        GetPolynom(nodeRef, range);

        unsigned nonZero = 0;
        for (size_t i = 0; i < nodes[nodeRef].polynom.size(); ++i)
            if (nodes[nodeRef].polynom[i].k != 0)
                ++nonZero;
         
        //maximize ComputeProbability()
        double max_x = nodes[nodeRef].timestamp;
        double max_p = ComputeProbability(max_x, nodeRef, true);
//		cerr << max_p << endl;
        if (nonZero != 0)
            for (size_t i = 0; i < nodes[nodeRef].polynom.size() - 1; ++i)
            {
                double a = nodes[nodeRef].polynom[i].t;
                double b = nodes[nodeRef].polynom[i+1].t;;
                double new_p = 0;
                double x = 0;
                Root root;
                root = RootBisection(a, b, nodeRef);
                if (!root.success){
                    // Probability == NA ?
                    continue;
                }
                x = root.root;
                new_p = ComputeProbability(x, nodeRef, true);
//				cerr << new_p << endl;
                if (max_p < new_p && x >= 0)
                {
                    max_p = new_p;
                    max_x = x;
                }
            }
        
        for (size_t i = 0; i < nodes[nodeRef].polynom.size(); ++i){
            if (nodes[nodeRef].polynom[i].k != 0)
                continue;
            double x = nodes[nodeRef].polynom[i].t;
            double new_p = ComputeProbability(x, nodeRef, true);
//			cerr << new_p << endl;
            if (max_p < new_p && x >= 0)
            {
                max_p = new_p;
                max_x = x;
            }
        }

        if (nonZero != 0)
        {
	        // FIRST
/*	        Root root = RootNewton(nodeRef, nodes[nodeRef].polynom[0].t, true);//true for -inf, false for +inf
	        if (root.success){
	            double x = root.root;
	            double new_p = ComputeProbability(x, nodeRef, true);
	            if (max_p < new_p && x >= 0)
	            {
	                max_p = new_p;
	                max_x = x;
	            }
	        }*/
	        // LAST 
			Root root = RootNewton(nodeRef, nodes[nodeRef].polynom.back().t, false);//true for -inf, false for +inf
	        if (root.success){
	            double x = root.root;
	            double new_p = ComputeProbability(x, nodeRef, true);
	            if (max_p < new_p && x >= 0)
	            {
	                max_p = new_p;
	                max_x = x;
	            }
	        }
		}
		double abs_change = abs(nodes[nodeRef].timestamp - max_x);
		double rel_change = 0;
		if (nodes[nodeRef].timestamp != 0)
			rel_change = abs(nodes[nodeRef].timestamp - max_x)/nodes[nodeRef].timestamp;
		mean_abs_change += abs_change;
		mean_rel_change += rel_change;
		if ( abs_change > it_norm_abs )
			it_norm_abs = abs_change;
		if ( rel_change > it_norm_rel )
			it_norm_rel = rel_change;
        nodes[nodeRef].timestamp = max_x;
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
				return j;
//                return nodes[j].timestamp;
            j = nodes[j].getParent(i); // Move upwards in the tree
        }
		return 0;
//        return nodes[0].timestamp;
    }


    struct compareEvents : std::binary_function<pair<NodeId,bool>,pair<NodeId,bool>,bool>
    {
        compareEvents(ARGraph * p)
            : instance(p)
            { }
        bool operator() (const pair<NodeId,bool>& o1, const pair<NodeId,bool>& o2) {
            return instance->compareEvents_(o1, o2);
        }
        ARGraph * instance;
    };
    bool compareEvents_(const pair<NodeId,bool> &e1, const pair<NodeId,bool> &e2)
    {
        double ts1 = nodes[e1.first].timestamp;
        double ts2 = nodes[e2.first].timestamp;
        return ts1 < ts2;
    }
    
    
    std::vector<ARNode> nodes;
    std::vector<unsigned> knuthShuffle;
    unsigned nleaves; // Number of leaves
    Position rRangeMax; // Largest position value encountered
    bool ok_;
    double mu, rho; //mutation and recombination rates
	double it_norm_abs, it_norm_rel;
	double mean_abs_change, mean_rel_change;
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
    srand ( time(NULL) );
//	srand(0);
    if (argc != 6)
    {
        cerr << "usage: " << argv[0] << " [pairs.txt] [pop1] [pop2] [max_iter] [dis_out] [n_edges] [output_prefix] < input > output" << endl;
        cerr << "  where" <<endl;
        cerr << "     pops_map.txt  - text file that lists pairs of <node id, pop id>" << endl;
        cerr << "     pop1          - Population to compare against" << endl;
        cerr << "     pop2          - another population/same population." << endl;
        cerr << "     max_iter      - How many iterations to make." << endl;
		cerr << "     dis_out       - Disable(0)/enable(1) the output of population distances." << endl;
        cerr << "     n_edges       - Output nodes with exactly <n_edges> edges." << endl;
        cerr << "     output_prefix - Output file prefix." << endl;
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
    unsigned max_iter = atoi_min(argv[4], 0);
	unsigned dis_out = atoi_min(argv[5], 0);
	if (dis_out != 0 && dis_out != 1)
		assert(false);
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

    arg.initializeEdges();
    cerr << "Updating times..." << endl;

    // Output debug information
    // Debug function is disabled for now...
/*    unsigned nEdges = atoi_min(argv[5], 0);
    const char * outputPrefix = argv[6];
    arg.outputDebug(nEdges, outputPrefix);
    return 0;*/

    // Iterate time updates max_iter times
    unsigned iter = 0;
    while (iter < max_iter)
    {
        iter ++;
		if (iter % 1 == 0 ){
			cerr << "Iterating times (" << iter << "/" << max_iter << ")" << endl;
			arg.iterateTimes(true);
		}
		else
			//arg.updateTimes(); // Linear traversal
			arg.iterateTimes();  // Random traversal
    }
	if (dis_out == 1){
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
	}
    
    return 0;
}
