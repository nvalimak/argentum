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
#include <limits>
using namespace std;

// Note: disable this for large inputs
#define VALIDATE_STRUCTURE

// Enable this to output range distributions (per node and total sum)
//#define OUTPUT_RANGE_DISTRIBUTION

typedef unsigned NodeId;
typedef unsigned ChildId;
typedef unsigned Position;

// Example class structure to be constructed
class ARGraph
{
public:
	void PrintARG(){
		for (unsigned i = 0; i < nodes.size(); i++){
			cerr << "id = " << i;
			for (unsigned j = 0; j < nodes[i].child.size(); j++){
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
        Position lbp;    // Range over bp's
        Position rbp; 
        bool include;
        bool recomb; // Was cut at position rRange due to recomb 
        ARGchild()
            : id(0), lRange(0), rRange(0), include(true), recomb(false)
        { }
    };
	
	struct Edge
	{
		unsigned childRef;
		double weight;
		double nevents;
	};
	
    struct Mutation
    {
        ChildId id;   // id points to the child array
        Position pos; 
        Position bp;  // Position as bp's
        Mutation(ChildId id_, Position pos_, Position bp_)
            : id(id_), pos(pos_), bp(bp_)
        { }
    };
    
    class ARNode
    {
    public:
        enum event_t { mutation_event = 0, insert_child_event, delete_child_event }; 
        vector<struct ARGchild> child;
        map<Position,NodeId> parent;
        vector<struct Mutation> mutation; 
#ifdef OUTPUT_RANGE_DISTRIBUTION
        map<unsigned,unsigned> rangeDist;
#endif

		bool edgesKnown;
		map<NodeId, pair<double, double> > edgesCh;
		map<NodeId, pair<double, double> > edgesP;
		unsigned edgesChNum, edgesPNum;

        unsigned childToCheck;
        double timestamp;
        bool inStack;
		bool timeAssigned;
        int popl_count;
        int popr_count;
        unsigned long total_npairs;
		//TODO set node range at the construction lNodeRange = min rangeMode?child.lbp:child.lRange over all child nodes
		unsigned lNodeRange, rNodeRange;
		unsigned sliceDegree;
		bool inComponent;
		unsigned idInSlice;
		unsigned clustId;
		bool reached, reset;
        ARNode()
            : child(), mutation(), edgesKnown(false), edgesChNum(0), edgesPNum(0), childToCheck(0), timestamp(-1.0), inStack(false), timeAssigned(false), popl_count(0), popr_count(0), total_npairs(0), sliceDegree(0), inComponent(false), clustId(0), reached(false), reset(true)
        {   }

        // Return parent node at step i and next rRange 
        pair<NodeId,Position> getParent(Position i)
        {
            map<Position,NodeId>::iterator z = parent.lower_bound(i);
            assert (z != parent.end());
            return make_pair(z->second, z->first);
        }
        
        //Return value: first of each pair is event type: mutation - 0, inserte child - 1, delete child - 2; second of each pair is a pointer to event: if mutation - index in ARNode.mutation, otherwise index in ARNode.child
        vector<pair<event_t,unsigned> > getEvents()
        {
            vector<pair<event_t,unsigned> > events;
            for( size_t i = 0; i < mutation.size(); i++ )
                if ( child[mutation[i].id].include )//TODO - prevents from compiling due to referring to "ARGraph.nodes" instead of "ARNode.child"
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
                return rangeMode?mutation[ev.second].bp:mutation[ev.second].pos;
            case insert_child_event:
                return rangeMode?child[ev.second].lbp:child[ev.second].lRange; 
            case delete_child_event:
                return ( rangeMode?child[ev.second].rbp:child[ev.second].rRange ) + 1;
            default:
                assert (0);
            }
            return 0;
        }
        unsigned getID(const pair<event_t,unsigned> &ev)
        {
            switch (ev.first)
            {
            case mutation_event:
                return child[mutation[ev.second].id].id;  
            case insert_child_event:
                return child[ev.second].id; 
            case delete_child_event:
                return child[ev.second].id;
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
            pos2 = getPosition(e2);
            return pos1 < pos2;
        }
    };
    

    
    ARGraph()
        : nodes(), knuthShuffle(), nleaves(0), rRangeMax(0), rbpMax(0), ok_(true), mu(0.00000001), rho(0.00000001)
    {
        ok_ = construct();
    }

	void SetRangeMode(bool mode){
		rangeMode = mode;
		cerr << "Ranges are measured in ";
		if (mode)
			cerr << "BP (basepairs).";
		else
			cerr << "SNP (positions).";
		cerr << endl;
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
        cerr << "warning: VALIDATE_STRUCTURE was defined in enumerate_example.cpp; validating data structure now, which can take large amount of time on large inputs." << endl;
        for (Position i = 0; i <= rRangeMax; i++)
            if (!traverseCol(i))
                return false;
#else
        cerr << "warning: VALIDATE_STRUCTURE was not defined; no validation performed" << endl;
#endif
        return true;
    }
	
	void SetSlice(Position lSlice, Position rSlice, double mint = -1.0, double maxt = -1.0){
		lSliceRange = lSlice;
		rSliceRange = rSlice;
		min_time = mint;
		max_time = maxt;
		if (max_time != -1.0 && min_time >= max_time){
			cerr << "Error in ARGraph::SetSlice(): max_time < min_time" << endl;
			exit(0);
		}
		setNodeRanges();
		for (vector< ARNode >::iterator it = nodes.begin(); it != nodes.end(); ++it){
			NodeId nodeRef = it - nodes.begin();
			if ( isInSlice( nodeRef ) ){
				nodes[nodeRef].idInSlice = SliceNodes.size();
				SliceNodes.push_back(nodeRef);
			}
		}
	}
	
	bool isInSlice(NodeId nodeRef){
		unsigned leftR = nodes[nodeRef].lNodeRange;
		unsigned rightR = nodes[nodeRef].rNodeRange;
		if (rightR < lSliceRange)
			return false;
		if (rSliceRange < leftR)
			return false;
		if (nodes[nodeRef].timestamp <= min_time)
			return false;
		if (max_time != -1.0 && max_time < nodes[nodeRef].timestamp)
			return false;
		return true;
	}

	void setNodeRanges(){
		for (vector< ARNode >::iterator it = nodes.begin(); it != nodes.end(); ++it){
			it->lNodeRange = rangeMode?rbpMax:rRangeMax + 2;
			it->rNodeRange = 0;
			if (it->child.size() == 0){
				it->lNodeRange = 0;
				it->rNodeRange = rangeMode?rbpMax:rRangeMax + 1;
			}
			for (vector< ARGchild >::iterator itt = it->child.begin(); itt != it->child.end(); ++itt){
				if ( it->lNodeRange > rangeMode?itt->lbp:itt->lRange ){
					it->lNodeRange = rangeMode?itt->lbp:itt->lRange;
				}
				if ( it->rNodeRange < rangeMode?itt->rbp:itt->rRange ){
					it->rNodeRange = rangeMode?itt->rbp:itt->rRange;
				}
			}
		}
	}
	
	unsigned VisitComponent(NodeId nodeRef, bool output = false){
		static unsigned nodesInComponents = 0;
		unsigned nodesInOtherComponents = nodesInComponents;
		nodesInComponents++;
		
		nodes[nodeRef].inComponent = true;
		for (map<NodeId, pair<double, double> >::iterator it = nodes[ nodeRef ].edgesCh.begin(); it != nodes[ nodeRef ].edgesCh.end(); ++it){
			if ( !isInSlice(it->first) )
				continue;
			if (output)
				cout << nodes[nodeRef].idInSlice << "\t" << nodes[it->first].idInSlice << "\t1\n";
			if ( !nodes[it->first].inComponent )
				VisitComponent(it->first, output);
		}
	
		for (map<NodeId, pair<double, double> >::iterator it = nodes[ nodeRef ].edgesP.begin(); it != nodes[ nodeRef ].edgesP.end(); ++it){
			if ( !isInSlice(it->first) )
				continue;
			if ( !nodes[it->first].inComponent )
				VisitComponent(it->first, output);
		}
		return nodesInComponents - nodesInOtherComponents;
	}

	NodeId CheckConnectedness(bool output = false){
		unsigned NumComponents = 0;
		unsigned sliceCurId;
		unsigned maxSize = 0;
		NodeId nodeSeed = 0;
		for (std::vector< NodeId >::iterator it = SliceNodes.begin(); it != SliceNodes.end(); ++it){
			if ( nodes[*it].inComponent )
				continue;
			if (output)
				cout << "Component " << NumComponents << endl;
			NumComponents++;
			sliceCurId = SliceNodes.size();
			unsigned componentSize = VisitComponent( *it, output );
			if (maxSize < componentSize){
				maxSize = componentSize;
				nodeSeed = *it;
			}
			cerr << "Component contains " << componentSize << " nodes." << endl;
		}
		cerr << "Number of components found: " << NumComponents << endl;
		cerr << "Total number of nodes in the slice " << SliceNodes.size() << endl;
		return nodeSeed;
	}
	
	void ResetComponents(){
		for (std::vector< NodeId >::iterator it = SliceNodes.begin(); it != SliceNodes.end(); ++it)
			nodes[*it].inComponent = false;
	}
	
	void OutputSlice(){
		cout << "Slice nodes" << endl;
		for (std::vector< NodeId >::iterator it = SliceNodes.begin(); it != SliceNodes.end(); ++it)
			cout << *it << "\t" << nodes[*it].idInSlice << "\t" << nodes[*it].timestamp << "\t" << nodes[*it].lNodeRange << "\t" << nodes[*it].rNodeRange << endl;
	}
	
	void ReadClust(){
		int nentries;
		NodeId nId;
		unsigned clust;
		cin >> nentries;
		cerr << "Reading " << nentries << " cluster entries." << endl;
        while (nentries--)
        {
            // Read parent header
            cin >> nId;
			cin >> clust;
			nodes[nId].clustId = clust;
		}	
	}
	
	class PaintChunk{
	public:
			NodeId hap;
			Position lRange, rRange;
			NodeId nodeRef;
			
			PaintChunk(NodeId h){
				hap = h;
			}
			
			void Set(unsigned lr, unsigned rr, NodeId nid){
				lRange = lr;
				rRange = rr;
				nodeRef = nid;
			}
	};
	
	void PaintHaplotype(NodeId hap, vector< PaintChunk > &chunks){
		vector< PaintChunk > chunksStack;
		PaintChunk chu(hap), curChu(hap), chuTmp(hap);
		curChu.Set(lSliceRange, rSliceRange, hap);
		chunksStack.push_back(curChu);
		while( chunksStack.size() ){
			curChu = chunksStack.back();
			chunksStack.pop_back();
			for (std::map<NodeId, std::pair<double, double> >::iterator it = nodes[curChu.nodeRef].edgesP.begin(); it != nodes[curChu.nodeRef].edgesP.end(); ++it){
				Position lr = rSliceRange+1, rr = rSliceRange+1;
				for (vector< ARGchild >::iterator itt = nodes[it->first].child.begin(); itt != nodes[it->first].child.end(); ++itt){
					if(itt->id == curChu.nodeRef){
						lr = itt->lbp;
						rr = itt->rbp;
						break;
					}
				}
				if (rr < curChu.lRange || lr > curChu.rRange)
					continue;
				else{
					lr = lr>curChu.lRange?lr:curChu.lRange;
					rr = rr<curChu.rRange?rr:curChu.rRange;
				}
				chuTmp.Set(lr, rr, it->first);
				if (nodes[it->first].clustId == 0){
					chunksStack.push_back(chuTmp);
				}
				else{
					chunks.push_back(chuTmp);
				}
			}
		}
	}
	
	void PaintHaps(){
		vector< PaintChunk > chunks;
		cerr << "Painting haplotypes..." << endl;
		unsigned inClustNodes = 0;
		for (vector< ARNode >::iterator it = nodes.begin(); it != nodes.end(); ++it){
			if (it->clustId != 0)
				inClustNodes++;
		}
		cout << "Number of nodes in clusters: " << inClustNodes << endl;
		for (unsigned i = 1; i <= nleaves; i++)
			PaintHaplotype(i, chunks);
		cout << "Total number of chunks: " << chunks.size() << endl;
		double fin[4], gbr[4], sar[4];
		for (int j = 0; j < 4; j++){
			fin[j] = 0; gbr[j] = 0; sar[j] = 0;
		}
		for (vector< PaintChunk >::iterator it = chunks.begin(); it != chunks.end(); ++it){
			if (it->hap <= 2000)
				fin[ nodes[it->nodeRef].clustId ] += (it->rRange - it->lRange);
			if (2001 <= it->hap && it->hap <= 4000)
				gbr[ nodes[it->nodeRef].clustId ] += (it->rRange - it->lRange);
			if (4001 <= it->hap)
				sar[ nodes[it->nodeRef].clustId ] += (it->rRange - it->lRange);
		}
		double sum = fin[0] + fin[1] + fin[2] + fin[3];
		sum = sum==0?1:sum;
		for (int i = 0; i < 4; i++)
			cout << fin[i]/sum << "\t";
		cout << endl;
		sum = gbr[0] + gbr[1] + gbr[2] + gbr[3];
		sum = sum==0?1:sum;
		for (int i = 0; i < 4; i++)
			cout << gbr[i]/sum << "\t";
		cout << endl;
		sum = sar[0] + sar[1] + sar[2] + sar[3];
		sum = sum==0?1:sum;
		for (int i = 0; i < 4; i++)
			cout << sar[i]/sum << "\t";
		cout << endl;
	}
	
	void NodeImpactDistribution(string filename, unsigned npop, double t_min, double t_max, bool allLeaves = true, double sampleRate = 0.01){
		vector<bool> leaves;
		vector<unsigned> pops;
		unsigned selectedNodes = 0;
		unsigned eligibleNodes = 0;
		ofstream fh;
		fh.open(filename.c_str());
		cerr << "NodeImpactDistribution() called" << endl;
		if (nleaves % npop != 0){
			cerr << "nleaves is not devisible by npop" << endl;
			exit(0);
		}
		for (unsigned i = 0; i < nleaves; i++)
			leaves.push_back(false);
		DebugReset(leaves);
		for (unsigned i = 0; i < npop; i++)
			pops.push_back(0);
		for (vector< ARNode >::iterator it = nodes.begin() + nleaves + 1; it != nodes.end(); ++it){
			if (it->timestamp < t_min || it->timestamp > t_max)
				continue;
			if (it->lNodeRange > rSliceRange || it->rNodeRange < lSliceRange)
				continue;
			eligibleNodes++;
			if (rand() % 1 > sampleRate)
				continue;
			selectedNodes++;
			FindReachableLeaves( it - nodes.begin(), leaves, allLeaves );
			unsigned sum = ComputePopImpact( leaves, npop, pops );
			fh << sum ;
			for (unsigned i = 0; i < npop; i++){
				fh << "\t" << double(pops[i])/double(sum);
			}
			fh << "\n";
			ResetFlagsReachableLeaves( it - nodes.begin(), leaves );
			DebugReset(leaves);
			for (unsigned i = 0; i < npop; i++)
				pops[i] = 0;
		}
		fh.close();
		cerr << selectedNodes << " nodes sampled for node impact distribution out of " << eligibleNodes << "within time period " << t_min << " to " << t_max << "." << endl;
	}
	
//	gunzip -c data/hrc_chr20_ARG_subset.txt.gz | ./enumerate-example data/pop_map.txt 1 1 200 1 3 1 0 5 1000 > tmp.txt
	
	void DebugReset(vector<bool>& leaves){
		for (vector< ARNode >::iterator it = nodes.begin(); it != nodes.end(); ++it){
			if ( it->reached ){
				cerr << "reached not reset" << endl ;
				exit(0);
			}
			if ( !it->reset ){
				cerr << "reset not reset" << endl ;
				exit(0);
			}	
		}
		for (unsigned i = 0; i < nleaves; i++){
			if ( leaves[i] ){
				cerr << i << " id, " << leaves[i] << endl;
				cerr << "leaf not reset" << endl ;
				exit(0);
			}
		}
	}
	
	unsigned ComputePopImpact( vector<bool>& leaves, unsigned npop, vector<unsigned>& pops ){
		unsigned sum = 0;
		for (unsigned i = 0; i < nleaves; i++){
			if (!leaves[i])
				continue;
			unsigned popSize = nleaves/npop;
			unsigned j = i/popSize;
			if (j >= npop ){
				cerr << "ComputePopImpact() segmentation fault." << endl;
				cerr << "leaf id-1 = " << i << "\tpopulation id " << j << endl;
				exit(0);
			}
			pops[j]++;
			sum++;
		}
		return sum;
	}
	
	void FindReachableLeaves( NodeId nodeRef, vector<bool>& leaves, bool allLeaves){
		vector<NodeId> stack;
		stack.push_back(nodeRef);
		nodes[nodeRef].reached = true;
		nodes[nodeRef].reset = false;
		unsigned lr = nodes[nodeRef].lNodeRange;
		unsigned rr = nodes[nodeRef].rNodeRange;
		while( stack.size() ){
			NodeId curNode = stack.back();
			stack.pop_back();
			for (vector< ARGchild >::iterator itt = nodes[curNode].child.begin(); itt != nodes[curNode].child.end(); ++itt){
				if (itt->id <= nleaves && itt->id != 0){
					leaves[itt->id - 1] = true;
					continue;
				}
				if (!allLeaves){
					if ( nodes[itt->id].lNodeRange > rr || nodes[itt->id].rNodeRange < lr)
						continue;
				}
				if (nodes[itt->id].reached)
					continue;
				stack.push_back(itt->id);
				nodes[itt->id].reached = true;
				nodes[itt->id].reset = false;
			}
		}
	}
	
	void ResetFlagsReachableLeaves( NodeId nodeRef, vector<bool>& leaves ){
		vector<NodeId> stack;
		stack.push_back(nodeRef);
		nodes[nodeRef].reached = false;
		nodes[nodeRef].reset = true;
		while( stack.size() ){
			NodeId curNode = stack.back();
			stack.pop_back();
			for (vector< ARGchild >::iterator itt = nodes[curNode].child.begin(); itt != nodes[curNode].child.end(); ++itt){
				if (itt->id <= nleaves && itt->id != 0){
					leaves[itt->id - 1] = false;
					continue;
				}
				if (nodes[itt->id].reset)
					continue;
				stack.push_back(itt->id);
				nodes[itt->id].reached = false;
				nodes[itt->id].reset = true;
			}
		}
	}

    void assignTimes(int method)
    {
        std::vector<NodeId> nodeStack;
        nodeStack.push_back(0); //root node
        while( nodeStack.size() )
        {
            NodeId nodeRef = nodeStack.back();
            if( nodes[ nodeRef ].childToCheck == nodes[ nodeRef ].child.size() )
            {
				switch(method){
					case 1:
						assignTime(nodeRef);
						break;
					case 2:
						assignTime2(nodeRef);
						break;
					default:
						cerr << "Unknown assign time method." << endl;
						exit(0);
						break;
				}
				nodes[ nodeRef ].timeAssigned = true;
                nodes[ nodeRef ].inStack = false;
                nodeStack.pop_back();
            }
            else
            {
				NodeId childRef = nodes[ nodeRef ].child[ nodes[nodeRef].childToCheck ].id;
				nodes[nodeRef].childToCheck++;
				if (nodes[childRef].timeAssigned)
					continue;
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

    void iterateTimes(bool output = false)//1 for f1, 2 for f2, 3 for f3
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
        for (unsigned i = 1; i < nodes.size(); ++i)
        {
            unsigned j = rand() % (nodes.size()-1) + 1; // Not to Skip root
            unsigned tmp = knuthShuffle[i];
            knuthShuffle[i] = knuthShuffle[j]; // Swap values
            knuthShuffle[j] = tmp;
        }

		it_norm_abs = 0.0;
		NodeId it_norm_abs_node = 0;
		it_norm_rel = 0.0;
		NodeId it_norm_rel_node = 0;
		mean_abs_change = 0.0;
		mean_rel_change = 0.0;


		for (unsigned i = 1; i < nodes.size(); ++i){
			NodeId nodeRef = knuthShuffle[i];
			double old_time = nodes[nodeRef].timestamp;
			UpdateTime(nodeRef);
			double abs_change = nodes[nodeRef].timestamp - old_time;
			double rel_change = 0;
			if (nodes[nodeRef].timestamp != 0)
				rel_change = abs(nodes[nodeRef].timestamp - old_time)/old_time;
			mean_abs_change += abs(abs_change);
			mean_rel_change += rel_change;
			if ( abs(abs_change) > abs(it_norm_abs) ){
				it_norm_abs = abs_change;
				it_norm_abs_node = i;
			}
			if ( rel_change > it_norm_rel ){
				it_norm_rel = rel_change;
				it_norm_rel_node = i;
			}
		}

		if (output){
			cerr << "Iteration max  norm: absolute = " << it_norm_abs << " for node = " << it_norm_abs_node << endl;
			cerr << "Iteration max  norm: relative = " << it_norm_rel << " for node = " << it_norm_rel_node << endl;
			cerr << "Iteration mean norm: absolute = " << mean_abs_change/nodes.size() << "\trelative = " << mean_rel_change/nodes.size() << endl;
		}
    }
	
	void timeUpdateReset(){ //Can be done in timeAssign()? TODO
		for (std::vector<ARNode>::iterator it = nodes.begin(); it != nodes.end(); it++){
			assert(!it->inStack);
			it->childToCheck = 0;
		}
	}
	
	void getMutAndRecombNumber(){
		nrecomb = 0;
		nmutation = 0;
		for (vector< ARNode >::iterator it = nodes.begin() + 1; it != nodes.end(); ++it){
			for (vector< ARGchild >::iterator itt = it->child.begin(); itt != it->child.end(); ++itt)
				nrecomb += itt->recomb;
			nmutation += it->mutation.size();
		}
	}
	
	void initializeEdges(bool exclude = false){
		double recombWeight, penalty = 1.0;
		getMutAndRecombNumber();
		recombWeight = float(nmutation)/float(nrecomb)*rho/mu*penalty;
		PDrate = penalty*rho+mu;
//		recombWeight = 1;
//		rate = mu;
		cerr << "rate = " << PDrate << "\t recombWeight = " << recombWeight << endl;
		for (vector< ARNode >::iterator it = nodes.begin() + 1; it != nodes.end(); ++it){
			for (vector< ARGchild >::iterator itt = it->child.begin(); itt != it->child.end(); ++itt){
				unsigned leftR = rangeMode?itt->lbp:itt->lRange;
				unsigned rightR = rangeMode?itt->rbp:itt->rRange;
				if (exclude && !itt->include)
					continue;
				if (rightR - leftR + 1 == 0)
					continue;
				else
					it->edgesCh[ itt->id ].second += PDrate*(rightR - leftR + 1);
				if (itt->recomb)
					it->edgesCh[ itt->id ].first += recombWeight;
			}
			it->edgesChNum = it->edgesCh.size();
			for (std::map<NodeId, std::pair<double, double> >::iterator itt = it->edgesCh.begin(); itt != it->edgesCh.end(); ++itt)
			{
				if (itt->second.second == 0){
					it->edgesCh.erase(itt);
					assert(false);
				}
			}
	//		extract mutations
			for (vector<struct Mutation>::iterator itt = it->mutation.begin(); itt != it->mutation.end(); ++itt){
				if (it->edgesCh.find(it->child[itt->id].id) == it->edgesCh.end())
					continue;
				it->edgesCh[ it->child[itt->id].id ].first++;
			}
	//		add	information to child
			for (std::map<NodeId, std::pair<double, double> >::iterator itt = it->edgesCh.begin(); itt != it->edgesCh.end(); ++itt){
				nodes[ itt->first ].edgesP[ it - nodes.begin() ].first = itt->second.first;
				nodes[ itt->first ].edgesP[ it - nodes.begin() ].second = itt->second.second;
			}
		}
	}
	
	void CountConflictEdges(){
		unsigned counter = 0, nedges = 0;
		for (vector< ARNode >::iterator it = nodes.begin() + 1; it != nodes.end(); ++it)
			for (std::map<NodeId, std::pair<double, double> >::iterator itt = it->edgesCh.begin(); itt != it->edgesCh.end(); ++itt){
				nedges++;
				if (it->timestamp < nodes[ itt->first ].timestamp)
					counter++;
			}
		cerr << "Number of conflicts is " << counter << " of " << nedges << endl;
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

    Position updatePopCounters(NodeId leaf, Position i, int direction, bool popl, bool popr)
    {
        Position min_rRange = ~0u;
        NodeId j = leaf;        
        while (j != 0)
        {            
            pair<NodeId,Position> pinfo = nodes[j].getParent(i); // Returns pointer and rRange for the parent at pos i
            min_rRange = min(min_rRange, pinfo.second);
            if (popl)
                nodes[pinfo.first].popl_count += direction;
            if (popr)
                nodes[pinfo.first].popr_count += direction;

            j = pinfo.first; // Move upwards in the tree
        }
        assert(min_rRange != ~0u);
        return min_rRange;
    }

    // Returns for leaf node their next update position (next after position 0)
    map<Position,vector<NodeId> > initPopInfo(unsigned popl, unsigned popr, map<unsigned,unsigned> const &popmap)
    {
        map<Position,vector<NodeId> > updateNext;
        for (unsigned i = 1; i <= nleaves; ++i)
        {
            if (popmap.count(i) == 0)
                continue; // Skip leaves that are not in the given population map
            if (popl != popmap.at(i) && popr != popmap.at(i))
                continue; // Skip leaves that are not in popl and popr
            
            // Find the parent's rRange value
            Position rRange = updatePopCounters(i, 0, +1, popl == popmap.at(i), popr == popmap.at(i));
            updateNext[rRange].push_back(i);
        }
        return updateNext;
    }
        
    // Collect population popl vs popr information for the leaf nodes in the given popmap
    void outputPopInfo(unsigned popl, unsigned popr, map<unsigned,unsigned> const &popmap)
    {
        assert (popl > 0);
        assert (popr > 0);

        // Reset counts
        for (std::vector<ARNode>::iterator it = nodes.begin(); it != nodes.end(); ++it)
        {
            it->popl_count = 0;
            it->popr_count = 0;
            it->total_npairs = 0;
        }

        // Initialize counts and collect next update information (nextUpdate tells, for each leaf, the next update position)
        map<Position,vector<NodeId> > nextUpdate = initPopInfo(popl, popr, popmap);
        // Iterate over all positions
        Position i = 0;
        while (i <= rRangeMax)
        {
            if (i % 1000 == 0)
                cerr << "at step i = " << i << endl;

            // Assert: Sum of counts at root must be == number of leaves (in the population map)
            if (nodes[0].popl_count + nodes[0].popr_count != (int)popmap.size() && nodes[0].popl_count + nodes[0].popr_count != (int)popmap.size()*2)
            {
                cerr << "At step i = " << i << endl;
                cerr << "popl = " << nodes[0].popl_count << ", popr = " << nodes[0].popr_count << ", popmapsize = " << (int)popmap.size() << endl;
            }
            assert (nodes[0].popl_count + nodes[0].popr_count == (int)popmap.size() || (popl == popr && nodes[0].popl_count + nodes[0].popr_count == (int)popmap.size()*2));
            assert (nextUpdate.size() > 0);
            Position lRange = i;
            Position rRange = nextUpdate.begin()->first;

            /**
             * Process the range [lRange, rRange] and corresponding popl_count & popr_count values here
             *
             * Increments the total_npairs value of each node by (rRange - lRange + 1) * (number of pairs).
             */
            for (unsigned j = 0; j < nodes.size(); ++j)
                if (nodes[j].popl_count > 0 && nodes[j].popr_count > 0)
                {
                    unsigned long npairs = 0;
                    for (vector<struct ARGchild>::iterator it = nodes[j].child.begin(); it != nodes[j].child.end(); ++it)
                        if (it->lRange <= lRange && rRange <= it->rRange)
                            for (vector<struct ARGchild>::iterator jt = it + 1; jt != nodes[j].child.end(); ++jt)
                                if (jt->lRange <= lRange && rRange <= jt->rRange)
                                    npairs += (nodes[it->id].popl_count * nodes[jt->id].popr_count) + (nodes[it->id].popr_count * nodes[jt->id].popl_count);

                    // (number of base-pairs) * (number of pairs)
                    nodes[j].total_npairs += (unsigned long)(mapToBp[rRange] - mapToBp[lRange] + 1) * npairs;
                }
            
            // Update popl_count and popr_count values 
            for (vector<NodeId>::iterator it = nextUpdate[rRange].begin(); it != nextUpdate[rRange].end(); ++it)
            {
                if (popmap.count(*it) == 0)
                    cerr << "error at i = " << i << ", it = " << *it << " not in popmap" << endl;
                Position rr = updatePopCounters(*it, rRange, -1, popl == popmap.at(*it), popr == popmap.at(*it));
                assert (rr == rRange);
                if (rRange + 1 <= rRangeMax)
                {
                    rr = updatePopCounters(*it, rRange + 1, +1, popl == popmap.at(*it), popr == popmap.at(*it));
                    assert (rr > rRange);
                    nextUpdate[rr].push_back(*it);
                }
            }
            nextUpdate.erase(rRange); // Processing complete up to the position rRange.
            i = rRange + 1;
        }

        cout << "nodeid\ttimestamp\ttotal_npairs\n";
        for (unsigned j = 0; j < nodes.size(); ++j)
            if (nodes[j].total_npairs > 0)
                cout << j << '\t' << nodes[j].timestamp << '\t' << nodes[j].total_npairs << '\n';
    }
    
    void outputTimes(NodeId x, NodeId y)
    {
        assert (x > 0);
        assert (y > 0);
        assert (x <= nleaves);
        assert (y <= nleaves);

        NodeId prev_lca = ~0u;
        Position prev_l = ~0u; 
        Position prev_r = ~0u; 
        Position i = 0;
        while (i <= rRangeMax)
        {
            pair<NodeId,Position> lca = getLCATime(i, x, y);
            // Flush value if the new LCA NodeId is different than the previous LCA
            if (prev_lca != ~0u && prev_lca != lca.first) 
                cout << x << '\t' << y << '\t' << prev_l << '\t' << prev_r << '\t' << mapToBp[prev_l] << '\t' << mapToBp[prev_r] << '\t' << prev_lca << '\t' << nodes[prev_lca].timestamp << '\n';
            if (prev_lca != lca.first)
                prev_l = i;  // Previous LCA was the same as current LCA; skip output and share the left boundary
            prev_lca = lca.first;
            prev_r = lca.second;
            i = lca.second + 1;
        }
        // Flush the last value
        cout << x << '\t' << y << '\t' << prev_l << '\t' << prev_r << '\t' << mapToBp[prev_l] << '\t' << mapToBp[prev_r] << '\t' << prev_lca << '\t' << nodes[prev_lca].timestamp << '\n';
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
        mapToBp.reserve(1024);
        
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
                cin >> argchild.lbp;
                cin >> argchild.rbp;
                cin >> argchild.recomb;
                // Keeps track of maximum position value
                rRangeMax = rRangeMax > argchild.rRange ? rRangeMax : argchild.rRange;
				rbpMax = rbpMax > argchild.rbp ? rbpMax : argchild.rbp;

                if (mapToBp.size() < rRangeMax)
                    mapToBp.resize(rRangeMax+1); // Fixme; make it amortized constant time or better
                assert (mapToBp[argchild.lRange] == 0 || mapToBp[argchild.lRange] == argchild.lbp);
                assert (mapToBp[argchild.rRange] == 0 || mapToBp[argchild.rRange] == argchild.rbp);
                mapToBp[argchild.lRange] = argchild.lbp;
                mapToBp[argchild.rRange] = argchild.rbp;
                                
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
                Position pbp = 0;
                cin >> pbp;
                map<Position,pair<size_t,Position> >::iterator it = pos[cid].lower_bound(p);
                assert (pos.count(cid) > 0);
                if (it == pos[cid].end())
                    it = pos[cid].begin();
                pair<size_t,Position> tmp = it->second;
                if (tmp.second <= p && p <= it->first)
                    arnode.mutation.push_back(Mutation(tmp.first, p, pbp));
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
            pair<NodeId,Position> tmp = arnode.getParent(i);
            assert(active_parent == tmp.first);
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


    void assignTime2(NodeId nodeRef)
    {
        if (nodeRef == 0)
        {
            nodes[ nodeRef ].timestamp = -1;
            return;
        }
        if (nodeRef <= nleaves && nodeRef != 0)
        {
            nodes[ nodeRef ].timestamp = 0;
            return;
        }  
        nodes[nodeRef].timestamp = 0;
		
		double ts = -2.0;
		
		for (unsigned i = 0; i < nodes[nodeRef].child.size(); i++){
			if (!nodes[nodeRef].child[i].include)
				continue;
			if (ts < nodes[ nodes[nodeRef].child[i].id ].timestamp)
				ts = nodes[ nodes[nodeRef].child[i].id ].timestamp;
		}	
        nodes[ nodeRef ].timestamp = ts+1.0;
    }

    void assignTime(NodeId nodeRef)
    {
        if (nodeRef == 0)
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
				if (rRange == rangeMode?rbpMax:rRangeMax + 1){
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
	
/*	void ExcludeNodes(){
		for (vector< ARNode >::iterator it = nodes.begin(); it != nodes.end(); ++it){
			NodeId nodeRef = it - nodes.begin();
			double A1 = 0.0, B1 = 0.0, A2 = 0.0, B2 = 0.0;
	        double C1 = 0.0, C2 = 0.0;
		
	        for (map<NodeId, pair<double, double> >::iterator it = nodes[ nodeRef ].edgesCh.begin(); it != nodes[ nodeRef ].edgesCh.end(); ++it){
				A1 += it->second.second;
				B1 += it->second.second*nodes[it->first].timestamp;
				C1 += it->second.first;
			}
		
	        for (map<NodeId, pair<double, double> >::iterator it = nodes[ nodeRef ].edgesP.begin(); it != nodes[ nodeRef ].edgesP.end(); ++it){
				if (it->first == 0)
					continue;
				A2 += it->second.second;
				B2 += it->second.second*nodes[it->first].timestamp;
				C2 += it->second.first;
			}
			if (C1 + C2 == 0)
				continue;
	        
			map<NodeId, pair<double, double> > edgesChNew, edgesPNew;
			assert(edgesChNew.size() == 0 && edgesPNew.size() == 0);
			for (map<NodeId, pair<double, double> >::iterator et = nodes[ nodeRef ].edgesP.begin(); et != nodes[ nodeRef ].edgesP.end(); ++et){
				for (vector<struct ARGchild>::iterator ct = nodes[ et->first ].сhild.begin(); ct != nodes[ et->first ].end(); ++ct){
					if (ct->id != et->first)
						continue;
					for (vector<struct ARGchild>::iterator cсt = nodes[ nodeRef ].сhild.begin(); cct != nodes[ nodeRef ].end(); ++cct){
						if (ct.lbp < cct.rbp || ct.rbp > cct.lbp){
							continue;
						}
						Position leftR, rightR;
						if (rangeMode){
							leftR = ct.lbp>cct.lbp?ct.lbp:cct.lbp;
							rightR = ct.rbp<cct.rbp?ct.rbp:cct.rbp;
					 	}
						else{
							leftR = ct.lRange>cct.lRange?ct.lRange:cct.lRange;
							rightR = ct.rRange<cct.rRange?ct.rRange:cct.rRange;
						}
						edgesChNew[cct->id].first = 0.0;
						edgesChNew[cct->id].second += (rightR-leftR+1)*PDrate ;
						edgesPNew[et->first].first = 0.0;
						edgesPNew[et->first].second += (rightR-leftR+1)*PDrate;
					}
				}
			}
			for (map<NodeId, pair<double, double> >::iterator et = nodes[ nodeRef ].edgesP.begin(); et != nodes[ nodeRef ].edgesP.end(); ++et){
				nodes[ et->first ].edgesCh
			}
		}
		
	}*/
	
	double ComputeConfigProb(double x, double A1, double B1, double C1, double A2, double B2, double C2){
		double prob;
		if (B1/A1 < B2/A2){
			prob = pow(A1*x - B1, C1) * pow(-A2*x + B2, C2);
			prob *= exp(-(A1*x - B1 - A2*x + B2) );
		}
		else{
			prob = pow(-A1*x + B1, C1) * pow(A2*x - B2, C2);
			prob *= exp(-(-A1*x + B1 + A2*x - B2) );
		}
		return prob;
	}
	
	bool CheckDoubleEqulity(double x, double y){
		double diff = abs(x - y);
		double epsilon = max(nextafter(x, numeric_limits<double>::infinity() ) - x, nextafter(y, numeric_limits<double>::infinity() ) - y);
		if (diff < 10*epsilon)
			return true;
		else
			return false;
	}
	
    void UpdateTime1(NodeId nodeRef) //f3 - "children vs parents"
    {		
        if (nodeRef == 0)
        {
            nodes[ nodeRef ].timestamp = -1;
            return;
        }
        if (nodeRef <= nleaves && nodeRef != 0)
        {
            nodes[ nodeRef ].timestamp = 0;
            return;
        }

		double A1 = 0.0, B1 = 0.0, A2 = 0.0, B2 = 0.0;
        double C1 = 0.0, C2 = 0.0;
		
        for (map<NodeId, pair<double, double> >::iterator it = nodes[ nodeRef ].edgesCh.begin(); it != nodes[ nodeRef ].edgesCh.end(); ++it){
			A1 += it->second.second;
			B1 += it->second.second*nodes[it->first].timestamp;
			C1 += it->second.first;
		}
		
        for (map<NodeId, pair<double, double> >::iterator it = nodes[ nodeRef ].edgesP.begin(); it != nodes[ nodeRef ].edgesP.end(); ++it){
			if (it->first == 0)
				continue;
			A2 += it->second.second;
			B2 += it->second.second*nodes[it->first].timestamp;
			C2 += it->second.first;
		}

		if (A2 == 0.0){
			nodes[nodeRef].timestamp = (B1+C1)/A1;
			assert(!isnan(nodes[nodeRef].timestamp));
			return;
		}
		if (C1 == 0)//FIXME?
			C1 += 0.0001;
		if (C2 == 0)//FIXME?
			C2 += 0.0001;
		double new_time = nodes[nodeRef].timestamp;

		
		double a = A1-A2;
		if (a == 0){
			A1 = 2*A1;
			B1 = 2*B1;
			C1 = 2*C1;
			a = A1 - A2;
		}
		if (B1/A1 > B2/A2)
			a = -a;
		double b = -a*(B1/A1 + B2/A2) - (C1 + C2);
		double c = C2*B1/A1 + C1*B2/A2 + a*B1*B2/(A1*A2);
		double D = b * b - 4*a*c;
		if (D < 0.0){
			cerr << "node " << nodeRef << endl;
			cerr << "t = " << nodes[nodeRef].timestamp << endl;
			cerr << "A1 = " << A1 << "\tB1 = " << B1 << "\tC1 = " << C1 << endl;
			cerr << "A2 = " << A2 << "\tB2 = " << B2 << "\tC2 = " << C2 << endl;
			cerr << "B1/A1 = " << B1/A1 << endl;
			cerr << "B2/A2 = " << B2/A2 << endl;
			cerr << "a = " << a << "\tb = " << b << "\tc = " << c << endl;
			cerr << "D = " << D << endl;
			assert(false);
		}
		double x1 = (-b - sqrt(D))/(2*a);
		double x2 = (-b + sqrt(D))/(2*a);
		double lEnd = min(B1/A1, B2/A2);
		double rEnd = max(B1/A1, B2/A2);
		if (lEnd <= x1 && x1 <= rEnd)
			new_time = x1;
		else if (lEnd <= x2 && x2 <= rEnd)
			new_time = x2;
		else{
			cerr << "node " << nodeRef << endl;
			cerr << "t = " << nodes[nodeRef].timestamp << endl;
			cerr << "A1 = " << A1 << "\tB1 = " << B1 << "\tC1 = " << C1 << endl;
			cerr << "A2 = " << A2 << "\tB2 = " << B2 << "\tC2 = " << C2 << endl;
			cerr << "B1/A1 = " << B1/A1 << endl;
			cerr << "B2/A2 = " << B2/A2 << endl;
			cerr << "a = " << a << "\tb = " << b << "\tc = " << c << endl;
			cerr << "D = " << D << endl;
			cerr << "x1 = " << x1 << endl;
			cerr << "x2 = " << x2 << endl;
			assert(false);
			new_time = (C2*B1 + C1*B2)/(C1*A2+C2*A1);
		}
		assert( !isnan(new_time) );
        nodes[nodeRef].timestamp = new_time;
    }
	
    void UpdateTime(NodeId nodeRef) //f3 - "children vs parents"
    {		
        if (nodeRef == 0)
        {
            nodes[ nodeRef ].timestamp = -1;
            return;
        }
        if (nodeRef <= nleaves && nodeRef != 0)
        {
            nodes[ nodeRef ].timestamp = 0;
            return;
        }

		double A1 = 0.0, B1 = 0.0, A2 = 0.0, B2 = 0.0;
        double C1 = 0.0001, C2 = 0.0001;
		
        for (map<NodeId, pair<double, double> >::iterator it = nodes[ nodeRef ].edgesCh.begin(); it != nodes[ nodeRef ].edgesCh.end(); ++it){
			if (nodes[it->first].timestamp < 0)
				continue;
			A1 += it->second.second;
			B1 += it->second.second*nodes[it->first].timestamp;
			C1 += it->second.first;
		}
		
        for (map<NodeId, pair<double, double> >::iterator it = nodes[ nodeRef ].edgesP.begin(); it != nodes[ nodeRef ].edgesP.end(); ++it){
			if (it->first == 0)
				continue;
			if (nodes[it->first].timestamp < 0)
				continue;
			A2 += it->second.second;
			B2 += it->second.second*nodes[it->first].timestamp;
			C2 += it->second.first;
		}

		if (A1 + A2 == 0){
			nodes[nodeRef].timestamp = -2;
			return;
		}
		
		double prob = ComputeConfigProb(nodes[ nodeRef ].timestamp, A1, B1, C1, A2, B2, C2);
		double new_time;
		
		if (A1 == 0.0){//Probably, it'd be better to remove such vertices
			new_time = (B2-C2)/A2;
			assert( !isnan(new_time) );
			assert( !isinf(new_time) );
//			nodes[nodeRef].timestamp = -3;
//			return;
		}
		else if (A2 == 0.0){
			new_time = (B1+C1)/A1;
			assert( !isnan(new_time) );
			assert( !isinf(new_time) );
		}
		else if (C1 == 0.0 && C2 == 0.0 && CheckDoubleEqulity(A1, A2) ){//FIXME? may cause problems!
			new_time = (B1/A1+B2/A2)/2.0;
			assert( !isnan(new_time) );
			assert( !isinf(new_time) );
		}
		else if (CheckDoubleEqulity(A1, A2) && C1*C2 == 0){
			if (C1 == 0)
				new_time = B2/A2;
			else
				new_time = B1/A1;
			assert( !isnan(new_time) );
			assert( !isinf(new_time) );
		}
		else if ( CheckDoubleEqulity(A1, A2) ){
			double numer = A1*C1*B2 + A2*C2*B1;
			double denum = A1*A2*C1 + A1*A2*C2;
			new_time = numer/denum;
			assert( !isnan(new_time) );
			assert( !isinf(new_time) );
		}
		else if (C1 + C2 == 0){
			new_time = B1/A1;
			assert( !isnan(new_time) );
			assert( !isinf(new_time) );
		}
		else if (C1 == 0){
			new_time = (A2*B2-A1*B2-C2*A2)/(A2*(A2-A1) );
			assert( !isnan(new_time) );
			assert( !isinf(new_time) );
		}
		else if (C2 == 0){
			new_time = (A1*B1-A2*B1+C1*A1)/(A1*(A1-A2) );
			assert( !isnan(new_time) );
			assert( !isinf(new_time) );
		}
		else{
			assert(C1*C2 > 0);

			double a = A1-A2;
			if (a == 0){
				A1 = 2*A1;
				B1 = 2*B1;
				C1 = 2*C1;
				a = A1 - A2;
			}
			if (B1/A1 > B2/A2)
				a = -a;
			double b = -a*(B1/A1 + B2/A2) - (C1 + C2);
			double c = C2*B1/A1 + C1*B2/A2 + a*B1*B2/(A1*A2);
			double D = b * b - 4*a*c;
			if (D < 0.0){
				cerr << "node " << nodeRef << endl;
				cerr << "t = " << nodes[nodeRef].timestamp << endl;
				cerr << "A1 = " << A1 << "\tB1 = " << B1 << "\tC1 = " << C1 << endl;
				cerr << "A2 = " << A2 << "\tB2 = " << B2 << "\tC2 = " << C2 << endl;
				cerr << "B1/A1 = " << B1/A1 << endl;
				cerr << "B2/A2 = " << B2/A2 << endl;
				cerr << "a = " << a << "\tb = " << b << "\tc = " << c << endl;
				cerr << "D = " << D << endl;
				assert(false);
			}
			double x1 = (-b - sqrt(D))/(2*a);
			double x2 = (-b + sqrt(D))/(2*a);
			double lEnd = min(B1/A1, B2/A2);
			double rEnd = max(B1/A1, B2/A2);
			if (lEnd <= x1 && x1 <= rEnd)
				new_time = x1;
			else if (lEnd <= x2 && x2 <= rEnd)
				new_time = x2;
			else{
				cerr << "node " << nodeRef << endl;
				cerr << "t = " << nodes[nodeRef].timestamp << endl;
				cerr << "A1 = " << A1 << "\tB1 = " << B1 << "\tC1 = " << C1 << endl;
				cerr << "A2 = " << A2 << "\tB2 = " << B2 << "\tC2 = " << C2 << endl;
				cerr << "A1-A2 = " << A1 - A2 << endl;
				cerr << "B1/A1 = " << B1/A1 << endl;
				cerr << "B2/A2 = " << B2/A2 << endl;
				cerr << "a = " << a << "\tb = " << b << "\tc = " << c << endl;
				cerr << "D = " << D << endl;
				cerr << "x1 = " << x1 << endl;
				cerr << "x2 = " << x2 << endl;
				assert(false);
				new_time = (C2*B1 + C1*B2)/(C1*A2+C2*A1);
			}
			assert( !isnan(new_time) );
			assert( !isinf(new_time) );
		}
		if ( (B1/A1 < B2/A2 && B1/A1 < new_time && new_time < B2/A2) || (B2/A2 < B1/A1 && B2/A2 < new_time && new_time < B1/A1) )
			if ( ComputeConfigProb(new_time, A1, B1, C1, A2, B2, C2) > prob ){
				nodes[nodeRef].timestamp = new_time;
				prob = ComputeConfigProb(new_time, A1, B1, C1, A2, B2, C2);
			}
		if (ComputeConfigProb(B1/A1, A1, B1, C1, A2, B2, C2) > prob){
			nodes[nodeRef].timestamp = B1/A1;
			prob = ComputeConfigProb(B1/A1, A1, B1, C1, A2, B2, C2);
		}
		if (ComputeConfigProb(B2/A2, A1, B1, C1, A2, B2, C2) > prob){
			nodes[nodeRef].timestamp = B2/A2;
			prob = ComputeConfigProb(B2/A2, A1, B1, C1, A2, B2, C2);
		}
    }
    
    pair<NodeId,Position> getLCATime(Position i, NodeId x, NodeId y)
    {
        // Collect all nodes from x to root at position i
        assert (x > 0);
        assert (y > 0);
        map<NodeId,Position> xParents;
        Position min_rRange = ~0u;        
        NodeId j = x;
        while (j != 0)
        {
            pair<NodeId,Position> pinfo = nodes[j].getParent(i); // Returns pointer and rRange for the parent at pos i
            min_rRange = min(min_rRange, pinfo.second);
            xParents[pinfo.first] = min_rRange;
            j = pinfo.first; // Move upwards in the tree
        }
        
        // Similarly, traverse from y to root at position i, and check if LCA is found
        j = y;
        min_rRange = ~0u;        
        while (j != 0)
        {
            pair<NodeId,Position> pinfo = nodes[j].getParent(i); // Returns pointer and rRange for the parent at pos i
            min_rRange = min(min_rRange, pinfo.second);
            if (xParents.count(pinfo.first))
                return make_pair(pinfo.first, min(min_rRange, xParents[pinfo.first]));
            j = pinfo.first; // Move upwards in the tree
        }

        assert (0); // Should not happen (root is always included in xParents)
        return make_pair(0,0);
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
    std::vector<Position> mapToBp; // Map
    unsigned nleaves; // Number of leaves
    Position rRangeMax; // Largest position value encountered
	Position rbpMax;
    bool ok_;
    double mu, rho; //mutation and recombination rates
	unsigned nmutation, nrecomb;
	double it_norm_abs, it_norm_rel;
	double mean_abs_change, mean_rel_change;
	double PDrate;
//	unsigned nedges;
	Position lSliceRange, rSliceRange;
	double min_time, max_time;
	std::vector<NodeId> SliceNodes;
public:
	static bool rangeMode; //false = SNP, true = BP
};
bool ARGraph::rangeMode = true;

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

map<unsigned,unsigned> init_pop_map(const char *fn, unsigned popl, unsigned popr)
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
        if (p != popl && p != popr)
            continue; // Discard other than popl and popr
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
	Position sliceL, sliceR;
	unsigned min_time, max_time;
	unsigned popl = 1, popr = 1;
	map<unsigned,unsigned> popmap;
    srand ( 85871701 );
//	srand(0);
    if (argc < 7)
    {
        cerr << "usage: " << argv[0] << " [rangeMode] [initMethod] [max_iter] [exCycles] [counter] [output_mode] < input > output" << endl;
        cerr << "  where" <<endl;
		cerr << "     rangeMode     - Measuring ranges in SNP(0) or BP(1)." << endl;
		cerr << "     initMethod    - Which time assign method to use (1-2)." << endl;
        cerr << "     max_iter      - How many iterations to make." << endl;
		cerr << "     exCycles      - Exclude cycles for time update (include = 0, exclude = 1)." << endl;
		cerr << "     counter       - Step for output information while iteration process.." << endl;
		cerr << "     output_mode   - Disable(0)/enable(1) the output of population distances." << endl;
        cerr << "     pops_map.txt  - text file that lists pairs of <node id, pop id>" << endl;
        cerr << "     pop1          - Population to compare against" << endl;
        cerr << "     pop2          - another population/same population." << endl;
        cerr << "     input         - standard input (pipe from `./main --enumerate`)" << endl;
        cerr << "     output        - standard output" << endl;
        return 1;
    }

	unsigned rm = atoi_min(argv[1], 0);
	bool rangeMode = false;
	if (rm == 1)
		rangeMode = true;
	unsigned initMethod = atoi_min(argv[2], 1);
    unsigned max_iter = atoi_min(argv[3], 0);
	unsigned ec = atoi_min(argv[4], 0);
	bool excludeCycles = false;
	if (ec == 1)
		excludeCycles = true;
	unsigned counter = atoi_min(argv[5], 1);
	
    unsigned output_mode = atoi_min(argv[6], 0);
	
    if (output_mode > 5)
    {
        cerr << "usage error: [output_mode] must be in the range 0-5." << endl;
        return 1;
    }
	if (output_mode == 1 || output_mode == 2){
		if (argc != 10){
	        cerr << "usage error: [output_mode] = 1, 2 need population information. [pops_map.txt] [pop1] [pop2]" << endl;
	        return 1;
		}
		popl = atoi_min(argv[8], 1);
		popr = atoi_min(argv[9], 1);
		popmap = init_pop_map(argv[7], popl, popr);
		if (popmap.empty())
		{
			cerr << "error: unable to read pop map input" << endl;
        	return 1;
    	}
	}
    if (output_mode == 1)
        cerr << "comparing pairs from pop " << popl << " vs " << popr << endl;
    if (output_mode == 2)
        cerr << "comparing population " << popl << " vs " << popr << endl;
    if (output_mode == 3)
        cerr << "Searching for graph clustering within a slice." << endl;
    if (output_mode == 4)
        cerr << "Painting..." << endl;
    if (output_mode == 5)
        cerr << "Computing node impact..." << endl;
    // Read data from standard input
    ARGraph arg;
    if (!arg.ok())
    {
        cerr << "enumerate-example error: unable to read standard input and construct the ARGraph class!" << endl;
        return 1;
    }
    
    // Bidirectional tree
    arg.assignParentPtrs();
    
    cerr << "ARGraph class constructed OK" << endl;

	if (output_mode == 3 || output_mode == 4 || output_mode == 5){
		if (argc != 11){
	        cerr << "usage error: [output_mode] = 3, 4, 5 need slice parameters. [slice_left] [slice_right] [min_time] [max_time]" << endl;
	        return 1;
		}
		sliceL = atoi_min(argv[7], 0);//35000000;
		sliceR = atoi_min(argv[8], 0);//55000000;
		min_time = atoi_min(argv[9], 0);
		max_time = atoi_min(argv[10], 0);
		if (sliceL > sliceR){
	        cerr << "usage error: [sliceL] should be less or equal to [sliceR] " << endl;
	        return 1;
		}
		arg.SetSlice(sliceL, sliceR, min_time, max_time);
	}

#ifdef OUTPUT_RANGE_DISTRIBUTION
    arg.outputRangeDistributions();
#endif
	arg.SetRangeMode(rangeMode);
//	cerr << "UNCOMMENT TIME ASSIGNMENT + CONFLICT EDGES!!!" << endl;
	cerr << "Assigning times..." << endl;
	arg.assignTimes(initMethod);

    {
        pair<unsigned,unsigned> tmp = arg.numberOfExcludedNodes();
        cerr << "Number of edges = " << tmp.first << ", number of excluded edges = " << tmp.second << endl;
    }
    arg.initializeEdges(excludeCycles);
	arg.CountConflictEdges();
    cerr << "Updating times..." << endl;

    // Iterate time updates max_iter times
    unsigned iter = 0;
    while (iter < max_iter)
    {
        iter ++;
		if (iter % counter == 0 || iter < 3 ){
			cerr << "Iterating times (" << iter << "/" << max_iter << ")" << endl;
			arg.iterateTimes(true);
			arg.CountConflictEdges();
		}
		else
			arg.iterateTimes(false);  // Random traversal
    }

    if (output_mode == 1)
    {
        cerr << "Extracting leaf-pair times..." << endl;
        cout << "x\ty\tlRange\trRange\tlbp\trbp\tnode\ttimestamp\tweight\tprobability\n";
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
    if (output_mode == 2)
    {
        cerr << "Extracting population vs population times..." << endl;
        arg.outputPopInfo(popl, popr, popmap);
    }
	if (output_mode == 3){
		cerr << "Getting slice..." << endl;
		NodeId nodeSeed = arg.CheckConnectedness();
		arg.ResetComponents();
		int compSize = arg.VisitComponent( nodeSeed, true );
		cerr << "Printing component with " << compSize << " nodes seeded at node " << nodeSeed << "." << endl;
		arg.OutputSlice();
	}
	if (output_mode == 4){
		arg.ReadClust();
		arg.PaintHaps();
		return 1;
	}
	if (output_mode == 5){
		cerr << "Sampling nodes for impact distribution..." << endl;
		arg.NodeImpactDistribution("node_impact/output/hrc_chr20_impact_time0_2_pos35_55_pos.txt", 3, 0, 2, false);
	}
    return 0;
}
