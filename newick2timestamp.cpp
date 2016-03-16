/**
 * Helper script to parse Newick trees
 *
 * Usage:
 *         cat newicktree.txt | ./newick2dot dot-prefix
 */
#include "NewickTree.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>
#include <cassert>
#include <cstdlib>
using namespace std;

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
    if (argc != 5)
    {
        cerr << "usage: cat newicktree.txt | " << argv[0] << " <pop_map.txt> <pop1> <pop2> <max rows>" << endl;
        cerr << "     pops_map.txt  - text file that lists pairs of <node id, pop id>" << endl;
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

    unsigned maxr = atoi_min(argv[4], 1);
    
    // Init tree input from <stdin>
    NewickTree predt("-");
    unsigned n = 0;
    unsigned base = 0;
    while (predt.good() && n < maxr)
    {        
        base += predt.validForBases(); // bp position of the tree

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
                
                if (it->first != itt->first && (popl != popr || it->first < itt->first)){
                    // Output format:         leaf X               leaf Y                column       bp position    timestamp of pair X,Y at position n
                    cout << "TIME" << '\t' << it->first << '\t' << itt->first << '\t' << n << '\t' << base << '\t' << predt.getLCATime(it->first, itt->first) << '\n';
				}
                do ++itt;
                while (itt->second != popr && itt != popmap.end());
            }
        }
        
        do
        {
            // Find next inferred tree
            predt.next();
        } while (predt.good() && predt.validForBases() == 0); // Skip if valid only for 0 bases
        n++;
		if (n%1000 == 0)
			cerr << n << " lines processed" << endl;
    }
    cerr << n << " sites, " << base << " predt_bases" << endl;
    return 0;
}
