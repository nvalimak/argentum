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
    unsigned npop = 0;
    unsigned first_pop_size = 0;
    ifstream ifs(fn);
    map<unsigned,unsigned> popm;
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
        popm[lid] = p-1;
        npop = max(npop, p);
        if (p == 1)
            first_pop_size ++;
    }
    assert (npop == 2); // default assumption for now
    assert (first_pop_size > 0);
    return popm;
}

int main(int argc, char ** argv)
{
    if (argc != 4)
    {
        cerr << "usage: cat newicktree.txt | " << argv[0] << " <pop_map.txt> <output prefix> <max rows>" << endl;
        cerr << "     pops_map.txt  - text file that lists pairs of <node id, pop id>" << endl;
        return 1;
    }
    map<unsigned,unsigned> popmap = init_pop_map(argv[1]);
    if (popmap.empty())
    {
        cerr << "error: unable to read pop map input" << endl;
        return 1;
    }

    string fnprefix = string(argv[2]);
    unsigned maxr = atoi_min(argv[3], 1);
    
    // Init tree input from <stdin>
    NewickTree predt("-");
    unsigned n = 0;
    unsigned base = 0;
    while (predt.good() && n < maxr)
    {        
        base += predt.validForBases(); // bp position of the tree

        for(map<unsigned,unsigned>::iterator it = popmap.begin(); it != popmap.end(); ++it)
        {
            if (it->second != 0)
                continue;
            map<unsigned,unsigned>::iterator itt = popmap.begin();
            while (itt->second == it->second && itt != popmap.end())
                ++itt;
            while (itt != popmap.end())
            {
                assert(it->second == 0);  // Compare two different pops
                assert(itt->second == 1);
                assert(it->first != itt->first); // and two different leaf nodes

                // Output format:         leaf X               leaf Y                column       bp position    timestamp of pair X,Y at position n
                cout << "TIME" << '\t' << it->first << '\t' << itt->first << '\t' << n << '\t' << base << '\t' << predt.getLCATime(it->first, itt->first) << '\n';
                
                do ++itt;
                while (itt->second == it->second && itt != popmap.end());
            }
        }
        
        do
        {
            // Find next inferred tree
            predt.next();
        } while (predt.good() && predt.validForBases() == 0); // Skip if valid only for 0 bases
        n++;
    }
    cerr << n << " sites, " << base << " predt_bases" << endl;
    return 0;
}
