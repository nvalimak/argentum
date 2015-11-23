/**
 * Helper script to parse Newick trees
 *
 * Usage:
 *         ./main [options] | ./newick2clustersize
 *
 * Input: Predicted Newick trees from 'main' (as standard input).
 *
 * Example:
 *        ./main --scrm --input trees100.txt --verbose --newick --nrow 1000 | ./newick2clustersize
 */
#include "NewickTree.h"
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
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
    
int main(int argc, char ** argv)
{
    if (argc != 1)
    {
        cerr << "usage: ./main [options] | " << argv[0] << endl;
        return 1;
    }

    // Init inferred tree input from <stdin>
    NewickTree predt("-");
    assert(predt.nleaves() == nleaves);
    unsigned predt_base = 0;
    unsigned n = 0;
    while (predt.good())
    {        
        predt_base += predt.validForBases();
        vector<unsigned> cd(nleaves*nleaves, 0);
        predt.subtreeDistance(cd);

        // Output <base> <leaf id> <leaf id> <dist. in SCRM> <dist. in stdin>
        for (unsigned i = 0; i < nleaves-1; ++i)
            for (unsigned j = i+1; j < nleaves; ++j)
                cout << predt_base << '\t' << i+1 << '\t' << j+1 << '\t' << cd[i + nleaves * j] << endl;

        // Find next inferred tree
        predt.next();        
        n++;
    }
    cerr << n << " sites, " << predt_base << " predt_bases and " << scrm_base << " scrm_bases done." << endl;
    return 0;
}
