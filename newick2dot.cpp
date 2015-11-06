/**
 * Helper script to parse Newick trees
 *
 * Usage:
 *         cat newicktree.txt | ./newick2dot dot-prefix
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
    if (argc != 3)
    {
        cerr << "usage: cat newicktree.txt | " << argv[0] << " <output prefix> <max rows>" << endl;
        return 1;
    }

    string fnprefix = string(argv[1]);
    unsigned maxr = atoi_min(argv[2], 1);
    
    // Init tree input from <stdin>
    NewickTree predt("-");
    unsigned n = 0;
    unsigned base = 0;
    while (predt.good() && n < maxr)
    {        
        base += predt.validForBases();
        predt.outputDOT(fnprefix, n, base);

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
