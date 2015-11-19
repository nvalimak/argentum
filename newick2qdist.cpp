/**
 * Helper script to parse Newick trees
 *
 * Usage:
 *         ./main [options] | ./newick2qdist <SCRM>
 *
 * Input: Predicted Newick trees from 'main' (as standard input).
 *        <SCRM> SCRM output file name.
 *
 * Example:
 *        ./main --scrm --input trees100.txt --verbose --newick --nrow 1000 | ./newick2qdist trees100.txt 
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
    if (argc != 4)
    {
        cerr << "usage: ./main [options] | " << argv[0] << " <SCRM-input-file> <min-leaf1> <max-leaf>" << endl;
        return 1;
    }

    // Initialize SCRM input tree
    NewickTree scrm(argv[1]);
    unsigned nleaves = scrm.nleaves();
    unsigned scrm_base = scrm.validForBases();

    // Init inferred tree input from <stdin>
    NewickTree predt("-");
    assert(predt.nleaves() == nleaves);
    unsigned predt_base = 0;
    unsigned n = 0;
    while (scrm.good() && predt.good())
    {        
        predt_base += predt.validForBases();
        while (scrm_base < predt_base)
        {
            // Find the corresponding SCRM tree for the site at 'predt_base'
            scrm.next();
            if (!scrm.good())
                break;
            scrm_base += scrm.validForBases();
            assert (scrm.nleaves() == nleaves);
        }        
        if (!scrm.good())
            break;

        {
            string ptree(predt.originalRepresentation());
            string stree(scrm.originalRepresentation());        
            cout << ptree.substr(ptree.find_first_of("(")) << '\t'
                 << stree.substr(stree.find_first_of("(")) << endl;
        }
        
        do
        {
            // Find next inferred tree
            predt.next();
            if (!predt.good())
                break;
            assert (predt.nleaves() == nleaves);
        } while (predt.good() && predt.validForBases() == 0); // Skip if valid only for 0 bases
        n++;
    }
    cerr << n << " sites, " << predt_base << " predt_bases and " << scrm_base << " scrm_bases done." << endl;
    return 0;
}
