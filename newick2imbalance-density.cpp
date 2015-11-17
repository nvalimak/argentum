/**
 * Helper script to parse Newick trees
 *
 * Usage:
 *         ./main [options] | ./newick2imbalance-density <SCRM>
 *
 * Input: Predicted Newick trees from 'main' (as standard input).
 *        <SCRM> SCRM output file name.
 *
 * Output: tab-separated values (as standard output).
 *
 * Example:
 *        ./main --scrm --input trees100.txt --verbose --newick --nrow 1000 | ./newick2imbalance-density trees100.txt pops_map.txt > R.tsv
 */
#include "NewickTree.h"
#include <iostream>
#include <sstream>
#include <fstream>
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
        cerr << "usage: ./main [options] | " << argv[0] << " <SCRM-input-file> <population file>" << endl;
        return 1;
    }

    // Init population map
    unsigned npop = 0;
    unsigned first_pop_size = 0;
    ifstream ifs(argv[2]);
    vector<unsigned> popv;
    if (!ifs.good())
    {
        cerr << "error: unable to read file " << argv[2] << endl;
        return 1;
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
        assert (lid-1 == popv.size()); // Must be ordered by leaf id
        popv.push_back(p-1);
        npop = max(npop, p);
        if (p == 1)
            first_pop_size ++;
    }
    assert (npop == 2); // default assumption for now
    assert (first_pop_size > 0);
    
    // Initialize SCRM input tree
    NewickTree scrm(argv[1]);
    unsigned nleaves = scrm.nleaves();
    unsigned scrm_base = scrm.validForBases();

    assert (nleaves == popv.size()); // Population map must have the same size as the tree.
    
    // Init inferred tree input from <stdin>
    NewickTree predt("-");
    assert(predt.nleaves() == nleaves);
    unsigned predt_base = 0;

    // Init the imbalance density vectors
    vector<double> scrm_density(nleaves, 0);
    vector<double> predt_density(nleaves, 0);
    vector<unsigned> scrm_count(nleaves, 0);
    vector<unsigned> predt_count(nleaves, 0);

    unsigned n = 0;
    while (scrm.good() && predt.good())
    {        
        predt_base += predt.validForBases();
        predt.updateDensity(popv, first_pop_size, predt_density, predt_count);

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
        scrm.updateDensity(popv, first_pop_size, scrm_density, scrm_count);
                
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

    cout << "\tscrm_avg\tpred_avg" << endl;
    for (unsigned i = 0; i < nleaves; ++i)
    {
        scrm_density[i] /= scrm_count[i];
        predt_density[i] /= predt_count[i];
        cout << i << '\t' << scrm_density[i] << '\t' << predt_density[i] << endl; 
    }    
    return 0;
}
