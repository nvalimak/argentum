/**
 * Transpose SCRM output to text matrix of '0's and '1's
 *
 * Input: standard input
 * Output: standard output
 */
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
    if (argc != 1)
    {
        cerr << "usage: cat SCRM.txt | " << argv[0] << " > output.txt" << endl;
        return 1;
    }

    std::istream *fp = &cin;

    vector<string> rows;
    string row;
    unsigned i = 0;
    while (i++ < 6)
        if (!std::getline(*fp, row).good())
        {
            cerr << "error: unable to read input" << endl;
            return 1;
        }

    while (std::getline(*fp, row).good())
    {
        rows.push_back(row);
        if (rows[0].size() != rows.back().size())
        {
            cerr << "error: Check input format (SCRM tree output is not supported)" << endl;
            abort();
        }
    }

    cerr << "read " << rows.size() << " rows of length " << rows[0].size() << endl;

    cout << "transposeddata" << endl;
    for (i = 0; i < rows[0].size(); ++i)
    {
        for (unsigned j = 0; j < rows.size(); ++j)
            cout << (unsigned)(rows[j][i] == '1');
        cout << endl;
    }

    return 0;
}
