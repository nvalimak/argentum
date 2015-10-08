#include "InputReader.h"
#include <fstream>
#include <cstdlib> // abort()
#include <cassert>
using namespace std;

InputReader * InputReader::build(input_format_t input, string file)
{
    switch (input)
    {
    case input_plaintext:
        return new PlainTextInputReader(file);
        break;
    case input_vcf:
        // TODO
        //return VCFInputReader(file);
        //break;
    default:
        cerr << "Error: invalid file type at InputReader::build()" << endl; 
        abort();
    }
}

InputReader::InputReader(string file)
    : fp(0)
{
    // Open file handle, "-" uses standard input
    if (file == "-")
        fp = &std::cin;
    else
        fp = new ifstream(file.c_str());
}

PlainTextInputReader::PlainTextInputReader(string file)
    : InputReader(file)
{ }

bool PlainTextInputReader::next(InputColumn &ic)
{
    string row;
    if (!std::getline(*fp, row).good())
        return false;
    if (row.empty() || (ic.size() && row.size() != ic.size())) 
    {   
        std::cerr << "unable to read plaintext input: expecting equal length rows of input" << std::endl;
        return false;
    }
    
    // first row of input?
    if (ic.empty())
        ic.resize(row.size());
    
    for (unsigned i = 0; i < row.size(); ++i)
        ic[i] = (row[i] == '1');
    return true;
}
