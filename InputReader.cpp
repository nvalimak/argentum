#include "InputReader.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib> // abort()
#include <cassert>
using namespace std;


char SimpleVCFInputReader::upcasedna[] = {0};


InputReader * InputReader::build(input_format_t input, string file, unsigned nrows)
{
    switch (input)
    {
    case input_vcf:
        return new SimpleVCFInputReader(file, nrows);
        break;
    case input_scrm:
        return new SimpleSCRMInputReader(file, nrows);
        break;
    case input_plaintext:
        return new PlainTextInputReader(file, nrows);
        break;
    default:
        cerr << "Error: invalid file type at InputReader::build()" << endl; 
        abort();
    }
}

InputReader::InputReader(string file)
    : fp(0), cols(), good_(true)
{
    cols.reserve(1024);
    // Open file handle, "-" uses standard input
    if (file == "-")
        fp = &std::cin;
    else
        fp = new ifstream(file.c_str());

    if (fp == 0 || !fp->good())
        good_ = false;
}

/**
 * Plain text input
 */
PlainTextInputReader::PlainTextInputReader(string file, unsigned nrows)
    : InputReader(file)
{
    if (!good_)
        return;

    InputColumn ic;
    while (nrows && next(ic))
        --nrows;
}

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
    cols.push_back(ic);
    return true;
}

/**
 * Simplistic VCF input
 */
SimpleVCFInputReader::SimpleVCFInputReader(string file, unsigned nrows)
    : InputReader(file), positions()
{
    if (!good_)
        return;
    for (unsigned i = 0; i < 256; ++i)
        SimpleVCFInputReader::upcasedna[i] = i;
    SimpleVCFInputReader::upcasedna[(unsigned)'a'] = 'A'; // Oh my.. FIXME
    SimpleVCFInputReader::upcasedna[(unsigned)'c'] = 'C';
    SimpleVCFInputReader::upcasedna[(unsigned)'g'] = 'G';
    SimpleVCFInputReader::upcasedna[(unsigned)'t'] = 'T';

    InputColumn ic;
    while (nrows && next(ic))
        --nrows;
}

bool SimpleVCFInputReader::next(InputColumn &ic)
{
    string row;
    if (!std::getline(*fp, row).good())
        return false;

    // Assumes all '#' rows are at the beginning of the file
    while (row[0] == '#')
        if (!std::getline(*fp, row).good())
            return false;

    std::size_t AApos = 0, pos = 0;
    int gtaa = -1;
    string AA = "A";
    string ref = "B";
    string alt = "C";
    do
    {
        AApos = row.find("AA=");
        while (AApos == std::string::npos || row[AApos+3]=='|' || row[AApos+3]=='.')
        {
            if (!std::getline(*fp, row).good())
                return false;
            AApos = row.find("AA=");
        }
        
        // Assert: row now corresponds to an ancerstral allele        
        for (size_t i = 0; i<= AApos+10; ++i)
            row[i] = SimpleVCFInputReader::upcasedna[(unsigned)row[i]];
        size_t len = 1;
        while (row[AApos+3+len] != '|' && row[AApos+3+len] != ';')
            ++len;
        AA = row.substr(AApos+3,len);
        
        // Extract REF and ALT
        pos = 0;
        unsigned col = 0;
        while (col < 3)
            if (row[pos++] == '\t')
                ++col;
        ref = row.substr(pos, row.find_first_of("\t", pos + 1) - pos);
        pos = row.find_first_of("\t", pos + 1);
        alt = row.substr(pos+1, row.find_first_of("\t,", pos + 1) - pos - 1); // No support for multiallelic sites

        if (row[row.find_first_of("\t,", pos + 1)] == ',')
        {
            if (!std::getline(*fp, row).good())
                return false;        
            continue; // Skip multiallelic sites
        } 

        if (ref == AA)
            gtaa = 1;
        else if (alt == AA)
            gtaa = 0;

        if (gtaa == -1)
            if (!std::getline(*fp, row).good())
                return false;        
    } while(gtaa == -1);

    
    // Parse the genotypes
    pos = row.find_first_of("\t", AApos + 1);
    pos = row.find_first_of("\t", pos + 1);
    size_t n = (row.size() - pos)/2; // Assuming 4 bytes per individual
    
    if (row.empty() || (ic.size() && n != ic.size()))
    {   
        std::cerr << "unable to read VCF input: expecting equal length rows of input" << std::endl;
        abort();
    }
    
    // first row of input?
    if (ic.empty())
        ic.resize(n);
    
    for (unsigned i = 0; i < n; ++i, pos+=2)
    {
        assert (row[pos] == '\t' || row[pos] == '|'); // FIXME Remove
        ic[i] = (row[pos+1] == (char)(48+gtaa));
    }
    cols.push_back(ic);

    unsigned vcf_pos = 0;
    {
        pos = row.find_first_of("\t") + 1;
        istringstream iss(row.substr(pos));
        iss >> vcf_pos;
        assert (pos > 0);
    }
    positions.push_back(vcf_pos);
    return true;
}

/**
 * Simplistic SCRM input
 */
SimpleSCRMInputReader::SimpleSCRMInputReader(string file, unsigned nrows)
    : InputReader(file), positions()
{
    if (!good_)
        return;
    vector<string> rows;
    string row;
    unsigned i = 0;
    double totaln = 0;
    while (row.size() < 10 || row.substr(0,9) != "positions")
    {
        if (!std::getline(*fp, row).good())
        {
            good_ = false;
            return;
        }
        if (row[0] == '[')
        {
            istringstream iss(row.substr(1));
            unsigned z = 0;
            iss>>z;
            totaln += z;
        }
    }
    positions.reserve(1024);
    assert (row.substr(0,9) == "positions");
    istringstream iss(row.substr(11));
    double d = 0;
    unsigned prev_pos = 0;
    while ((iss >> d))
    {
        assert (d>0.0);
        unsigned pos = round(d*totaln);
        positions.push_back(pos - prev_pos);
        prev_pos = pos;
    }
    while (std::getline(*fp, row).good())
    {
        rows.push_back(row);
        if (rows[0].size() != rows.back().size())
        {
            cerr << "error at SimpleSCRMInputReader: Check input format (tree input is mandatory)" << endl;
            abort();
        }
    }

    cols.resize(rows[0].size() < nrows ? rows[0].size() : nrows ); // transpose
    for (i = 0; i < rows[0].size() && i < nrows; ++i)
        for (unsigned j = 0; j < rows.size(); ++j)
            cols[i].push_back((rows[j][i] == '1'));
}

