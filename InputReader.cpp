#include "InputReader.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cstdlib> // abort()
#include <cassert>
using namespace std;


char SimpleVCFInputReader::upcasedna[] = {0};


InputReader * InputReader::build(input_format_t input, string file)
{
    switch (input)
    {
    case input_plaintext:
        return new PlainTextInputReader(file);
        break;
    case input_vcf:
        return new SimpleVCFInputReader(file);
        break;
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

/**
 * Plain text input
 */
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

/**
 * Simplistic VCF input
 */
SimpleVCFInputReader::SimpleVCFInputReader(string file)
    : InputReader(file)
{
    for (unsigned i = 0; i < 256; ++i)
        SimpleVCFInputReader::upcasedna[i] = i;
    SimpleVCFInputReader::upcasedna[(unsigned)'a'] = 'A'; // Oh my.. horrible. FIXME
    SimpleVCFInputReader::upcasedna[(unsigned)'c'] = 'C';
    SimpleVCFInputReader::upcasedna[(unsigned)'g'] = 'G';
    SimpleVCFInputReader::upcasedna[(unsigned)'t'] = 'T';
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
        //cerr << "next row = " << row.substr(0,200)<<endl;
        //cerr << "found AA = " << AA << "'" << endl;
        
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
            continue;
        } 
            
        //cerr << "found ref = " << ref << "' and alt = " << alt << "'" << endl;
        
/*        int gtaa = -1;
        if (ref == AA)
            gtaa = 0;       // Ancestral allele is the REF
        else if (alt == AA)
            gtaa = 1;       // Ancestral allele is the first ALT
        else
        {
            // Ancestral allele is the second, or third, or ... ALT
            pos = row.find_first_of("\t,", pos + 1);
            assert(row[pos] != string::npos);
            int tmp = 2;
            while (row.substr(pos+1, row.find_first_of("\t,", pos + 1) - pos - 1) != AA)
            {
                ++tmp;
                pos = row.find_first_of("\t,", pos + 1);
            }
            if (row.substr(pos+1, row.find_first_of("\t,", pos + 1) - pos - 1) == AA)
                gtaa = tmp;
            assert(row.substr(pos+1, row.find_first_of("\t,", pos + 1) - pos -1) == AA);
            }
        
    if (gtaa == -1)
    {
        cerr << "error at SimpleVCFInputReader::SimpleVCFInputReader(): problem while deciding the ancestral allele with REF = "
             << ref << ", alt " << alt << ", aa = " << AA << endl;
        cerr << "error at SimpleVCFInputReader::SimpleVCFInputReader(): at row: " << row.substr(0,100) << endl;
        abort();
        }*/

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
    return true;
}
