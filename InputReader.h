#ifndef _InputReader_h_
#define _InputReader_h_
#include "default.h"
#include <string>
#include <vector>
#include <cassert>
#include <iostream>

/**
 * Base class for input
 *
 * Use InputReader::build() to construct an instance.
 */
class InputReader
{
public:
    enum input_format_t { input_unset, input_vcf, input_scrm, input_plaintext };

    // Constructor
    static InputReader* build(input_format_t, std::string, unsigned);
    
    // Misc helper functions
    virtual bool good() const
    { return good_; }
    InputColumn const & col(std::size_t i) const
    { return cols[i]; }
    InputColumn const & front() const
    { return cols.front(); }
    InputColumn const & back() const
    { return cols.back(); }
    std::size_t size() const
    { return cols.size(); }
    virtual unsigned position(size_t i)
    { return 0; } // SCRM+VCF input only
    
    virtual ~InputReader()
    { if (fp != &std::cin) delete fp; fp = 0; }
    
protected:
    explicit InputReader(std::string);
    std::istream *fp;
    std::vector<InputColumn> cols;
    bool good_;
    
private:
    InputReader();
    // No copy constructor or assignment
    InputReader(InputReader const&);
    InputReader& operator = (InputReader const&);
};


/**
 * Simple input reader
 */
class PlainTextInputReader : public InputReader
{
public:
    PlainTextInputReader(std::string, unsigned);
    virtual ~PlainTextInputReader()
    { }
private:
    bool next(InputColumn &);
    // No copy constructor or assignment
    PlainTextInputReader(PlainTextInputReader const&);
    PlainTextInputReader& operator = (PlainTextInputReader const&);
};

/**
 * Simple VCF input reader
 */
class SimpleVCFInputReader : public InputReader
{
public:
    SimpleVCFInputReader(std::string, unsigned);
    virtual unsigned position(size_t i)
    { assert(i<positions.size()); return positions[i]; }

    virtual ~SimpleVCFInputReader()
    { }
private:
    bool next(InputColumn &);
    static char upcasedna[256];
    // No copy constructor or assignment
    SimpleVCFInputReader(SimpleVCFInputReader const&);
    SimpleVCFInputReader& operator = (SimpleVCFInputReader const&);
    std::vector<unsigned> positions;
};

/**
 * Simple SCRM input reader
 */
class SimpleSCRMInputReader : public InputReader
{
public:
    SimpleSCRMInputReader(std::string, unsigned);
    virtual unsigned position(size_t i)
    { assert(i<positions.size()); return positions[i]; }
    virtual ~SimpleSCRMInputReader()
    { }
private:    
    // No copy constructor or assignment
    SimpleSCRMInputReader(SimpleSCRMInputReader const&);
    SimpleSCRMInputReader& operator = (SimpleSCRMInputReader const&);
    std::vector<unsigned> positions;
};
#endif
