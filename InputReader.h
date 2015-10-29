#ifndef _InputReader_h_
#define _InputReader_h_
#include "default.h"
#include <string>
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
    
    // Get next input row (or return False if EOF)
    virtual bool next(InputColumn &) = 0;
    
    // Misc helper functions
    virtual bool good() const
    { return (fp != 0 && fp->good()); }
    InputColumn const & col(std::size_t i) const
    { return cols[i]; }
    std::size_t size() const
    { return cols.size(); }

    
    virtual ~InputReader()
    { delete fp; fp = 0; }
    
protected:
    explicit InputReader(std::string);
    std::istream *fp;
    std::vector<InputColumn> cols;
    
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
    bool next(InputColumn &);
    virtual ~PlainTextInputReader()
    { }
    
    explicit PlainTextInputReader(std::string);
    
private:
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
    bool next(InputColumn &);
    virtual ~SimpleVCFInputReader()
    { }
    
    explicit SimpleVCFInputReader(std::string);
    
private:
    static char upcasedna[256];
    // No copy constructor or assignment
    SimpleVCFInputReader(SimpleVCFInputReader const&);
    SimpleVCFInputReader& operator = (SimpleVCFInputReader const&);
};

/**
 * Simple SCRM input reader
 */
class SimpleSCRMInputReader : public InputReader
{
public:
    bool next(InputColumn &);
    virtual ~SimpleSCRMInputReader()
    { }
    
    SimpleSCRMInputReader(std::string, unsigned);
    virtual bool good() const
    { return good_; }
    
private:    
    // No copy constructor or assignment
    SimpleSCRMInputReader(SimpleSCRMInputReader const&);
    SimpleSCRMInputReader& operator = (SimpleSCRMInputReader const&);
    std::size_t pos;
    bool good_;
};
#endif
