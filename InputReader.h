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
    enum input_format_t { input_unset, input_plaintext, input_vcf };

    // Constructor
    static InputReader* build(input_format_t, std::string);
    
    // Get next input row (or return False if EOF)
    virtual bool next(InputColumn &) = 0;
    
    // Misc helper functions
    bool good() const
    { return (fp != 0 && fp->good()); }
    
    virtual ~InputReader()
    { delete fp; fp = 0; }
    
protected:
    explicit InputReader(std::string);
    std::istream *fp;
    
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
 * VCF input reader (TODO)
 */

#endif
