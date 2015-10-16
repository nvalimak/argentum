#ifndef _default_h_
#define _default_h_
#include <vector>

// Default storage format
typedef unsigned char InputLabel;
typedef unsigned NodeId;
typedef unsigned LeafId;
typedef std::vector<InputLabel> InputColumn;
typedef double TreeDepth;

// Initial memory reservation
#define HISTORY_INIT_SIZE 1024*1024
#define BUFFER_INIT_SIZE 20*1024

#endif
