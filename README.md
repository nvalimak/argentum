argentum
====

To get started:

0) `git clone https://github.com/nvalimak/argentum.git`

1) Compile the software by issuing the command `make`.

2) Issue the test by `./main --input test.input --plaintext --verbose --debug --dot test`

You will need the /Graphwiz DOT/ software package to compile the *.dot files with:

3) `make dot`

Modify the `Makefile` to turn on compiler optimizations, if needed.

Input file format
---

VCF files can be processed as, for example:

    ./main --vcf --input <(gunzip -c input.vcf.gz) --verbose

SCRM files (without the SCRM tree) can be processed as, for example:

    ./main --scrm --input input.scrm --verbose

Time estimation
---
Use `--enumerate` for `./main` to extract an ARG as a single graph with node ranges.
`./enumerate-example` reads the standard output in the enumerate format. Please see guide\_time\_estimate.pdf for instructions.

Graph slice
---
Set `[dis_out]` parameter of `./enumerate-example` equal to 3.
There is no control of slice parameters from the command line (TODO).
The output consists of connected components in the wpairs format `[node 1]\t[node 2]\t[edge weight]\n` (native input format for Fast Community http://www.cs.unm.edu/~aaron/research/fastmodularity.htm) followed by output of node description in the format `[ID in ARG]\t[ID in slice]\t[timestamp]\n`.

TODO
----

* Optimize the tree representation
* Make the compiler optimizations on by default
