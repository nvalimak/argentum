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

TODO
----

* Optimize the tree representation
* Make the compiler optimizations on by default
