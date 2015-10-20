argentum
====

To get started:

0) git clone https://github.com/nvalimak/argentum.git

1) Compile the software by issuing the command `make'.

2) Issue the test by './main --input test.input --plaintext --verbose --debug --dot test'

You will need the /Graphwiz DOT/ software package to compile the *.dot files with:

3) 'make dot'

Modify the 'Makefile' to turn on compiler optimizations, if needed.

Input file format
---

VCF files can be processed as, for example:

    './main --vcf --input <(gunzip -c input.vcf.gz) --verbose'

SCRM files (without the SCRM tree) can be processed as, for example:

    './main --scrm --input input.scrm --verbose'

TODO
----

* Optimize the tree representation
* Make the compiler optimizations on by default