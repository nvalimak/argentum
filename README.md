arg
====

To get started:

0) git clone https://github.com/nvalimak/arg.git

1) Compile the software by issuing the command `make'.

2) Issue the test by './main --input test.input --plaintext --verbose --debug --dot test'

You will need the Graphwiz DOT software package to compile the *.dot files with:

3) 'make dot'

VCF files can be processed as, for example:

4) './main --vcf --input <(gunzip -c input.vcf.gz) --verbose'

Modify the 'Makefile' to turn on compiler optimizations, if needed.

TODO
----

* Choose target by subtree size (number of leaves)
* Choose targets by 0-1-genotype count values
