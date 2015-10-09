arg
====

To get started:

0) git clone https://github.com/nvalimak/arg.git

1) Compile the software by issuing the command `make'.

2) Issue the test by './main --input test.input --plaintext --verbose --debug --dot test'

You will need the Graphwiz DOT software package to compile the *.dot files with:

3) 'make dot'

TODO
----

* Optimize the second reduce() call in TreeController::process().
* Choose target by subtree size (number of leaves)
* VCF format
* 
