CC=g++
CPPFLAGS=-std=c++11 -Wall -ansi -pedantic -g -O2 -DNDEBUG
DISBALEDOPTIMIZATIONFLAGS = -O2 -DNDEBUG
OBJ = Configuration.o InputReader.o Tree.o TreeController.o TreeControllerSimple.o

all: main main-simple newickcomp rand-matrix scrm2bin

main: main.o $(OBJ)
	$(CC) $(CPPFLAGS) -o main main.o $(OBJ)

main-simple: main-simple.o $(OBJ)
	$(CC) $(CPPFLAGS) -o main-simple main-simple.o $(OBJ)

newickcomp: newickcomp.o
	$(CC) $(CPPFLAGS) -o newickcomp newickcomp.o

rand-matrix: rand-matrix.o
	$(CC) $(CPPFLAGS) -o rand-matrix rand-matrix.o

scrm2bin: scrm2bin.o
	$(CC) $(CPPFLAGS) -o scrm2bin scrm2bin.o

test: main
	./main --input test.input --plaintext --verbose --debug 1000 --dot test --rewind
	./main --input test.input2 --plaintext --verbose --debug 1000 --dot test --rewind
	./main --input test.input3 --plaintext --verbose --debug 1000 --dot test --rewind
	./main --input test.input4 --plaintext --verbose --debug 1000 --dot test --rewind

dot: main
	for i in *.dot; do dot -Tpng $$i > $${i%.dot}.png; done

clean:
	rm -f main main-simple newickcomp rand-matrix scrm2bin *.o *~ *.dot *.png

depend:
	g++ -MM -std=c++11 *.cpp > dependencies.mk

include dependencies.mk
