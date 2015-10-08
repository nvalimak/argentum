CC = g++
CPPFLAGS = -std=c++0x -Wall -ansi -pedantic -g
DISBALEDOPTIMIZATIONFLAGS = -O2 -DNDEBUG
OBJ = InputReader.o Tree.o TreeController.o

all: main test dot

main: main.o $(OBJ)
	$(CC) $(CPPFLAGS) -o main main.o $(OBJ)

test: main
	./main --input test.input --plaintext --verbose --debug --dot test

dot: main
	for i in *.dot; do dot -Tpng $$i > $${i%.dot}.png; done

clean:
	rm -f main *.o *~ *.dot *.png

depend:
	g++ -MM -std=c++0x *.cpp > dependencies.mk

include dependencies.mk
