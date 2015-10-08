CC = g++
CPPFLAGS = -std=c++0x -Wall -ansi -pedantic -g
DISBALEDOPTIMIZATIONFLAGS = -O2 -DNDEBUG
OBJ = InputReader.o Tree.o TreeController.o

main: main.o $(OBJ)
	$(CC) $(CPPFLAGS) -o main main.o $(OBJ)

dot:
	for i in *.dot; do dot -Tpng $$i > $${i%.dot}.png; done


clean:
	rm -f main *.o *~

depend:
	g++ -MM -std=c++0x *.cpp > dependencies.mk

include dependencies.mk
