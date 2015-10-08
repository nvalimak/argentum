InputReader.o: InputReader.cpp InputReader.h default.h
Tree.o: Tree.cpp Tree.h default.h
TreeController.o: TreeController.cpp TreeController.h default.h Tree.h
main.o: main.cpp InputReader.h default.h Tree.h TreeController.h
