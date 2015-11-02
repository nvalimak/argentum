Configuration.o: Configuration.cpp Configuration.h default.h \
  InputReader.h
InputReader.o: InputReader.cpp InputReader.h default.h
Tree.o: Tree.cpp Tree.h default.h
TreeController.o: TreeController.cpp TreeController.h default.h Tree.h
TreeControllerSimple.o: TreeControllerSimple.cpp TreeControllerSimple.h \
  default.h Tree.h TreeController.h
main-simple.o: main-simple.cpp InputReader.h default.h Tree.h \
  TreeControllerSimple.h
main.o: main.cpp Configuration.h default.h InputReader.h Tree.h \
  TreeController.h
rand-matrix.o: rand-matrix.cpp
scrm2bin.o: scrm2bin.cpp
