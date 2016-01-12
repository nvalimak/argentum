Configuration.o: Configuration.cpp Configuration.h default.h \
  InputReader.h TreeDistance.h Tree.h NewickTree.h
InputReader.o: InputReader.cpp InputReader.h default.h
NewickTree.o: NewickTree.cpp NewickTree.h default.h
Tree.o: Tree.cpp Tree.h default.h TreeDistance.h NewickTree.h \
  TreeEnumerator.h
TreeController.o: TreeController.cpp TreeController.h default.h Tree.h \
  TreeDistance.h NewickTree.h TreeEnumerator.h
TreeControllerSimple.o: TreeControllerSimple.cpp TreeControllerSimple.h \
  default.h Tree.h TreeDistance.h NewickTree.h TreeEnumerator.h \
  TreeController.h
TreeDistance.o: TreeDistance.cpp TreeDistance.h default.h Tree.h \
  NewickTree.h
enumerate-example.o: enumerate-example.cpp
main-simple.o: main-simple.cpp Configuration.h default.h InputReader.h \
  TreeDistance.h Tree.h NewickTree.h TreeControllerSimple.h \
  TreeEnumerator.h
main.o: main.cpp Configuration.h default.h InputReader.h TreeDistance.h \
  Tree.h NewickTree.h TreeController.h TreeEnumerator.h
newick2clustersize.o: newick2clustersize.cpp NewickTree.h default.h
newick2dot.o: newick2dot.cpp NewickTree.h default.h
newick2imbalance-density.o: newick2imbalance-density.cpp NewickTree.h \
  default.h
newick2qdist.o: newick2qdist.cpp NewickTree.h default.h
newickcomp.o: newickcomp.cpp NewickTree.h default.h
rand-matrix.o: rand-matrix.cpp
scrm2bin.o: scrm2bin.cpp
