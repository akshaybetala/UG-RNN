OBJECTST = ./obj/General.o ./obj/Graph.o ./obj/Molecule.o ./obj/Layer.o ./obj/NNt.o ./obj/NNperc.o ./obj/NN.o ./obj/NNr.o ./obj/main.o ./obj/UGRNN.o ./obj/myOptions.o ./obj/nodeIO.o
HEADERST = ./sources/General.h ./sources/Node.h ./sources/Graph.h ./sources/DataSet.h ./sources/Molecule.h ./sources/General.h ./sources/Layer.h ./sources/NNt.h ./sources/NNperc.h ./sources/NN.h ./sources/NNr.h ./sources/UGRNN.h ./sources/myOptions.h ./sources/nodeIO.h

CC	= g++
#CFLAGS  = -O3 -Q -v -da -g
CFLAGS = -O3

ugrnnTrain :$(OBJECTST)
	$(CC) $(CFLAGS) -o $@ $(OBJECTST) $(LLIBS)

./obj/main.o : ./sources/main.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/main.cxx
	mv main.o ./obj
./obj/General.o : ./sources/General.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/General.cxx
	mv General.o ./obj
./obj/Layer.o : ./sources/Layer.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/Layer.cxx
	mv Layer.o ./obj
./obj/NNt.o : ./sources/NNt.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/NNt.cxx
	mv NNt.o ./obj 
./obj/NNperc.o : ./sources/NNperc.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/NNperc.cxx
	mv NNperc.o ./obj
./obj/NNr.o : ./sources/NNr.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/NNr.cxx
	mv NNr.o ./obj
./obj/NN.o : ./sources/NN.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/NN.cxx
	mv NN.o ./obj
./obj/Molecule.o: ./sources/Molecule.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/Molecule.cxx
	mv Molecule.o ./obj
./obj/Graph.o : ./sources/Graph.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/Graph.cxx
	mv Graph.o ./obj
./obj/UGRNN.o : ./sources/UGRNN.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/UGRNN.cxx
	mv UGRNN.o ./obj
./obj/myOptions.o : ./sources/myOptions.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/myOptions.cxx
	mv myOptions.o ./obj
./obj/nodeIO.o : ./sources/nodeIO.cxx $(HEADERST)
	$(CC) -c $(CFLAGS) ./sources/nodeIO.cxx
	mv nodeIO.o ./obj
clean:
	rm ./obj/*.o
