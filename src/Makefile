IDIR=../include
CC=g++
CFLAGS=-I$(IDIR) -O3

ODIR=../obj
EXEC=../main


_DEPS = binarytree.h driftFunctions.h fastmarch.h fastMarchEJM.h fastmarchOLIM.h  linalg.h updaterBase.h updaterLinear.h updaterHermite.h utility.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=main.o binarytree.o driftFunctions.o fastmarch.o fastmarchEJM.o fastmarchOLIM.o updaterHermite.o updaterLinear.o utility.o
OBJ= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f $(OBJ) $(EXEC)