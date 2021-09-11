IDIR=../include
CC=g++
CFLAGS=-I$(IDIR) -O3

_DEPS = binarytree.h fastmarch.h fastmarchASR.h fastmarchOLIM.h \
       structs.h updater.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

ODIR=../obj
EXEC=../main

_OBJ=binarytree.o fastmarch.o fastmarchASR.o fastmarchOLIM.o main.o updater.o
OBJ= $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	rm -f $(OBJ) $(EXEC)
