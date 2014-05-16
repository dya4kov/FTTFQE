IDIR=./header
LDIR=/usr/local/lib
CC=g++
CFLAGS=-I$(IDIR)
ODIR=./obj
SDIR=./src
LIBS=-lgsl -lgslcblas -lm

_DEPS=Interface.h \
      Model.h \
      ODESolver.h \
      ODESystem.h \
      ODESystemList.h \
      Solution.h \
      SpecialFunctions.h \
      ThermodynamicFunctions.h \
      ThomasFermi.h \
      ThomasFermiCorrection.h \
      ThomasFermiPotential.h \
      ThomasFermiPotentialCorrection.h \
      Units.h

DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ=Interface.o \
     main.o \
     Model.o \
     ODESystemList.o \
     Solution.o \
     SpecialFunctions.o \
     ThermodynamicFunctions.o \
     ThomasFermi.o \
     ThomasFermiCorrection.o \
     ThomasFermiPotential.o \
     ThomasFermiPotentialCorrection.o

OBJ=$(patsubst %,$(ODIR)/%,$(_OBJ))

all: makeobjdir $(OBJ) fttfqe

makeobjdir:
	mkdir $(ODIR)

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS) 
	$(CC) -c -o $@ $< $(CFLAGS)

fttfqe: $(OBJ) 
	$(CC) -o $@ $^ $(LIBS)

.PHONY: clean

clean:
	rm -rf $(ODIR) *~ $(IDIR)/*~ $(SDIR)/*~ fttfqe

