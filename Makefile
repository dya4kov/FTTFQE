IDIR=./hdr
LDIR=./lib
ODIR=./obj
BDIR=./bin
SDIR=./src
LOGDIR=./log
CC=g++
CFLAGS=-I$(IDIR)
LIBS=-lm

TF_POT_OBJ=$(ODIR)/FTTFpotential.o
TFQE_POT_OBJ=$(ODIR)/FTTFQEpotential.o
TF_MOD_OBJ=$(ODIR)/FTTFmodel.o
TFQE_MOD_OBJ=$(ODIR)/FTTFQEmodel.o
Y_OBJ=$(ODIR)/Yfunction.o
PRINT_OBJ=$(ODIR)/printer.o
TIMER_OBJ=$(ODIR)/Timer.o
PER_TABL_OBJ=$(ODIR)/periodicTableTest.o

objdir:
	mkdir -p $(ODIR) 

bindir:
	mkdir -p $(BDIR) 

logdir:
	mkdir -p $(LOGDIR)

$(ODIR)/%.o: $(SDIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< -I$(IDIR)

FTTFpot: objdir bindir logdir $(TF_POT_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TF_POT_OBJ) $(LIBS)

FTTFpot: DEPS=$(TF_POT_DEPS)

FTTFQEpot: objdir bindir logdir $(TFQE_POT_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TFQE_POT_OBJ) $(LIBS)

FTTFQEpot: DEPS=$(TFQE_POT_DEPS)	

FTTFmodel: objdir bindir logdir $(TF_MOD_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TF_MOD_OBJ) $(LIBS)

FTTFmodel: DEPS=$(TF_MOD_DEPS)

FTTFQEmodel: objdir bindir logdir $(TFQE_MOD_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TFQE_MOD_OBJ) $(LIBS)

FTTFQEmodel: DEPS=$(TFQE_MOD_DEPS)

Ytest: objdir bindir $(Y_OBJ) 
	$(CC) -o $(BDIR)/$@ $(Y_OBJ) $(LIBS)

Ytest: DEPS=$(Y_DEPS)

printTest: objdir bindir $(PRINT_OBJ) 
	$(CC) -o $(BDIR)/$@ $(PRINT_OBJ) $(LIBS)

timer: objdir bindir $(TIMER_OBJ) 
	$(CC) -o $(BDIR)/$@ $(TIMER_OBJ) $(LIBS)

periodicTableTest: objdir bindir $(PER_TABL_OBJ) 
	$(CC) -o $(BDIR)/$@ $(PER_TABL_OBJ) $(LIBS)

.PHONY: clean

clean:
	rm -rf $(ODIR) *~ $(IDIR)/*~ $(SDIR)/*~ $(BDIR) $(LOGDIR)
