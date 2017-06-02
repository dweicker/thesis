CC=gcc
CFLAGS= -I /Users/davidweicker/Desktop/Thesis/p4est-1.1/local/include/       
LDFLAGS= -L /Users/davidweicker/Desktop/Thesis/p4est-1.1/local/lib/ -lp4est -lsc -lz
GLLIBS = 
SRCDIR = 
HEADDIR = 
LIBDIR = build/
BINDIR = 

EXEC=unit
# SRC= unit.c problemDef.c p4estFunc.c geometry.c sem.c
SRC= $(wildcard *.c)
OBJ= $(SRC:%.c=$(LIBDIR)%.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $(BINDIR)$@ $^ $(LDFLAGS) $(GLLIBS)

$(LIBDIR)unit.o: problemDef.h geometry.h p4estFunc.h sem.h conjGrad.h
$(LIBDIR)sem.o:  problemDef.h geometry.h p4estFunc.h
$(LIBDIR)conjGrad.o : sem.h

$(LIBDIR)%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS) $(GLLIBS)

print-%  : ; @echo $* = $($*)

.PHONY: clean mrproper

clean:
	rm -rf $(LIBDIR)*.o
	rm -rf Out/*.vtu
	rm -rf Out/*.pvtu
	rm -rf Out/*.visit

mrproper: clean
	rm -rf $(BINDIR)$(EXEC)