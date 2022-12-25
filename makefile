# Executable
EXE = program.exe

# Compiler
CC = gcc
FC = gfortran

# Compiler Library
CLIB = -lm -lgfortran -lquadmath

# Compiler Flags
CFLG = -Wall -Wextra -pedantic -Ofast -march=native
F90FLG = -Wall -Ofast -static -march=native
F77FLG = -Ofast -static -march=native

# Directories
BDIR = bin
HDIR = include
ODIR = obj
SDIR = src

# Files
CSRC = $(wildcard $(SDIR)/*.c)
COBJ = $(patsubst $(SDIR)/%.c, $(ODIR)/%.o, $(CSRC))
F90SRC = $(wildcard $(SDIR)/*.f90)
F90OBJ = $(patsubst $(SDIR)/%.f90, $(ODIR)/%.o, $(F90SRC))
F77SRC = $(wildcard $(SDIR)/*.f)
F77OBJ = $(patsubst $(SDIR)/%.f, $(ODIR)/%.o, $(F77SRC))

# Targets
all: $(BDIR)/$(EXE)

$(BDIR)/$(EXE): $(COBJ) $(F90OBJ) $(F77OBJ)
	$(CC) -o $@ $^ -L./ $(CLIB)

$(ODIR)/%.o: $(SDIR)/%.c
	$(CC) $(CFLG) -c $^ -o $@ -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.f90
	$(FC) $(F90FLG) -c $^ -o $@ -I $(HDIR)

$(ODIR)/%.o: $(SDIR)/%.f
	$(FC) $(F77FLG) -c $^ -o $@ -I $(HDIR)

# Clean
.PHONY: clean
clean:
	rm -f $(BDIR)/$(EXE) $(ODIR)/*.o
