# makefile for Forward_Acoustic

FC = ifort
OPTS = -O3 -parallel -qopt-matmul -qopenmp -mkl=parallel -fpp -Ofast
LINKOPT =
OBJ = o

# Debug
#OPTS = -O0 -WB -Winline -traceback -heap-arrays -fp-stack-check -warn -CB -CU -fpp -qopenmp

#LIB =


UNAME := $(shell uname -a)
ifeq ($(findstring CYGWIN_NT, $(UNAME)), CYGWIN_NT)
   OBJ = obj
   OPTS = /nologo /O3 /Qparallel /Qopt-matmul /Qopenmp /Qmkl:parallel /fpp /QxHost
   LINKOPT = /F150994944
endif

EXE = Acoustic_Forward

OBJS = \
	$ butter.$(OBJ) \
	$ driver.$(OBJ)

#all: build install

#build: $(OBJS)
all: $(OBJS)
	rm -rf $(EXE)
	$(FC) $(OPTS) $(LINKOPT) -c $< -o $@
%.$(OBJ): %.F90
	$(FC) $(OPTS) $(LINKOPT) -c $< -o $@
%.$(OBJ): %.f90
	$(FC) $(OPTS) $(LINKOPT) -c $< -o $@
%.$(OBJ): %.for
	$(FC) $(OPTS) $(LINKOPT) -c $< -o $@
%.$(OBJ): %.For
	$(FC) $(OPTS) $(LINKOPT) -c $< -o $@

install:
	$(FC) $(OPTS) $(LINKOPT) $(OBJS) -o $(EXE)

clean:
	rm -rf *.mod *__genmod.f90 *.$(OBJ) *.dwf *.pdb
