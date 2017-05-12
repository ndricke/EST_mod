#makefile: makes the fci program


FC = ifort
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -O3  -heap-arrays
LDFLAGS = -Wl,--no-as-needed -L/opt/intel/mkl/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm  -m64 -I${MKLROOT}/include 


all: hred

debug: FCFLAGS = -O0 -g -heap-arrays
debug: hred

hred: FCIman.o fciutils.mod Hred.o util90.o
	$(FC) $(FCFLAGS) -o hred FCIman.o fciutils.o Hred.o util90.o $(LDFLAGS)

fciutils.mod: fciutils.f
	$(FC) $(FCFLAGS) -c fciutils.f

FCIman.o: FCIman.f fciutils.mod
	$(FC) $(FCFLAGS) -c FCIman.f

Hred.o: Hred.f90 fciutils.mod
	$(FC) $(FCFLAGS) -c Hred.f90

util90.o: util90.f
	$(FC) $(FCFLAGS) -c util90.f

clean:
	rm hred *.o *.mod



