#makefile: makes the fci program


FC = ifort
# flags for debugging or for maximum performance, comment as necessary
FCFLAGS = -O3  -heap-arrays
LDFLAGS = -Wl,--no-as-needed -L/opt/intel/mkl/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread -lm  -m64 -I${MKLROOT}/include 


all: fci

debug: FCFLAGS = -O0 -g -heap-arrays
debug: fci

fci: FCIman.o fciutils.mod FCIWrap.o util90.o
	$(FC) $(FCFLAGS) -o fci FCIman.o fciutils.o FCIWrap.o util90.o $(LDFLAGS)

fciutils.mod: fciutils.f
	$(FC) $(FCFLAGS) -c fciutils.f

FCIman.o: FCIman.f fciutils.mod
	$(FC) $(FCFLAGS) -c FCIman.f

FCIWrap.o: FCIWrap.f90 fciutils.mod
	$(FC) $(FCFLAGS) -c FCIWrap.f90

util90.o: util90.f
	$(FC) $(FCFLAGS) -c util90.f

clean:
	rm fci *.o *.mod



