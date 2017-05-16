FC=ifort
FC=gfortran
subjects1=cut_events.o sacio.o taper.o
subjects2=cut_events_2016_12_21.o sacio.o taper.o
all:sacio.mod cut_events  cut_events_2016_12_21
sacio.mod:sacio.f90
	$(FC) $^ -c
cut_events:$(subjects1)
	$(FC) $^ -o $@ 
cut_events_2016_12_21:$(subjects2)
	$(FC) $^ -o $@ 
%.o:%.f90
	$(FC) $^ -c
install:
	cp cut_events cut_events_2016_12_21 ../bin
clean:
	rm *.o *.mod
