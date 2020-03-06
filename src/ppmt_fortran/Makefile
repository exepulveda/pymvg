OBJECTS = backtrmod.o solve.o eig.o acorni.o quicksort.o normaldist.o ppmt_module.o ppmt_ext.o
FC := gfortran
LINKFLAGS := -llapack
FCFLAGS := -Wall #unused-variable

all: $(OBJECTS)

python_extension: ppmt.pyf $(OBJECTS)
	python -m numpy.f2py -c ppmt.pyf $(OBJECTS) $(LINKFLAGS)

ppmt.pyf: ppmt_ext.f90
	python -m numpy.f2py ppmt_ext.f90 -m ppmt_interface -h ppmt.pyf --overwrite-signature

acorni.o: acorni.f90
	$(FC) $(FCFLAGS) -c acorni.f90

backtrmod.o: backtrmod.f90
	$(FC) $(FCFLAGS) -c backtrmod.f90

solve.o: solve.f90
	$(FC) $(FCFLAGS) -c solve.f90

eig.o: eig.f90
	$(FC) $(FCFLAGS) -c eig.f90

quicksort.o: quicksort.f90
	$(FC) $(FCFLAGS) -c quicksort.f90

normaldist.o: normaldist.f90
	$(FC) $(FCFLAGS) -c normaldist.f90

ppmt_module.o: ppmt_module.f90
	$(FC) $(FCFLAGS) -c ppmt_module.f90

ppmt_ext.o: ppmt_ext.f90
	$(FC) $(FCFLAGS) -c ppmt_ext.f90

clean:
	rm -fr *.o *.mod *.pyf *.so
