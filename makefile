# compilador
FF = gfortran
FFLAGSDEB = -Wall -g -cpp -fbounds-check -Wmaybe-uninitialized -Wuninitialized
FFLAGSOP = -O3 -cpp
# dependencias para la mac:
DEPS = -I/usr/local/include -L/usr/local/lib -lfftw3 -framework Accelerate

# dependencias para la mac con graficacion con DISLIN
DEPSGRAF = -I/usr/local/include -L/usr/local/lib -lfftw3 -I/usr/local/dislin/gf -L/opt/X11/lib -lXt -lXm -L/usr/OpenMotif/lib -L/usr/local/dislin -ldislin -lm -framework Accelerate

ARCHSUELTOS = guardarResultados.f90 preProcessAtThisfrec.f90 pontertable.f90 predimension.f90 setupmodel.f90 usefftw.a fitting.a matrizGlobal.a variables.a

GreenDT : main.F90 $(ARCHSUELTOS)
	$(FF) -o $@ $< $(ARCHSUELTOS) $(FFLAGSDEB) $(DEPS)

matrizGlobal.a: matrizGlobal.F90 variables.a
	$(FF) -O3 -ffast-math -c matrizGlobal.F90 $(FFLAGSDEB) -O3 $(DEPS)
	ar cr matrizGlobal.a matrizGlobal.o
	rm *.o

usefftw.a: useFFTW.f90 variables.a
	$(FF) -O3 -ffast-math -c useFFTW.f90 $(FFLAGSOP) $(DEPS)
	ar cr usefftw.a usefftw.o
	rm *.o

fitting.a: fitting.f90 variables.a
	$(FF) -O3 -ffast-math -c fitting.f90 $(FFLAGSOP) $(DEPS)
	ar cr fitting.a fitting.o
	rm *.o

variables.a: variables.f90
	$(FF) -O3 -ffast-math -c variables.f90 $(FFLAGSOP) $(DEPS)
	ar cr variables.a variables.o
	rm *.o

clean:
	rm *.mod
	rm *.a
	rm GreenDT

