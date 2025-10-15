#make main

# -Wno-c++11-extensions —— turns off warnings about recent c++ features

CPP = clang++
CPPFLAGS = -c -O2 -Wno-c++11-extensions -Xclang -fopenmp
#CPPFLAGS = -c -O2 -Wno-c++11-extensions 

main: main.o stokes.o surface.o utilities.o
	$(CPP) main.o stokes.o surface.o utilities.o -lfftw3 -lm -lomp -o main.out
	# $(CPP) main.o stokes.o surface.o utilities.o -lfftw3 -lm -o main.out
	rm *.o

main.o: main.cpp
	$(CPP) $(CPPFLAGS) main.cpp

stokes.o: stokes.cpp stokes.h
	$(CPP) $(CPPFLAGS) stokes.cpp

surface.o: surface.cpp surface.h
	$(CPP) $(CPPFLAGS) surface.cpp

utilities.o: utilities.cpp utilities.h
	$(CPP) $(CPPFLAGS) utilities.cpp

clean:
	rm *.o *~

#when name is different from makefile:
#make -f makefile main
