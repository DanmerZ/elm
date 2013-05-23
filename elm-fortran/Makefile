OBJECTS = Profiles.o main.o Flux.o ProfilesArrays.o ProfilesToFiles.o Q.o Equation.o

main.exe: $(OBJECTS)
	gfortran $(OBJECTS) -o main.exe
	./main.exe
	gnuplot graphics.graph
	ristretto graphics/Q.jpeg

%.o : %.f90
	gfortran -c $<

clean:
	rm -f $(OBJECTS)

