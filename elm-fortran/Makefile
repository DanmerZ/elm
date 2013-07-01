OBJECTS = Profiles.o ArDef.o main.o Flux.o ProfilesArrays.o ProfilesToFiles.o \
          Q.o Bootstrap.o Evans60.o compare.o#Equation.o

main.exe: $(OBJECTS)
	gfortran $(OBJECTS) -o main.exe
	./main.exe

%.o : %.f90
	gfortran -c $<

clean:
	rm -f $(OBJECTS)

