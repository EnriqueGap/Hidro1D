#METHOD = LAX
METHOD = MACCORMACK
include methods/Makefile

hidro : globals.o \
		general.o \
			methods.o \
				main.o
	@echo ""
	@echo ""
	@echo "...Compiling files..."
	@echo ""
	@-gfortran -o hidro main.o globals.o general.o method.o
	@echo ""
	@echo "################### Makefile parameters ###################"
	@echo "" 
	@echo " METHOD = "$(METHOD)
	@echo ""
	@echo "###########################################################"
	@echo "" 
	
globals.o : globals.f90
	gfortran -c globals.f90

general.o : general.f90
	gfortran -c general.f90
	
methods.o : $(SOLVE)
	gfortran -c $(SOLVE)

main.o : main.f90
	gfortran -c main.f90

clean:
	rm -f *.o *.mod *.dat *png *.gif hidro
