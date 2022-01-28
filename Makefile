METHOD = LAX
#METHOD = MACCORMACK
MODEL = SN
#MODEL = SHOCK
include methods/Makefile
include initial/Makefile

hidro : globals.o \
		initial.o \
			physics.o \
				methods.o \
					general.o \
						main.o
	@echo ""
	@echo ""
	@echo "...Compiling files..."
	@echo ""
	@-gfortran -o hidro main.o globals.o general.o method.o initial.o physics.o
	@echo ""
	@echo "################### Makefile parameters ###################"
	@echo "" 
	@echo " METHOD = "$(METHOD)
	@echo " PROBLEM = "$(MODEL)
	@echo ""
	@echo "###########################################################"
	@echo "" 
	
globals.o : globals.f90
	gfortran -c -k8 globals.f90

initial.o : $(INITIAL)
	gfortran -c -k8 $(INITIAL)

physics.o : physics/physics.f90
	gfortran -c -k8 physics/physics.f90

methods.o : $(SOLVE)
	gfortran -c -k8 $(SOLVE)

general.o : general.f90
	gfortran -c -k8 general.f90

main.o : main.f90
	gfortran -c -k8 main.f90

clean:
	rm -f *.o *.mod *.dat *png *.gif hidro
