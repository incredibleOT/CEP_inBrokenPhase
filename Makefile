CXX=g++-4.7
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic -O2

# # CXX=icpc
# # CXXFLAGS= -Wall -w1 -O2

# ALLOFILES = ./src/constrainedEffectivePotential_inBrokenPhase.o \
#          ./src/CEPscan_inBrokenPhase_helper.o  \
#          ./src/CEPscan_inBrokenPhase.o  \
#          ./src/generate_list_of_fermionic_contribution.o  \
#          ./src/test_CEP_inBrokenPhase.o  \
#          ./src/plotPotential_inBrokenPhase.o  \
#          ./speedTest_CEP_inBrokenPhase.o  

LIBOFILES = ./src/constrainedEffectivePotential_inBrokenPhase.o \
            ./src/CEPscan_inBrokenPhase_helper.o  \
            ./src/plotPotential_inBrokenPhase_helper.o


INCL=-I/opt/products/gsl/1.15/include

LIBSLOC=-L/opt/products/gsl/1.15/lib64

LIBS= -lgsl -lgslcblas -lm -static


# ##########################################
###   mainprogs   ###

CEPscan_inBrokenPhase: ./src/CEPscan_inBrokenPhase.o ${LIBOFILES}
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}


plotPotential_inBrokenPhase: ./src/plotPotential_inBrokenPhase.o ${LIBOFILES} ./src/plotPotential_inBrokenPhase_helper.o
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}


test_CEP_inBrokenPhase: ./src/test_CEP_inBrokenPhase.o ${LIBOFILES}
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}

speedTest_CEP_inBrokenPhase: ./src/speedTest_CEP_inBrokenPhase.o ${LIBOFILES}
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}

generate_list_of_fermionic_contribution: ./src/generate_list_of_fermionic_contribution.o ${LIBOFILES}
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}

###########################################3
###   OFILES   ###

src/%.o:
	cd src; make $(@F); cd ..

# constrainedEffectivePotential_inBrokenPhase.o: constrainedEffectivePotential_inBrokenPhase.cc constrainedEffectivePotential_inBrokenPhase.h
# 	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $< 
# 
# CEPscan_inBrokenPhase_helper.o: CEPscan_inBrokenPhase_helper.cc CEPscan_inBrokenPhase_helper.h
# 	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $< 
# 
# plotPotential_inBrokenPhase_helper.o: plotPotential_inBrokenPhase_helper.cc plotPotential_inBrokenPhase_helper.h
# 	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $< 

clean:
	rm ./src/*.o ./src/*~


