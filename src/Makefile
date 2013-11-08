CXX=g++-4.7
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic -O2

# # CXX=icpc
# # CXXFLAGS= -Wall -w1 -O2

LIBOFILES = constrainedEffectivePotential_inBrokenPhase.o CEPscan_inBrokenPhase_helper.o

EXECOFILES=CEPscan_inBrokenPhase.o  generate_list_of_fermionic_contribution.o  \
            test_CEP_inBrokenPhase.o  plotPotential_inBrokenPhase.o  speedTest_CEP_inBrokenPhase.o  


INCL=-I/opt/products/gsl/1.15/include

LIBSLOC=-L/opt/products/gsl/1.15/lib64

LIBS= -lgsl -lgslcblas -lm -static

libobjects: ${LIBOFILES}

# ##########################################
###   mainprogs   ###

# CEPscan_inBrokenPhase: CEPscan_inBrokenPhase.o ${OFILES}
# 	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}


# plotPotential_inBrokenPhase: plotPotential_inBrokenPhase.o constrainedEffectivePotential_inBrokenPhase.o plotPotential_inBrokenPhase_helper.o
# 	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}


# test_CEP_inBrokenPhase: test_CEP_inBrokenPhase.o ${OFILES}
# 	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}


# speedTest_CEP_inBrokenPhase: speedTest_CEP_inBrokenPhase.o ${OFILES}
# 	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}


# generate_list_of_fermionic_contribution: generate_list_of_fermionic_contribution.o ${OFILES}
# 	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}


###########################################3
###   OFILES   ###
##for libs
constrainedEffectivePotential_inBrokenPhase.o: constrainedEffectivePotential_inBrokenPhase.cc constrainedEffectivePotential_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $< 

CEPscan_inBrokenPhase_helper.o: CEPscan_inBrokenPhase_helper.cc CEPscan_inBrokenPhase_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $< 

plotPotential_inBrokenPhase_helper.o: plotPotential_inBrokenPhase_helper.cc plotPotential_inBrokenPhase_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $< 

##for executables
CEPscan_inBrokenPhase.o: CEPscan_inBrokenPhase.cc 
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<

generate_list_of_fermionic_contribution.o: generate_list_of_fermionic_contribution.cc 
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<

test_CEP_inBrokenPhase.o: test_CEP_inBrokenPhase.cc
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<

speedTest_CEP_inBrokenPhase.o: speedTest_CEP_inBrokenPhase.cc 
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<

plotPotential_inBrokenPhase.o: plotPotential_inBrokenPhase.cc
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<


clean:
	rm *.o *~


