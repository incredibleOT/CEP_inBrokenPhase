CXX=g++-4.7
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic -O2

# # CXX=icpc
# # CXXFLAGS= -Wall -w1 -O2


##for zeuthen
# INCL=-I/opt/products/gsl/1.15/include

LIBS= -lgsl -lgslcblas -lm -static



###########################################
###   OFILES   ###
###########################################
##for libs
constrainedEffectivePotential_inBrokenPhase.o: constrainedEffectivePotential_inBrokenPhase.cc \
                                               constrainedEffectivePotential_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

CEPscan_inBrokenPhase_helper.o: CEPscan_inBrokenPhase_helper.cc \
                                CEPscan_inBrokenPhase_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

###########################################
##for executables
CEPscan_inBrokenPhase.o: CEPscan_inBrokenPhase.cc \
                         constrainedEffectivePotential_inBrokenPhase.h \
                         CEPscan_inBrokenPhase_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

generate_list_of_fermionic_contribution.o: generate_list_of_fermionic_contribution.cc \
                                           constrainedEffectivePotential_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

test_CEP_inBrokenPhase.o: test_CEP_inBrokenPhase.cc \
                          constrainedEffectivePotential_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

speedTest_CEP_inBrokenPhase.o: speedTest_CEP_inBrokenPhase.cc \
                               constrainedEffectivePotential_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

plotPotential_inBrokenPhase.o: plotPotential_inBrokenPhase.cc \
                               constrainedEffectivePotential_inBrokenPhase.h \
                               CEPscan_inBrokenPhase_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo


clean:
	rm *.o


