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
CEP_inBrokenPhase.o: CEP_inBrokenPhase.cc \
                                               CEP_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

CEPscan_helper.o: CEPscan_helper.cc \
                                CEPscan_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

CEP_withFullBosDet.o: CEP_withFullBosDet.cc \
                      CEP_withFullBosDet.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo
###########################################
##for executables
CEPscan_inBrokenPhase.o: CEPscan_inBrokenPhase.cc \
                         CEP_inBrokenPhase.h \
                         CEPscan_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

generate_list_of_fermionic_contribution.o: generate_list_of_fermionic_contribution.cc \
                                           CEP_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

test_CEP_inBrokenPhase.o: test_CEP_inBrokenPhase.cc \
                          CEP_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

speedTest_CEP_inBrokenPhase.o: speedTest_CEP_inBrokenPhase.cc \
                               CEP_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

plotPotential_inBrokenPhase.o: plotPotential_inBrokenPhase.cc \
                               CEP_inBrokenPhase.h \
                               CEPscan_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

CEPscan_withFullBosDet.o: CEPscan_withFullBosDet.cc \
                          CEP_withFullBosDet.h \
                          CEPscan_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

example_minimizer_CEP_withFullBosDet.o: example_minimizer_CEP_withFullBosDet.cc \
                                        CEP_withFullBosDet.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

plotPotential_withFullBosDet.o: plotPotential_withFullBosDet.cc \
                               CEP_withFullBosDet.h \
                               CEPscan_helper.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<; echo

	
clean:
	rm *.o


