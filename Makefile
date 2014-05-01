CXX=g++-4.7
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic -O2

# # CXX=icpc
# # CXXFLAGS= -Wall -w1 -O2



##for zeuthen
# LIBSLOC=-L/opt/products/gsl/1.15/lib64

LIBS= -lgsl -lgslcblas -lm -static

##all executables
EXECUTABLES = CEPscan_inBrokenPhase \
              plotPotential_inBrokenPhase \
              test_CEP_inBrokenPhase \
              speedTest_CEP_inBrokenPhase \
              generate_list_of_fermionic_contribution \
              example_minimizer_CEP_withFullBosDet

all: ${EXECUTABLES}

# ##########################################
###   executables  ###

CEPscan_inBrokenPhase: ./src/CEPscan_inBrokenPhase.o   \
                       ./src/constrainedEffectivePotential_inBrokenPhase.o  \
                       ./src/CEPscan_inBrokenPhase_helper.o 
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}; echo


plotPotential_inBrokenPhase: ./src/plotPotential_inBrokenPhase.o \
                             ./src/constrainedEffectivePotential_inBrokenPhase.o   \
                             ./src/CEPscan_inBrokenPhase_helper.o 
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}; echo


test_CEP_inBrokenPhase: ./src/test_CEP_inBrokenPhase.o \
                             ./src/constrainedEffectivePotential_inBrokenPhase.o
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}; echo

speedTest_CEP_inBrokenPhase: ./src/speedTest_CEP_inBrokenPhase.o \
                             ./src/constrainedEffectivePotential_inBrokenPhase.o
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}; echo

generate_list_of_fermionic_contribution: ./src/generate_list_of_fermionic_contribution.o \
                             ./src/constrainedEffectivePotential_inBrokenPhase.o
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}; echo

example_minimizer_CEP_withFullBosDet: ./src/example_minimizer_CEP_withFullBosDet.o \
                                      ./src/CEP_withFullBosDet.o
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}; echo

###########################################3
###   OFILES for executables   ###

./src/CEPscan_inBrokenPhase.o: ./src/CEPscan_inBrokenPhase.cc \
                               ./src/constrainedEffectivePotential_inBrokenPhase.h \
                               ./src/CEPscan_inBrokenPhase_helper.h
	cd ./src; make CEPscan_inBrokenPhase.o; echo

./src/plotPotential_inBrokenPhase.o: ./src/plotPotential_inBrokenPhase.cc \
                                     ./src/constrainedEffectivePotential_inBrokenPhase.h \
                                     ./src/CEPscan_inBrokenPhase_helper.h
	cd ./src; make plotPotential_inBrokenPhase.o

./src/test_CEP_inBrokenPhase.o: ./src/test_CEP_inBrokenPhase.cc \
                                ./src/constrainedEffectivePotential_inBrokenPhase.h
	cd ./src; make test_CEP_inBrokenPhase.o; echo

./src/speedTest_CEP_inBrokenPhase.o: ./src/speedTest_CEP_inBrokenPhase.cc \
                                     ./src/constrainedEffectivePotential_inBrokenPhase.h
	cd ./src; make speedTest_CEP_inBrokenPhase.o; echo

./src/generate_list_of_fermionic_contribution.o: ./src/generate_list_of_fermionic_contribution.cc \
                                                 ./src/constrainedEffectivePotential_inBrokenPhase.h
	cd ./src; make generate_list_of_fermionic_contribution.o; echo


./src/example_minimizer_CEP_withFullBosDet.o: ./src/example_minimizer_CEP_withFullBosDet.cc \
                                              ./src/CEP_withFullBosDet.h
	cd ./src; make example_minimizer_CEP_withFullBosDet.o; echo
	
#############################################
### OFILES for libs
./src/constrainedEffectivePotential_inBrokenPhase.o: ./src/constrainedEffectivePotential_inBrokenPhase.cc \
                                                     ./src/constrainedEffectivePotential_inBrokenPhase.h
	cd ./src; make constrainedEffectivePotential_inBrokenPhase.o; echo

./src/CEPscan_inBrokenPhase_helper.o : ./src/CEPscan_inBrokenPhase_helper.cc\
                                        ./src/CEPscan_inBrokenPhase_helper.h
	cd ./src; make CEPscan_inBrokenPhase_helper.o; echo

./src/CEP_withFullBosDet.o: ./src/CEP_withFullBosDet.cc \
                            ./src/CEP_withFullBosDet.h
	cd ./src; make CEP_withFullBosDet.o; echo
clean:
	rm ./src/*.o; rm ${EXECUTABLES}
	


