CXX=g++-4.7

# CXXFLAGS=-DF_
CXXFLAGS= -Wall -Wextra -Wno-long-long -pedantic 


OFILES = constrainedEffectivePotential_inBrokenPhase.o



# CEPFiles= constrainedEffectivePotential_computeFunctionOnly.cc constrainedEffectivePotential_computeGradientOnly.cc constrainedEffectivePotential_computeFunctionAndGradient.cc constrainedEffectivePotential_computeBosonicLoop.cc constrainedEffectivePotential_computeSecondDerivative.cc

INCL=-I/opt/products/gsl/1.15/include

LIBSLOC=-L/opt/products/gsl/1.15/lib64

LIBS= -lgsl -lgslcblas -lm -static
# ##########################################
###   mainprogs   ###

test_CEP_inBrokenPhase: test_CEP_inBrokenPhase.o ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} ${LIBSLOC}  -o $@ $^ ${LIBS}

test_CEP_inBrokenPhase.o: test_CEP_inBrokenPhase.cc ${OFILES}
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $<

###########################################3
###   OFILES   ###

constrainedEffectivePotential_inBrokenPhase.o: constrainedEffectivePotential_inBrokenPhase.cc constrainedEffectivePotential_inBrokenPhase.h
	${CXX} ${CXXFLAGS} ${INCL} -c -o $@ $< 

clean:
	rm *.o *~

clean_exec:
	rm test_CEP_inBrokenPhase

