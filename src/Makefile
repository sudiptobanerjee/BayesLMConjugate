JAVA_HOME = ..
SRC = $(JAVA_HOME)/src

# For package BayesLMConjugate
JC = javac
JFLAGS = -d
JAVA_SRC = $(SRC)/BayesLMC.java $(SRC)/fileinput.java
PACKAGE = $(JAVA_HOME)/class/BayesLMConjugate
SUBPACKAGES_CLASS_PATH = $(JAVA_HOME)/class

all: subpackage1 subpackage2 $(PACKAGE)

subpackage1:
	@mkdir -p $(@D)
	$(MAKE) -C ../JAMAJniLite/src

subpackage2:
	@mkdir -p $(@D)
	$(MAKE) -C ../java_Rmath/

$(PACKAGE): subpackage1 subpackage2 $(JAVA_SRC)
	@mkdir -p $(@D)
	$(JC) $(JFLAGS) $(JAVA_HOME)/class -classpath $(SUBPACKAGES_CLASS_PATH):. $(JAVA_SRC)

clean: 
	rm -fr $(JAVA_HOME)/class
	rm -fr $(JAVA_HOME)/lib
