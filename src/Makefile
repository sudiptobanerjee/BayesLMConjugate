JAVA_HOME = ..
SRC = $(JAVA_HOME)/src

# For package BayesLMConjugate
JC = javac
JFLAGS = -d
JAVA_SRC = $(SRC)/BayesLMC.java $(SRC)/fileinput.java
PACKAGE = $(JAVA_HOME)/class/BayesLMConjugate
JAMA_RMATH_CLASS_PATH = $(JAVA_HOME)/class

all: $(PACKAGE)

$(PACKAGE): $(JAVA_SRC)
	@mkdir -p $(@D)
	$(JC) $(JFLAGS) $(JAVA_HOME)/class -classpath $(JAMA_RMATH_CLASS_PATH):. $(JAVA_SRC)

clean: 
	rm -fr $(JAVA_HOME)/class/BayesLMConjugate
	rm -f $(JAVA_HOME)/test/*.class
	rm -f $(JAVA_HOME)/test/*.txt
	rm -f $(JAVA_HOME)/test/*.out
	rm -f $(JAVA_HOME)/test/*.serial
