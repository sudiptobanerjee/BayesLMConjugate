JAVA_HOME = ..

TEST_CLASS = Test.class
TEST_RUN = test_result.txt
JAMA_RMATH_CLASS_PATH = $(JAVA_HOME)/class

all:  $(TEST_CLASS) $(TEST_RUN) 

$(TEST_CLASS):
	javac -classpath $(JAVA_HOME)/class:$(JAMA_RMATH_CLASS_PATH):. Test.java

$(TEST_RUN):
	java -Djava.library.path=$(JAVA_HOME)/lib -classpath $(JAVA_HOME)/class:$(JAMA_RMATH_CLASS_PATH):. Test > test_result.txt

clean: 
	rm -f *.txt *.class *.out *.serial
