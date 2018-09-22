DEST = .
CC = gcc
CFLAGS = $(INCLUDE) $(LIB)

JAVA_HOME = ..
JAVA_PACKAGES = $(JAVA_HOME)/class
JAVA_LIB = $(JAVA_HOME)/lib

JNI_SRC = $(JAVA_HOME)/java_Rmath

INCLUDE = -I$(JAVA_HOME)/java_Rmath/include
LIB = -L/usr/local/lib -L$(JAVA_LIB)
C_SRC = $(JNI_SRC)/c/jniRmath.c
DYLIB = $(JAVA_HOME)/lib/libjni_java_Rmath.so
	

# For package java_Rmath
JC = javac
JFLAGS = -d
JAVA_SRC = $(JNI_SRC)/java/jniRmath.java
PACKAGE = $(JAVA_PACKAGES)/java_Rmath

# header file jniRmath.h


all: $(PACKAGE) $(DYLIB)

$(PACKAGE): 
	$(JC) $(JFLAGS) $(JAVA_PACKAGES) $(JAVA_SRC)

$(DYLIB): $(C_SRC)
	$(CC) $(CFLAGS) -shared -fPIC $(C_SRC) -lRmath -o $(DYLIB)

clean: 
	rm -fr $(JAVA_PACKAGES)/java_Rmath
	rm -f $(JAVA_HOME)/lib/libjni_java_Rmath.so
