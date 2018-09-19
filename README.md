# BayesLMConjugate
BayesLMConjugate is a Java package providing methods for Bayesian Conjugate Linear Modeling which is widely used in Bayesian Geostatisitics. This package uses Java classes created by two other packages, JAMAJniLite and java_Rmath. JAMAJniLite is a highly computing efficient and user-friendly numerical linear algebra package, which helps to perform matrix operations in BayesLMConjugate. java_Rmath is a random number generator which helps BayesLMConjugate to do sampling. 


Build Instructions
------------------

* To compile the package, enter src directory and execute "make".

* To clean generated files, type “make clean” on the command line. 

Data Format Preprocessing
-----------------

* To preprocess data, you need to install R software in advance.

* Enter data directory. Then open R environment and type > source("DataFormat.R"). Then you will find two text files named X and Y generated under data directory.

* Notice: in your raw data, you should put Y in the first column and X in the rest columns.

Running the test
-----------------
* For testing, enter test directory and execute “make” . If you want to clean testing results and all class files, type "make clean". 

* Notice: before running test, you have to preprocess data and execute make in the src directory.

Source Repository
-----------------
BayesLMConjugate's source-code repository is hosted here on GitHub: 
https://github.com/pkuxchen/BayesLMConjugate


Authors
---------

| Name   | Email       |              |
|:------ |:----------- | :----------- |
| Xiang Chen | pkuchenxiang@pku.edu.cn   | Visiting student, Department of Biostatistics  UCLA|                         
| Sudipto Banerjee | sudipto@ucla.edu   | Professor, Department of Biostatistics  UCLA |
<!--- --->
                             


Licensing
---------
BayesLMConjugate is licensed under the Creative Commons Attribution License. 



