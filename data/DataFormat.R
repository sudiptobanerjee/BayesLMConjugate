## R Program to format data for input in ../test/Test.java

data.input = read.table("FORMGMT.dat.tab", header=F)

Y = data.input[,1]
write.table(Y, "Y.txt", row.names=F, col.names=F)

X = data.input[,2:7]
X = cbind(1, X)
write.table(X, "X.txt", row.names=F, col.names=F)

