DataQA = read.table("D:\\covidDATA.txt")
modulo11=function(x) {x- floor(x/11)*11}
studentnumber=772586 # fill in your student number here, this is an example!
mycol=modulo11(studentnumber)
mydataA=DataQA[,c((1+mycol*4):(4+mycol*4))]

