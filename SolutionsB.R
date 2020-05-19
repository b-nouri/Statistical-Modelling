DataQB = read.table("D:\\BankDefaultData.txt",header=T)
studentnumber = 772586 # fill in your student number here, this is an example!
set.seed(studentnumber)
rownumbers = sample(1:6436,size=1000)
mydataB = DataQB[rownumbers,]

