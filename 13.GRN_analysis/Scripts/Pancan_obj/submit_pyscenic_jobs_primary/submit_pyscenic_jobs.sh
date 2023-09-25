###Run over 10 iterations
for iter in iter_1 iter_2 iter_3 iter_4 iter_5 iter_6 iter_7 iter_8 iter_9 iter_10

do

bash ../Pyscenic.mainSteps.cli.v4.sh -i $iter &>>pyscenic.err.log.p1.txt

done

