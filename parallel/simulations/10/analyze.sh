fname1="PT_10"
symb="_"
ext=".ods"
source=""
fname3="avg"
dest=""
for k in `seq 0 9`
do
	dest="$fname3$symb$k$ext"
	for i in `seq 0 49`
	do 
		Temp=`expr 10 - $k`
		#echo $Temp
		source="$fname1$symb$i$ext"
		cat $source |tail -$Temp | head -1 >> $dest 
	done
	echo '\t=AVERAGE(B1:B50)\t=AVERAGE(C1:C50)\t=AVERAGE(D1:D50)\t=AVERAGE(E1:E50)\t=AVERAGE(F1:F50)\t=AVERAGE(G1:G50)\t=AVERAGE(H1:H50)\t=AVERAGE(I1:I50)\t=AVERAGE(J1:J50)\t=AVERAGE(K1:K50)\t' >> $dest
done 
#num=`expr 10 - 2`;
#tail -$num "$fname1$symb$num$ext" >> dest.ods