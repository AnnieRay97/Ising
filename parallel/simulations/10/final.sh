fname="avg_"
ext=".ods"
fname2=""
num=0
for k in `seq 0 9`
do
	fname2="$fname$k$ext"
	tail -1 $fname2 >> final.ods
done	