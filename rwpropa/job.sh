#!/bin/sh
file_name="/home/patrick/Documents/rwpropa/test/test_"
nr=5
for i in `seq 1 $nr`; 
do python3 server_script.py $file_name ${i} ${nr}& 
done
wait
python3 server_merge_script.py $file_name $nr
exit