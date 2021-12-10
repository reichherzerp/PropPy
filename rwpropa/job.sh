#!/bin/sh
file_name="test_new_"
nr=10
for i in `seq 1 $nr`; 
do python3 server_script.py "$file_name${i}"& 
done
wait
python3 server_merge_script.py $file_name $nr
exit