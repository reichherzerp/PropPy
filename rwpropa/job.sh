#!/bin/sh
for i in `seq 1 10`; 
do python3 python_script.py ${i}& 
done
wait
exit