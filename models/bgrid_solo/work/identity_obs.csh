#!/bin/csh

set model_size = 28200
set fout = bob.out

rm $fout
set num_obs = 48000

echo "$num_obs"  >> $fout # upper bound on number of observations in sequence
echo "0"         >> $fout # number of copies of data
echo "0"         >> $fout # quality control

set iobs = 0
#while ($i <= $model_size)
set i = 1
while ($i <= $model_size)
   echo "0"   >> $fout
   echo "-$i" >> $fout
   echo "0 0" >> $fout
   echo "1.0" >> $fout

   @ i = $i + 3
   @ iobs++
end

echo "-1" >> $fout
echo "set_def.out" >> $fout
echo ""            >> $fout

cat $fout | ./create_obs_sequence
