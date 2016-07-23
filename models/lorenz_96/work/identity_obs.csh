#!/bin/csh

set model_size = 40
set fout = bob.out 

rm $fout
set num_obs = 40

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

   @ i = $i + 1
   @ iobs++
end

#echo "-1" >> $fout
echo "set_def.out" >> $fout
echo ""            >> $fout

cat $fout | ./create_obs_sequence


# exit(0)

set fout2


set num_obs = 4000

echo "set_def.out" >> $fout2
echo "num_obs"     >> $fout2
echo "1"           >> $fout2 # regularly repeated time sequence
echo "0 0"         >> $fout2 # initial time
echo "0 3600"      >> $fout2 # time step days seconds

echo "obs_seq.in"  >> $fout2 

cat $fout2 | ./create_fixed_network_sequence

# run perfect model obs to generate obs_seq.out

./mkmf_perfect_model_obs
make
./perfect_model_obs


