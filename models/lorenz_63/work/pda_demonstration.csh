#!/bin/csh


./mkmf_preprocess
make
./mkmf_create_obs_sequence 
make

set model_size = 3
set fout = t1_63.out 
rm $fout
set num_obs = 3

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

./mkmf_create_fixed_network_seq 
make

set fout2=t2_63.out

rm $fout2

set num_obs = 500

echo "set_def.out" >> $fout2
echo "1"           >> $fout2 # regularly repeated time sequence
echo "num_obs"     >> $fout2
echo "0 0"         >> $fout2 # initial time
echo "0 3600"      >> $fout2 # time step days seconds

echo "obs_seq.in"  >> $fout2 

cat $fout2 | ./create_fixed_network_sequence

# run perfect model obs to generate obs_seq.out

./mkmf_perfect_model_obs
make
./perfect_model_obs

# exit(0)

csh ../shell_scripts/quickbuild.csh

cp obs_seq.out ../shell_scripts/

##### edit split_obs_seq.csh, make modification on start and end day and cycle_interval if necessary

cd ../shell_scripts

csh split_obs_seq.csh 


