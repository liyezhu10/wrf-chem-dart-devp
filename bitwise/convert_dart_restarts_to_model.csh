#!/bin/tcsh

if ($#argv >= 5) then 
   set dart_to_model  = "$1"
   set dart_restart   = "$2"
   set model_restart  = "$3"
   set out_stub       = "$4"
   set filter_restart = "$5"
else
   echo " expecting -dart_to_model -dart_restart -model_restart -out_stub -filter_restart"
   exit(1)
endif

set stop = 1000
if ($#argv == 6) then 
   set stop = "$6"
endif

if ($#argv == 7) then 
   set wrf2dom = "$7"
endif

set n = 1
foreach file (`ls $filter_restart.*`)


   cp $file $dart_restart
   ./$dart_to_model > out.txt

   if ("$wrf2dom" == "true") then
     set num     = `echo $file | tail -c 5`
     set d01 = "_d01"
     set d02 = "_d02"
     set outfile1 = `printf "%s%4.4d%s" $out_stub$d01. $n .nc`
     set outfile2 = `printf "%s%4.4d%s" $out_stub$d02. $n .nc`

     echo "converting $file to $outfile1, $outfile2"
     cp $model_restart$d01 $outfile1
     cp $model_restart$d02 $outfile2
   else
     set num     = `echo $file | tail -c 5`
     set outfile1 = "$out_stub.$num.nc"

     echo "converting $file to $outfile1"
     cp $model_restart $outfile1
   endif

   if ($n == $stop) exit(0)

   @ n = $n + 1
   
end

rm out.txt
