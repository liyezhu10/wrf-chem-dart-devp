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

set n = 1
foreach file (`ls $filter_restart.*`)

   set num     = `echo $file | tail -c 5`
   set outfile = "$out_stub.$num.nc"

   echo "converting $file to $outfile"

   cp $file $dart_restart
   ./$dart_to_model > out.txt

   cp $model_restart $outfile

   if ($n == $stop) exit(0)

   @ n = $n + 1
   
end

rm out.txt
