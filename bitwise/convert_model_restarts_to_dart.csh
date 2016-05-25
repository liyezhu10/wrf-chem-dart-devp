#!/bin/tcsh

# convert model restart to dart restart
if ($#argv >= 5) then 
   set model_to_dart  = "$1"
   set dart_restart   = "$2"
   set model_restart  = "$3"
   set out_stub       = "$4"
   set filter_restart = "$5"
else
   echo " expecting -model_to_dart -dart_restart -model_restart -out_stub -filter_restart"
   exit(1)
endif

set stop = 1000
if ($#argv == 6) then 
   set stop = $6
endif

set n = 1
foreach file (../restarts/*)

   #set num     = `echo $file | tail -c 5`
   set fileout = `printf "%s%4.4d" $filter_restart. $n`

   if (-d $file) then
      echo "restarts directory contains directory of restarts"
      echo "converting $file/$model_restart to $fileout"

      cp $file/$model_restart $model_restart
      ./$model_to_dart >& model_to_dart.out

      cp $dart_restart $fileout
   else
      echo "converting $file to $fileout"

      cp $file $model_restart
      ./$model_to_dart >& model_to_dart.out

      cp $dart_restart $fileout
   endif

   if ($n == stop) exit(0)
 
   @ n = $n + 1

end
