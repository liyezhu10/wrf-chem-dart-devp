#!/bin/csh

if ($#argv >= 3) then
   set restart_stub = $1
   set out_stub     = $2
   set test_dir     = $3
else
   echo "expecting stage_restarts.csh -restart_stub -out_stub -test_dir"
endif

# optional argument
set stop = 1000 # something really large

if ($#argv == 4) then
   set stop = $4
endif

set n = 1
foreach file (../restarts/*)

   set fileout = `printf "%s%4.4d%s" $out_stub. $n .nc`
   
   if (-d $file) then
      echo "restarts are contained in directory $file"
      echo "linking $file/$restart_stub to $test_dir run directory"
      ln -sf $file . 

      echo "copying $file/$restart_stub to $fileout in $test_dir directory"
      cp $file/$restart_stub $fileout 
   else
      echo "linking $file to $test_dir run directory"
      ln -sf $file .

      echo "copying $file to $fileout in $test_dir directory"
      cp $file $fileout 
   endif 
 
   if ($n == $stop) exit(0)

   @ n = $n + 1

end


