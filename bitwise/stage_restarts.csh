#!/bin/csh

if ($#argv >= 3) then
   set restart_stub = $1
   set out_stub     = $2
   set test_dir     = $3
else
   echo "expecting stage_restarts.csh -restart_stub -out_stub -test_dir -stop "
endif

# optional argument
set stop = 1000 # something really large

if ($#argv == 4) then
   set stop = $4
endif

set n = 1
foreach file (../restarts/*.r.*)

   set dom = "_d01"
   set fileout1 = `printf "%s%4.4d%s" $out_stub$dom. $n .nc`

   echo "linking $file to $test_dir run directory"
   ln -sf $file .

   touch     $fileout1
   chmod u+w $fileout1

   echo "copying $file to $fileout1 in $test_dir directory"
   cp $file $fileout1 
 
   if ($n == $stop) break

   @ n = $n + 1

end


set n = 1
foreach file (../restarts/*.h0.*)

   set dom = "_d02"
   set fileout1 = `printf "%s%4.4d%s" $out_stub$dom. $n .nc`

   echo "linking $file to $test_dir run directory"
   ln -sf $file .

   touch     $fileout1
   chmod u+w $fileout1

   echo "copying $file to $fileout1 in $test_dir directory"
   cp $file $fileout1 
 
   if ($n == $stop) break

   @ n = $n + 1

end

set n = 1
foreach file (../restarts/*.h1.*)

   set dom = "_d03"
   set fileout1 = `printf "%s%4.4d%s" $out_stub$dom. $n .nc`

   echo "linking $file to $test_dir run directory"
   ln -sf $file .

   touch     $fileout1
   chmod u+w $fileout1

   echo "copying $file to $fileout1 in $test_dir directory"
   cp $file $fileout1 
 
   if ($n == $stop) break

   @ n = $n + 1

end
