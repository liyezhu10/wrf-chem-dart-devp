#!/bin/csh

if ($#argv >= 3) then
   set restart_stub = $1
   set out_stub     = $2
   set test_dir     = $3
else
   echo "expecting stage_restarts.csh -restart_stub -out_stub -test_dir -stop -wrf2dom"
endif

# optional argument
set stop = 1000 # something really large

if ($#argv == 4) then
   set stop = $4
endif

if ($#argv == 5) then
   set wrf2dom = $5
else
   set wrf2dom = "false" 
endif

set n = 1
foreach file (../restarts/*)

   if ("$wrf2dom" == "true") then
     set d01 = "_d01"
     set d02 = "_d02"
     set fileout1 = `printf "%s%4.4d%s" $out_stub$d01. $n .nc`
     set fileout2 = `printf "%s%4.4d%s" $out_stub$d02. $n .nc`
     set restart_stub1 = "$restart_stub$d01"
     set restart_stub2 = "$restart_stub$d02"
   else 
     set fileout1 = `printf "%s%4.4d%s" $out_stub. $n .nc`
   endif

   if (-d $file) then
      echo "restarts are contained in directory $file"
      echo "linking $file/$restart_stub to $test_dir run directory"
      ln -sf $file . 

      touch $fileout1
      chmod u+w $fileout1

      if ("$wrf2dom" == "true") then
         echo "copying $file/$restart_stub1 to $fileout1 in $test_dir directory"
         echo "copying $file/$restart_stub2 to $fileout2 in $test_dir directory"
         touch     $fileout2
         chmod u+w $fileout2
         cp $file/$restart_stub1 $fileout1
         cp $file/$restart_stub2 $fileout2
      else
         echo "copying $file/$restart_stub to $fileout1 in $test_dir directory"
         cp $file/$restart_stub $fileout1 
      endif
   else
      echo "linking $file to $test_dir run directory"
      ln -sf $file .

      touch     $fileout1
      chmod u+w $fileout1

      if ("$wrf2dom" == "true") then
         echo "copying $file/$restart_stub1 to $fileout1 in $test_dir directory"
         echo "copying $file/$restart_stub2 to $fileout2 in $test_dir directory"
         cp $file/$restart_stub1 $fileout1
         cp $file/$restart_stub2 $fileout2
      else
         echo "copying $file to $fileout1 in $test_dir directory"
         cp $file $fileout1 
      endif
   endif 
 
   if ($n == $stop) exit(0)

   @ n = $n + 1

end


