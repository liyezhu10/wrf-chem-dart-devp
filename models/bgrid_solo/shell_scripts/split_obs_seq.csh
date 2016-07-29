#!/bin/csh 
#
# DART software - Copyright 2004 - 2013 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id$
#
# A shell script to split a series of obs_sequence files
# into a series of multiple smaller files.  This version assumes no calendar
# is being used.  (For gregorian calendar calcs see the version that uses
# advance_time to advance the time.)
#

# set start and end day
@ startd = 0
@   endd = 1 

# set time window for each output file.  in this version of the program
# the cycle_interval MUST be an even divisor of 1 day (86400 seconds)
@ cycle_interval = 3600                        # total interval duration in seconds
@ last_half      = $cycle_interval / 2         # exact half interval
@ first_half     = $cycle_interval / 2 - 1     # half interval - 1 second

echo $cycle_interval, $last_half, $first_half


# set once and should be able to leave as-is
set output_dir   = .
set input_file   = "obs_seq.out"
set nml_template = ../work/input.nml

# loop from start to end time.
@ dtd = $startd
while ( $dtd <= $endd )
  
  echo processing observation files for day $dtd
  @ dts = 0
  @ endwindow = $dts + $last_half
echo dts, endwindow = $dts $endwindow
  while ( $endwindow <= 86400 )
     echo processing window centered on day $dtd  seconds $dts
     @ start_day = $dtd
     @ start_sec = $dts - $first_half
     if ( $start_sec < 0 ) then
        # special case for day 0 only - first window is second half only
        if ( $start_day == 0 ) then
           @ start_sec = 0
        else
          @ start_day --
          @ start_sec += 86400
        endif
     endif
     @ end_day = $dtd
     @ end_sec = $dts + $last_half
     if ( $end_sec > 86400 ) then
         @ end_day ++
         @ end_sec -= 86400
     endif

echo start, end = $start_day, $start_sec, $end_day, $end_sec

    # use time to construct unique output file names
    set dat=`printf %03d_%05d $dtd $dts`
    set filn_out = ${output_dir}/obs_seq.${dat}.out
  
    echo output will go into file $filn_out
  
    # change the values in the template file and overwrite the
    # old input.nml - this controls what the execution of the
    # obs_sequence_tool will do.

    sed -e "s;^ *first_obs_days *=.*;first_obs_days = $start_day;" \
        -e "s;^ *first_obs_seconds *=.*;first_obs_seconds = $start_sec;" \
        -e "s;^ *last_obs_days *=.*;last_obs_days = $end_day;" \
        -e "s;^ *last_obs_seconds *=.*;last_obs_seconds = $end_sec;" \
        -e "s;^ *filename_seq *=.*;filename_seq='${input_file}';" \
        -e "s;^ *filename_out *=.*;filename_out='${filn_out}';" \
        $nml_template >! input.nml 
  
    # do the splitting here
    # If it errors out, just stop. Usually means we have run out of files.
    ../work/obs_sequence_tool  || exit 1
  
    # advance to next time on this day
    @ dts = $dts + $cycle_interval
    @ endwindow = $dts + $last_half
  end

  # advance to the next day
  @ dtd ++

end

echo 'Finished'

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

