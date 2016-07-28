#!/bin/csh
#
# DART software - Copyright 2004 - 2016 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# DART $Id: $

# From Charlie Zender, NCO grand-meister:
# use ncwa to average-over, and thus remove the dimensions you do not want.
# if you don't want to average-over a dimension before plotting, then hyperslab the region you do want, e.g.,
# 
# ncwa -a Band in.nc out.nc
# 
# will remove the Band dimension by averaging over it whereas
# 
# ncwa -a Band,3 in.nc out.nc 
# 
# selects only band 3 data and then removes the band dimension, 
# which is now degenerate, from it.

rm filter_ics.*
rm truth.nc

set fname1 = True_state.nc
set fname2 = truth.nc

ncwa -a copy,0 $fname1 $fname2

ncwa -O -a time $fname2 $fname2
ncwa -O -a metadatalength $fname2 $fname2
ncwa -O -a NMLnlines $fname2 $fname2
ncwa -O -a NMLlinelen $fname2 $fname2
ncrename -d TmpI,lon -d TmpJ,lat -d VelI,slon -d VelJ,slat $fname2
ncrename -v ps,PS -v t,T -v u,U -v v,V -v TmpI,lon -v TmpJ,lat -v VelI,slon -v VelJ,slat $fname2
ncks -C -O -x -v time $fname2 $fname2 
# ncks -C -O -x -v TmpI $fname2 $fname2
# ncks -C -O -x -v TmpJ $fname2 $fname2
# ncks -C -O -x -v VelI $fname2 $fname2
# ncks -C -O -x -v VelJ $fname2 $fname2
ncks -C -O -x -v copy $fname2 $fname2
ncks -C -O -x -v lev $fname2 $fname2
ncks -C -O -x -v CopyMetaData $fname2 $fname2
ncks -C -O -x -v inputnml $fname2 $fname2
ncks -A time.nc $fname2

@ trip = 2
while ( $trip < 82 )

  @ member = $trip - 1
  set fname = `printf filter_ics.%04d.nc $member`

  echo "making trip $trip to create file $fname"

  #ncwa -a copy,$member Prior_Diag.nc $fname
  ncks -d  copy,$member,$member Prior_Diag.nc $fname


  ncwa -O -a time $fname $fname
  ncwa -O -a metadatalength $fname $fname
  ncwa -O -a NMLnlines $fname $fname
  ncwa -O -a NMLlinelen $fname $fname

  ncrename -d TmpI,lon -d TmpJ,lat -d VelI,slon -d VelJ,slat $fname
  ncrename -v ps,PS -v t,T -v u,U -v v,V -v TmpI,lon -v TmpJ,lat -v VelI,slon -v VelJ,slat $fname
  ncks -C -O -x -v time $fname $fname 
  # ncks -C -O -x -v TmpI $fname $fname
  # ncks -C -O -x -v TmpJ $fname $fname
  # ncks -C -O -x -v VelI $fname $fname
  # ncks -C -O -x -v VelJ $fname $fname
  ncks -C -O -x -v copy $fname $fname
  ncks -C -O -x -v lev $fname $fname
  ncks -C -O -x -v CopyMetaData $fname $fname
  ncks -C -O -x -v inputnml $fname $fname
  ncks -A time.nc $fname

  @ trip ++

end

exit 0

# <next few lines under version control, do not edit>
# $URL: $
# $Revision: $
# $Date: $

