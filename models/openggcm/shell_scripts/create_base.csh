#!/bin/csh

setenv OPENGGCMDIR $HOME/openggcm-2016-04-18
set baserundir=$PWD/baserun

# ----------------------------------------------------------------------
# create_base.csh
#
# This script creates an openggcm run and copies necessary files to 
#   create similar runs with identical source code.
#
# Inputs:
#   baserundir :: location to create run
#   runme :: normal runme file
#   solar wind input files ::
#     - SWFILE :: usual solar wind input file
#     - (SWFILE=auto,minvar) SWMON.xxx files 
#   TARGET :: edited target file (in conf/target) that 
#               - prevents submission 
#               - sets PPN
#   ORBITFILE :: satellite orbit data file (optional)
#
# To do/problems:
#   - what causes all the STOP messages?
# ----------------------------------------------------------------------

alias die 'echo FATAL: \!*; exit 1'

if ( ! -r runme ) then
  die "runme missing"
endif

if ( -d $baserundir ) then
  rm -rf $baserundir
  echo "previous base run removed"
endif
mkdir $baserundir

#copy runme
cp runme $baserundir

#get solar wind, orbit data file info from runme
$OPENGGCMDIR/bin/script.includes $OPENGGCMDIR runme $OPENGGCMDIR/include/input.defines 0
if ( $status != 0 ) then
  die "script.includes failed."
endif
if ( ! -r tmp.csh ) then
  die "tmp.csh not found"
endif
source tmp.csh

#copy solar wind, orbit data files
if ( ( $SWFILE == "auto" ) || ( $SWFILE == "minvar" ) ) then
  # if relative path, copy to base run dir
  if ("$IDATA" !~ /*) then
    if ( ! -d $baserundir/$IDATA ) then
      mkdir $baserundir/$IDATA
    endif
    cp $IDATA/$SWMON.* $baserundir/$IDATA
  endif
else
  cp $SWFILE $baserundir
endif
if ( $ORBITFILE != "none" ) then
  cp $ORBITFILE $baserundir
endif

# build base run
cd $baserundir
./runme  #use special target file to prevent submission and set PPN

#re-get runme params (necessary?)
if ( ! -r $baserundir/tmp.csh ) then
  die "tmp.csh not found"
endif
source $baserundir/tmp.csh

#if not a swdata file run (auto or minvar), create generated swdata file
if ( ! -r $baserundir/swdata ) then
expand $baserundir/tmp.swdata1 |  awk '{ \
    printf "%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f ", \
    $1,$2,$3,$4,$5,$6,$7,$8,$9;  \
    printf "%8.5f %8.5f %8.5f\n",$11,$12,$13; }'  >! $baserundir/swdata
endif

#save files needed for "copied" runs
set GGCM_SOURCES = ( ggcm-cliches.for \
    ggcm-util.for \
    ggcm-vec.for \
    mhd-cliches.for \
    ggcm.for \
    new_cotr.for \
    mhd-modules.for \
    mhd-scon.for \
    mhd-corea.for \
    mhd-coreb.for \
    mhd-bnd.for \
    mhd-cotr.for \
    glob-util.for \
    mhd-machspec.for \
    mhd-diag.for \
    mhd-dipole.for \
    mhd-iono-map.for \
    mhd-iono.for \
    mhd-fltrace.for \
    mhd-misig.for \
    io-psol.for \
    io-post.for \
    mhd-commu1.for \
    mhd-bndsw.for \
    mhd-ini.for \
    mhd-rcm.for \
    mhd-diag2.for )

if ( -r $GGCMDIR/src/mhd-debug.for ) then
    set GGCM_SOURCES = ( $GGCM_SOURCES mhd-debug.for )
endif

if ( "$SATOUT" == "true" ) then
    set GGCM_SOURCES = ( $GGCM_SOURCES mhd-satout.for )
endif

#save openggcm.f, ctim.f source files
mkdir $baserundir/target/src
foreach src ( $GGCM_SOURCES )
    cp -p $GGCMDIR/src/$src $baserundir/target/src
end
cp -p $GGCMDIR/src/ctim-core.for $baserundir/target/src

#save input parameter definitions
mkdir $baserundir/target/include
cp -p $GGCMDIR/include/input.defines $baserundir/target/include

#save current target-build, ctim-indata, dktable, enchan.dat
mv $baserundir/target/target-build $baserundir/target/target-build-link
cp -prL $baserundir/target/target-build-link $baserundir/target/target-build
mv $baserundir/target/ctim-indata $baserundir/target/ctim-indata-link
cp -prL $baserundir/target/ctim-indata-link $baserundir/target/ctim-indata
mv $baserundir/target/dktable $baserundir/target/dktable-link
cp -prL $baserundir/target/dktable-link $baserundir/target/dktable
mv $baserundir/target/enchan.dat $baserundir/target/enchan.dat-link
cp -prL $baserundir/target/enchan.dat-link $baserundir/target/enchan.dat

#create openggcm.f template that points to saved source files
set HOST=`hostname`
set UNAM=`uname -a`
setenv DATE `date`
cat >! tmp.pre <<eeeeee
c++++++++++++++++++++++++++++++++++++++
c+++++++>>> creation date: $DATE
c+++++++>>> creation script: $0
c+++++++>>> creation hostname: $HOST
c+++++++>>> creation uname: $UNAM
c++++++++++++++++++++++++++++++++++++++
c++ the following sources were used: ++
c++++++++++++++++++++++++++++++++++++++
eeeeee
foreach src ( $GGCM_SOURCES )
    echo "c+++++++>>> $baserundir/target/src/$src" >> tmp.pre
end
cat >> tmp.pre <<eeeeee
c++++++++++++++++++++++++++++++++++++++
.include tmp.h
.include $baserundir/grid_include.$RUN
eeeeee
foreach src ( $GGCM_SOURCES )
    echo ".include $baserundir/target/src/$src" >> tmp.pre
end
/bin/mv tmp.$RUN.1.pre tmp.$RUN.1.pre.save #save old version
/bin/mv tmp.pre tmp.$RUN.1.pre

