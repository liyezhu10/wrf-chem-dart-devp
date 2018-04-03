#!/bin/csh -ef

setenv OPENGGCMDIR $HOME/openggcm-2016-04-18
#set baserundir=/home/space/dcramer/lustre/dart/cir07_test
set baserundir=/mnt/lustre/lus0/home/dart/dart/rma_openggcm/models/openggcm/test/baserun
set runprefix="target"
set ensemble_size=3

# ----------------------------------------------------------------------
# create_runs.csh
#
# This script will create multiple runs based on an existing run. The
#   runs will be created in the current directory.
#
# Inputs:
#   baserundir :: path to existing run
#   ensemble_size :: number of "copies"
#   runprefix :: what to name "copies" ("target_", etc.)
#   (targetfile) :: special targetfile that sets PPN and doesn't submit
#   runme (optional) :: for different parameter values (BE CAREFUL!)
#   swdata (optional) :: to replace base run's version with different values
# 
# Created directory structure (below current directory):
#   ./build_all.sh (recreates and builds all runs)
#   ./submit_all.sh (submits all runs)
#   ./recreate.csh (rebuilds in.$RUN, openggcm.f, ctim.f from runme, swdata)
#   ./target_001 (runme[opt], swdata)
#   ./target_001/target (qsub, sh)
#   ...
#   ./target_NNN (")
#   ./target_NNN/target (")
#
# To do/problems:
#   - in.$RUN rebuilt even if swdata, runme not changed
#   - RUN and SWFILE parameters are updated by adding to the top and
#     bottom of tmp.csh, tmp.h, and tmp.i files - should change existing
#     value instead
# ----------------------------------------------------------------------

alias die 'echo FATAL: \!*; exit 1'

if ( ! -d $baserundir ) then
    die "base run directory $baserundir not found!"
endif

set BASERUN=`basename $baserundir`

# ----------------------------------------------------------------------
# Read all the variables

if ( ! -r $baserundir/tmp.csh ) then
    die "tmp.csh not found"
endif
source $baserundir/tmp.csh

# ----------------------------------------------------------------------
# Load per-target defaults and use them if they were not overridden 
# previously

set targetfile=$GGCMDIR/conf/target/$TARGET
if ( ! -r $targetfile ) then
    set targetfile=$HOME/.openggcm/conf/target/$TARGET
endif

if ( ! -r $targetfile ) then
    die "Cannot read $GGCMDIR/conf/target/$TARGET nor $targetfile."
endif

rm -f tmp.submit
source $targetfile

# ----------------------------------------------------------------------

# FIXME: Use proper MPI autoconf macros?

set configure_args
set configure_args=( $configure_args CC='"'$MPICC'"' )
if ( $?MPICXX ) then 
  set configure_args=( $configure_args CXX='"'$MPICXX'"' )
endif
set configure_args=( $configure_args F77='"'$MPIF90'"' )
set configure_args=( $configure_args FC='"'$MPIF90'"' )
if ( "$MPICC_FLAGS" != "none" ) then
    set configure_args=( $configure_args CFLAGS='"'$MPICC_FLAGS'"' )
endif    
if ( $?MPICXX_FLAGS ) then 
  if ( "$MPICXX_FLAGS" != "none" ) then
      set configure_args=( $configure_args CXXFLAGS='"'$MPICXX_FLAGS'"' )
  endif    
endif    
if ( "$MPIF77_FLAGS" != "none" ) then
    set configure_args=( $configure_args FFLAGS='"'$MPIF77_FLAGS'"' )
endif    
if ( "$MPIF90_FLAGS" != "none" ) then
    set configure_args=( $configure_args FCFLAGS='"'$MPIF90_FLAGS'"' )
endif

if ( "$OUTPUTMODE" == "hdf5" || "$OUTPUTMODE" == "xdmf" || "$OUTPUTMODE" == "xdmf_serial" ) then
    set configure_args=( $configure_args --with-hdf5='"'$HDF5_DIR'"' )
endif

if ( "$OUTPUTMODE" == "cdf" ) then
    set configure_args=( $configure_args --with-cdf='"'$CDF_DIR'"' )
endif

set configure_args=( $configure_args $CONF_ARGS )

# ----------------------------------------------------------------------
# create each run
@ NODES = ( $TOTNODES + $PPN - 1 ) / $PPN
@ CORES = $NODES * $PPN

set INST_NUM=0
while ($INST_NUM < $ensemble_size)
  @ INST_NUM++
  set INST=`printf "%03d" $INST_NUM`
  set RUN=$runprefix"_"$INST
  echo $RUN
  set RUNDIR=$PWD/$RUN

  # ----------------------------------------------------------------------
  # Put all the files needed on the target into the "target" subdir

  if ( -e $RUNDIR ) then
      rm -rf $RUNDIR
  endif

  mkdir $RUNDIR
  mkdir $RUNDIR/target
  mkdir $RUNDIR/target/build
  mkdir $RUNDIR/target/$RUN     #workaround for Jimmy's special dart code

  #get input parameters
  cp $baserundir/tmp.csh $RUNDIR
  cp $baserundir/tmp.h $RUNDIR
  cp $baserundir/tmp.i $RUNDIR
  cp $baserundir/BASETIME $RUNDIR

  #solar wind, in.$RUN data
  cp $baserundir/swdata $RUNDIR
  cp $baserundir/in.$BASERUN $RUNDIR/in.$RUN
  ln -s $RUNDIR/in.$RUN $RUNDIR/target
  #switch swfile to swdata, if necessary (in.$RUN must be recreated)
  if ( ( $SWFILE == "auto" ) || ( $SWFILE == "minvar" ) ) then
    sed -i '1s/^/setenv swfile "./swdata"\n/' tmp.csh
    sed -i '1s/^/setenv SWFILE "./swdata"\n/' tmp.csh
    echo 'setenv swfile "./swdata"' >> $RUNDIR/tmp.csh #use swdata file from base run for openggcm.f creation
    echo 'setenv SWFILE "./swdata"' >> $RUNDIR/tmp.csh
  endif

  #switch run name from base run name
  #  (for recreation of openggcm.f, ctim.f, $RUN.f, in.$RUN files)
  sed -i '1s/^/setenv run "'$RUN'"\n/' $RUNDIR/tmp.csh
  sed -i '1s/^/setenv RUN "'$RUN'"\n/' $RUNDIR/tmp.csh
  sed -i "1s/^/#run=$RUN\n/" $RUNDIR/tmp.h
  sed -i "1s/^/#RUN=$RUN\n/" $RUNDIR/tmp.h
  sed -i "1s/^/run=    $RUN\n/" $RUNDIR/tmp.i
  echo 'setenv run "'$RUN'"' >> $RUNDIR/tmp.csh
  echo 'setenv RUN "'$RUN'"' >> $RUNDIR/tmp.csh
  echo "#run=$RUN" >> $RUNDIR/tmp.h
  echo "#RUN=$RUN" >> $RUNDIR/tmp.h
  echo "run=    $RUN" >> $RUNDIR/tmp.i

  #recreate in.$RUN file in include correct run name (and swfile)
  if ( $DOINP ) then
    cd $RUNDIR
    setenv INPF $baserundir/runme
    $GGCMDIR/bin/input.csh 0
    if ( $status != 0 ) then
        die "input.csh failed."
    endif
    cd -
  endif

  #copy openggcm.f fortran template file
  cp $baserundir/tmp.$BASERUN.1.pre $RUNDIR/tmp.$RUN.1.pre

  #grid info
  cp $baserundir/$BASERUN.smf $RUNDIR/$RUN.smf
  ln -s $RUNDIR/$RUN.smf $RUNDIR/target

  #satellite orbits
  if ( -f $ORBITFILE ) then
      ln -s $baserundir/$ORBITFILE $RUNDIR
      ln -s $RUNDIR/$ORBITFILE $RUNDIR/target
  endif

  #ctim init data
  ln -s $baserundir/target/ctim-indata $RUNDIR/target

  #rcm input files
  if ( "$RCMCODE" == "rice" ) then
      ln -s $baserundir/target/enchan.dat $RUNDIR/target
      ln -s $baserundir/target/dktable $RUNDIR/target
  endif

  #source and build files from base run
  ln -s $baserundir/target/target-build $RUNDIR/target
  ln -s $baserundir/target/build/* $baserundir/target/build/.deps $RUNDIR/target/build
  rm $RUNDIR/target/build/openggcm* $RUNDIR/target/build/ctim*

  # ----------------------------------------------------------------------
  # Create a script to build and submit the code on the target

  if ( "$TARGETPREP" == "none" ) then
      set TARGETPREP=""
  endif
  
cat > $RUNDIR/target/$RUN.sh <<EOF
#! /bin/bash -l
set -e
export RUN="$RUN"
export TOTNODES="$TOTNODES"
export RUNTIME="$RUNTIME"
$TARGETPREP
cd $RUNDIR
../recreate.csh
cd target
#build code.  If it fails to build, bail out right away.
cd build 
make || exit 1
cd -
ln -sf build/openggcm .
set +e
EOF
  rm -f tmp.submit
  source $targetfile
  if ( -r tmp.submit ) then
      cat tmp.submit >> $RUNDIR/target/$RUN.sh
      echo "cd "'$PWD'" && qsub $RUN.qsub" >> $RUNDIR/target/$RUN.sh  #re-add submit directive
  endif
  
  chmod 755 $RUNDIR/target/$RUN.sh

set PPNSTR = :ppn=$PPN  # avoid "Bad modifier" error
cat > $RUNDIR/target/$RUN.qsub <<XEOF
#! /bin/bash
#PBS -l nodes=$NODES$PPNSTR
#PBS -l walltime=$RUNTIME

cd $RUNDIR/target

export ATP_ENABLED=1
export GFORTRAN_UNBUFFERED_ALL=1
aprun -n $TOTNODES ./$RUN.exe &> $RUN.log

XEOF

end  #instance loop

  # ----------------------------------------------------------------------
  # create script to recreate in.$RUN, openggcm.f, ctim.f, $RUN.f to:
  #   - change run name in code
  #   - handle new runme (updated parameters)
  #   - include SWFILE that may have changes from the base run

cat > recreate.csh <<EOF
#!/bin/csh -ef

source tmp.csh

# recreate tmp.csh, tmp.h, tmp.i if new runme provided
if ( -r runme ) then
  cp tmp.csh tmp.csh.old
  cp tmp.h tmp.h.old
  cp tmp.i tmp.i.old
  \$GGCMDIR/bin/script.includes \$GGCMDIR runme $baserundir/target/include/input.defines 0 #GGCMDIR isn't used in script.includes, thankfully
  if ( \$status != 0 ) then
      die "script.includes failed."
  endif
  sed -i '1s/^/setenv swfile "./swdata"\n/' tmp.csh
  sed -i '1s/^/setenv SWFILE "./swdata"\n/' tmp.csh
  sed -i '1s/^/setenv run "'\$RUN'"\n/' tmp.csh
  sed -i '1s/^/setenv RUN "'\$RUN'"\n/' tmp.csh
  sed -i "1s/^/#run=\$RUN\n/" tmp.h
  sed -i "1s/^/#RUN=\$RUN\n/" tmp.h
  sed -i "1s/^/run=    \$RUN\n/" tmp.i
  echo 'setenv swfile "./swdata"' >> tmp.csh #use swdata file
  echo 'setenv SWFILE "./swdata"' >> tmp.csh
  echo 'setenv run "'\$RUN'"' >> tmp.csh
  echo 'setenv RUN "'\$RUN'"' >> tmp.csh
  echo "#run=\$RUN" >> tmp.h #switch run name from base run name for openggcm.f creation
  echo "#RUN=\$RUN" >> tmp.h
  echo "run=    \$RUN" >> tmp.i
  echo "$STARTTIME" >! BASETIME
endif

# recreate in.\$RUN due to possible change in swdata
if ( \$DOINP ) then
  cp in.\$RUN in.\$RUN.old
  setenv INPF runme
  \$GGCMDIR/bin/input.csh 0
  if ( \$status != 0 ) then
      die "input.csh failed."
  endif
endif

# recreate openggcm.f
\$GGCMDIR/bin/script.precomp  \$RUN \$PRECISION 0
if ( \$status != 0 ) then
   die "Pre-compilation code generation failed."
endif
cat \$RUN.f > openggcm.f

# recreate ctim.f
cat tmp.h $baserundir/target/src/ctim-core.for >! tmp.pre
set FPPN="fppn -q"
\$FPPN -c -- tmp.pre >! ctim.f
cat ctim.f >> \$RUN.f

ln -sf ../../openggcm.f target/build
ln -sf ../../ctim.f target/build

if( ! -f \$RUN.f ) then
    die "No \$RUN.f!"
endif
EOF
  chmod 755 recreate.csh

# ----------------------------------------------------------------------
# create build_all and submit_all scripts

cat > build_all.sh <<EOF
#! /bin/bash -l
set -e
export RUN="$RUN"
export TOTNODES="$TOTNODES"
export RUNTIME="$RUNTIME"
$TARGETPREP
EOF

rm -f submit_all.sh

set INST_NUM=0
while ($INST_NUM < $ensemble_size)
  @ INST_NUM++
  set INST=`printf "%03d" $INST_NUM`
  set RUN=$runprefix"_"$INST
  set RUNDIR=$PWD/$RUN

cat >> build_all.sh <<EOF
cd $RUNDIR
../recreate.csh
cd target 
#build code.  If it fails to build, bail out right away.
cd build 
make || exit 1
cd -
ln -sf build/openggcm $RUN.exe
cd ../..
EOF

cat >> submit_all.sh <<EOF
echo === Submitting job $RUN.qsub
cd $RUNDIR/target && qsub $RUN.qsub
EOF

end

cat >> submit_all.sh <<EOF
cd $PWD
qsub ../shell_scripts/dart_test.qsub
EOF

chmod 755 build_all.sh
chmod 755 submit_all.sh

