echo '*********************************************'
echo
echo Testing dart restarts. 
echo Assuming you have both wrfHydro_to_dart and dart_to_wrfHydro built. 
echo
## must work from here. 
cd /home/jamesmcc/DART/lanai/models/wrfHydro/work/tests
date
pwd

## load my nco Functions
source ~/ncScripts/ncFilters.sh

## test Noah with and without the hydro component
echo
echo '*********************************************'
echo Testing NOAH configurations
echo Setting Noah restart files.
cp noah.RESTART.2013091102_DOMAIN3 ../restart.nc
cp noah.HYDRO_RST.2013-09-11_02:00_DOMAIN3 ../restart.hydro.nc 

## noah only
echo
echo Passing and retrieving noah restarts without the hydro component
echo Setting the input.nml
cp input.nml.noahOnlyPerfectRestart ../input.nml
echo wrfHydro_to_dart
./wrfHydro_to_dart &>/dev/null
echo dart_to_wrfHydro
./dart_to_wrfHydro &>/dev/null
echo Diffing the noah restart  \(All differences printed immediately below\)
ncVarDiff ../restart.nc noah.RESTART.2013091102_DOMAIN3 

## noah with hydro component
echo 
echo Passing and retrieving noah restarts WITH the hydro component
echo Setting the input.nml
cp input.nml.noahWrfHydroPerfectRestart ../input.nml
echo wrfHydro_to_dart
./wrfHydro_to_dart &>/dev/null
echo dart_to_wrfHydro
./dart_to_wrfHydro &>/dev/null
echo Diffing the noah restart \(All differences printed immediately below\)
ncVarDiff ../restart.nc noah.RESTART.2013091102_DOMAIN3 
echo Diffing the HYDRO restart \(All differences printed immediately below\)
ncVarDiff ../restart.hydro.nc noah.HYDRO_RST.2013-09-11_02:00_DOMAIN3 
## Just for testing that differences actually do get printed to screen
#ncVarDiff ../restart.hydro.nc /d2/weiyu/wrf_release/EXE/Noah_mp_test/HYDRO_RST.2013-09-11_04:00_DOMAIN3 

###########################
## test Noah with and without the hydro component
## set the proper restarts
echo
echo '*********************************************'
echo Testing NoahMP configurations.
echo Setting NoahMP restart files.
cp noahMP.RESTART.2013091200_DOMAIN3 ../restart.nc 
cp noahMP.HYDRO_RST.2013-09-12_00:00_DOMAIN3 ../restart.hydro.nc 

## noahMP only
echo 
echo Passing and retrieving noahMP restarts without the hydro component
echo Setting the input.nml
cp input.nml.noahMPOnlyPerfectRestart ../input.nml
echo wrfHydro_to_dart
./wrfHydro_to_dart &>/dev/null
echo dart_to_wrfHydro
./dart_to_wrfHydro &>/dev/null
echo Diffing the noahMP restart \(All differences printed immediately below\)
ncVarDiff ../restart.nc noahMP.RESTART.2013091200_DOMAIN3 

## noahMP with hydro component
echo 
echo Passing and retrieving noahMP restarts WITH the hydro component
echo Setting the input.nml
cp input.nml.noahMPWrfHydroPerfectRestart ../input.nml
echo wrfHydro_to_dart
./wrfHydro_to_dart &>/dev/null
echo dart_to_wrfHydro
./dart_to_wrfHydro &>/dev/null
echo Diffing the noahMP restart \(All differences printed immediately below\)
ncVarDiff ../restart.nc noahMP.RESTART.2013091200_DOMAIN3 
echo Diffing the HYDRO restart \(All differences printed immediately below\)
ncVarDiff ../restart.hydro.nc noahMP.HYDRO_RST.2013-09-12_00:00_DOMAIN3 

