#!/bin/csh

# Purpose: 
#  Apply forcing adjustments, and parameters adjustments (etc) supplied in 
#  restart.assimOnly.nc restart files to the model run in current directory. 

# Assumptions: 
#  You are in a directory with 
#  1) one restart.assimOnly.nc
#  2) all the other information necessary to advance the model are properly specified in 
#     the namelist.hrldas and hydro.namelist which also reside in this dir. A summary of 
#     file information used:
#     namelist.hrldas
#         KHOUR: the number of forcing timesteps to which forcing adjustments are applied. 
#         INDIR: where the forcing files are found
#     hydro.namelist
#         GEO_FINEGRID_FLNM: where certain model parameters are set

# Notes: 
#  1) Previously, when the model was to be advanced for long periods, I would 
#     consider the specified precipMult as a mean and I would generate a timeseries
#     about this mean with some variance. It's not clear that's appropriate. So it
#     been abandoned. It seems like it could simply be added in at a given timestep
#     by introducing zero biased noise (rather than generating a full timeseries 
#     upfront).
#  2) Following on 1, eliminating the addtional perturbation eliminates any need to 
#     keep perturbed timeseries because they can be recreated from restart.assimOnly.nc
#     and the forcing. If the additional perturbation is desired, then the numbers/noise
#     should be just written to a separate file. 

#===============================================================================
## Setup
set restartTime = `ncdump -h restart.hydro.nc | grep Restart_ | cut -d'=' -f2 | tr -d ' ";'`
set startYyyy = `echo $restartTime | cut -d- -f1`
set startMm   = `echo $restartTime | cut -d- -f2`
set startDd   = `echo $restartTime | cut -d- -f3 | cut -d_ -f1`
set startHh   = `echo $restartTime | cut -d_ -f2 | cut -d: -f1`

set assimOnlyVars = `ncks -m restart.assimOnly.nc | grep ': type' | cut -f1 -d':'`

#===============================================================================
# Forcings
# Edit the path to the forcing dir in the namelist.hrldas to use perturbed forcing
# in the current dir

set forcingVars = ( precipMult otherVarsHere )
set forcingGrep = `echo $forcingVars | tr ' ' '|'`
set forcingGrep = `echo "($forcingGrep)"`
set forcingActive =  `echo $assimOnlyVars | egrep $forcingGrep | wc -l`

# Original forcing location . Note that there no required name for this.
# We need this even if forcing is not active because the path is one dir deeper
set inDirIn = `grep INDIR namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
set inDirIn = `echo $inDirIn | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
set inDirIn = `echo $inDirIn | sed -e 's#"# #g'`
set inDirIn = `echo $inDirIn[$#inDirIn]`

if ( $forcingActive ) then 

    set forcDir = FORCING.assimOnly
    ## break it down to build it back up (make sure it's clean)
    ## cant delete it in this routine since the model is called after exit.
    \rm -rf $forcDir
    \mkdir $forcDir

#make the namelist aware of the new forcings
ex namelist.hrldas <<ex_end
g;INDIR ;s;= .*;= "${forcDir}/";
wq
ex_end

    # How many forcing times need adjusted?
    set khourIn = `grep KHOUR namelist.hrldas | tr -d ' '| egrep -v '^[!\]'`
    set khourIn = `echo $khourIn | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
    set khourIn = `echo $khourIn | sed -e 's#"# #g'`
    set khourIn = `echo $khourIn[$#khourIn]`

    set endYyyy = `date -ud "UTC $startYyyy-$startMm-$startDd $startHh hour + $khourIn hours" +%Y`

    ## This is not easy to do/optimize when there are lots of files. This is my best shot.
    ## It is relatively fast. 
    set yearSeq = `seq -s, $startYyyy $endYyyy | tr ',' '|'`
    ## the annchor (^) supposedly gives a speed up
    set forcFilesSuper = \
        `find ../$inDirIn/ -regextype posix-extended \
            -regex "^../$inDirIn.*/(${yearSeq}).*LDASIN.*" | sort`
    set forcFilesPosn = `echo $forcFilesSuper | tr ' ' '\n' | \
                          grep -n "${startYyyy}${startMm}${startDd}${startHh}" | cut -d':' -f1`
    set forcFilesIn = `echo $forcFilesSuper | tr ' ' '\n' | tail -n+$forcFilesPosn | head -$khourIn`
    \cp $forcFilesIn $forcDir/.
    set forcFiles = `ls $forcDir/*`

    #-----------------------------------------------------------
    # precipMult
    # See if it is among the assimOnlyVars
    if ( `echo $assimOnlyVars | grep precipMult | wc -l` ) then 
        set precipMultDims = `ncks -m restart.assimOnly.nc | grep 'precipMult dimension' | wc -l`

        if ( $precipMultDims == 1 ) then
            set thePrecipMult = `ncdump -v precipMult restart.assimOnly.nc | tail -n2 | head -1 | \
                                    cut -d'=' -f2 | tr -d ' ;'`
            foreach iForc ( $forcFiles )
                ncap2 -O -s "RAINRATE=RAINRATE*${thePrecipMult}" $iForc $iForc
            end
        endif 

        if ( $precipMultDims == 2 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # precipMult

    #-----------------------------------------------------------
    # more vars... to come.

else  ## forcingActive else not
    ## we are one dir level deeper than we should be for the ensembles
    cd ../   
    ln -sf $inDirIn .
    cd -
endif

#===============================================================================
# Parameters
set paramVars = ( OVROUGHRTFAC RETDEPRTFAC gwCoeff gwExpon ksatMult slope maxSmcMult )
set paramGrep = `echo $paramVars | tr ' ' '|'`
set paramGrep = `echo "($paramGrep)"`
set paramActive =  `echo $assimOnlyVars | egrep $paramGrep | wc -l`

if ( $paramActive ) then 

    #---------------------------------------------------------------------
    #---------------------------------------------------------------------
    # GEO_FINEGRID: params which live in this file need set there 
    # for each ensemble member individually. copy to local and change the 
    # location of the file in use.
    set geoFineVars = ( OVROUGHRTFAC RETDEPRTFAC )
    set geoFineGrep = `echo $geoFineVars | tr ' ' '|'`
    set geoFineGrep = `echo "($geoFineGrep)"`
    set geoFineActive =  `echo $assimOnlyVars | egrep $geoFineGrep | wc -l`

    if ($geoFineActive) then 

        # original file
        set geoFineFileIn = `grep GEO_FINEGRID_FLNM hydro.namelist | tr -d ' '| egrep -v '^[!\]'`
        set geoFineFileIn = `echo $geoFineFileIn | cut -d'!' -f1  | sed -e "s#[=,']# #g"`
        set geoFineFileIn = `echo $geoFineFileIn | sed -e 's#"# #g'`
        set geoFineFileIn = `echo $geoFineFileIn[$#geoFineFileIn]`
        # ensemble member specific file
        set geoFineFile   = `echo $geoFineFileIn | tr '/' ' '`
        set geoFineFile   = `echo $geoFineFile[$#geoFineFile]`
        ## copy the file  locally 
        \cp $geoFineFileIn $geoFineFile
        ## change location
ex hydro.namelist <<ex_end
g;GEO_FINEGRID_FLNM;s;=.*;= "$geoFineFile";
wq
ex_end

        # OVROUGHRTFAC -------
        # See if it is among the assimOnlyVars
        if ( `echo $assimOnlyVars | grep OVROUGHRTFAC | wc -l` ) then 

            set ovRoughDims = `ncks -m restart.assimOnly.nc | grep 'OVROUGHRTFAC dimension' | wc -l`

            if ( $ovRoughDims == 1 ) then
                set theOvRough = `ncdump -v OVROUGHRTFAC restart.assimOnly.nc | tail -n2 | head -1 | \
                                cut -d'=' -f2 | tr -d ' ;'`
                ncap2 -O -s "OVROUGHRTFAC=OVROUGHRTFAC*0+${theOvRough}f" $geoFineFile $geoFineFile
            endif 
                
            if ( $ovRoughDims == 2 ) then
                echo 'Not yet configured'
                exit 3
            endif 
        endif  # OVROUGHRTFAC
   

        # RETDEPRTFAC --------
        # See if it is among the assimOnlyVars
        if ( `echo $assimOnlyVars | grep RETDEPRTFAC | wc -l` ) then 

            set retDepDims = `ncks -m restart.assimOnly.nc | grep 'RETDEPRTFAC dimension' | wc -l`

            if ( $retDepDims == 1 ) then
                set theRetDep = `ncdump -v RETDEPRTFAC restart.assimOnly.nc | tail -n2 | head -1 | \
                                cut -d'=' -f2 | tr -d ' ;'`
                ncap2 -O -s "RETDEPRTFAC=RETDEPRTFAC*0+${theRetDep}f" $geoFineFile $geoFineFile
            endif 
                
            if ( $retDepDims == 2 ) then
                echo 'Not yet configured'
                exit 3
            endif 
        endif  # RETDEPRTFAC

    endif # geoFineActive
 


    # PARAMETER FILES - *.TBL
    # If running filter, these WERE ALREADY copied locally by advance_model.csh if the assimOnly is active.
    ## OTherwise, it dosent seem to matter as adjustments are to a tmp file which is *moved* to the original.
    # NOTE: the files are format sensitive.

    #--------------------------------------------    
    #--------------------------------------------    
    ## GENPARM.TBL
    
    ## SLOPE -----------------
    if ( `echo $assimOnlyVars | grep slope | wc -l` ) then 
        set slopeDims = `ncks -m restart.assimOnly.nc | grep 'slope dimension' | wc -l`
        if ( $slopeDims == 1 ) then
            set theSlope = `ncdump -v slope restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            set theSlope = `echo $theSlope | cut -c1-5`
            set lineNumSlope = `grep -n SLOPE_DATA GENPARM.TBL | cut -d: -f1`
            @ setLineNum = $lineNumSlope + 2
            sed -i "${setLineNum}s/.*/$theSlope/" GENPARM.TBL
        endif 
        if ( $slopeDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif


    #--------------------------------------------    
    #--------------------------------------------    
    ## GWBUCKPARM.TBL

    ## gwCoeff ---------------
    if ( `echo $assimOnlyVars | grep gwCoeff | wc -l` ) then 
        set gwCoeffDims = `ncks -m restart.assimOnly.nc | grep 'gwCoeff dimension' | wc -l`
        if ( $gwCoeffDims == 1 ) then
            set theGwCoeff = `ncdump -v gwCoeff restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            set theGwCoeff = `echo $theGwCoeff | cut -c1-6`
            cat GWBUCKPARM.TBL | \
                awk -v theGwCoeff=$theGwCoeff 'BEGIN{FS=",";OFS=","};{$2=theGwCoeff;print}' \
                > GWBUCKPARM.TBL.NEW
            mv GWBUCKPARM.TBL.NEW GWBUCKPARM.TBL
        endif 
        if ( $gwCoeffDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 

    ## gwExpon ---------------
    if ( `echo $assimOnlyVars | grep gwExpon | wc -l` ) then 
        set gwExponDims = `ncks -m restart.assimOnly.nc | grep 'gwExpon dimension' | wc -l`
        if ( $gwExponDims == 1 ) then
            set theGwExpon = `ncdump -v gwExpon restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            set theGwExpon = `echo $theGwExpon | cut -c1-6`
            \rm -rf GWBUCKPARM.TBL.NEW
            touch GWBUCKPARM.TBL.NEW
            sed -n 1p GWBUCKPARM.TBL > GWBUCKPARM.TBL.NEW
            sed -n 2p GWBUCKPARM.TBL | \
                awk -v theGwExpon=$theGwExpon 'BEGIN{FS=",";OFS=","};{$3=theGwExpon;print}' \
                >> GWBUCKPARM.TBL.NEW
            mv GWBUCKPARM.TBL.NEW GWBUCKPARM.TBL
        endif 
        if ( $gwExponDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 

    #--------------------------------------------    
    #--------------------------------------------    
    ## SOILPARM.TBL and HYDRO.TBL
    ## Oh joy, these are in two places!

    ## ksatMult --------------
    if ( `echo $assimOnlyVars | grep ksatMult | wc -l` ) then 

        set ksatMultDims = `ncks -m restart.assimOnly.nc | grep 'ksatMult dimension' | wc -l`
        if ( $ksatMultDims == 1 ) then
            set theKsatMult = `ncdump -v ksatMult restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## SOILPARM.TBL
            set nlines = `cat SOILPARM.TBL | wc -l`
            \rm -rf SOILPARM.TBL.NEW
            touch SOILPARM.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3|23|24|25|45)$' | wc -l` ) then
                    sed -n ${ll}p SOILPARM.TBL >> SOILPARM.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f8`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theKsatMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p SOILPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$8="  "newValue;print}' \
                    >> SOILPARM.TBL.NEW 
            end

            mv SOILPARM.TBL.NEW   SOILPARM.TBL

            ## HYDRO.TBL
            set nlines = `cat HYDRO.TBL | wc -l`
            \rm -rf HYDRO.TBL.NEW
            touch HYDRO.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( $ll <= 31 ) then
                    sed -n ${ll}p HYDRO.TBL >> HYDRO.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p HYDRO.TBL | cut -d',' -f1`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theKsatMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p HYDRO.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$1=newValue;print}' \
                    >> HYDRO.TBL.NEW 
            end

            mv HYDRO.TBL.NEW   HYDRO.TBL

        endif ## 1D

        if ( $ksatMultDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 

    ## maxSmcMult ------------
    if ( `echo $assimOnlyVars | grep maxSmcMult | wc -l` ) then 

        set maxSmcMultDims = `ncks -m restart.assimOnly.nc | grep 'maxSmcMult dimension' | wc -l`
        if ( $maxSmcMultDims == 1 ) then
            set theMaxSmcMult = `ncdump -v maxSmcMult restart.assimOnly.nc | tail -n2 | head -1 | \
                              cut -d'=' -f2 | tr -d ' ;'`
            ## SOILPARM.TBL
            set nlines = `cat SOILPARM.TBL | wc -l`
            \rm -rf SOILPARM.TBL.NEW
            touch SOILPARM.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( `echo $ll | egrep '^(1|2|3|23|24|25|45)$' | wc -l` ) then
                    sed -n ${ll}p SOILPARM.TBL >> SOILPARM.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p SOILPARM.TBL | cut -d',' -f5`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theMaxSmcMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p SOILPARM.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$5="  "newValue;print}' \
                    >> SOILPARM.TBL.NEW 
            end

            mv SOILPARM.TBL.NEW   SOILPARM.TBL

            ## HYDRO.TBL
            set nlines = `cat HYDRO.TBL | wc -l`
            \rm -rf HYDRO.TBL.NEW
            touch HYDRO.TBL.NEW
            foreach ll (`seq 1 $nlines`)
                if ( $ll <= 31 ) then
                    sed -n ${ll}p HYDRO.TBL >> HYDRO.TBL.NEW || exit 7
                    continue
                endif 
                set origValue = `sed -n ${ll}p HYDRO.TBL | cut -d',' -f2`
                set origValue = `printf '%.16f' $origValue`
                set newValue = `echo "$origValue * $theMaxSmcMult" | bc`
                set newValue = `printf '%.2e' $newValue`
                set newValue = `echo $newValue | sed 's/-0/-/' | sed 's/+0/+/'`
                sed -n ${ll}p HYDRO.TBL | \
                    awk -v newValue=$newValue 'BEGIN{FS=",";OFS=","};{$2=newValue;print}' \
                    >> HYDRO.TBL.NEW 
            end

            mv HYDRO.TBL.NEW   HYDRO.TBL

        endif ## 1D

        if ( $maxSmcMultDims > 1 ) then
            echo 'Not yet configured'
            exit 3
        endif 
    endif  # 

endif  # paramActive


exit 0
