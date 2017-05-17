; DART software - Copyright 2004 - 2011 UCAR. This open source software is
; provided by UCAR, "as is", without charge, subject to all terms of use at
; http://www.image.ucar.edu/DAReS/DART/DART_download
;***********************************************************************************
; Tarkeshwar Singh
; PhD, IIT Delhi
; Email: tarkphysics87@gmail.com
;************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;************************************************
begin
;************************************************
; file handling
;************************************************
   fn  = "input.nc" ; define filename
   in  = addfile(fn,"r")                         ; open netcdf file
;************************************************
; read needed variables from file
;************************************************
   T      =  in->T(:,0,:,:,:)     
   U      =  in->U(:,0,:,:,:) 
   V      =  in->V(:,0,:,:,:)
   Q      =  in->Q(:,0,:,:,:)
   CLDLIQ =  in->CLDLIQ(:,0,:,:,:)
   PS     =  in->PS(:,0,:,:)

   time   = in->time  
   lat    = in->lat   
   lon    = in->lon  
   slat   = in->slat   
   slon   = in->slon  
   sigs   = in->sigs 
   hyam   = in->apm 
   hybm   = in->bpm

      
;________________
; Due to periodicity in longitude in LMDZ, two values in lon and slon are same
; but 'linint2' function require monotonic increasing/decreasing lon values.
; following do loop identify those same values in lon and slon using where function
; and then decrease one of them by a small value ( e.g. 0.05)

   i = 0
   do while(i.le.dimsizes(lon)-2)

    ; following logic is applied by where function to make monotonic increasing
    ; if lon(i) = lon(i+1) then  lon(i) = lon(i) -0.05 ; else lon(i) = lon(i)

     lon(i)  =  where ( lon(i)  .eq. lon(i+1), lon(i) -0.05, lon(i))
     slon(i) =  where ( slon(i) .eq. slon(i+1), slon(i) -0.05, slon(i))

     i=i+1

   end do
;________________

 
;********************************************************************************************
;Interpolate  U and V at scalar point 

Unew = linint2(slon,lat,U(:,:,:,:),True,lon,lat,0)
Vnew = linint2(lon,slat,V(:,:,:,:),True,lon,lat,0)


;*********************************************************************************************
; Write ouput in .nc file
;*********************************************************************************************

    diro = "./"                     ; Output directory
    filo = "State_vector_ens_mean_ml.nc"             ; Output file
    system("/bin/rm -f " + diro + filo)    ; remove if exists
    fout  = addfile (diro + filo, "c")  ; open output file
    setfileoption(fout,"DefineMode",True)

    dims = dimsizes(T)
    print(dims)

    TEMP      = new((/dims(0),dims(1),dims(2),dims(3)/),float)
    UWIND     = new((/dims(0),dims(1),dims(2),dims(3)/),float)
    VWIND     = new((/dims(0),dims(1),dims(2),dims(3)/),float)
    Qsh       = new((/dims(0),dims(1),dims(2),dims(3)/),float)
    CDLIQ     = new((/dims(0),dims(1),dims(2),dims(3)/),float)
    Psurf     = new((/dims(0),dims(2),dims(3)/),float)

    dimNames = (/"time", "plev","lat", "lon"/)   ; coordinate variables
    dimSizes = (/  -1  ,  dims(1),  dims(2), dims(3)    /)   ; coordinate dimensions
    dimUnlim = (/ True , False, False, False /)   ; define unlimited d

    filedimdef(fout, dimNames  , dimSizes,  dimUnlim )

;put var names you want to define as "varNames3D"
     varNames3D        = (/ "TEMP" /)
     varNames3D        = (/ "UWIND" /)
     varNames3D        = (/ "VWIND" /)
     varNames3D        = (/ "Qsh" /)
     varNames3D        = (/ "CDLIQ" /)
     varTypes3D        = (/ "float" /)

     varNames2D        = (/ "Psurf" /)
     varTypes2D        = (/ "float" /)

     filevardef(fout, varNames3D, varTypes3D, (/"time", "plev","lat" , "lon"/) )
     filevardef(fout, varNames2D, varTypes2D, (/"time", "lat" , "lon"/) )

     ;define variables
     TEMP  = T
     UWIND = Unew
     VWIND = Vnew
     Qsh   = Q
     CDLIQ = CLDLIQ
     Psurf = PS
     plev  = sigs

     TEMP!0 = "time"
     TEMP&time = time
     TEMP!1 = "plev"
     TEMP&plev = plev
     TEMP!2 = "lat"
     TEMP&lat = lat
     TEMP!3 = "lon"
     TEMP&lon = lon

     UWIND!0 = "time"
     UWIND&time = time
     UWIND!1 = "plev"
     UWIND&plev = plev
     UWIND!2 = "lat"
     UWIND&lat = lat
     UWIND!3 = "lon"
     UWIND&lon = lon

     VWIND!0 = "time"
     VWIND&time = time
     VWIND!1 = "plev"
     VWIND&plev = plev
     VWIND!2 = "lat"
     VWIND&lat = lat
     VWIND!3 = "lon"
     VWIND&lon = lon

     Qsh!0 = "time"
     Qsh&time = time
     Qsh!1 = "plev"
     Qsh&plev = plev
     Qsh!2 = "lat"
     Qsh&lat = lat
     Qsh!3 = "lon"
     Qsh&lon = lon

     CDLIQ!0 = "time"
     CDLIQ&time = time
     CDLIQ!1 = "plev"
     CDLIQ&plev = plev
     CDLIQ!2 = "lat"
     CDLIQ&lat = lat
     CDLIQ!3 = "lon"
     CDLIQ&lon = lon

     Psurf!0 = "time"
     Psurf&time = time
     Psurf!1 = "lat"
     Psurf&lat = lat
     Psurf!2 = "lon"
     Psurf&lon = lon
; write metadata
    fout->plev  = plev
    fout->lat   = lat
    fout->lon   = lon
    fout->time  = time
    fout->TEMP  = TEMP
    fout->UWIND = UWIND
    fout->VWIND = VWIND
    fout->Qsh   = Qsh
    fout->CDLIQ   = CDLIQ
    fout->Psurf   = Psurf
end
