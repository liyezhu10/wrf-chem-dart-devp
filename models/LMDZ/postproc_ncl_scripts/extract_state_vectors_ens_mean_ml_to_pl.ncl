; DART software - Copyright UCAR. This open source software is provided
; by UCAR, "as is", without charge, subject to all terms of use at
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

   T      =  in->T
   U      =  in->U
   V      =  in->V
   Q      =  in->Q
   CLDLIQ =  in->CLDLIQ
   PS     =  in->PS
   PSmean =  in->PS(:,0,:,:)

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

Unew = linint2(slon,lat,U(:,0,:,:,:),True,lon,lat,0)
Vnew = linint2(lon,slat,V(:,0,:,:,:),True,lon,lat,0)

;*********************************************************************************************
; Model Hyprid level to pressure level conversion
;*********************************************************************************************
   P0mb =1000.
   hyam = hyam /100000   ; Divide by 100000  to mak similar coefficient like CAM

   pnew = (/ 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50,30, 20, 10/)
   ;pnew =(/10,20,30,50,70,100,150,200,250,300,400,500,600,700,850,925,1000/)
; type of interpolation: 1 = linear, 2 = log, 3 = loglog
   interp = 2
; is extrapolation desired if data is outside the range of PS
   extrap = False
; create an array of desired pressure levels:

; calculate ens mean var  (copy = 0) on pressure levels
; note, the 7th argument is not used, and so is set to 1.

   TonP      = vinth2p( T(:,0,:,:,:),hyam,hybm,pnew,PS(:,0,:,:),interp,P0mb,1,extrap)
   UonP      = vinth2p(Unew(:,:,:,:),hyam,hybm,pnew,PS(:,0,:,:),interp,P0mb,1,extrap)
   VonP      = vinth2p(Vnew(:,:,:,:),hyam,hybm,pnew,PS(:,0,:,:),interp,P0mb,1,extrap)
   QonP      = vinth2p( Q(:,0,:,:,:),hyam,hybm,pnew,PS(:,0,:,:),interp,P0mb,1,extrap)
   CLDLIQonP = vinth2p( CLDLIQ(:,0,:,:,:),hyam,hybm,pnew,PS(:,0,:,:),interp,P0mb,1,extrap)

;*********************************************************************************************
; Write ouput in .nc file
;*********************************************************************************************

    diro = "./"                     ; Output directory
    filo = "State_vector_ens_mean_pl.nc"             ; Output file
    system("/bin/rm -f " + diro + filo)    ; remove if exists
    fout  = addfile (diro + filo, "c")  ; open output file
    setfileoption(fout,"DefineMode",True)

    dims = dimsizes(TonP)
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


  ; coordinate var lon lat and time are already defined in input netcdf file
  ; but plev is new and need to define here

    plev_name  = "Atmospheric Pressure levels"                   ; define lev coordinates
    plev_type  = "float"
    plev_units = "hPa"



;put var names you want to define as "varNames3D"
     varNames3D        = (/ "TEMP" /)
     varNames3D        = (/ "UWIND" /)
     varNames3D        = (/ "VWIND" /)
     varNames3D        = (/ "Qsh" /)
     varNames3D        = (/ "CDLIQ" /)
     varTypes3D        = (/ "float" /)

     varNames2D        = (/ "Psurf" /)
     varTypes2D        = (/ "float" /)

     varNames1D        = (/ "plev" /)
     varTypes1D        = (/ "float" /)


     filevardef(fout, varNames3D, varTypes3D, (/"time", "plev","lat" , "lon"/) )
     filevardef(fout, varNames2D, varTypes2D, (/"time", "lat" , "lon"/) )
     filevardef(fout, varNames1D, varTypes1D, (/"plev"/) )

 ;define variables
     TEMP  = TonP
     UWIND = UonP
     VWIND = VonP
     Qsh   = QonP
     CDLIQ = CLDLIQonP
     Psurf = PSmean
     plev  = pnew

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
