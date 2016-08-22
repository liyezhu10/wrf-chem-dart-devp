
idart states:

0 --> not in DA cycle, OpenGGCM iterates until next DA time
1 --> begin DA cycle, OpenGGCM writes 
2 --> OpenGGCM finished write, wait for DART
3 --> startup, only read id_dem

Main
  call iono_srv() loops infinite, synchronizes with MPI   file:mhd-iono.for
  idart=3
  Loop:
    .
    call ionos2()  file:mhd-iono.for
      .
      call mhd_ctim()  file:mhd-misig.for
         return if tim < #ctim_start_time (300 sec)
         return if < 60 sec elapsed since last call (CTIM time step)
         .
         call ctim_interface()  file:ctim-core.for
            .
            if idart >= 2:
               Wait for `from_dart.semaphore'
               read id_dem (seconds until next assim)
               if id_dem > 0 (must be < 0 first time around, next assim will be uttime+abs(id_dem) )
                 Open dart_posterior.nc
                 Insert fields from DART
                 dart_nextwrite = uttime+abs(id_dem)
                 idart=0
               endif
            endif
            .
            if uttime > dart_nextwrite :
               open DATA.ionos2.nc
               write fields
               idart=1  signals ionos2 to write potential etc.
            endif
      back to mhd_ctim()
    back to ionos2()  file:mhd-iono.for
       if idart > 0:
          write potential, dB etc.  all diagnostic variables
          close DATA.ionos2.nc
          write 'to_dart.semaphore'
          idart=2
       endif
       
  Return to Loop:

hint:  to check, run OpenGGCM with DART=true and grep 'DART' from the log file

