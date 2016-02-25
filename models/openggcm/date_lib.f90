!This module contains (pure) f90 code in fixed-format with modules.  
!There is no fppn necessary to compile it. For the sake of future 
!projects, it might be nice to keep it that way :-) -- MLG
! 
!As such, it needs to be
!compiled BEFORE any other code which uses it.  It has no dependencies,
!so there is no harm in compiling it first.

!doc -
!doc - These modules contains stuff originally in mhd-cotr.for
!doc -
!doc - Original AUTHOR:
!doc -           Joachim Raeder, IGPP/UCLA, June 1994
!doc - The new modules are updated to minimize global data
!doc - and to use some of the neater features of fortran90
!doc -
!doc - Updating AUTHOR:
!doc -           Matthew L Gilson, SSC/UNH, January 2013

      module date_lib

      integer,private,parameter,dimension(12,2) :: mdays = RESHAPE( &
           (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, &
             0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335/), &
           (/12,2/))

      contains

!----------------------------------------------------------------
      subroutine epoch1966(dsecs,iy,mo,id,ih,mi,sec,iimod)
!----------------------------------------------------------------
!doc-
!doc-FUNCTION:
!doc-  converts seconds since JAN 1, 1966 0UT to year,month,day,hour,
!doc-  minute and second or vice versa
!doc-
!doc-PARAMETERS:
!doc-  dsec (in/out): seconds of epoch 1966
!doc-  iy,mo,id,ih,mi,se (in):  time in UT
!doc-  iimod (in,integer): flag
!doc-                     iimod=0: convert iy,... to dsec
!doc-                     iimod!=0: convert dsec to iy,...
!doc-
!doc-NOTE:
!doc-  works only for times after JAN 1, 1966, i.e. dsec positive
!doc-
      implicit none

      integer iy,mo,id,ih,mi,iimod
      real*8 dsecs
      real sec

      integer i,ij,jj
      real*8 ss,ds,dso
      integer,parameter,dimension(0:11) ::  &
           rdays = (/ 31,28,31,30,31,30,31,31,30,31,30,31 /)
      integer,parameter,dimension(0:11) ::  &
           jdays = (/ 31,29,31,30,31,30,31,31,30,31,30,31 /)
!
      if(iimod.eq.0) then
      jj=iy
      if(jj.lt.100)jj=jj+1900
      ss=0.0
      do i=1966,jj-1
         if(mod(i,4).eq.0) ss=ss+366.0d0*24.0d0*3600.0d0
         if(mod(i,4).ne.0) ss=ss+365.0d0*24.0d0*3600.0d0
      end do
      do i=1,mo-1
         if( mod(jj,4).eq.0 ) ss=ss+dble(jdays(i-1))*24.0d0*3600.0d0
         if( mod(jj,4).ne.0 ) ss=ss+dble(rdays(i-1))*24.0d0*3600.0d0
      end do
      ss=ss+dble(id-1)*3600.0d0*24.0d0
      ss=ss+dble(ih)*3600.0d0
      ss=ss+dble(mi)*60.0d0
      ss=ss+dble(sec)
      dsecs=ss
!
      else
      ss=dsecs
      ds=0.0d0
      dso=0.0d0
      jj=1966
      do i=0,100
         if(mod(jj,4).eq.0) then
            ds=ds+(366.0d0*24.0d0*3600.0d0)
         else
            ds=ds+(365.0d0*24.0d0*3600.0d0)
         endif
         if(ds.gt.ss) exit
         dso=ds
         jj=jj+1
      end do

      ss=ss-dso
      iy=jj-1900
      ij=0
      if(mod(jj,4).eq.0)ij=1
      ds=0.0d0
      dso=0.0d0
      jj=0
      do i=0,11
         if(ij.eq.0) then
            ds=ds+dble(rdays(i))*24.0d0*3600.0d0
         else
            ds=ds+dble(jdays(i))*24.0d0*3600.0d0
         endif
         if(ds.gt.ss) exit
         dso=ds
         jj=jj+1
      end do

      ss=ss-dso
      mo=jj+1
      ds=ss/(24.0d0*3600.0d0)
      id=int(ds)
      ss=ss-dble(id)*24.0d0*3600.0d0
      id=id+1
      ds=ss/3600.0d0
      ih=int(ds)
      ss=ss-dble(ih)*3600.0d0
      ds=ss/60.0d0
      mi=int(ds)
      ss=ss-dble(mi)*60.0d0
      sec=ss
      endif
!
      return
      end subroutine epoch1966
      
!---------------------------------------------------------------
      real*8 function djul(iy,mo,id,ih,mi,se)
!---------------------------------------------------------------
!doc-
!doc-FUNCTION:
!doc-  converts time to fractional julian days
!doc-
!doc-PARAMETERS:
!doc-  iy,mo,id,ih,mi,se (in):  time in UT
!doc-  djul (out,real*8):   fractional julian day
!doc-
!doc-NOTE:
!doc-  expected to work for any time after 4500 BC
!doc-
      implicit none

      integer iy,mo,id,ih,mi
      real*4 se

      integer jy
      integer i1,i2,i3,i4,ier
      real*8 d1,zero,one,xjd,uthour
      parameter(zero=0.0d0,one=1.0d0)
      jy=iy
      if(jy.lt.100)jy=jy+1900
      ier=0
      if(jy.lt.1901) ier=1
      if(jy.gt.2049) ier=7
      if(mo.lt.1.or.mo.gt.12) ier=2
      if(id.lt.1.or.id.gt.31) ier=3
      if(ih.lt.0.or.ih.gt.23) ier=4
      if(mi.lt.0.or.mi.gt.63) ier=5
!      if(se.lt.0.0.or.se.gt.60.0) ier=6
      if(ier.ne.0) then
         write(0,*)'error djul, ier= ',ier
!      stop
      endif
      uthour=dble(float(ih))+(dble(float(mi))/60.0d0)+ &
          dble(se)/3600.0d0
      i1=367*jy
      i2=(mo+9)/12
      i2=(i2+jy)*7
      i2=i2/4
      i3=(275*mo)/9
      i4=100*jy+mo
      d1=one
      if((dble(i4)-190002.5d0).lt.zero)d1=(-d1)
      xjd=dble(i1) - dble(i2) + dble(i3) + dble(id) + 1721013.5d0 &
           + uthour/24.0d0 - 0.5d0*d1 + 0.5
      djul=xjd
      return
      end function djul

!----------------------------------------------------
      integer function julday(iy,mo,id)
!----------------------------------------------------
!doc-
!doc-FUNCTION:
!doc-  convert date from day - month system to julian day system
!doc-
!doc-PARAMETERS:
!doc-  iy (in,integer):  year (nn or nnnn notation)
!doc-  mo (in,integer): month
!doc-  id (in,integer): day of month
!doc-  julday (out,integer): julian day of year
!doc-
      implicit none

      integer iy,mo,id
      integer leap,jmo

      jmo=max0(1,min0(12,mo))
      leap=1
      if(mod(iy,4).eq.0)leap=2
      julday=mdays(jmo,leap)+id
      return
      end function julday


!----------------------------------------------------
      subroutine daymon(iy,jd,mo,id)
!----------------------------------------------------
!doc-
!doc-FUNCTION:
!doc-  convert date from julian day system to day - month system
!doc-
!doc-PARAMETERS:
!doc-  iy (in,integer):  year (nn or nnnn notation)
!doc-  jd (in,integer):  julian day (1<=jd<=366), january 1 <> jd=1
!doc-  mo (out,integer): month
!doc-  id (out,integer): day of month
!doc-
      implicit none

      integer iy,jd,mo,id
      integer leap,i,k

      leap=1
      if(mod(iy,4).eq.0)leap=2
      k=1
      do i=1,12
         if(jd.le.mdays(i,leap)) exit
         k=i
      end do

      mo=k
      id=jd-mdays(k,leap)
      return
      end subroutine daymon

      end module date_lib

