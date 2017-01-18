! DART software - Copyright © 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
program iasi_ascii_to_obs
!
! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   iasi_ascii_to_obs - a program that only needs minor customization to read
!      in a iasi_ascii-based dataset - either white-space separated values or
!      fixed-width column data.
!      
!     this is work in progress. IASI dataset are in HDF format. I do not
!     have HDF libraries for now, so Gabi Pfister reads the hdf file in IDL and
!     did some processing before she dumped the data in ascii. what you are
!     reading here is a 'processed' dataset of IASI O3
!
!     created 29 Mar 2010   nancy collins NCAR/IMAGe
!     modified 16 Mar 2012  ave arellano UArizona
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   use         types_mod, only : r8, PI, DEG2RAD
   use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                                 open_file, close_file
   use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                                 operator(>=), increment_time, get_time, &
                                 operator(-), GREGORIAN, operator(+), print_date
   use      location_mod, only : VERTISUNDEF, VERTISHEIGHT,VERTISPRESSURE
   use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                                 static_init_obs_sequence, init_obs, write_obs_seq, & 
                                 init_obs_sequence, get_num_obs, & 
                                 set_copy_meta_data, set_qc_meta_data
   use      obs_kind_mod, only : IASI_O3_RETRIEVAL
!
   implicit none
!
   character(len=64), parameter :: iasi_ascii_input_file = 'iasi_asciidata.input'
   character(len=64), parameter :: obs_out_file          = 'iasi_obs_seq.out'
   character (len=84) :: input_line
   character (len=10) :: otype_char
!
!
   logical, parameter :: debug = .false.  ! set to .true. to print info
   logical  :: first_obs,file_exist
!
   integer, parameter :: nlevels = 40
   integer, parameter :: nlevelsp = nlevels+1
   integer, parameter :: lwrk=5*nlevels
   integer  :: iunit, max_obs, num_copies, num_qc, qc_count
   integer  :: rcio, year, month, day, hour, minute, second
   integer  :: osec, oday, otype
   integer  :: i, j, jj, k, nlev_use, nlevp_use, klev, nrows
   integer  :: apm_accept = 0
   integer  :: apm_reject =0
!
   real(r8) :: fac,seconds, sec_sav, lat, lon, sza, cloud
   real(r8) :: zlev1, znlev1, qstatus, avg_tmp
   real(r8) :: alt_ag(nlevelsp) = 0.0_r8
   real(r8) :: press(nlevelsp) = 0.0_r8
   real(r8) :: press_out(nlevelsp) = 0.0_r8
   real(r8) :: AK(nlevels,nlevels) = 0.0_r8
   real(r8) :: AK_vmr(nlevels,nlevels) = 0.0_r8
   real(r8) :: AK_col(nlevels,nlevels) = 0.0_r8
   real(r8) :: ya_col(nlevels) = 0.0_r8
   real(r8) :: ya_vmr(nlevels) = 0.0_r8
   real(r8) :: yr_col(nlevels) = 0.0_r8
   real(r8) :: yr_vmr(nlevels) = 0.0_r8
   real(r8) :: err2_r(nlevels) = 0.0_r8
   real(r8) :: Ca_col(nlevels,nlevels) = 0.0_r8
   real(r8) :: Ca_vmr(nlevels,nlevels) = 0.0_r8
   real(r8) :: Cm_col(nlevels,nlevels) = 0.0_r8
   real(r8) :: Cm_vmr(nlevels,nlevels) = 0.0_r8
   real(r8) :: Cr_col(nlevels,nlevels) = 0.0_r8
   real(r8) :: Cr_vmr_rev(nlevels,nlevels) = 0.0_r8
   real(r8) :: Cr_vmr(nlevels,nlevels) = 0.0_r8
   real(r8) :: Cv_use(nlevels,nlevels) = 0.0_r8
   real(r8) :: Cv_use_col(nlevels,nlevels) = 0.0_r8
   real(r8) :: Cv_use_vmr(nlevels,nlevels) = 0.0_r8
!
   real(r8) :: qc, trm, rnlev_use, vert, retlev2
   real(r8) :: ImAya_col(nlevels) = 0.0_r8
   real(r8) :: ImAya_vmr(nlevels) = 0.0_r8
   real(r8) :: aircol(nlevels) = 0.0_r8
   real(r8) :: yqr_col(nlevels) = 0.0_r8
   real(r8) :: yqr_vmr(nlevels) = 0.0_r8
!
   double precision,allocatable,dimension(:,:)    :: Z
   double precision,allocatable,dimension(:,:)    :: U,U_T,V,V_T,SV_m
   double precision,allocatable,dimension(:)      :: SV,ZC
   double precision,dimension(lwrk)               :: wrk
   integer                                        :: info,nlev_trc,nlev_out
   real,dimension(nlevels,nlevels)                :: AK_T,ImAk,AKmI,AKmI_T,IpAKmI_T
   real,dimension(nlevels,nlevels)                :: AK_p,Ca_p,Cm_p,Cr_p 
   real(r8),dimension(nlevels,nlevels)            :: r_AK,r_Cv_use
   real(r8),dimension(nlevels)                    :: r_yr_col,r_ya_col,r_yqr_col,r_ImAya_col
   real(r8),dimension(nlevels)                    :: r_yr_vmr,r_ya_vmr,r_yqr_vmr,r_ImAya_vmr
   real                                           :: eps
!
   type(obs_type)          :: obs, prev_obs
   type(obs_sequence_type) :: obs_seq
   type(time_type)         :: time_obs, prev_time
!
! Start of executable code
   call initialize_utilities('iasi_ascii_to_obs')
!
! Time setup
   call set_calendar_type(GREGORIAN)
! 
! Open IASI ascii input file
   iunit = open_file(iasi_ascii_input_file, 'formatted', 'read')
   if (debug) print *, 'opened input file ' // trim(iasi_ascii_input_file)
!
! Each observation in this series will have a single observation value 
! and a quality control flag.  the max possible number of obs needs to
! be specified but it will only write out the actual number created.
   max_obs    = 10000000
   num_copies = 1
   num_qc     = 1
!
! Call the initialization code, and initialize two empty observation types
   call static_init_obs_sequence()
   call init_obs(obs,      num_copies, num_qc)
   call init_obs(prev_obs, num_copies, num_qc)
   first_obs = .true.
!
! Create a new, empty obs_seq file.  you must give a max limit
! on number of obs.  increase the size if too small.
   call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)
!
! The first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
   call set_copy_meta_data(obs_seq, 1, 'IASI O3 observation')
   call set_qc_meta_data(obs_seq, 1, 'IASI O3 QC index')
!
! If you want to append to existing files (e.g. you have a lot of
! small iasi_ascii files you want to combine), you can do it this way,
! or you can use the obs_sequence_tool to merge a list of files 
! once they are in DART obs_seq format.
!
! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
   qc = 0.0_r8
   qc_count = 0
!
! Loop for reading IASI ascii data file
   obsloop: do    ! No end limit - have the loop break when input ends
!
! Read in a line from the iasi_ascii file.   What you need to create an obs:
!  location: lat, lon, and height in pressure or meters
!  time: when the observation was taken
!  type: from the DART list of obs types
!  error: very important - the instrument error plus representativeness error
!        (see html file for more info)
!  averaging kernel and a priori profile
!  for now, we chose 2 'retrieval levels' corresponding to highest sensitivity
!  assume to be independent from each other
!  assume here a line is a type (1/2), location, time, value, obs error
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! READ IASI ASCII CODE BLOCK
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! READ RECORD ONE
!      print *, 'read record one'
      read(iunit, "(A)", iostat=rcio) input_line
      if (debug) print *, input_line
      if (rcio /= 0) then 
         if (debug) print *, 'got bad read code from input file, rcio = ', rcio
         exit obsloop
      endif
!
! Convert to date/time data
      read(input_line(1:10), *, iostat=rcio) otype_char
      if (rcio /= 0) then 
         if (debug) print *, 'got bad read code trying to get obs type, rcio = ', rcio
         exit obsloop
      endif
!   
! otype is fixed to 1 for IASI O3
      if (debug) print *, 'next observation type = ', otype_char
      read(input_line(12:15), *, iostat=rcio) year
      read(input_line(16:17), *, iostat=rcio) month
      read(input_line(18:19), *, iostat=rcio) day
!
! READ RECORD TWO
!      print *, 'read record two'
      read(iunit, *, iostat=rcio) seconds, lat, lon, sza, cloud, zlev1, &
      znlev1
      sec_sav=seconds
      if (debug) print *, 'zlev1, znlev1 ',seconds, lat, lon, sza, cloud, zlev1, znlev1
!
! Check lon and convert to 0-360
      if ( lon > 180.0_r8 .or. lon < -180.0_r8 ) cycle obsloop
      if ( lon < 0.0_r8 )  lon = lon + 360.0_r8
!
! Assign date/time data
      nlev_use = nint(znlev1)
      nlevp_use = nlev_use+1
      klev = nint(zlev1) 
      hour =  seconds/3600
      minute = (seconds - hour*3600)/60
      second = (seconds - hour*3600 - minute*60)
      if (debug) print *, "APM(2): yy, mm, dd, hh, mn, ss = ",year,month,day,hour,minute,second
      if (debug) print *, "APM(2): klev, nlev_use  = ",klev, nlev_use
!
! Put date/time data into DART format
      time_obs = set_date(year, month, day, hour, minute, second)
      if (debug) call print_date(time_obs, 'next obs time is')
!
! Put day and time into Gregorian format
      call get_time(time_obs, osec, oday)
!
! READ RECORD THREE
!      print *, 'read record three'
      read(iunit, *, iostat=rcio) (alt_ag(i),i=1,nlevp_use)
      if (debug) print *, "APM(3): alt_ag = ",(alt_ag(i),i=1,nlevp_use)
!
! READ RECORD FOUR
!      print *, 'read record four'
      read(iunit, *, iostat=rcio) (press(i),i=1,nlevp_use)
      if (debug) print *, "APM(4): press = ",(press(i),i=1,nlevp_use)
!
! READ RECORD FIVE (AK in partial columns)
!      print *, 'read record five'
      read(iunit,*, iostat=rcio) ((AK(j,i),i=1,nlev_use),j=1,nlev_use)
      if (debug .and. (int(sec_sav).eq.44755)) then
         do i=1,nlev_use
            print *, "APM(5): AK = ",i,(AK(i,j),j=1,nlev_use)
         enddo
         stop
      endif
!
! READ RECORD SIX
!      print *, 'read record six'
      read(iunit,*, iostat=rcio) (ya_col(i),i=1,nlev_use)
      if (debug) print *, "APM(6): ya_col = ", (ya_col(i),i=1,nlev_use)
!
! READ RECORD SEVEN
!      print *, 'read record seven'
      read(iunit,*, iostat=rcio) (ya_vmr(i),i=1,nlev_use)
      if (debug) print *, "APM(7): ya_vmr = ", (ya_vmr(i),i=1,nlev_use)
!
! READ RECORD EIGHT
!      print *, 'read record eight'
      read(iunit,*, iostat=rcio) (yr_col(i),i=1,nlev_use)
      if (debug) print *, "APM(8): yr_col = ", (yr_col(i),i=1,nlev_use)
!
! READ RECORD NINE
!      print *, 'read record nine'
      read(iunit,*, iostat=rcio) (yr_vmr(i),i=1,nlev_use)
      if (debug) print *, "APM(9): yr_vmr = ", (yr_vmr(i),i=1,nlev_use)
!
! READ RECORD TEN
!      print *, 'read record ten'
      read(iunit,*, iostat=rcio) ((Ca_col(j,i),i=1,nlev_use),j=1,nlev_use)
      if (debug) then
         do i=1,nlev_use
            print *, "APM(10): Ca_col = ",i,(Ca_vmr(i,j),j=1,nlev_use)
         enddo
      endif
!
! READ RECORD ELEVEN
!      print *, 'read record eleven'
      read(iunit,*, iostat=rcio) ((Cm_vmr(j,i),i=1,nlev_use),j=1,nlev_use)
      if (debug) then
         do i=1,nlev_use
            print *, "APM(11): Cm_col = ",i,(Cm_col(i,j),j=1,nlev_use)
         enddo
      endif
!
! READ RECORD TWELVE
!      print *, 'read record twelve'
      read(iunit,*, iostat=rcio) ((Cr_col(j,i),i=1,nlev_use),j=1,nlev_use)
      if (rcio /= 0) then 
         if (debug) print *, 'APM: Exit obs read loop, rcio = ', rcio
         exit obsloop
      endif
      if (debug) then
         do i=1,nlev_use
            print *, "APM(12): Cr_col = ",i,(Cr_col(i,j),j=1,nlev_use)
         enddo
      endif
!
! Check Ca, Cm, and Cr.  The following relationshipship should hold:
!    Cm=(I-AK) Ca (I+(AK-I)^T)
!    Cx=(I-AK) Ca
!    Ca, Cm, and Cr should have non-negative values along the diagonal
!    Ca, Cm, and Cr should be symmetric
!
! Calculate air column
      do i = 1, nlev_use
            aircol(i) = yr_col(i)/yr_vmr(i)
      enddo
!
! Convert error covariance from column to vmr
      do i = 1, nlev_use
         do j = 1, nlev_use
            Cr_vmr(i,j)=Cr_col(i,j)/aircol(i)/aircol(j)
            Ca_vmr(i,j)=Ca_col(i,j)/aircol(i)/aircol(j)
            Cm_vmr(i,j)=Cm_col(i,j)/aircol(i)/aircol(j)
         enddo
      enddo
!
! Convert AK_col to AK_vmr 
      AK_col(:,:)=AK(:,:)
      do i = 1, nlev_use
         do j = 1, nlev_use
            AK_vmr(i,j)=AK(i,j)/aircol(i)*aircol(j)
         enddo
      enddo
!
! Recalculate Cr_vmr from Ca_vmr
!      Cr_vmr_rev(:,:)=Cr_vmr(:,:)
!      Cr_vmr(:,:)=0.
!      do i = 1, nlev_use
!         do j = 1, nlev_use
!            do jj=1,nlev_use
!               if(i.eq.jj) then
!                  trm = (1. - AK_vmr(i,jj)) * Ca_vmr(jj,j) 
!               else
!                  trm = -1. * AK_vmr(i,jj) * Ca_vmr(jj,j)
!               endif
!               Cr_vmr(i,j)= Cr_vmr(i,j) + trm
!            enddo
!         enddo
!      enddo
!      do i=1,nlev_use
!         print *, 'Row ',i
!         print *, 'Cr_vmr_rev ',(Cr_vmr_rev(i,j),j=1,nlev_use)
!         print *, 'Cr_vmr ',(Cr_vmr(i,j),j=1,nlev_use)
!      enddo
!
! Assign error covariance to use
     fac=1.0
!     Cv_use_col(:,:)=fac*Cm_col(:,:)+(1.-fac)*Cr_col(:,:)
!     Cv_use_vmr(:,:)=fac*Cm_vmr(:,:)+(1.-fac)*Cr_vmr(:,:)
!     Cv_use_col(:,:)=fac*Cr_col(:,:)
     Cv_use_vmr(:,:)=fac*Cr_vmr(:,:)
!
     Cv_use(:,:)=Cv_use_vmr(:,:)
!
! Calculate (I-A)xa
!      ImAya_col(:) = 0.0_r8
!      do i = 1, nlev_use
!         do j = 1, nlev_use
!            if(i.eq.j) then
!               trm = (1. - AK_col(i,j)) * ya_col(j)
!            else
!               trm = -1. * AK_col(i,j) * ya_col(j)
!               ImAya_col(i) = ImAya_col(i) + trm
!            endif
!         enddo
!      enddo
      ImAya_vmr(:) = 0.0_r8
      do i = 1, nlev_use
         do j = 1, nlev_use
            if(i.eq.j) then
               trm = (1. - AK_vmr(i,j)) * ya_vmr(j)
            else
               trm = -1. * AK_vmr(i,j) * ya_vmr(j)
            endif
            ImAya_vmr(i) = ImAya_vmr(i) + trm
         enddo
      enddo
!
! Calculate "quasi-optimal" retriaval
!      yqr_col(:) = yr_col(:) - ImAya_col(:)
      yqr_vmr(:) = yr_vmr(:) - ImAya_vmr(:)
!
! Do SVD rotation based on Cr
     allocate(U(nlev_use,nlev_use),U_T(nlev_use,nlev_use))
     allocate(V(nlev_use,nlev_use),V_T(nlev_use,nlev_use))
     allocate(SV(nlev_use),ZC(nlev_use),SV_m(nlev_use,nlev_use),Z(nlev_use,nlev_use))
     Z(1:nlev_use,1:nlev_use)=dble(Cv_use(1:nlev_use,1:nlev_use))
     call dgesvd('A','A',nlev_use,nlev_use,Z,nlev_use,SV,U,nlev_use,V_T,nlev_use,wrk,lwrk,info)
! 
!     if(nlev_use.lt.40) then
!        SV_m(:,:)=dble(0.)
!        do i=1,nlev_use
!           SV_m(i,i)=SV(i)
!        enddo
!        call mat_tri_prd(U,SV_m,V_T,Z,nlev_use,nlev_use,nlev_use,nlev_use,nlev_use,nlev_use)
!        do i=1,nlev_use
!           print *,'Cv_use_orig ',(Cv_use(i,j),j=1,nlev_use)
!           print *,'Cv_use_tran ',(real(Z(i,j)),j=1,nlev_use)
!        enddo
!     endif
!     do i=1,nlev_trc
!        print *, 'r_Cv_use(i,j) ',(r_Cv_use(i,j),j=1,nlev_trc)
!     enddo
!
!      if (debug) then
!         do i=1,nlev_use
!            print *, "APM: U = ",i,(U(i,j),j=1,nlev_use)
!         enddo
!      endif
!
!     eps=1.e-5
!     do j=1,nlev_use
!        if(SV(j).ge.eps) then
!           nlev_trc=j
!        else
!           SV(j)=0.
!           U(:,j)=0.
!           V_T(j,:)=0.
!        endif
!     enddo
      nlev_trc=nlev_use
!
! Scale the singular values to adjust error (APM need to base this on pressure
! because number of MOP and IAS levels not the same)
!      fac=1.0
!      do j=1,nlev_trc
!         if(j.eq.1) then
!            fac=0.7
!         else if(j.eq.3.or.j.eq.4.or.j.eq.5.or.j.eq.6.or.j.eq.7.or.j.eq.8.or.j.eq.9.or. &
!         j.eq.10.or.j.eq.11.or.j.eq.12.or.j.eq.13.or.j.eq.14.or.j.eq.15.or.j.eq.16.or. &
!         j.eq.17.or.j.eq.18) then
!            fac=.6
!         endif
!         SV(j)=SV(j)*fac
!     enddo
!
! Scale the singular vectors
!     do j=1,nlev_trc
!        U(:,j)=U(:,j)/sqrt(SV(j))
!     enddo
!     call mat_transpose(U,Z,nlev_use,nlev_use)
!     U_T(:,:)=Z(:,:)
!     call vec_to_mat(SV,SV_m,nlev_use)
!     do i=1,nlev_use
!        print *, 'U_T(i,j) ',(real(U_T(i,j)),j=1,nlev_use)
!     enddo
!
! Do the rotation
!     call mat_prd(U_T,dble(AK_vmr),Z,nlev_use,nlev_use,nlev_use,nlev_use)
!     r_AK(:,:)=real(Z(:,:))
!     do i=1,nlev_trc
!        print *, 'r_AK(i,j) ',(r_AK(i,j),j=1,nlev_use)
!     enddo
!     call mat_tri_prd(U_T,dble(Cv_use),U,Z,nlev_use,nlev_use,nlev_use,nlev_use,nlev_use,nlev_use)
!     r_Cv_use(:,:)=real(Z(:,:))
!     do i=1,nlev_trc
!        print *, 'r_Cv_use(i,j) ',(r_Cv_use(i,j),j=1,nlev_trc)
!     enddo
!
!     call lh_mat_vec_prd(U_T,dble(yr_col),ZC,nlev_use)
!     r_yr_col(:)=real(ZC(:))
!     call lh_mat_vec_prd(U_T,dble(yr_vmr),ZC,nlev_use)
!     r_yr_vmr(:)=real(ZC(:))
!     print *, 'r_yr_vmr(:) ',r_yr_vmr(:)
!     call lh_mat_vec_prd(U_T,dble(ya_col),ZC,nlev_use)
!     r_ya_col(:)=real(ZC(:))
!     call lh_mat_vec_prd(U_T,dble(ya_vmr),ZC,nlev_use)
!     r_ya_vmr(:)=real(ZC(:))
!    print *, 'r_ya_vmr(:) ',r_ya_vmr(:)
!     call lh_mat_vec_prd(U_T,dble(yqr_col),ZC,nlev_use)
!     r_yqr_col(:)=real(ZC(:))
!     call lh_mat_vec_prd(U_T,dble(yqr_vmr),ZC,nlev_use)
!     r_yqr_vmr(:)=real(ZC(:))
!    print *, 'r_yqr_vmr(:) ',r_yqr_vmr(:)
!     call lh_mat_vec_prd(U_T,dble(ImAya_col),ZC,nlev_use)
!     r_ImAya_col(:)=real(ZC(:))
!     call lh_mat_vec_prd(U_T,dble(ImAya_vmr),ZC,nlev_use)
!     r_ImAya_vmr(:)=real(ZC(:))
!     print *, 'r_ImAya_vmr(:) ',r_ImAya_vmr(:)
!
! Get new errors (check if err2_r < 0 the qstatus=1)
     qstatus=0.0
     err2_r(:)=0.
!     do i=1,nlev_use
!        if (nint(r_Cv_use(i,i)) .ne. 1) then
!           qstatus = 1.0
!!           print *, 'APM:Obs error problem: ',nlev_use,i,r_Cv_use(i,i)
!        else
!           err2_r(i)=sqrt(r_Cv_use(i,i))
!        endif
!     enddo
     do i=1,nlev_use
         err2_r(i)=sqrt(Cv_use(i,1))
     enddo
!
     nrows=nlev_use
!     do i=1,nlev_use
!        if(abs(SV(i)).ge.eps) then
!           nrows=nrows+1
!        endif
!     enddo 
     if(nint(qstatus).eq.1) then
        apm_reject=apm_reject+1
        print *, 'nlev_use ',nlev_use
        print *, 'SV ',SV(1:nlev_use)
        print *, 'err2 ',err2_r(1:nlev_use)
     else
        apm_accept=apm_accept+1
     endif
     deallocate(U,U_T,V,V_T,SV,ZC,SV_m,Z)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! WRITE DATA TO OBS_SEQ FILE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Loop over all IASI vertical levels
      if (qstatus.eq.0) then
         nlev_out=nrows
         press_out(1:nlev_out)=press(1:nlev_out)*100.
         do i = 1, nlev_out
            qc_count = qc_count + 1
!
! This example assumes there is an obs type, where otype=1 is
! a temperature measured in height, and if otype=2, there's a wind
! speed and direction and height is pressure.  any kind of observation
! can use any of the vertical types; this is just an example.
!
! Fixed otype=1 for IASI O3
! No height since it is a column integrated quantity
!
! Make an obs derived type, and then add it to the sequence
!            call create_3d_obs(lat, lon, press_out(i), VERTISPRESSURE, r_yqr_vmr(i), &
!                            IASI_O3_RETRIEVAL, err2_r(i), oday, osec, qc, obs, &
!                            r_AK(i,1:nlev_use), r_ya_vmr(i), alt_ag(1:nlev_use), &
!                            press_out(1:nlev_use),aircol(1:nlev_use), qc_count, nlev_use)
            call create_3d_obs(lat, lon, press_out(i), VERTISPRESSURE, yqr_vmr(i), &
                            IASI_O3_RETRIEVAL, err2_r(i), oday, osec, qc, obs, &
                            AK_vmr(i,1:nlev_use), ImAya_vmr(i), alt_ag(1:nlev_use), &
                            press_out(1:nlev_use),aircol(1:nlev_use), qc_count, nlev_use)
            call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
            if (debug) print *, 'added iasi obs to output seq'
         enddo 
      endif  
   end do obsloop
!
! If we added obs to the sequence, write it out to a file
!   print *,"APM: obsloop exit num_obs ",get_num_obs(obs_seq)
   if ( get_num_obs(obs_seq) > 0 ) then
      if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
      call write_obs_seq(obs_seq, obs_out_file)
   endif
   print *,'APM: accept ',apm_accept
   print *,'APM: reject ',apm_reject
!
! End of main program
   call finalize_utilities()
!
   contains
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.  
!
!       NOTE: assumes the code is using the threed_sphere locations module, 
!             that the observation has a single data value and a single
!             qc value, and that this obs type has no additional required
!             data (e.g. gps and radar obs need additional data per obs)
!
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    okind - observation kind
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
!    obs   - observation type
!
!    metadata (AFAJ)
!    akcol - averaging kernel
!    apcol_val - a priori sub-column (I-A)xa
!    nlevels - number of levels
!    aircol_val - air sub column
!    qc_count - obs count (key)
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
   subroutine create_3d_obs(lat, lon, vval, vkind, obsv, okind, oerr, day, sec, qc, &
                            obs, akcol, apcol_val, altretlev, press, aircol_val, qc_count, nlev_use)
      use        types_mod, only : r8
      use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_kind, &
                             set_obs_def_error_variance, set_obs_def_location, set_obs_def_key
      use obs_def_iasi_O3_mod, only : set_obs_def_iasi_o3
      use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
      use time_manager_mod, only : time_type, set_time
      use     location_mod, only : set_location
!
      integer, parameter            :: nlevels = 40
      integer,        intent(in)    :: nlev_use
      integer,        intent(in)    :: okind, vkind, day, sec
      integer,        intent(in)    :: qc_count
      real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
      real(r8),       intent(in)    :: akcol(nlevels)
      real(r8),       intent(in)    :: altretlev(nlevels)
      real(r8),       intent(in)    :: press(nlevels)
      real(r8),       intent(in)    :: aircol_val(nlevels)
      real(r8),       intent(in)    :: apcol_val
      type(obs_type), intent(inout) :: obs
!
      real(r8)           :: obs_val(1), qc_val(1)
      type(obs_def_type) :: obs_def
!
      call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
      call set_obs_def_kind(obs_def, okind)
      call set_obs_def_time(obs_def, set_time(sec, day))
      call set_obs_def_error_variance(obs_def, oerr * oerr)
      call set_obs_def_iasi_o3(qc_count, apcol_val, altretlev(1:nlev_use), press(1:nlev_use), &
      akcol(1:nlev_use), aircol_val(1:nlev_use), nlev_use)
      call set_obs_def_key(obs_def, qc_count)
!
      obs_val(1) = obsv
      call set_obs_values(obs, obs_val)
      qc_val(1)  = qc
      call set_qc(obs, qc_val)
!
      call set_obs_def(obs, obs_def)
   end subroutine create_3d_obs
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   add_obs_to_seq -- adds an observation to a sequence.  inserts if first
!           obs, inserts with a prev obs to save searching if that's possible.
!
!     seq - observation sequence to add obs to
!     obs - observation, already filled in, ready to add
!     obs_time - time of this observation, in dart time_type format
!     prev_obs - the previous observation that was added to this sequence
!                (will be updated by this routine)
!     prev_time - the time of the previously added observation (will also
!                be updated by this routine)
!     first_obs - should be initialized to be .true., and then will be
!                updated by this routine to be .false. after the first obs
!                has been added to this sequence.
!
!     created Mar 8, 2010   nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine add_obs_to_seq(seq, obs, obs_time, prev_obs, prev_time, first_obs)
      use types_mod, only        : r8
      use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq
      use time_manager_mod, only : time_type, operator(>=)
!
      type(obs_sequence_type), intent(inout) :: seq
      type(obs_type),          intent(inout) :: obs, prev_obs
      type(time_type),         intent(in)    :: obs_time
      type(time_type),         intent(inout) :: prev_time
      logical,                 intent(inout) :: first_obs
!
! insert(seq,obs) always works (i.e. it inserts the obs in
! proper time format) but it can be slow with a long file.
! supplying a previous observation that is older (or the same
! time) as the new one speeds up the searching a lot.
!
      if(first_obs) then
         call insert_obs_in_seq(seq, obs)
         first_obs = .false.
      else               
         if(obs_time >= prev_time) then
            call insert_obs_in_seq(seq, obs, prev_obs)
         else
            call insert_obs_in_seq(seq, obs)
         endif
      endif
!
! update for next time
      prev_obs = obs
      prev_time = obs_time
   end subroutine add_obs_to_seq
!
    subroutine mat_prd(A_mat,B_mat,C_mat,na,ma,nb,mb)
!
! compute dot product of two matrics
    integer :: ma,na,mb,nb,i,j,k
    double precision :: A_mat(na,ma),B_mat(nb,mb),C_mat(na,mb)
!
! check that na=mb
    if(ma .ne. nb) then
       print *, 'Error in matrix dimension ma (cols) must equal nb (rows) ',ma,' ',nb
       stop
    endif
!
! initialze the product array
    C_mat(:,:)=0.
!
! calculate inner product
    do i=1,na
       do j=1,mb
          do k=1,mb
             C_mat(i,j)=C_mat(i,j)+A_mat(i,k)*B_mat(k,j) 
          enddo
       enddo
    enddo
    return
    end subroutine mat_prd
!
    subroutine mat_tri_prd(A_mat,B_mat,C_mat,D_mat,na,ma,nb,mb,nc,mc)
!
! compute dot product of three matrics D=A*B*C
    integer :: na,ma,nb,mb,nc,mc,i,j,k
    double precision :: A_mat(na,ma),B_mat(nb,mb),C_mat(nc,mc),D_mat(na,mc)
    double precision :: Z_mat(nb,mc)
!
! check that na=mb
    if(ma .ne. nb) then
       print *, 'Error in matrix dimension ma (cols) must equal nb (rows) ',ma,' ',nb
       stop
    endif
    if(mb .ne. nc) then
       print *, 'Error in matrix dimension mb (cols) must equal nc (rows) ',mb,' ',nc
       stop
    endif
!
! initialze the product array
    Z_mat(:,:)=0.
    D_mat(:,:)=0.
!
! calculate first inner product Z=B*C
    do i=1,nb
       do j=1,mc
          do k=1,mb
             Z_mat(i,j)=Z_mat(i,j)+B_mat(i,k)*C_mat(k,j) 
          enddo
       enddo
    enddo
!
! calculate second inner product D=A*Z
    do i=1,na
       do j=1,mc
          do k=1,ma
             D_mat(i,j)=D_mat(i,j)+A_mat(i,k)*Z_mat(k,j) 
          enddo
       enddo
    enddo
    return
    end subroutine mat_tri_prd
!
    subroutine vec_to_mat(a_vec,A_mat,n)
!
! compute dot product of two matrics
    integer :: n,i
    double precision :: a_vec(n),A_mat(n,n)
!
! initialze the product array
    A_mat(:,:)=0.
!
! calculate inner product
    do i=1,n
       A_mat(i,i)=a_vec(i) 
    enddo
    return
    end subroutine vec_to_mat
!
    subroutine diag_inv_sqrt(A_mat,n)
!
! calculate inverse square root of diagonal elements
    integer :: n,i
    double precision :: A_mat(n,n)
    do i=1,n
       if(A_mat(i,i).le.0.) then
          print *, 'Error in Subroutine vec_to_mat arg<=0 ',i,' ',A_mat(i,i)
          call abort
       endif
       A_mat(i,i)=1./sqrt(A_mat(i,i)) 
    enddo
    return
    end subroutine diag_inv_sqrt
!
    subroutine diag_sqrt(A_mat,n)
!
! calculate square root of diagonal elements
    integer :: n,i
    double precision :: A_mat(n,n)
    do i=1,n
       if(A_mat(i,i).lt.0.) then
          print *, 'Error in Subroutine vec_to_mat arg<0 ',i,' ',A_mat(i,i)
          call abort
       endif
       A_mat(i,i)=sqrt(A_mat(i,i)) 
    enddo
    return
    end subroutine diag_sqrt
!
    subroutine lh_mat_vec_prd(SCL_mat,a_vec,s_a_vec,n)
!
! calculate left hand side scaling of column vector
    integer :: n,i,j
    double precision :: SCL_mat(n,n),a_vec(n),s_a_vec(n)
!
! initialize s_a_vec
    s_a_vec(:)=0.
!
! conduct scaling
    do i=1,n
       do j=1,n
          s_a_vec(i)=s_a_vec(i)+SCL_mat(i,j)*a_vec(j)
       enddo 
    enddo
    return
    end subroutine lh_mat_vec_prd
!
    subroutine rh_vec_mat_prd(SCL_mat,a_vec,s_a_vec,n)
!
! calculate right hand side scaling of a row vector
    integer :: n,i,j
    double precision :: SCL_mat(n,n),a_vec(n),s_a_vec(n)
!
! initialize s_a_vec
    s_a_vec(:)=0.
!
! conduct scaling
    do i=1,n
       do j=1,n
          s_a_vec(i)=s_a_vec(i)+a_vec(j)*SCL_mat(j,i) 
       enddo
    enddo
    return
    end subroutine rh_vec_mat_prd
!
    subroutine mat_transpose(A_mat,AT_mat,n,m)
!
! calculate matrix transpose
    integer :: n,m,i,j
    double precision :: A_mat(n,m),AT_mat(m,n)
    do i=1,n
       do j=1,m
          AT_mat(j,i)=A_mat(i,j) 
       enddo
    enddo
    return
    end subroutine mat_transpose
!
    subroutine diag_vec(A_mat,a_vec,n)
!
! calculate square root of diagonal elements
    integer :: n,i
    double precision :: A_mat(n,n),a_vec(n)
    do i=1,n
       a_vec(i)=A_mat(i,i) 
    enddo
    return
    end subroutine diag_vec
!
end program iasi_ascii_to_obs


