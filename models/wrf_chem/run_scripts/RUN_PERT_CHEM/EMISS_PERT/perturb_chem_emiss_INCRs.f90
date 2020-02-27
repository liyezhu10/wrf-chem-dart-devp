! Copyright 2019 NCAR/ACOM
! 
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
! 
!     http://www.apache.org/licenses/LICENSE-2.0
! 
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

! code to ensemble mean and spread adjustemnts to the wrfchem emission files

! ifort -CB -C perturb_chem_emiss_INCRs.f90 -o perturb_chem_emiss_INCRs.exe -L${NETCDF}/lib -lnetcdff -lnetcdf -I${NETCDF}/include -L${CURC_MKL_LIB} -lmkl_rt -lmkl_lapack95_lp64 -lmkl_blas95_lp64

program main

implicit none

! version controlled file description for error handling, do not edit
character(len=*), parameter :: source   = 'perturb_chem_emiss_INCRs.f90'
character(len=*), parameter :: revision = ''
character(len=*), parameter :: revdate  = ''

             integer                                  :: unit,nx,ny,nz,nzp,nz_chem,nz_fire,nz_biog
             integer                                  :: nchem_spcs,nfire_spcs,nbiog_spcs
             integer                                  :: i,j,k,isp,num_mem,imem
             real                                     :: pi,grav,nnum_mem,mean,std
             real                                     :: sprd_chem,sprd_fire,sprd_biog
             real                                     :: corr_lngth_hz
             real                                     :: corr_lngth_vt
             real                                     :: corr_lngth_tm
             real                                     :: corr_tm_delt
             real,allocatable,dimension(:,:,:,:)      :: chem_data3d_sav
             real,allocatable,dimension(:,:,:)        :: chem_mean_prior,chem_mean_post,chem_mean_incr
             real,allocatable,dimension(:,:,:)        :: chem_sprd_prior,chem_sprd_post,chem_sprd_incr
             real,allocatable,dimension(:,:,:)        :: fire_mean_prior,fire_mean_post,fire_mean_incr
             real,allocatable,dimension(:,:,:)        :: fire_sprd_prior,fire_sprd_post,fire_sprd_incr
             real,allocatable,dimension(:,:,:)        :: biog_mean_prior,biog_mean_post,biog_mean_incr
             real,allocatable,dimension(:,:,:)        :: biog_sprd_prior,biog_sprd_post,biog_sprd_incr
             real,allocatable,dimension(:)            :: mems,pers
             character(len=150)                       :: pert_path_pr,pert_path_po
             character(len=150)                       :: wrfchemi,wrffirechemi,wrfbiogchemi
             character(len=20)                        :: cmem
             character(len=150)                       :: wrfchem_file,wrffire_file,wrfbiog_file
             character(len=150),allocatable,dimension(:) :: ch_chem_spc 
             character(len=150),allocatable,dimension(:) :: ch_fire_spc 
             character(len=150),allocatable,dimension(:) :: ch_biog_spc 
             logical                                  :: sw_corr_tm,sw_seed,sw_chem,sw_fire,sw_biog
             namelist /perturb_chem_emiss_corr_nml/nx,ny,nz,nz_chem,nchem_spcs,nfire_spcs,nbiog_spcs, &
             pert_path_pr,pert_path_po,nnum_mem,wrfchemi,wrffirechemi,wrfbiogchemi,sprd_chem,sprd_fire,sprd_biog, &
             sw_corr_tm,sw_seed,sw_chem,sw_fire,sw_biog,corr_lngth_hz,corr_lngth_vt,corr_lngth_tm,corr_tm_delt
             namelist /perturb_chem_emiss_spec_nml/ch_chem_spc,ch_fire_spc,ch_biog_spc
!
! Assign constants
             pi=4.*atan(1.)
             grav=9.8
             nz_fire=1
             nz_biog=1
!
! Read control namelist
             unit=20
             open(unit=unit,file='perturb_chem_emiss_corr_nml.nl',form='formatted', &
             status='old',action='read')
             read(unit,perturb_chem_emiss_corr_nml)
             close(unit)
             print *, 'nx                 ',nx
             print *, 'ny                 ',ny
             print *, 'nz                 ',nz
             print *, 'nz_chem            ',nz_chem
             print *, 'nchem_spcs         ',nchem_spcs
             print *, 'nfire_spcs         ',nfire_spcs
             print *, 'nbiog_spcs         ',nbiog_spcs
             print *, 'pert_path_pr       ',trim(pert_path_pr)
             print *, 'pert_path_po       ',trim(pert_path_po)
             print *, 'num_mem            ',nnum_mem
             print *, 'wrfchemi           ',trim(wrfchemi)
             print *, 'wrffirechemi       ',trim(wrffirechemi)
             print *, 'wrfbiogchemi       ',trim(wrfbiogchemi)
             print *, 'sprd_chem          ',sprd_chem
             print *, 'sprd_fire          ',sprd_fire
             print *, 'sprd_biog          ',sprd_biog
             print *, 'sw_corr_tm         ',sw_corr_tm
             print *, 'sw_seed            ',sw_seed
             print *, 'sw_chem            ',sw_chem
             print *, 'sw_fire            ',sw_fire
             print *, 'sw_biog            ',sw_biog
             print *, 'corr_lngth_hz      ',corr_lngth_hz
             print *, 'corr_lngth_vt      ',corr_lngth_vt
             print *, 'corr_lngth_tm      ',corr_lngth_tm
             print *, 'corr_tm_delt       ',corr_tm_delt
             nzp=nz+1
             num_mem=nint(nnum_mem)
!
! Allocate arrays
             allocate(ch_chem_spc(nchem_spcs),ch_fire_spc(nfire_spcs),ch_biog_spc(nbiog_spcs))
!
! Read the species namelist
             unit=20
             open(unit=unit,file='perturb_emiss_chem_spec_nml.nl',form='formatted', &
             status='old',action='read')
             read(unit,perturb_chem_emiss_spec_nml)
             close(unit)
!
! Calcuate ensemble mean adjustments
             if(sw_chem) then 
                allocate(chem_data3d_sav(nx,ny,nz_chem,num_mem))
                allocate(chem_mean_prior(nx,ny,nz_chem))
                allocate(chem_sprd_prior(nx,ny,nz_chem))
                allocate(chem_mean_post(nx,ny,nz_chem))
                allocate(chem_sprd_post(nx,ny,nz_chem))
                allocate(chem_mean_incr(nx,ny,nz_chem))
                allocate(chem_sprd_incr(nx,ny,nz_chem))
                allocate(mems(num_mem),pers(num_mem))
!
! Chemi Emissions
                do isp=1,nchem_spcs
                   print *, 'process the chemi EMISSs ',trim(ch_chem_spc(isp))
!
! Calcuate ensemble mean adjustments prior
                   do imem=1,num_mem
                      if(imem.ge.0.and.imem.lt.10) write(cmem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(cmem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('e',i3)"),imem
                      wrfchem_file='wrk_wrf_'//trim(cmem)//'/wrfchemi_d01'
                      call get_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp), &
                      chem_data3d_sav(:,:,:,imem),nx,ny,nz_chem)
                   enddo
                   do i=1,nx
                      do j=1,ny
                         do k=1,nz_chem
                            mems(:)=chem_data3d_sav(i,j,k,1:num_mem)
                            mean=sum(mems)/real(num_mem)
                            pers(:)=(mems(:)-mean)*(mems(:)-mean)
                            std=sqrt(sum(pers)/real(num_mem-1))
                            chem_mean_prior(i,j,k)=mean
                            chem_sprd_prior(i,j,k)=std
                         enddo
                      enddo
                   enddo
!
! Calcuate ensemble mean adjustments posterior
                   do imem=1,num_mem
                      if(imem.ge.0.and.imem.lt.10) write(cmem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(cmem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('e',i3)"),imem
                      wrfchem_file='wrk_dart_'//trim(cmem)//'/wrfchemi_d01'
                      call get_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp), &
                      chem_data3d_sav(:,:,:,imem),nx,ny,nz_chem)
                   enddo
                   do i=1,nx
                      do j=1,ny
                         do k=1,nz_chem
                            mems(:)=chem_data3d_sav(i,j,k,1:num_mem)
                            mean=sum(mems)/real(num_mem)
                            pers(:)=(mems(:)-mean)*(mems(:)-mean)
                            std=sqrt(sum(pers)/real(num_mem-1))
                            chem_mean_post(i,j,k)=mean
                            chem_sprd_post(i,j,k)=std
                         enddo
                      enddo
                   enddo
!
! Calcuate ensemble mean increment
                   chem_mean_incr(:,:,:)=chem_mean_post(:,:,:)-chem_mean_prior(:,:,:)
                   chem_sprd_incr(:,:,:)=chem_sprd_post(:,:,:)-chem_sprd_prior(:,:,:)
                   wrfchem_file='wrfchemi_d01_mean_prior'
                   call put_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_mean_prior,nx,ny,nz_chem)           
                   wrfchem_file='wrfchemi_d01_sprd_prior'
                   call put_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_sprd_prior,nx,ny,nz_chem)
                   wrfchem_file='wrfchemi_d01_mean_post'
                   call put_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_mean_post,nx,ny,nz_chem)
                   wrfchem_file='wrfchemi_d01_sprd_post'
                   call put_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_sprd_post,nx,ny,nz_chem)
                   wrfchem_file='wrfchemi_d01_mean_incr'
                   call put_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_mean_incr,nx,ny,nz_chem)
                   wrfchem_file='wrfchemi_d01_sprd_incr'
                   call put_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_sprd_incr,nx,ny,nz_chem)
                enddo                              ! species loop
                deallocate(chem_data3d_sav)
                deallocate(chem_mean_prior)
                deallocate(chem_sprd_prior)
                deallocate(chem_mean_post)
                deallocate(chem_sprd_post)
                deallocate(chem_mean_incr)
                deallocate(chem_sprd_incr)
                deallocate(mems,pers)                   
             endif
!
! Fire emissions
             if(sw_fire) then 
                allocate(chem_data3d_sav(nx,ny,nz_fire,num_mem))
                allocate(fire_mean_prior(nx,ny,nz_fire))
                allocate(fire_sprd_prior(nx,ny,nz_fire))
                allocate(fire_mean_post(nx,ny,nz_fire))
                allocate(fire_sprd_post(nx,ny,nz_fire))
                allocate(fire_mean_incr(nx,ny,nz_fire))
                allocate(fire_sprd_incr(nx,ny,nz_fire))
                allocate(mems(num_mem),pers(num_mem))
                do isp=1,nfire_spcs
                   print *, 'process the firei EMISSs ',trim(ch_fire_spc(isp))
!
! Calcuate ensemble mean adjustments prior
                   do imem=1,num_mem
                      if(imem.ge.0.and.imem.lt.10) write(cmem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(cmem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('e',i3)"),imem
                      wrffire_file='wrk_wrf_'//trim(cmem)//'/wrffirechemi_d01'
                      call get_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp), &
                      chem_data3d_sav(:,:,:,imem),nx,ny,nz_fire)
                   enddo
                   do i=1,nx
                      do j=1,ny
                         do k=1,nz_fire
                            mems(:)=chem_data3d_sav(i,j,k,1:num_mem)
                            mean=sum(mems)/real(num_mem)
                            pers(:)=(mems(:)-mean)*(mems(:)-mean)
                            std=sqrt(sum(pers)/real(num_mem-1))
                            fire_mean_prior(i,j,k)=mean
                            fire_sprd_prior(i,j,k)=std
                         enddo
                      enddo
                   enddo
!
! Calcuate ensemble mean adjustments posterior
                   do imem=1,num_mem
                      if(imem.ge.0.and.imem.lt.10) write(cmem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(cmem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('e',i3)"),imem
                      wrffire_file='wrk_dart_'//trim(cmem)//'/wrffirechemi_d01'
                      call get_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp), &
                      chem_data3d_sav(:,:,:,imem),nx,ny,nz_fire)
                   enddo
                   do i=1,nx
                      do j=1,ny
                         do k=1,nz_fire
                            mems(:)=chem_data3d_sav(i,j,k,1:num_mem)
                            mean=sum(mems)/real(num_mem)
                            pers(:)=(mems(:)-mean)*(mems(:)-mean)
                            std=sqrt(sum(pers)/real(num_mem-1))
                            fire_mean_post(i,j,k)=mean
                            fire_sprd_post(i,j,k)=std
                         enddo
                      enddo
                   enddo
!
! Calcuate ensemble mean increment
                   fire_mean_incr(:,:,:)=fire_mean_post(:,:,:)-fire_mean_prior(:,:,:)
                   fire_sprd_incr(:,:,:)=fire_sprd_post(:,:,:)-fire_sprd_prior(:,:,:)
                   wrffire_file='wrffirechemi_d01_mean_prior'
                   call put_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),fire_mean_prior,nx,ny,nz_fire)
                   wrffire_file='wrffirechemi_d01_sprd_prior'
                   call put_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),fire_sprd_prior,nx,ny,nz_fire)
                   wrffire_file='wrffirechemi_d01_mean_post'
                   call put_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),fire_mean_post,nx,ny,nz_fire)
                   wrffire_file='wrffirechemi_d01_sprd_post'
                   call put_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),fire_sprd_post,nx,ny,nz_fire)
                   wrffire_file='wrffirechemi_d01_mean_incr'
                   call put_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),fire_mean_incr,nx,ny,nz_fire)
                   wrffire_file='wrffirechemi_d01_sprd_incr'
                   call put_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),fire_sprd_incr,nx,ny,nz_fire)
                enddo                              ! species loop
                deallocate(chem_data3d_sav)
                deallocate(fire_mean_prior)
                deallocate(fire_sprd_prior)
                deallocate(fire_mean_post)
                deallocate(fire_sprd_post)
                deallocate(fire_mean_incr)
                deallocate(fire_sprd_incr)
                deallocate(mems,pers)                   
             endif
!
! Bio emissions
             if(sw_biog) then 
                allocate(chem_data3d_sav(nx,ny,nz_biog,num_mem))
                allocate(biog_mean_prior(nx,ny,nz_biog))
                allocate(biog_sprd_prior(nx,ny,nz_biog))
                allocate(biog_mean_post(nx,ny,nz_biog))
                allocate(biog_sprd_post(nx,ny,nz_biog))
                allocate(biog_mean_incr(nx,ny,nz_biog))
                allocate(biog_sprd_incr(nx,ny,nz_biog))
                allocate(mems(num_mem),pers(num_mem))
                do isp=1,nbiog_spcs
                   print *, 'process the bio EMISSs ',trim(ch_biog_spc(isp))
!
! Calcuate ensemble mean adjustments prior
                   do imem=1,num_mem
                      if(imem.ge.0.and.imem.lt.10) write(cmem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(cmem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('e',i3)"),imem
                      wrfbiog_file='wrk_wrf_'//trim(cmem)//'/wrfbiochemi_d01'
                      call get_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp), &
                      chem_data3d_sav(:,:,:,imem),nx,ny,nz_biog)
                   enddo
                   do i=1,nx
                      do j=1,ny
                         do k=1,nz_biog
                            mems(:)=chem_data3d_sav(i,j,k,1:num_mem)
                            mean=sum(mems)/real(num_mem)
                            pers(:)=(mems(:)-mean)*(mems(:)-mean)
                            std=sqrt(sum(pers)/real(num_mem-1))
                            biog_mean_prior(i,j,k)=mean
                            biog_sprd_prior(i,j,k)=std
                         enddo
                      enddo
                   enddo
!
! Calcuate ensemble mean adjustments posterior
                   do imem=1,num_mem
                      if(imem.ge.0.and.imem.lt.10) write(cmem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(cmem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('e',i3)"),imem
                      wrfbiog_file='wrk_dart_'//trim(cmem)//'/wrfbiochemi_d01'
                      call get_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp), &
                      chem_data3d_sav(:,:,:,imem),nx,ny,nz_biog)
                   enddo
                   do i=1,nx
                      do j=1,ny
                         do k=1,nz_biog
                            mems(:)=chem_data3d_sav(i,j,k,1:num_mem)
                            mean=sum(mems)/real(num_mem)
                            pers(:)=(mems(:)-mean)*(mems(:)-mean)
                            std=sqrt(sum(pers)/real(num_mem-1))
                            biog_mean_post(i,j,k)=mean
                            biog_sprd_post(i,j,k)=std
                         enddo
                      enddo
                   enddo
!
! Calcuate ensemble mean increment
                   biog_mean_incr(:,:,:)=biog_mean_post(:,:,:)-biog_mean_prior(:,:,:)
                   biog_sprd_incr(:,:,:)=biog_sprd_post(:,:,:)-biog_sprd_prior(:,:,:)
                   wrfbiog_file='wrfbiochemi_d01_mean_prior'
                   call put_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp),biog_mean_prior,nx,ny,nz_biog)
                   wrfbiog_file='wrfbiochemi_d01_sprd_prior'
                   call put_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp),biog_sprd_prior,nx,ny,nz_biog)
                   wrfbiog_file='wrfbiochemi_d01_mean_post'
                   call put_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp),biog_mean_post,nx,ny,nz_biog)
                   wrfbiog_file='wrfbiochemi_d01_sprd_post'
                   call put_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp),biog_sprd_post,nx,ny,nz_biog)
                   wrfbiog_file='wrfbiochemi_d01_mean_incr'
                   call put_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp),biog_mean_incr,nx,ny,nz_biog)
                   wrfbiog_file='wrfbiochemi_d01_sprd_incr'
                   call put_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp),biog_sprd_incr,nx,ny,nz_biog)
                enddo                              ! species loop
                deallocate(chem_data3d_sav)
                deallocate(biog_mean_prior)
                deallocate(biog_sprd_prior)
                deallocate(biog_mean_post)
                deallocate(biog_sprd_post)
                deallocate(biog_mean_incr)
                deallocate(biog_sprd_incr)
                deallocate(mems,pers)                   
             endif
             stop
          end program main
!
          subroutine get_WRFCHEM_emiss_data(file,name,data,nx,ny,nz_chem)
             implicit none
             include 'netcdf.inc'
             integer, parameter                    :: maxdim=6
             integer                               :: nx,ny,nz_chem
             integer                               :: i,rc
             integer                               :: f_id
             integer                               :: v_id,v_ndim,typ,natts
             integer,dimension(maxdim)             :: one
             integer,dimension(maxdim)             :: v_dimid
             integer,dimension(maxdim)             :: v_dim
             real,dimension(nx,ny,nz_chem)         :: data
             character(len=150)                    :: v_nam
             character*(*)                         :: name
             character*(*)                         :: file
!
! open netcdf file
             rc = nf_open(trim(file),NF_SHARE,f_id)
!             print *, trim(file)
             if(rc.ne.0) then
                print *, 'nf_open error in get ',rc, trim(file)
                stop
             endif
!
! get variables identifiers
             rc = nf_inq_varid(f_id,trim(name),v_id)
!             print *, v_id
             if(rc.ne.0) then
                print *, 'nf_inq_varid error ', v_id
                stop
             endif
!
! get dimension identifiers
             v_dimid=0
             rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!             print *, v_dimid
             if(rc.ne.0) then
                print *, 'nf_inq_var error ', v_dimid
                stop
             endif
!
! get dimensions
             v_dim(:)=1
             do i=1,v_ndim

                rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
             enddo
!             print *, v_dim
             if(rc.ne.0) then
                print *, 'nf_inq_dimlen error ', v_dim
                stop
             endif
!
! check dimensions
             if(nx.ne.v_dim(1)) then
                print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
                stop
             else if(ny.ne.v_dim(2)) then
                print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
                stop
             else if(nz_chem.ne.v_dim(3)) then             
                print *, 'ERROR: nz_chem dimension conflict ',nz_chem,v_dim(3)
                stop
             else if(1.ne.v_dim(4)) then             
                print *, 'ERROR: time dimension conflict ',1,v_dim(4)
                stop
             endif
!
! get data
             one(:)=1
             rc = nf_get_vara_real(f_id,v_id,one,v_dim,data)
             if(rc.ne.0) then
                print *, 'nf_get_vara_real ', data(1,1,1)
                stop
             endif
             rc = nf_close(f_id)
             return
          end subroutine get_WRFCHEM_emiss_data
!
          subroutine put_WRFCHEM_emiss_data(file,name,data,nx,ny,nz_chem)
             implicit none
             include 'netcdf.inc'
             integer, parameter                    :: maxdim=6
             integer                               :: nx,ny,nz_chem
             integer                               :: i,rc
             integer                               :: f_id
             integer                               :: v_id,v_ndim,typ,natts
             integer,dimension(maxdim)             :: one
             integer,dimension(maxdim)             :: v_dimid
             integer,dimension(maxdim)             :: v_dim
             real,dimension(nx,ny,nz_chem)         :: data
             character(len=150)                    :: v_nam
             character*(*)                         :: name
             character*(*)                         :: file
!
! open netcdf file
             rc = nf_open(trim(file),NF_WRITE,f_id)
             if(rc.ne.0) then
                print *, 'nf_open error in put ',rc, trim(file)
                stop
             endif
!             print *, 'f_id ',f_id,trim(file)
!
! get variables identifiers
             rc = nf_inq_varid(f_id,trim(name),v_id)
!             print *, v_id
             if(rc.ne.0) then
                print *, 'nf_inq_varid error ', v_id
                stop
             endif
!             print *, 'v_id ',v_id
!
! get dimension identifiers
             v_dimid=0
             rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!             print *, v_dimid
             if(rc.ne.0) then
                print *, 'nf_inq_var error ', v_dimid
                stop
             endif
!             print *, 'v_ndim, v_dimid ',v_ndim,v_dimid      
!
! get dimensions
             v_dim(:)=1
             do i=1,v_ndim
                rc = nf_inq_dimlen(f_id,v_dimid(i),v_dim(i))
             enddo
!             print *, v_dim
             if(rc.ne.0) then
                print *, 'nf_inq_dimlen error ', v_dim
                stop
             endif
!
! check dimensions
             if(nx.ne.v_dim(1)) then
                print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
                stop
             else if(ny.ne.v_dim(2)) then
                print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
                stop
             else if(nz_chem.ne.v_dim(3)) then             
                print *, 'ERROR: nz_chem dimension conflict ',nz_chem,v_dim(3)
                stop
             else if(1.ne.v_dim(4)) then             
                print *, 'ERROR: time dimension conflict ',1,v_dim(4)
                stop
             endif
!
! put data
             one(:)=1
!             rc = nf_close(f_id)
!             rc = nf_open(trim(file),NF_WRITE,f_id)
             rc = nf_put_vara_real(f_id,v_id,one(1:v_ndim),v_dim(1:v_ndim),data)
             if(rc.ne.0) then
                print *, 'nf_put_vara_real return code ',rc
                print *, 'f_id,v_id ',f_id,v_id
                print *, 'one ',one(1:v_ndim)
                print *, 'v_dim ',v_dim(1:v_ndim)
                stop
             endif
             rc = nf_close(f_id)
             return
          end subroutine put_WRFCHEM_emiss_data
