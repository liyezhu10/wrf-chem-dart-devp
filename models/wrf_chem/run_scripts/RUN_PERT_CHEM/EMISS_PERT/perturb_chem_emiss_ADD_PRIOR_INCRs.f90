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
!
! DART $Id: perturb_chem_emiss_INCRs.f90 13171 2019-05-09 16:42:36Z thoar@ucar.edu $

! code to ensemble mean and spread adjustemnts to the wrfchem emission files

! ifort -CB -C perturb_chem_emiss_ADD_PRIOR_INCRs.f90 -o perturb_chem_emiss_ADD_PRIOR_INCRs.exe -L${NETCDF}/lib -lnetcdff -lnetcdf -I${NETCDF}/include -L${CURC_MKL_LIB} -lmkl_rt -lmkl_lapack95_lp64 -lmkl_blas95_lp64


          program main
             implicit none
             integer                                  :: unit,nx,ny,nz,nzp,nz_chem,nz_fire,nz_biog
             integer                                  :: nchem_spcs,nfire_spcs,nbiog_spcs
             integer                                  :: i,j,k,isp,num_mem,imem
             real                                     :: pi,grav,nnum_mem,mean,std
             real                                     :: sprd_chem,sprd_fire,sprd_biog
             real                                     :: corr_lngth_hz
             real                                     :: corr_lngth_vt
             real                                     :: corr_lngth_tm
             real                                     :: corr_tm_delt
             real,allocatable,dimension(:,:,:)        :: chem_data3d
             real,allocatable,dimension(:,:,:)        :: chem_mean_incr
             real,allocatable,dimension(:,:,:)        :: chem_sprd_incr
             real,allocatable,dimension(:,:,:)        :: fire_mean_incr
             real,allocatable,dimension(:,:,:)        :: fire_sprd_incr
             real,allocatable,dimension(:,:,:)        :: biog_mean_incr
             real,allocatable,dimension(:,:,:)        :: biog_sprd_incr
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
! Chemi Emissions
             if(sw_chem) then 
                allocate(chem_data3d(nx,ny,nz_chem))
                allocate(chem_mean_incr(nx,ny,nz_chem))
                do isp=1,nchem_spcs
                   print *, 'process the chemi EMISSs ',trim(ch_chem_spc(isp))
                   do imem=1,num_mem
                      if(imem.ge.0.and.imem.lt.10) write(cmem,"('.e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(cmem,"('.e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('.e',i3)"),imem
! read updated
                      if(imem.eq.1) then
                         wrfchem_file='wrfchemi_d01_mean_incr'
                         call get_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp), &
                         chem_mean_incr,nx,ny,nz_chem)
                      endif
! read unadjusted emission
                      wrfchem_file=trim(wrfchemi)//trim(cmem)
                      call get_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp), &
                      chem_data3d,nx,ny,nz_chem)
! apply update
                      chem_data3d(:,:,:)=chem_data3d(:,:,:)+chem_mean_incr(:,:,:)
! Write update file
                      call put_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_data3d,nx,ny,nz_chem)
                   enddo
                enddo
                deallocate(chem_data3d)
                deallocate(chem_mean_incr)
             endif
!
! Fire Emissions
             if(sw_fire) then 
                allocate(chem_data3d(nx,ny,nz_fire))
                allocate(fire_mean_incr(nx,ny,nz_fire))
                do isp=1,nfire_spcs
                   print *, 'process the firei EMISSs ',trim(ch_fire_spc(isp))
                   do imem=1,num_mem
                      if(imem.ge.0.and.imem.lt.10) write(cmem,"('.e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(cmem,"('.e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('.e',i3)"),imem
! read update
                      if(imem.eq.1) then
                         wrffire_file='wrffirechemi_d01_mean_incr'
                         call get_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp), &
                         fire_mean_incr,nx,ny,nz_fire)
                      endif
! read unadjusted emission
                      wrffire_file=trim(wrffirechemi)//trim(cmem)
                      call get_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp), &
                      chem_data3d,nx,ny,nz_fire)
! apply update
                      chem_data3d(:,:,:)=chem_data3d(:,:,:)+fire_mean_incr(:,:,:)
! Write updated file
                      call put_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),chem_data3d,nx,ny,nz_fire)
                   enddo
                enddo
                deallocate(chem_data3d)
                deallocate(fire_mean_incr)
             endif 
!
! Biog Emissions
             if(sw_biog) then 
                allocate(chem_data3d(nx,ny,nz_biog))
                allocate(biog_mean_incr(nx,ny,nz_biog))
                do isp=1,nbiog_spcs
                   print *, 'process the biog EMISSs ',trim(ch_biog_spc(isp))
                   do imem=1,num_mem
                      if(imem.ge.0.and.imem.lt.10) write(cmem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(cmem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(cmem,"('e',i3)"),imem
! read update
                      if(imem.eq.1) then
                         wrfbiog_file='wrfbiochemi_d01_mean_incr'
                         call get_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp), &
                         biog_mean_incr,nx,ny,nz_biog)
                      endif
! read unadjusted emission
                      wrfbiog_file=trim(wrfbiogchemi)//trim(cmem)
                      call get_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp), &
                      chem_data3d,nx,ny,nz_biog)
! apply update
                      chem_data3d(:,:,:)=chem_data3d(:,:,:)+biog_mean_incr(:,:,:)
! Write updated file
                      call put_WRFCHEM_emiss_data(wrfbiog_file,ch_biog_spc(isp),chem_data3d,nx,ny,nz_biog)
                   enddo
                enddo
                deallocate(chem_data3d)
                deallocate(biog_mean_incr)
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
!             print *, 'f_id ',f_id
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
