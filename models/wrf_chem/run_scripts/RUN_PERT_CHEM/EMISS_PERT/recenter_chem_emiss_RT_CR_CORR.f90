!
! code to recenter the perturbed wrfchem emission files
!
! ifort -C recenter_chem_emiss_RT_CR_CORR.f90 -o recenter_chem_emiss_RT_CR_CORR.exe -lgfortran -lnetcdff -lnetcdf
!
          program main
             implicit none
             integer,parameter                        :: num_mem=20
             integer,parameter                        :: nx=179
             integer,parameter                        :: ny=139
             integer,parameter                        :: nz=36
             integer,parameter                        :: nz_chem=11
             integer,parameter                        :: nchem_spc=43
             integer,parameter                        :: nfire_spc=31
             integer,parameter                        :: nbio_spc=1
             integer                                  :: unit,isp,imem
             real,allocatable,dimension(:,:)          :: parent_2d,ens_mean_2d
             real,allocatable,dimension(:,:,:)        :: parent_3d,ens_mean_3d
             real,allocatable,dimension(:,:,:)        :: chem_data2d
             real,allocatable,dimension(:,:,:,:)      :: chem_data3d
             character(len=10)                        :: ch_mem
             character(len=150)                       :: wrfchemi_parent,wrffirechemi_parent
             character(len=150)                       :: wrfbiochemi_parent
             character(len=150)                       :: wrfchemi,wrffirechemi,wrfbiochemi
             character(len=150)                       :: wrfchemi_mem,wrffirechemi_mem,wrfbiochemi_mem
             character(len=150),dimension(nchem_spc)  :: ch_chem_spc 
             character(len=150),dimension(nfire_spc)  :: ch_fire_spc 
             character(len=150),dimension(nbio_spc)   :: ch_bio_spc 
             logical                                  :: pert_chem,pert_fire,pert_bio
             namelist /recenter_chem_emiss_nml/wrfchemi_parent,wrffirechemi_parent, &
             wrfbiochemi_parent,wrfchemi,wrffirechemi,wrfbiochemi,pert_chem, &
             pert_fire,pert_bio
!
! Assign emission species names
             ch_chem_spc(1)='E_CO'
             ch_chem_spc(2)='E_NO'
             ch_chem_spc(3)='E_NO2'
             ch_chem_spc(4)='E_BIGALK'
             ch_chem_spc(5)='E_BIGENE'
             ch_chem_spc(6)='E_C2H4'
             ch_chem_spc(7)='E_C2H5OH'
             ch_chem_spc(8)='E_C2H6'
             ch_chem_spc(9)='E_C3H6'
             ch_chem_spc(10)='E_C3H8'
             ch_chem_spc(11)='E_CH2O'
             ch_chem_spc(12)='E_CH3CHO'
             ch_chem_spc(13)='E_CH3COCH3'
             ch_chem_spc(14)='E_CH3OH'
             ch_chem_spc(15)='E_MEK'
             ch_chem_spc(16)='E_SO2'
             ch_chem_spc(17)='E_TOLUENE'
             ch_chem_spc(18)='E_NH3'
             ch_chem_spc(19)='E_ISOP'
             ch_chem_spc(20)='E_C10H16'
             ch_chem_spc(21)='E_sulf'
             ch_chem_spc(22)='E_CO_A'
             ch_chem_spc(23)='E_CO_BB'
             ch_chem_spc(24)='E_CO02'
             ch_chem_spc(25)='E_CO03'
             ch_chem_spc(26)='E_XNO'
             ch_chem_spc(27)='E_XNO2'
             ch_chem_spc(28)='E_PM25I'
             ch_chem_spc(29)='E_PM25J'
             ch_chem_spc(30)='E_PM_10'
             ch_chem_spc(31)='E_ECI'
             ch_chem_spc(32)='E_ECJ'
             ch_chem_spc(33)='E_ORGI'
             ch_chem_spc(34)='E_ORGJ'
             ch_chem_spc(35)='E_SO4I'
             ch_chem_spc(36)='E_SO4J'
             ch_chem_spc(37)='E_NO3I'
             ch_chem_spc(38)='E_NO3J'
             ch_chem_spc(39)='E_NH4I'
             ch_chem_spc(40)='E_NH4J'
             ch_chem_spc(41)='E_PM_25'
             ch_chem_spc(42)='E_OC'
             ch_chem_spc(43)='E_BC'
!
             ch_fire_spc(1)='ebu_in_co'
             ch_fire_spc(2)='ebu_in_no'
             ch_fire_spc(3)='ebu_in_so2'
             ch_fire_spc(4)='ebu_in_bigalk'
             ch_fire_spc(5)='ebu_in_bigene'
             ch_fire_spc(6)='ebu_in_c2h4'
             ch_fire_spc(7)='ebu_in_c2h5oh'
             ch_fire_spc(8)='ebu_in_c2h6'
             ch_fire_spc(9)='ebu_in_c3h8'
             ch_fire_spc(10)='ebu_in_c3h6'
             ch_fire_spc(11)='ebu_in_ch2o'
             ch_fire_spc(12)='ebu_in_ch3cho'
             ch_fire_spc(13)='ebu_in_ch3coch3'
             ch_fire_spc(14)='ebu_in_ch3oh'
             ch_fire_spc(15)='ebu_in_mek'
             ch_fire_spc(16)='ebu_in_toluene'
             ch_fire_spc(17)='ebu_in_nh3'
             ch_fire_spc(18)='ebu_in_no2'
             ch_fire_spc(19)='ebu_in_open'
             ch_fire_spc(20)='ebu_in_c10h16'
             ch_fire_spc(21)='ebu_in_ch3cooh'
             ch_fire_spc(22)='ebu_in_cres'
             ch_fire_spc(23)='ebu_in_glyald'
             ch_fire_spc(24)='ebu_in_mgly'
             ch_fire_spc(25)='ebu_in_gly'
             ch_fire_spc(26)='ebu_in_acetol'
             ch_fire_spc(27)='ebu_in_isop'
             ch_fire_spc(28)='ebu_in_macr'
             ch_fire_spc(29)='ebu_in_mvk'
             ch_fire_spc(30)='ebu_in_oc'
             ch_fire_spc(31)='ebu_in_bc'
!
             ch_bio_spc(1)='MSEBIO_ISOP'
!
! Allocate arrays
             allocate(chem_data3d(nx,ny,nz_chem,num_mem),chem_data2d(nx,ny,num_mem))
             allocate(parent_3d(nx,ny,nz_chem),parent_2d(nx,ny))
             allocate(ens_mean_3d(nx,ny,nz_chem),ens_mean_2d(nx,ny))
!
! Read namelist
             unit=20
             open(unit=unit,file='recenter_chem_emiss_nml.nl',form='formatted', &
             status='old',action='read')
             read(unit,recenter_chem_emiss_nml)
             close(unit)
!
             if(pert_chem) then
                do isp=1,nchem_spc
!                   print *, 'processing species ',isp
                   call get_WRFCHEM_emiss_data(wrfchemi_parent,ch_chem_spc(isp), &
                   parent_3d,nx,ny,nz_chem)
                   do imem=1,num_mem
                      if(imem.lt.10) write(ch_mem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(ch_mem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(ch_mem,"('e',i3)"),imem
                      wrfchemi_mem=trim(wrfchemi)//trim(ch_mem)
                      call get_WRFCHEM_emiss_data(wrfchemi_mem,ch_chem_spc(isp), &
                      chem_data3d(:,:,:,imem),nx,ny,nz_chem)
                      ens_mean_3d(:,:,:)=ens_mean_3d(:,:,:)+ &
                      chem_data3d(:,:,:,imem)/real(num_mem)
                   enddo
                   do imem=1,num_mem
                      if(imem.lt.10) write(ch_mem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(ch_mem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(ch_mem,"('e',i3)"),imem
                      wrfchemi_mem=trim(wrfchemi)//trim(ch_mem)
                      chem_data3d(:,:,:,imem)=chem_data3d(:,:,:,imem) &
                      -ens_mean_3d(:,:,:)+parent_3d(:,:,:)   
                      call put_WRFCHEM_emiss_data(wrfchemi_mem,ch_chem_spc(isp), &
                      chem_data3d(:,:,:,imem),nx,ny,nz_chem)
                   enddo
                enddo
             endif
             if(pert_fire) then
                do isp=1,nfire_spc
!                   print *, 'processing species ',isp
                   call get_WRFCHEM_emiss_data(wrffirechemi_parent,ch_fire_spc(isp), &
                   parent_2d,nx,ny,1)
                   do imem=1,num_mem
                      if(imem.lt.10) write(ch_mem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(ch_mem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(ch_mem,"('e',i3)"),imem
                      wrffirechemi_mem=trim(wrffirechemi)//trim(ch_mem)
                      call get_WRFCHEM_emiss_data(wrffirechemi_mem,ch_fire_spc(isp), &
                      chem_data2d(:,:,imem),nx,ny,1)
                      ens_mean_2d(:,:)=ens_mean_2d(:,:)+ &
                      chem_data2d(:,:,imem)/real(num_mem)
                   enddo
                   do imem=1,num_mem
                      if(imem.lt.10) write(ch_mem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(ch_mem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(ch_mem,"('e',i3)"),imem
                      wrffirechemi_mem=trim(wrffirechemi)//trim(ch_mem)
                      chem_data2d(:,:,imem)=chem_data2d(:,:,imem) &
                      -ens_mean_2d(:,:)+parent_2d(:,:)   
                      call put_WRFCHEM_emiss_data(wrffirechemi_mem,ch_fire_spc(isp), &
                      chem_data2d(:,:,imem),nx,ny,1)
                   enddo
                enddo
             endif
             if(pert_bio) then
                do isp=1,nbio_spc
!                   print *, 'processing species ',isp
                   call get_WRFCHEM_emiss_data(wrfbiochemi_parent,ch_bio_spc(isp), &
                   parent_2d,nx,ny,1)
                   do imem=1,num_mem
                      if(imem.lt.10) write(ch_mem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(ch_mem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(ch_mem,"('e',i3)"),imem
                      wrfbiochemi_mem=trim(wrfbiochemi)//trim(ch_mem)
                      call get_WRFCHEM_emiss_data(wrfbiochemi_mem,ch_bio_spc(isp), &
                      chem_data2d(:,:,imem),nx,ny,1)
                      ens_mean_2d(:,:)=ens_mean_2d(:,:)+ &
                      chem_data2d(:,:,imem)/real(num_mem)
                   enddo
                   do imem=1,num_mem
                      if(imem.lt.10) write(ch_mem,"('e00',i1)"),imem
                      if(imem.ge.10.and.imem.lt.100) write(ch_mem,"('e0',i2)"),imem
                      if(imem.ge.100.and.imem.lt.1000) write(ch_mem,"('e',i3)"),imem
                      wrfbiochemi_mem=trim(wrfbiochemi)//trim(ch_mem)
                      chem_data2d(:,:,imem)=chem_data2d(:,:,imem) &
                      -ens_mean_2d(:,:)+parent_2d(:,:)   
                      call put_WRFCHEM_emiss_data(wrfbiochemi_mem,ch_bio_spc(isp), &
                      chem_data2d(:,:,imem),nx,ny,1)
                   enddo
                enddo
             endif
!
! Deallocate arrays
             deallocate(chem_data3d,chem_data2d)
             deallocate(parent_3d,parent_2d)
             deallocate(ens_mean_3d,ens_mean_2d)
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
             rc = nf_open(trim(file),NF_NOWRITE,f_id)
!             print *, trim(file)
             if(rc.ne.0) then
                print *, 'nf_open error ',trim(file)
                call abort
             endif
!
! get variables identifiers
             rc = nf_inq_varid(f_id,trim(name),v_id)
!             print *, v_id
             if(rc.ne.0) then
                print *, 'nf_inq_varid error ', v_id
                call abort
             endif
!
! get dimension identifiers
             v_dimid=0
             rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!             print *, v_dimid
             if(rc.ne.0) then
                print *, 'nf_inq_var error ', v_dimid
                call abort
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
                call abort
             endif
!
! check dimensions
             if(nx.ne.v_dim(1)) then
                print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
                call abort
             else if(ny.ne.v_dim(2)) then
                print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
                call abort
             else if(nz_chem.ne.v_dim(3)) then             
                print *, 'ERROR: nz_chem dimension conflict ',nz_chem,v_dim(3)
                call abort
             else if(1.ne.v_dim(4)) then             
                print *, 'ERROR: time dimension conflict ',1,v_dim(4)
                call abort
             endif
!
! get data
             one(:)=1
             rc = nf_get_vara_real(f_id,v_id,one,v_dim,data)
             if(rc.ne.0) then
                print *, 'nf_get_vara_real ', data(1,1,1)
                call abort
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
                print *, 'nf_open error ',trim(file)
                call abort
             endif
!
! get variables identifiers
             rc = nf_inq_varid(f_id,trim(name),v_id)
!             print *, v_id
             if(rc.ne.0) then
                print *, 'nf_inq_varid error ', v_id
                call abort
             endif
!
! get dimension identifiers
             v_dimid=0
             rc = nf_inq_var(f_id,v_id,v_nam,typ,v_ndim,v_dimid,natts)
!             print *, v_dimid
             if(rc.ne.0) then
                print *, 'nf_inq_var error ', v_dimid
                call abort
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
                call abort
             endif
!
! check dimensions
             if(nx.ne.v_dim(1)) then
                print *, 'ERROR: nx dimension conflict ',nx,v_dim(1)
                call abort
             else if(ny.ne.v_dim(2)) then
                print *, 'ERROR: ny dimension conflict ',ny,v_dim(2)
                call abort
             else if(nz_chem.ne.v_dim(3)) then             
                print *, 'ERROR: nz_chem dimension conflict ',nz_chem,v_dim(3)
                call abort
             else if(1.ne.v_dim(4)) then             
                print *, 'ERROR: time dimension conflict ',1,v_dim(4)
                call abort
             endif
!
! put data
             one(:)=1
             rc = nf_put_vara_real(f_id,v_id,one,v_dim,data)
             if(rc.ne.0) then
                print *, 'nf_put_vara_real ', data(1,1,1)
                call abort
             endif
             rc = nf_close(f_id)
             return
          end subroutine put_WRFCHEM_emiss_data   




