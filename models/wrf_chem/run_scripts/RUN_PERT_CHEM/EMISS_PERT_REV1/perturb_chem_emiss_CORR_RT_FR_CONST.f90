
! code to perturb the wrfchem emission files
!
! ifort -C perturb_chem_emiss_CORR_RT_FR_CONST.f90 -o perturb_chem_emiss_CORR_RT_FR_CONST.exe -lgfortran -lnetcdff -lnetcdf
!
          program main
             implicit none
             integer                                  :: i,j,k,num_mem,imem,iimem
             integer                                  :: unit,unita,isp
             integer                                  :: nx,ny,nz,nz_chem,nchem_spc,nfire_spc,nbio_spc
             real                                     :: pi,u_ran_1,u_ran_2,nnum_mem
             real                                     :: sprd_chem,sprd_fire,sprd_biog
             real                                     :: pert_land_chem,pert_watr_chem
             real                                     :: pert_land_fire,pert_watr_fire
             real                                     :: pert_land_biog,pert_watr_biog
             real                                     :: pert_land_chem_d,pert_watr_chem_d
             real                                     :: pert_land_fire_d,pert_watr_fire_d
             real                                     :: pert_land_biog_d,pert_watr_biog_d
             real                                     :: pert_land_chem_m,pert_watr_chem_m
             real                                     :: pert_land_fire_m,pert_watr_fire_m
             real                                     :: pert_land_biog_m,pert_watr_biog_m
             real,allocatable,dimension(:,:)          :: xland
             real,allocatable,dimension(:,:)          :: chem_data2d
             real,allocatable,dimension(:,:,:)        :: chem_data3d
             character(len=20)                        :: cmem
             character(len=150)                       :: pert_path
             character(len=150)                       :: wrfchemi,wrffirechemi,wrfbiochemi
             character(len=150)                       :: wrfchem_file,wrffire_file,wrfbio_file
             character(len=150),allocatable,dimension(:) :: ch_chem_spc 
             character(len=150),allocatable,dimension(:) :: ch_fire_spc 
             character(len=150),allocatable,dimension(:) :: ch_bio_spc 
             logical                                  :: sw_gen,sw_chem,sw_fire,sw_biog
             namelist /perturb_chem_emiss_CORR_nml/nx,ny,nz,nz_chem,nchem_spc,nfire_spc,nbio_spc, &
             pert_path,nnum_mem,wrfchemi,wrffirechemi,wrfbiochemi,sprd_chem,sprd_fire,sprd_biog, &
             sw_gen,sw_chem,sw_fire,sw_biog
!
! Assign constants
             pi=4.*atan(1.)
             pert_land_chem=-9999
             pert_watr_chem=-9999
             pert_land_fire=-9999
             pert_watr_fire=-9999
             pert_land_biog=-9999
             pert_watr_biog=-9999
!
! Read namelist
             unit=20
             open(unit=unit,file='perturb_chem_emiss_CORR_nml.nl',form='formatted', &
             status='old',action='read')
             read(unit,perturb_chem_emiss_CORR_nml)
             close(unit)
!
             print *, 'nx                ',nx
             print *, 'ny                ',ny
             print *, 'nz                ',nz
             print *, 'nz_chem           ',nz_chem
             print *, 'nchem_spc         ',nchem_spc
             print *, 'nfire_spc         ',nfire_spc
             print *, 'nbio_spc          ',nbio_spc
             print *, 'pert_path         ',trim(pert_path)
             print *, 'num_mem           ',nnum_mem
             print *, 'wrfchemi          ',trim(wrfchemi)
             print *, 'wrffirechemi      ',trim(wrffirechemi)
             print *, 'wrfbiochemi       ',trim(wrfbiochemi)
             print *, 'sprd_chem         ',sprd_chem
             print *, 'sprd_fire         ',sprd_fire
             print *, 'sprd_biog         ',sprd_biog
             print *, 'sw_gen            ',sw_gen
             print *, 'sw_chem           ',sw_chem
             print *, 'sw_fire           ',sw_fire
             print *, 'sw_biog           ',sw_biog
             num_mem=nint(nnum_mem)
             allocate(ch_chem_spc(nchem_spc),ch_fire_spc(nfire_spc),ch_bio_spc(nbio_spc))
!
             unita=30
             open(unit=unita,file=trim(pert_path)//'/pert_file_emiss',form='unformatted', &
             status='unknown')
!
! Assign emission species names
!
! FRAPPE
!             ch_chem_spc(1)='E_CO'
!             ch_chem_spc(2)='E_NO'
!             ch_chem_spc(3)='E_NO2'
!             ch_chem_spc(4)='E_BIGALK'
!             ch_chem_spc(5)='E_BIGENE'
!             ch_chem_spc(6)='E_C2H4'
!             ch_chem_spc(7)='E_C2H5OH'
!             ch_chem_spc(8)='E_C2H6'
!             ch_chem_spc(9)='E_C3H6'
!             ch_chem_spc(10)='E_C3H8'
!             ch_chem_spc(11)='E_CH2O'
!             ch_chem_spc(12)='E_CH3CHO'
!             ch_chem_spc(13)='E_CH3COCH3'
!             ch_chem_spc(14)='E_CH3OH'
!             ch_chem_spc(15)='E_MEK'
!             ch_chem_spc(16)='E_SO2'
!             ch_chem_spc(17)='E_TOLUENE'
!             ch_chem_spc(18)='E_NH3'
!             ch_chem_spc(19)='E_ISOP'
!             ch_chem_spc(20)='E_C10H16'
!             ch_chem_spc(21)='E_sulf'
!             ch_chem_spc(22)='E_CO_A'
!             ch_chem_spc(23)='E_CO_BB'
!             ch_chem_spc(24)='E_CO02'
!             ch_chem_spc(25)='E_CO03'
!             ch_chem_spc(26)='E_XNO'
!             ch_chem_spc(27)='E_XNO2'
!             ch_chem_spc(28)='E_PM25I'
!             ch_chem_spc(29)='E_PM25J'
!             ch_chem_spc(30)='E_PM_10'
!             ch_chem_spc(31)='E_ECI'
!             ch_chem_spc(32)='E_ECJ'
!             ch_chem_spc(33)='E_ORGI'
!             ch_chem_spc(34)='E_ORGJ'
!             ch_chem_spc(35)='E_SO4I'
!             ch_chem_spc(36)='E_SO4J'
!             ch_chem_spc(37)='E_NO3I'
!             ch_chem_spc(38)='E_NO3J'
!             ch_chem_spc(39)='E_NH4I'
!             ch_chem_spc(40)='E_NH4J'
!             ch_chem_spc(41)='E_PM_25'
!             ch_chem_spc(42)='E_OC'
!             ch_chem_spc(43)='E_BC'
!             ch_chem_spc(44)='E_BALD'
!             ch_chem_spc(45)='E_C2H2'
!             ch_chem_spc(46)='E_BENZENE'
!             ch_chem_spc(47)='E_XYLENE'
!             ch_chem_spc(48)='E_CRES'
!             ch_chem_spc(49)='E_HONO'
!
! PANDA
             ch_chem_spc(1)='E_CO'
             ch_chem_spc(2)='E_NO'
             ch_chem_spc(3)='E_NO2'
             ch_chem_spc(4)='E_SO2'
             ch_chem_spc(5)='E_BIGALK'
             ch_chem_spc(6)='E_C2H4'
             ch_chem_spc(7)='E_C2H5OH'
             ch_chem_spc(8)='E_C2H6'
             ch_chem_spc(9)='E_C3H6'
             ch_chem_spc(10)='E_C3H8'
             ch_chem_spc(11)='E_CH2O'
             ch_chem_spc(12)='E_CH3CHO'
             ch_chem_spc(13)='E_BIGENE'
             ch_chem_spc(14)='E_CH3COCH3'
             ch_chem_spc(15)='E_CH3OH'
             ch_chem_spc(16)='E_MEK'
             ch_chem_spc(17)='E_TOLUENE'
             ch_chem_spc(18)='E_ISOP'
             ch_chem_spc(19)='E_C10H16'
             ch_chem_spc(20)='E_NH3'
             ch_chem_spc(21)='E_OC'
             ch_chem_spc(22)='E_BC'
             ch_chem_spc(23)='E_PM_10'
             ch_chem_spc(24)='E_PM_25'
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
! Get the land mask data
             print *, 'At read for xland'
             allocate(xland(nx,ny))
             call get_WRFINPUT_land_mask(xland,nx,ny)
!
! Allocate arrays
             allocate(chem_data3d(nx,ny,nz_chem),chem_data2d(nx,ny))
!
             do iimem=1,num_mem,2
                if(sw_gen) then
                   do while (pert_land_chem.le.-1 .or. pert_land_chem.ge.1)
                      call random_number(u_ran_1)
                      if(u_ran_1.eq.0.) call random_number(u_ran_1)
                      call random_number(u_ran_2)
                      if(u_ran_2.eq.0.) call random_number(u_ran_2)
                      pert_land_chem=sqrt(-2.*log(u_ran_1))*cos(2.*pi*u_ran_2)*sprd_chem
                   end do
                   do while (pert_watr_chem.le.-1 .or. pert_watr_chem.ge.1)
                      call random_number(u_ran_1)
                      if(u_ran_1.eq.0.) call random_number(u_ran_1)
                      call random_number(u_ran_2)
                      if(u_ran_2.eq.0.) call random_number(u_ran_2)
                      pert_watr_chem=sqrt(-2.*log(u_ran_1))*cos(2.*pi*u_ran_2)*sprd_chem
                   end do
                   do while (pert_land_fire.le.-1 .or. pert_land_fire.ge.1)
                      call random_number(u_ran_1)
                      if(u_ran_1.eq.0.) call random_number(u_ran_1)
                      call random_number(u_ran_2)
                      if(u_ran_2.eq.0.) call random_number(u_ran_2)
                      pert_land_fire=sqrt(-2.*log(u_ran_1))*cos(2.*pi*u_ran_2)*sprd_fire
                   end do
                   do while (pert_watr_fire.le.-1 .or. pert_watr_fire.ge.1)
                      call random_number(u_ran_1)
                      if(u_ran_1.eq.0.) call random_number(u_ran_1)
                      call random_number(u_ran_2)
                      if(u_ran_2.eq.0.) call random_number(u_ran_2)
                      pert_watr_fire=sqrt(-2.*log(u_ran_1))*cos(2.*pi*u_ran_2)*sprd_fire
                   end do
                   do while (pert_land_biog.le.-1 .or. pert_land_biog.ge.1)
                      call random_number(u_ran_1)
                      if(u_ran_1.eq.0.) call random_number(u_ran_1)
                      call random_number(u_ran_2)
                      if(u_ran_2.eq.0.) call random_number(u_ran_2)
                      pert_land_biog=sqrt(-2.*log(u_ran_1))*cos(2.*pi*u_ran_2)*sprd_biog
                   end do
                   do while (pert_watr_biog.le.-1 .or. pert_watr_biog.ge.1)
                      call random_number(u_ran_1)
                      if(u_ran_1.eq.0.) call random_number(u_ran_1)
                      call random_number(u_ran_2)
                      if(u_ran_2.eq.0.) call random_number(u_ran_2)
                      pert_watr_biog=sqrt(-2.*log(u_ran_1))*cos(2.*pi*u_ran_2)*sprd_biog
                   end do
                   print *, 'At emiss pert write ',iimem
                   pert_land_chem_d=pert_land_chem
                   pert_watr_chem_d=pert_watr_chem
                   pert_land_fire_d=pert_land_fire
                   pert_watr_fire_d=pert_watr_fire
                   pert_land_biog_d=pert_land_biog
                   pert_watr_biog_d=pert_watr_biog
                   pert_land_chem_m=-1.*pert_land_chem
                   pert_watr_chem_m=-1.*pert_watr_chem
                   pert_land_fire_m=-1.*pert_land_fire
                   pert_watr_fire_m=-1.*pert_watr_fire
                   pert_land_biog_m=-1.*pert_land_biog
                   pert_watr_biog_m=-1.*pert_watr_biog
                   write(unita) pert_land_chem_d,pert_watr_chem_d,pert_land_fire_d,pert_watr_fire_d, &
                   pert_land_biog_d,pert_watr_biog_d
                   write(unita) pert_land_chem_m,pert_watr_chem_m,pert_land_fire_m, &
                   pert_watr_fire_m,pert_land_biog_m,pert_watr_biog_m
                else
                   print *, 'At emiss pert read ',iimem
                   read(unita) pert_land_chem_d,pert_watr_chem_d,pert_land_fire_d,pert_watr_fire_d, &
                   pert_land_biog_d,pert_watr_biog_d
                   read(unita) pert_land_chem_m,pert_watr_chem_m,pert_land_fire_m, &
                   pert_watr_fire_m,pert_land_biog_m,pert_watr_biog_m
                endif 
                print *, pert_land_chem_d,pert_watr_chem_d
                print *, pert_land_chem_m,pert_watr_chem_m
                print *, pert_land_fire_d,pert_watr_fire_d
                print *, pert_land_fire_m,pert_watr_fire_m
                print *, pert_land_biog_d,pert_watr_biog_d
                print *, pert_land_biog_m,pert_watr_biog_m
!
! draw member
                if(sw_chem) then
                   imem=iimem
                   if(imem.ge.0.and.imem.lt.10) write(cmem,"('.e00',i1)"),imem
                   if(imem.ge.10.and.imem.lt.100) write(cmem,"('.e0',i2)"),imem
                   if(imem.ge.100.and.imem.lt.1000) write(cmem,"('.e',i3)"),imem
                   wrfchem_file=trim(wrfchemi)//trim(cmem)
                   do isp=1,nchem_spc
                      call get_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_data3d,nx,ny,nz_chem)
                      do i=1,nx
                         do j=1,ny
                            if(xland(i,j).eq.1) then
                               chem_data3d(i,j,1:nz_chem)=(1.+pert_land_chem_d)*chem_data3d(i,j,1:nz_chem)
                            else
                               chem_data3d(i,j,1:nz_chem)=(1.+pert_watr_chem_d)*chem_data3d(i,j,1:nz_chem)
                            endif
                         enddo
                      enddo  
                      call put_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_data3d,nx,ny,nz_chem)
                   enddo
!
! mirror member
                   imem=iimem+1
                   if(sw_gen) then
                      pert_land_chem=-1.*pert_land_chem
                      pert_watr_chem=-1.*pert_watr_chem
                   endif
                   if(imem.ge.0.and.imem.lt.10) write(cmem,"('.e00',i1)"),imem
                   if(imem.ge.10.and.imem.lt.100) write(cmem,"('.e0',i2)"),imem
                   if(imem.ge.100.and.imem.lt.1000) write(cmem,"('.e',i3)"),imem
                   wrfchem_file=trim(wrfchemi)//trim(cmem)
                   do isp=1,nchem_spc
                      call get_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_data3d,nx,ny,nz_chem)
                      do i=1,nx
                         do j=1,ny
                            if(xland(i,j).eq.1) then
                               chem_data3d(i,j,1:nz_chem)=(1.+pert_land_chem_m)*chem_data3d(i,j,1:nz_chem)
                            else
                               chem_data3d(i,j,1:nz_chem)=(1.+pert_watr_chem_m)*chem_data3d(i,j,1:nz_chem)
                            endif
                         enddo
                      enddo  
                      call put_WRFCHEM_emiss_data(wrfchem_file,ch_chem_spc(isp),chem_data3d,nx,ny,nz_chem)
                   enddo
                endif
!
                if(sw_fire) then
!
! draw member
                   imem=iimem
                   if(imem.ge.0.and.imem.lt.10) write(cmem,"('.e00',i1)"),imem
                   if(imem.ge.10.and.imem.lt.100) write(cmem,"('.e0',i2)"),imem
                   if(imem.ge.100.and.imem.lt.1000) write(cmem,"('.e',i3)"),imem
                   wrffire_file=trim(wrffirechemi)//trim(cmem)
                   do isp=1,nfire_spc
                      call get_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),chem_data2d,nx,ny,1)
                      do i=1,nx
                         do j=1,ny
                            if(xland(i,j).eq.1) then
                               chem_data2d(i,j)=(1.+pert_land_fire_d)*chem_data2d(i,j)
                            else
                               chem_data2d(i,j)=(1.+pert_watr_fire_d)*chem_data2d(i,j)
                            endif
                         enddo
                      enddo  
                      call put_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),chem_data2d,nx,ny,1)
                   enddo
!
! mirror member
                   imem=iimem+1
                   if(sw_gen) then
                      pert_land_fire=-1.*pert_land_fire
                      pert_watr_fire=-1.*pert_watr_fire
                   endif
                   if(imem.ge.0.and.imem.lt.10) write(cmem,"('.e00',i1)"),imem
                   if(imem.ge.10.and.imem.lt.100) write(cmem,"('.e0',i2)"),imem
                   if(imem.ge.100.and.imem.lt.1000) write(cmem,"('.e',i3)"),imem
                   wrffire_file=trim(wrffirechemi)//trim(cmem)
                   do isp=1,nfire_spc
                      call get_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),chem_data2d,nx,ny,1)
                      do i=1,nx
                         do j=1,ny
                            if(xland(i,j).eq.1) then
                               chem_data2d(i,j)=(1.+pert_land_fire_m)*chem_data2d(i,j)
                            else
                               chem_data2d(i,j)=(1.+pert_watr_fire_m)*chem_data2d(i,j)
                            endif
                         enddo
                      enddo  
                      call put_WRFCHEM_emiss_data(wrffire_file,ch_fire_spc(isp),chem_data2d,nx,ny,1)
                   enddo
                endif
!
                if(sw_biog) then
!
! draw member
                   imem=iimem
                   if(imem.ge.0.and.imem.lt.10) write(cmem,"('.e00',i1)"),imem
                   if(imem.ge.10.and.imem.lt.100) write(cmem,"('.e0',i2)"),imem
                   if(imem.ge.100.and.imem.lt.1000) write(cmem,"('.e',i3)"),imem
                   wrfbio_file=trim(wrfbiochemi)//trim(cmem)
                   do isp=1,nbio_spc
                      call get_WRFCHEM_emiss_data(wrfbio_file,ch_bio_spc(isp),chem_data2d,nx,ny,1)
                      do i=1,nx
                         do j=1,ny
                            if(xland(i,j).eq.1) then
                               chem_data2d(i,j)=(1.+pert_land_biog_d)*chem_data2d(i,j)
                            else
                               chem_data2d(i,j)=(1.+pert_watr_biog_d)*chem_data2d(i,j)
                            endif
                         enddo
                      enddo  
                      call put_WRFCHEM_emiss_data(wrfbio_file,ch_bio_spc(isp),chem_data2d,nx,ny,1)
                   enddo
!
! mirror member
                   imem=iimem+1
                   if(sw_gen) then
                      pert_land_biog=-1.*pert_land_biog
                      pert_watr_biog=-1.*pert_watr_biog
                   endif
                   if(imem.ge.0.and.imem.lt.10) write(cmem,"('.e00',i1)"),imem
                   if(imem.ge.10.and.imem.lt.100) write(cmem,"('.e0',i2)"),imem
                   if(imem.ge.100.and.imem.lt.1000) write(cmem,"('.e',i3)"),imem
                   wrfbio_file=trim(wrfbiochemi)//trim(cmem)
                   do isp=1,nbio_spc
                      call get_WRFCHEM_emiss_data(wrfbio_file,ch_bio_spc(isp),chem_data2d,nx,ny,1)
                      do i=1,nx
                         do j=1,ny
                            if(xland(i,j).eq.1) then
                               chem_data2d(i,j)=(1.+pert_land_biog_m)*chem_data2d(i,j)
                            else
                               chem_data2d(i,j)=(1.+pert_watr_biog_m)*chem_data2d(i,j)
                            endif
                         enddo
                      enddo  
                      call put_WRFCHEM_emiss_data(wrfbio_file,ch_bio_spc(isp),chem_data2d,nx,ny,1)
                   enddo
                endif
                pert_land_chem=-9999
                pert_watr_chem=-9999
                pert_land_fire=-9999
                pert_watr_fire=-9999
                pert_land_biog=-9999
                pert_watr_biog=-9999
             enddo
!
! Deallocate arrays
             deallocate(chem_data3d,chem_data2d,xland)
             deallocate(ch_chem_spc,ch_fire_spc,ch_bio_spc)
             close(unita)
          end program main
!
          subroutine get_WRFINPUT_land_mask(xland,nx,ny)
             implicit none
             include 'netcdf.inc'
             integer, parameter                    :: maxdim=6
             integer                               :: nx,ny
             integer                               :: i,rc
             integer                               :: f_id
             integer                               :: v_id,v_ndim,typ,natts
             integer,dimension(maxdim)             :: one
             integer,dimension(maxdim)             :: v_dimid
             integer,dimension(maxdim)             :: v_dim
             real,dimension(nx,ny)                 :: xland
             character(len=150)                    :: v_nam
             character*(80)                         :: name
             character*(80)                         :: file
!
! open netcdf file
             file='wrfinput_d02'
             name='XLAND'
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
             else if(1.ne.v_dim(3)) then             
                print *, 'ERROR: nz dimension conflict ','1',v_dim(3)
                call abort
!             else if(1.ne.v_dim(4)) then             
!                print *, 'ERROR: time dimension conflict ',1,v_dim(4)
!                call abort
             endif
!
! get data
             one(:)=1
             rc = nf_get_vara_real(f_id,v_id,one,v_dim,xland)
             if(rc.ne.0) then
                print *, 'nf_get_vara_real ', xland(1,1)
                call abort
             endif
             rc = nf_close(f_id)
             return
          end subroutine get_WRFINPUT_land_mask   
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




