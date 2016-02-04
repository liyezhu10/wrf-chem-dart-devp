! DART software - Copyright 2004 - 2016 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module cross_comp_localization_mod

! Computes information necessary to determine the alpha associated with a 
! cross interface incremental update in a three dimensional spherical shell. 


use      types_mod, only : r8, MISSING_R8, MISSING_I, PI, RAD2DEG, DEG2RAD
use  utilities_mod, only : register_module, error_handler, E_ERR, ascii_file_format, &
                           nc_check, E_MSG, open_file, close_file, set_output,       &
                           logfileunit, nmlfileunit, find_namelist_in_file,          &
                           check_namelist_read, do_output, do_nml_file,              &
                           do_nml_term, is_longitude_between
use random_seq_mod, only : random_seq_type, init_random_seq, random_uniform
use   obs_kind_mod, only : get_num_obs_kinds, get_obs_kind_name
use mpi_utilities_mod, only : my_task_id, task_count
use location_mod
use cov_cutoff_mod, only: comp_cov_factor

implicit none
private

public :: get_vert_alpha, get_vert_dist, get_xcomp_dist 

character(len=256), parameter :: source   = &
   "$URL$"
character(len=32 ), parameter :: revision = "$Revision$"
character(len=128), parameter :: revdate  = "$Date$"

logical, save         :: module_initialized = .false.

character(len = 512) :: errstring

! Global storage for vertical distance normalization factors
real(r8)              :: vert_normalization(4)

character(len = 512) :: metafile                    = "xcomp-sec-tatmtocn_metafile.txt"
character(len = 512) :: atm_lookup_table            = "xcomp-sec-tatmtocn_atmfile.txt"
character(len = 512) :: ocn_lookup_table            = "xcomp-sec-tatmtocn_ocnfile.txt"

logical  :: cross_component_updates           = .false.
real(r8)     :: cutoff_cross = 0.2 
real(r8)     :: cutoff = .1

!nlat, nlon,nlev are the number of grid points in the lookup table
!these will be set from the "metafile", along with lon_grid, lat_grid,lev_grid

!A(nlon, nlat,nlev) holds the atmospheric alpha at levels from 500-1000mb
!note that A is defined in *mb*
! At this point A and L must be on regular grids.
!L holds the latitude dependent ocean-depth alpha correction (called alpha2)
!it is a length scale such that alpha2 = exp(-(d/L)^2)

integer  :: nlat,nlon,nlev
real(r8), allocatable :: lat_grid(:),lon_grid(:),lev_grid(:)
real(r8), allocatable :: A(:,:,:), L(:)
real(r8)              :: lat_min, lon_min, lev_min, lat_dx,lon_dx,lev_dx

namelist /cross_comp_localization_nml/ cross_component_updates, metafile, &
                                       atm_lookup_table, ocn_lookup_table,&
                                       cutoff_cross               
 
contains

!----------------------------------------------------------------------------

subroutine initialize_module
 
! things which need doing exactly once.

integer :: iunit, io, i,j,k
character(len=129) :: str1

if (module_initialized) return

call register_module(source, revision, revdate)
module_initialized = .true.

! Read the namelist entry
call find_namelist_in_file("input.nml", "cross_comp_localization_nml", iunit)
read(iunit, nml = cross_comp_localization_nml, iostat = io)
call check_namelist_read(iunit, io, "cross_comp_localization_nml")

! Write the namelist values to the log file

if(do_nml_file()) write(nmlfileunit, nml=cross_comp_localization_nml)
if(do_nml_term()) write(     *     , nml=cross_comp_localization_nml)


! read in the metadata 
open (1, file = metafile, status = 'old')
READ (1, *) nlon,nlat,nlev
!print *,nlat,nlon,nlev

!allocate the arrays
allocate(lat_grid(nlat))
allocate(lon_grid(nlon))
allocate(lev_grid(nlev))

allocate(A(nlon,nlat,nlev))
allocate(L(nlat))

!read in the grids (at this point they must be regular!)
READ (1, *) lon_grid
READ (1, *) lat_grid
READ (1, *) lev_grid

!compute the min and delta from each grid
lat_min = minval(lat_grid)
lon_min = minval(lon_grid)
lev_min = minval(lev_grid)

lat_dx = lat_grid(2)-lat_grid(1)
lon_dx = lon_grid(2)-lon_grid(1)
lev_dx = lev_grid(2)-lev_grid(1)

!open the lookup tables 
open (2, file = atm_lookup_table, status = 'old')
open (3, file = ocn_lookup_table, status = 'old')


!read in A and L
do k=1,nlev
do j=1,nlat
do i=1,nlon
READ (2, *) A(i,j,k)
enddo
enddo
enddo

READ (3, *) L

end subroutine initialize_module



! ---------------------------------------------------------------------------


function get_vert_alpha(loc_obs,loc_state,kind_obs,kind_state)

!this function will return a cross-component (ocean/atm) localization coefficient for the vertical ONLY!
!this alpha must be combined with the horizontal distances to get the final localization (unless 1d mod!)
!this is not a sampling error correction... it is completely apriori specified
!based only on the location of the observation and state

! the algorithm works as follows:
!  use the observation to specify lat and lon
!  if observation is atm (query using 'vert') then
                !atm depth from obs (in mb) is z_atm
                ! ocean depth from state (in meters) is z_ocn
!else if observation is ocean  
                !atm depth from state
                !ocean depth from obs
!alpha = A(lon,lat,z_atm)*exp(-(z_ocn/L(lat).^2))
      

!there is an asymmetry in alpha from ocean-atm and atm-ocean since
!lat and lon are define always from the obs side.
!in the future we could imagine defining a mid-way point lat and lon
                  


type(location_type), intent(in) :: loc_obs, loc_state
integer, optional,   intent(in) :: kind_obs,kind_state
real(r8)                        :: get_vert_alpha
real(r8)                        :: this_loc_obs(3),this_vert_obs
real(r8)                        :: this_loc_state(3),this_vert_state
integer                         :: ip,ilat,ilon
real(r8)                        :: alpha, alpha2
real(r8)                        :: press_in_mb

if (module_initialized == .false.) then
call initialize_module()
endif

this_loc_obs = get_location(loc_obs)
this_vert_obs = query_location(loc_obs, 'which_vert')

this_loc_state = get_location(loc_state)
this_vert_state = query_location(loc_state, 'which_vert')
 
! tester stub
!print *, 'location obs  is ', get_location(loc_obs)
!print *, 'obs kind  is ', kind_obs
!print *, 'obs vert is' , this_vert_obs

!print *, 'location state is ', get_location(loc_state)
!print *, 'state kind is ', kind_state
!print *, 'location vert is' , this_vert_state


!LON
ilon=ceiling( (this_loc_obs(1) + lon_dx/2.0 - lon_min)/lon_dx)
ilon=min0(ilon,nlon);
ilon=max0(ilon,1);

!LAT
ilat=ceiling( (this_loc_obs(2) + lat_dx/2.0 - lat_min)/lat_dx)
ilat=min0(ilat,nlat);
ilat=max0(ilat,1);

!LE
if ( vert_is_pressure(loc_obs) ) then
press_in_mb = this_loc_obs(3)*.01 !convert Pa to mb
        ip=ceiling( (press_in_mb + lev_dx/2.0 - lev_min)/lev_dx)
        ip=min0(ip,nlev)
        ip=max0(ip,1)
        if (press_in_mb .lt.  500.0) then
                alpha = 0.0
        else
                alpha = A(ilon,ilat,ip)
        endif
!        print*,'OBS HAS THE PRESSURE COORDINATE'
elseif ( vert_is_pressure(loc_state) ) then
press_in_mb = this_loc_state(3)*.01 !convert Pa to mb
        ip=ceiling( (press_in_mb + lev_dx/2.0 - lev_min)/lev_dx)
        ip=min0(ip,nlev)
        ip=max0(ip,1)
        if (press_in_mb .lt.  500.0) then
                alpha = 0.0
        else
                alpha = A(ilon,ilat,ip)
        endif
!        print*,'STATE HAS THE PRESSURE COORDINATE'
else
        print*,'ERROR,EITHER OBS OR STATE MUST HAVE A PRESSURE COORDINATE'
        alpha = 0.0
endif 

!DEPTH
if ( vert_is_height(loc_obs) ) then
     alpha2 = exp(-(this_loc_obs(3)/L(ilat))**2)
!     print*,'OBS HAS THE DEPTH COORDINATE'
else if ( vert_is_height(loc_state) ) then
     alpha2 = exp(-(this_loc_state(3)/L(ilat))**2)
!     print*,'STATE HAS THE DEPTH COORDINATE'
else
print*,'ERROR,EITHER OBS OR STATE MUST HAVE A HEIGHT COORDINATE'
     alpha2= 0.0
endif


get_vert_alpha = alpha*alpha2


end function get_vert_alpha

!----------------------------------------------------------------------------



function get_vert_dist(alpha,cutoff_override)

!this function will return a distance (in radians)  that is implied by an alpha and a cutoff
!if not present, the cutoff is read from the namelist
!in the future, this cutoff should be read from the "right" namelist (assim)

real(r8),optional, intent(in)  :: cutoff_override
real(r8)                        :: alpha,cutoff
real(r8)                        :: get_vert_dist


!we are actually computing an approximate distance assumpint that the GC
!function is like:
!exp(-(d/(c/1.25))**2)

if (module_initialized == .false.) then
call initialize_module()
endif

get_vert_dist = 0.0

if (present(cutoff_override)) then 
cutoff = cutoff_override
else
cutoff = cutoff_cross
endif

get_vert_dist = cutoff/1.25*(-log(alpha))**.5

end function get_vert_dist


function get_xcomp_dist(loc1,loc2,kind1,kind2)

!this function will return a distance that sums the vertical and horizontal

type(location_type), intent(in) :: loc2, loc1
integer, optional,   intent(in) :: kind1,kind2
real(r8)                        :: alpha,d1,d2,d_vert,d_horiz
real(r8)                        :: get_xcomp_dist


if (module_initialized == .false.) then
call initialize_module()
endif


get_xcomp_dist = 0.0

!VERTICAL
!get the vertical alpha 1
alpha = get_vert_alpha(loc1,loc2,kind1,kind2)
!get corresponding distance
d1 = get_vert_dist(alpha)

!get the vertical alpha 2
alpha = get_vert_alpha(loc2,loc1,kind2,kind1)
!get corresponding distance
d2 = get_vert_dist(alpha)

!get average distance vertical
d_vert = 0.5*(d1 + d2)

!HORIZONTAL
d_horiz=get_dist(loc1,loc2,0,0,.true.)

!TOTAL
get_xcomp_dist  = (d_vert**2  + d_horiz**2)**0.5

end function get_xcomp_dist



end module cross_comp_localization_mod




! <next few lines under version control, do not edit>
! $URL: https://subversion.ucar.edu/DAReS/DART/trunk/location/threed_sphere/$
! $Id: cross_comp_localization_mod.f90 alicia K$
! $Revision$
! $Date: 2016-01-29$
