! DART software - Copyright � 2004 - 2010 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download

program create_obs_grid

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

use    utilities_mod, only : register_module, open_file, close_file, &
                             initialize_utilities 
use obs_sequence_mod, only : obs_sequence_type, interactive_obs, write_obs_seq, &
                             static_init_obs_sequence
use  assim_model_mod, only : static_init_assim_model

implicit none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(obs_sequence_type) :: seq
character(len = 129)    :: file_name

! Record the current time, date, etc. to the logfile
call initialize_utilities('create_obs_grid')
call register_module(source,revision,revdate)

! Initialize the assim_model module, need this to get model
! state meta data for locations of identity observations
call static_init_assim_model()

! Initialize the obs_sequence module
call static_init_obs_sequence()

! Create grid of obs
seq = create_grid()

! Write the sequence to a file
write(*, *) 'Input filename for sequence (  obs_seq.in   usually works well)'
read(*, *) file_name
call write_obs_seq(seq, file_name)

! Clean up
call finalize_utilities('create_obs_grid')

contains

function create_grid()
 type(obs_sequence_type) :: create_grid

type(obs_type)     :: obs, prev_obs
type(obs_def_type) :: obs_def
type(time_type)    :: obs_time, prev_time
type(location_type) :: loc
integer            :: max_num_grids, num_copies, num_qc, end_it_all
integer            :: num_dim, ni, nj, nk, i, j, k, l

! these things aren't prompted for - they're fixed in the code
num_copies = 1
num_qc     = 0
max_num_obs = 1000000   ! FIXME: made up

write(*, *) 'Input upper bound on number of grids of observations in sequence'
read(*, *) max_num_grids

! Initialize an obs_sequence structure
call init_obs_sequence(create_grid, num_copies, num_qc, max_num_obs)

do i = 1, num_copies
   create_grid%copy_meta_data(i) = 'observations'
end do

! Initialize the obs variable
call init_obs(obs, num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)

! Loop to initialize each observation in turn; terminate by -1
do l = 1, max_num_grids
   write(*, *) 'input a -1 if there are no more grids'
   read(*, *) end_it_all
   if(end_it_all == -1) exit

   ! FIXME: this is the corner of a grid, need to prompt for
   ! extents in each dim (2d or 3d) and number of points in same.
   ! then loop below using same type and error, just bumping
   ! location each time.

   write(*, *) 'the location of the next observation defines the corner of a box'

   ! Need to have key available for specialized observation modules
   call interactive_obs(num_copies, num_qc, obs, i)

   write(*, *) 'enter the location of the opposite corner of the box'
   call interactive_location(loc)

   num_dim = -1
   while (num_dim < 1 .and. num_dim > 3) do
      write(*, *) 'input 1, 2 or 3 for 1d, 2d or 3d grid'
      read(*, *) num_dim
   done

   do i=1, num_dim
   enddo

   do k=1, nk
      do j=1, nj
         do i=1, ni
            ! set an obs based on the corner one
            ! compute new location and set it into an obs

            if(i == 1) then
               call insert_obs_in_seq(create_grid, obs)
            else
               ! if this is not the first obs, make sure the time is larger
               ! than the previous observation.  if so, we can start the
               ! linked list search at the location of the previous obs.
               ! otherwise, we have to start at the beginning of the entire
               ! sequence to be sure the obs are ordered correctly in
               ! monotonically increasing times. 
               call get_obs_def(obs, obs_def)
               obs_time = get_obs_def_time(obs_def)
               call get_obs_def(prev_obs, obs_def)
               prev_time = get_obs_def_time(obs_def)
               if(prev_time > obs_time) then
                  call insert_obs_in_seq(create_grid, obs)
               else
                  call insert_obs_in_seq(create_grid, obs, prev_obs)
               endif
            endif
            prev_obs = obs
         end do    ! i
      end do    ! j
   end do    ! k
end do    ! max_grids

call destroy_obs(obs)
call destroy_obs(prev_obs)

end function create_grid

end program create_obs_grid
