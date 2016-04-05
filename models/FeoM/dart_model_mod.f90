! DART software - Copyright 2004 - 2013 UCAR. This open source software is
! provided by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

module dart_model_mod

use netcdf

implicit none
private

public :: get_model_nLons, get_model_nLats, get_model_nAlts, &
                             get_nSpecies, get_nSpeciesTotal, get_nIons,     &
                             get_nSpeciesAll, decode_model_indices

!======================================================================
contains
!======================================================================


function get_model_nLons()
!------------------------------------------------------------------
integer :: get_model_nLons

get_model_nLons = 70270

end function get_model_nLons


function get_model_nLats()
!------------------------------------------------------------------
integer :: get_model_nLats

get_model_nLats = 70270

end function get_model_nLats

function get_model_nAlts()
!------------------------------------------------------------------
integer :: get_model_nAlts

get_model_nAlts = 110

end function get_model_nAlts

function get_nSpecies()
!------------------------------------------------------------------
integer :: get_nSpecies

get_nSpecies = 1

end function get_nSpecies


function get_nSpeciesTotal()
!------------------------------------------------------------------
integer :: get_nSpeciesTotal

get_nSpeciesTotal = -1

end function get_nSpeciesTotal

function get_nIons()
!------------------------------------------------------------------
integer :: get_nIons

get_nIons = -1

end function get_nIons


function get_nSpeciesAll()
!------------------------------------------------------------------
integer :: get_nSpeciesAll

get_nSpeciesAll = -1

end function get_nSpeciesAll

subroutine decode_model_indices (varname, model_varname, model_dim, &
                             model_index, long_name, units)

   character(len=NF90_MAX_NAME) :: varname       ! crazy species name
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=NF90_MAX_NAME) :: model_varname  ! NDensityS, IDensityS, ...
   integer :: model_dim                           ! dimension defining species
   integer :: model_index                         ! 'iSpecies' or u,v,w ...

   !code goes here
   
end subroutine

end module dart_model_mod

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
