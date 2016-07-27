#!/usr/bin/perl
#
# DART software - Copyright 2004 - 2011 UCAR. This open source software is
# provided by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download
#
# $Id$
#
# AUTHOR:   T. R. Whitcomb
#           Naval Research Laboratory
#
# Take the file generated by the parser that contains information
# on the field name, the dimension type, the numerical position
# within that dimension type, and the type of I/O (single or 
# multiple processor) and generate a Fortran routine that will:
#  1. Take the input of a name
#  2. Search through the fields and return:
#      a) an integer dimension type
#      b) a two-dimensional integer array containing the record
#         number for both the single and the multiple processor
#         I/O cases
# 
# This just creates a separate subroutine file - use an INTERFACE
# block in the module for calling.  The command line argument
# gives the file containing the properly formatted field 
# information - the easiest way to do this is have the output of
# the parse_field_list.pl script.
######

$DIM_TYPE = 'dim_type';

$fields_file = $ARGV[0];

if (! -e $fields_file)
{
    die("File $fields_file does not exist!\n");
}

# First, just read in the data and populate our local array
open(FIELD_LIST,"<$fields_file");
while (<FIELD_LIST>)
{
        chomp;
    @field_data = split;

    $variable_name = $field_data[0];
    $dimension_type = $field_data[1];
    $variable_index = $field_data[2];
    $iotype         = $field_data[3];

    $fields{$variable_name}{$DIM_TYPE} = $dimension_type;
    $fields{$variable_name}{$iotype}    = $variable_index 
}

close(FIELD_LIST);

open(ROUTINE_FILE,">get_name_info.f90");

# Write function header
print ROUTINE_FILE <<END_HEADER;
! get_name_info
! -------------
! Given the value of module-specific constant values and a name of 
! a variable in the COAMPS restart file, returns an integer 
! representing the dimension type and a two-dimensional integer 
! array containing the position of that variable in that particular 
! dimension set in the COAMPS restart file written by either 
! single-processor or multi-processor I/O. 
!  PARAMETERS
!   IN  DIM_TYPE_2D         Constant for 2-D dimension
!   IN  DIM_TYPE_3D         Constant for 3-D dimension
!   IN  DIM_TYPE_3DW        Constant for 3-D (w level) dimension
!   IN  SINGLEIO            Constant for single-processor I/O
!   IN  MULTIIO             Constant for multi-processor I/O
!   IN  var_name            the name of the variable to look up
!   OUT var_dim_type        the dimension type (2d/3d/3dw)
!                           1: 2D  2: 3D  3: 3DW
!   OUT var_record_num      the position of the variable in its
!                           particular dimension
subroutine get_name_info(DIM_TYPE_2D, DIM_TYPE_3D, DIM_TYPE_3DW,   &
                         SINGLEIO, MULTIIO, var_name, var_dim_type,&
                         var_record_num)
  integer, intent(in)                :: DIM_TYPE_2D
  integer, intent(in)                :: DIM_TYPE_3D
  integer, intent(in)                :: DIM_TYPE_3DW
  integer, intent(in)                :: SINGLEIO
  integer, intent(in)                :: MULTIIO
  character(len=*), intent(in)       :: var_name
  integer, intent(out)               :: var_dim_type
  integer, dimension(2), intent(out) :: var_record_num

  select case (trim(var_name))
END_HEADER

# Get the actual data
foreach $field_name (sort keys %fields)
{
    $var_dim_type = $fields{$field_name}{$DIM_TYPE};
    $single_io_field = $fields{$field_name}{'SINGLEIO'};
    if (!$single_io_field) {$single_io_field = -1;}
    $multiple_io_field = $fields{$field_name}{'MULTIIO'};
    if (!$multiple_io_field) {$multiple_io_field = -1;}
    print "\t$field_name\t$var_dim_type\t$single_io_field";
    print "\t$multiple_io_field\n";

    print ROUTINE_FILE "  case('$field_name')\n";
    print ROUTINE_FILE "    var_dim_type = $var_dim_type\n"; 
    print ROUTINE_FILE "    var_record_num(SINGLEIO) = $single_io_field\n";
    print ROUTINE_FILE "    var_record_num(MULTIIO)  = $multiple_io_field\n";
}

print ROUTINE_FILE <<END_FOOTER;
  case default
    write (*,*) "Can't match name " // var_name // "in restart"
  end select
end subroutine get_name_info
END_FOOTER
close(ROUTINE_FILE);

exit 0

# <next few lines under version control, do not edit>
# $URL$
# $Revision$
# $Date$

