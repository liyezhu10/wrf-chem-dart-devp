MODULE cosmo_data_mod

! This module contains the routines to read or write
! COSMO grib files, convert the COSMO data to DART and
! vice versa
!
! Author: Jan D. Keller
!         Meteorological Institute, University of Bonn, Germany
!         2011-09-15
!
! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$

  use        types_mod, only : r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                               rad2deg, deg2rad, PI

  use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                               print_time, print_date, set_calendar_type,        &
                               operator(*),  operator(+), operator(-),           &
                               operator(>),  operator(<), operator(/),           &
                               operator(/=), operator(<=),                       &
                               julian_day,GREGORIAN,JULIAN


  use     obs_kind_mod, only : KIND_U_WIND_COMPONENT,                            &
                               KIND_V_WIND_COMPONENT,                            &
                               KIND_VERTICAL_VELOCITY,                           &
                               KIND_TEMPERATURE,                                 &
                               KIND_PRESSURE,                                    &
                               KIND_SPECIFIC_HUMIDITY,                           &
                               KIND_CLOUD_LIQUID_WATER,                          &
                               KIND_CLOUD_ICE,                                   &
                               KIND_SURFACE_ELEVATION,                           &
                               KIND_SURFACE_GEOPOTENTIAL

  USE          byte_mod, ONLY: concat_bytes1,concat_bytes1_sign,to_positive,&
                               byte_to_word_signed,to_float1,&
                               byte_to_word_data,get_word

  USE     grib_info_mod, ONLY: get_level,get_varname,get_dims,get_dart_kind

  IMPLICIT none

! version controlled file description for error handling, do not edit
character(len=128), parameter :: &
   source   = "$URL: ", &
   revision = "$Revision$", &
   revdate  = "$Date$"

  type cosmo_meta
    character(len=16)  :: varname_short
    character(len=256) :: varname_long
    character(len=32)  :: units
    integer            :: grib_file_position
    integer            :: grib_table
    integer            :: grib_code
    integer            :: grib_leveltype
    integer            :: grib_level(2)
    integer            :: grib_nlevels
    integer            :: ilevel
    integer            :: numdims
    integer            :: dims(3)
    integer            :: hcoord_type
    integer            :: dart_sindex
    integer            :: dart_eindex
    integer            :: dart_kind
    real(r8)           :: dart_level
  end type cosmo_meta

  type cosmo_non_state_data
    real(r8),allocatable :: surface_orography(:,:)
    logical              :: vertical_coords_set
    integer              :: nfl
    integer              :: nhl
    real(r8),allocatable :: vct_a(:)
    real(r8),allocatable :: vct_b(:)
    real(r8),allocatable :: hhl(:,:,:)
    real(r8),allocatable :: hfl(:,:,:)
    real(r8),allocatable :: phl(:,:,:)
    real(r8),allocatable :: pfl(:,:,:)
    real(r8),allocatable :: p0hl(:,:,:)
    real(r8),allocatable :: p0fl(:,:,:)
    real(r8)             :: p0sl
    real(r8)             :: t0sl
    real(r8)             :: dt0lp
    real(r8)             :: vcflat
  end type cosmo_non_state_data

  type cosmo_hcoord
    real(r8),allocatable :: lon(:)
    real(r8),allocatable :: lat(:)
  end type cosmo_hcoord

  type grib_header_type
    integer                     :: start_record
    integer                     :: start_offset
    integer                     :: data_record
    integer                     :: data_offset
    integer                     :: data_length
    integer(kind=1),allocatable :: pds(:)
    integer(kind=1),allocatable :: gds(:)
    integer(kind=1),allocatable :: bds(:)
  end type grib_header_type
  
  INTEGER :: nmaxrec,nmaxvar,nhcoord
  LOGICAL :: hcoord_read(3)

  PARAMETER(nmaxrec=1000000000)
  PARAMETER(nmaxvar=1000)
  PARAMETER(nhcoord=3)

  public  :: get_cosmo_info
  public  :: set_vertical_coords
  private :: read_grib_data
  private :: find_grib_records
  private :: read_grib_header
  private :: find_in_varlist
  private :: read_variable_info
  private :: read_hcoords
  private :: read_lonlat
  private :: read_time

CONTAINS

  SUBROUTINE get_cosmo_info(filename,v,h,header,allowed,tf)

    ! get_cosmo_info_and_data calls subroutines to read the data, extract the desired information
    ! and put the data into the state vector

    CHARACTER(len=256),INTENT(in)            :: filename
    LOGICAL,INTENT(in)                       :: allowed(:)
    TYPE(cosmo_meta),ALLOCATABLE,INTENT(out) :: v(:)
    TYPE(cosmo_hcoord),INTENT(out)           :: h(3)
    TYPE(grib_header_type),ALLOCATABLE,INTENT(out) :: header(:)

!    INTEGER(kind=1),ALLOCATABLE,INTENT(out)  :: bdata(:)
!    INTEGER,ALLOCATABLE,INTENT(out)          :: bytepos(:,:) ! byte positions for  (1,:) beginning of GRIB record
                                                             !                     (2,:) beginning of PDS section
                                                             !                     (3,:) beginning of GDS section
                                                             !                     (4,:) beginning of data section
!    INTEGER,ALLOCATABLE,INTENT(out)          :: bytelen(:,:) ! byte length for the (2,:) PDS section
                                                             !                     (3,:) GDS section
                                                             !                     (4,:) data section
    TYPE(time_type),INTENT(out)              :: tf

    INTEGER                                  :: nbyte,nvar

    CALL read_grib_header(filename,header,nvar)

!    CALL read_grib_data(filename,bdata,nbyte)
!    CALL find_grib_records(bdata,nbyte,nvar,bytepos,bytelen)

    ALLOCATE(v(1:nvar))

    CALL read_variable_info(header,nvar,v,allowed)
    CALL read_hcoords(header,nvar,v,h)
    CALL read_time(header(1),tf)

!    CALL read_variable_info(bdata,nvar,bytepos,bytelen,v,allowed)
!    CALL read_hcoords(bdata,nvar,bytepos,v,h)
!    CALL read_time(bdata,nvar,bytepos,tf)

  END SUBROUTINE get_cosmo_info

  SUBROUTINE read_grib_data(filename,bdata,nbyte)
    
    IMPLICIT NONE

    CHARACTER(len=256),INTENT(in)           :: filename
    INTEGER(kind=1),INTENT(out),ALLOCATABLE :: bdata(:)
    INTEGER,INTENT(out)                     :: nbyte

    INTEGER(kind=1)                         :: temp(1:nmaxrec)
    INTEGER                                 :: ibyte,irec,istat,nrec
    INTEGER(kind=4)                         :: word(4)
    INTEGER(kind=1)                         :: bin1(4)
    
    OPEN(10,FILE=TRIM(filename),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)

    ibyte=1
    irec=1
    istat=0
    DO WHILE(istat==0)
      READ(10,rec=irec,iostat=istat) bin1
      word=bin1
      irec=irec+1
      temp(ibyte:ibyte+3)=word
      ibyte=ibyte+4
    END DO
    CLOSE(10)

    nrec=irec-1
    ALLOCATE(bdata(1:nrec*4))
    bdata(1:nrec*4)=temp(1:nrec*4)
    nbyte=nrec*4

    RETURN
  END SUBROUTINE read_grib_data

  SUBROUTINE find_grib_records(bdata,nbyte,nvar,bytepos,bytelen)
    
    INTEGER,INTENT(in)              :: nbyte
    INTEGER(kind=1),INTENT(in)      :: bdata(1:nbyte)

    INTEGER,INTENT(out)             :: nvar
    INTEGER,ALLOCATABLE,INTENT(out) :: bytepos(:,:)
    INTEGER,ALLOCATABLE,INTENT(out) :: bytelen(:,:)

    INTEGER(kind=4)                 :: grib
    INTEGER(kind=1)                 :: gribword(4)
    INTEGER                         :: ibyte,ivar,len
    INTEGER                         :: temppos(1:4,1:nmaxvar),templen(1:4,1:nmaxvar)

    grib=ICHAR('G')*256**0+ICHAR('R')*256**1+ICHAR('I')*256**2+ICHAR('B')*256**3
    gribword(1)=ICHAR('G')
    gribword(2)=ICHAR('R')
    gribword(3)=ICHAR('I')
    gribword(4)=ICHAR('B')

    ivar=0
    ibyte=1
    DO WHILE (ibyte .LE. nbyte-4)
      IF (get_word(gribword).EQ.get_word(bdata(ibyte:ibyte+3))) THEN
        ivar=ivar+1
        temppos(1,ivar)=ibyte
        ibyte=ibyte+8
        temppos(2,ivar)=ibyte
        len=concat_bytes1(bdata(ibyte:ibyte+2),3,.TRUE.)
        templen(2,ivar)=len
        ibyte=ibyte+len
        temppos(3,ivar)=ibyte
        len=concat_bytes1(bdata(ibyte:ibyte+2),3,.TRUE.)
        templen(3,ivar)=len
        ibyte=ibyte+len
        temppos(4,ivar)=ibyte
        len=concat_bytes1(bdata(ibyte:ibyte+2),3,.TRUE.)
        templen(4,ivar)=len
        ibyte=ibyte+len
      ELSE
        ibyte=ibyte+1
      END IF
    END DO
    nvar=ivar

    ALLOCATE(bytepos(1:4,1:nvar))
    ALLOCATE(bytelen(1:4,1:nvar))
    bytepos=temppos(1:4,1:nvar)
    bytelen=templen(1:4,1:nvar)

    RETURN
  END SUBROUTINE find_grib_records

  SUBROUTINE read_grib_header(filename,header,nvar)
    
    CHARACTER(len=256),                INTENT(in)  :: filename
    TYPE(grib_header_type),ALLOCATABLE,INTENT(out) :: header(:)
    INTEGER,                           INTENT(out) :: nvar

    TYPE(grib_header_type)                         :: temp(1:nmaxvar)
    INTEGER                                           :: irec,istat,nrec,irec2,ibyte,m
    INTEGER(kind=4)                                   :: word(4)
    INTEGER(kind=1)                                   :: bin4(4),bin1,bin(4),bin8(8)
    INTEGER(kind=1),ALLOCATABLE                       :: pds(:),gds(:)

    INTEGER(kind=1)                 :: gribword(4)
    INTEGER(kind=1),ALLOCATABLE     :: bytearr(:)
    INTEGER                         :: ivar,mylen,pdslen,gdslen,datalen,bdslen

    INTEGER                         :: grib_start,data_start,ioff
    LOGICAL                         :: foundoff

    OPEN(10,FILE=TRIM(filename),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)

    gribword(1)=ICHAR('G')
    gribword(2)=ICHAR('R')
    gribword(3)=ICHAR('I')
    gribword(4)=ICHAR('B')

    ivar=1
    irec=1
    istat=0

    DO WHILE(istat==0)
      READ(10,rec=irec,iostat=istat) bin4
      bin8(1:4)=bin4
      READ(10,rec=irec+1,iostat=istat) bin4
      bin8(5:8)=bin4

      foundoff=.FALSE.
      findoff : DO ioff=0,3
        bin4(1:4)=bin8(1+ioff:4+ioff)
        if (bin4(1)==gribword(1) .AND. bin4(2)==gribword(2) .AND. &
            bin4(3)==gribword(3) .AND. bin4(4)==gribword(4)) THEN

          temp(ivar)%start_record=irec
          temp(ivar)%start_offset=ioff

          foundoff=.TRUE.
          exit findoff
        end if
      end DO findoff
      
      IF (foundoff) THEN
        grib_start=irec

        irec=irec+2

        READ(10,rec=irec,iostat=istat) bin4
        bin8(1:4)=bin4
        READ(10,rec=irec+1,iostat=istat) bin4
        bin8(5:8)=bin4
        
        pdslen=concat_bytes1(bin8(1+ioff:3+ioff),3,.TRUE.)

        m=MOD(pdslen+ioff,4)

        READ(10,rec=irec+(pdslen+ioff)/4,iostat=istat) bin4
        IF (m==0) THEN
          bin(1:4)=bin4(1:4)
        ELSE
          bin(1:m)=bin4(m+1:4)
        END IF
        READ(10,rec=irec+(pdslen+ioff)/4+1,iostat=istat) bin4
        if (m>0) THEN
          bin(m+1:4)=bin4(1:m)
        END if
        gdslen=concat_bytes1(bin(1:3),3,.TRUE.)

        bdslen=11

        ALLOCATE(bytearr(1:pdslen+gdslen+bdslen+ioff+8+4))

        ibyte=1
        DO irec2=irec,irec+CEILING((pdslen+gdslen+bdslen+ioff)/4.)+1
          READ(10,rec=irec2,iostat=istat) bin4
          bytearr(ibyte:ibyte+3)=bin4
          ibyte=ibyte+4
        END DO
        bytearr(1:pdslen+gdslen+bdslen+8-ioff)=bytearr(1+ioff:pdslen+gdslen+bdslen+8)

        ALLOCATE(temp(ivar)%pds(1:pdslen))
        ALLOCATE(temp(ivar)%gds(1:gdslen))
        ALLOCATE(temp(ivar)%bds(1:bdslen))

        temp(ivar)%pds(1:pdslen)=bytearr(1:pdslen)
        temp(ivar)%gds(1:gdslen)=bytearr(pdslen+1:pdslen+gdslen)
        temp(ivar)%bds(1:bdslen)=bytearr(pdslen+gdslen+1:pdslen+gdslen+bdslen)

        irec2=irec+FLOOR((pdslen+gdslen+ioff)/4.)
        temp(ivar)%data_record=irec2
        temp(ivar)%data_offset=MOD(pdslen+gdslen+ioff,4)
        
        mylen=concat_bytes1(bytearr(pdslen+gdslen+1:pdslen+gdslen+3),3,.TRUE.)
        temp(ivar)%data_length=mylen

        DEALLOCATE(bytearr)

        irec=irec+FLOOR((pdslen+gdslen+mylen)/4.)

        ivar=ivar+1
      else
        irec=irec+1
      end if
    END DO
    CLOSE(10)

    nvar=ivar-1

    ALLOCATE (header(1:nvar))
    DO ivar=1,nvar
      ALLOCATE(header(ivar)%pds(1:SIZE(temp(ivar)%pds)))
      ALLOCATE(header(ivar)%gds(1:SIZE(temp(ivar)%gds)))
      ALLOCATE(header(ivar)%bds(1:SIZE(temp(ivar)%bds)))
      header(ivar)%start_record=temp(ivar)%start_record
      header(ivar)%data_record=temp(ivar)%data_record
      header(ivar)%data_length=temp(ivar)%data_length
      header(ivar)%data_offset=temp(ivar)%data_offset
      header(ivar)%pds=temp(ivar)%pds
      header(ivar)%gds=temp(ivar)%gds
      header(ivar)%bds=temp(ivar)%bds
    END DO

    RETURN
  END SUBROUTINE read_grib_header

  FUNCTION find_in_varlist(tab,code,ltype,varlist,nlist)
    INTEGER            :: find_in_varlist

    INTEGER,INTENT(in) :: tab,code,ltype,nlist
    INTEGER,INTENT(in) :: varlist(1:nmaxvar,1:4)
    
    INTEGER            :: i
    
    find_in_varlist=nlist+1
    DO i=1,nlist
      IF (tab==varlist(i,1) .AND. code==varlist(i,2) .AND. ltype==varlist(i,3)) THEN
        find_in_varlist=i
      END IF
    END DO

    RETURN
  END FUNCTION find_in_varlist

  SUBROUTINE read_variable_info(header,nvar,v,allowed)
    
    TYPE(grib_header_type),INTENT(in)  :: header(1:nvar)
    INTEGER,               INTENT(in)  :: nvar
    LOGICAL,               INTENT(in)  :: allowed(:)
    TYPE(cosmo_meta),      INTENT(out) :: v(1:nvar)

    INTEGER                      :: ivar,ibyte,level(1:3),index,nlevels(200)
    CHARACTER(len=256)           :: varname(3)

    index=0
    nlevels(:)=0

    DO ivar=1,nvar
      v(ivar)%grib_file_position =ivar
      v(ivar)%grib_table         =to_positive(header(ivar)%pds(4))
      v(ivar)%grib_code          =to_positive(header(ivar)%pds(9))
      v(ivar)%grib_leveltype     =to_positive(header(ivar)%pds(10))

      level=get_level(v(ivar)%grib_leveltype,header(ivar)%pds(11:12))
      v(ivar)%grib_nlevels       =level(3)
      v(ivar)%grib_level         =level(1:2)
      IF (v(ivar)%grib_nlevels==1) THEN
        v(ivar)%dart_level       =level(1)
      ELSE
        v(ivar)%dart_level       =SUM(level(1:2))/2.
      END IF
      varname=get_varname(v(ivar)%grib_table,v(ivar)%grib_code,v(ivar)%grib_leveltype)
      v(ivar)%varname_short=varname(2)
      v(ivar)%varname_long=varname(1)
      v(ivar)%units=varname(3)
      v(ivar)%dart_kind=get_dart_kind(v(ivar)%grib_table,v(ivar)%grib_code,v(ivar)%grib_leveltype)
      v(ivar)%dims=get_dims(header(ivar)%gds(1:10))

      if (allowed(v(ivar)%dart_kind)) then
        v(ivar)%dart_sindex=index+1
        v(ivar)%dart_eindex=index+(v(ivar)%dims(1)*v(ivar)%dims(2))
        index=v(ivar)%dart_eindex
        v(ivar)%ilevel=nlevels(v(ivar)%dart_kind)+1
        nlevels(v(ivar)%dart_kind)=nlevels(v(ivar)%dart_kind)+1
      ELSE
        v(ivar)%dart_sindex=-1
        v(ivar)%dart_eindex=-1
        v(ivar)%ilevel=-1
      END if

      v(ivar)%hcoord_type=1
      IF (v(ivar)%dart_kind==KIND_U_WIND_COMPONENT) v(ivar)%hcoord_type=2
      IF (v(ivar)%dart_kind==KIND_V_WIND_COMPONENT) v(ivar)%hcoord_type=3

!      print*,TRIM(v(ivar)%varname_short)," - HCOORD_TYPE: ",v(ivar)%hcoord_type
!      print*,TRIM(v(ivar)%varname_short)," - DART_KIND: ",v(ivar)%dart_kind," - DART_SINDEX: ",v(ivar)%dart_sindex
!      print*,v(ivar)%grib_table,v(ivar)%grib_code,v(ivar)%grib_leveltype,v(ivar)%grib_level(1:v(ivar)%grib_nlevels),v(ivar)%dart_level
!      print*,v(ivar)%dims
    END DO

    RETURN

  END SUBROUTINE read_variable_info

  SUBROUTINE read_hcoords(header,nvar,v,h)
    
    TYPE(grib_header_type),INTENT(in)   :: header(1:nvar)
    INTEGER,INTENT(in)             :: nvar
    TYPE(cosmo_meta),INTENT(in)    :: v(1:nvar)
    TYPE(cosmo_hcoord),INTENT(out) :: h(1:3)

    INTEGER                        :: ivar,nx,ny
    REAL(r8),ALLOCATABLE           :: lons(:,:),lats(:,:)
    
    hcoord_read(:)=.FALSE.

    DO ivar=1,nvar
      IF (.NOT. hcoord_read(v(ivar)%hcoord_type)) THEN
        nx=v(ivar)%dims(1)
        ny=v(ivar)%dims(2)
        ALLOCATE(lons(1:nx,1:ny))
        ALLOCATE(lats(1:nx,1:ny))
        ALLOCATE(h(v(ivar)%hcoord_type)%lon(1:nx*ny))
        ALLOCATE(h(v(ivar)%hcoord_type)%lat(1:nx*ny))
        CALL read_lonlat(header(ivar)%gds,lons,lats,v(ivar)%dims(1),v(ivar)%dims(2))
        h(v(ivar)%hcoord_type)%lon(1:nx*ny)=RESHAPE( lons, (/ nx*ny /) )
        h(v(ivar)%hcoord_type)%lat(1:nx*ny)=RESHAPE( lats, (/ nx*ny /) )
        DEALLOCATE(lons)
        DEALLOCATE(lats)
        hcoord_read(v(ivar)%hcoord_type)=.TRUE.
      END IF
    END DO

    RETURN

  END SUBROUTINE read_hcoords

  SUBROUTINE read_lonlat(gds,lons,lats,nx,ny)

    INTEGER(kind=1),INTENT(in) :: gds(:)
    INTEGER,INTENT(in)         :: nx,ny
    REAL(r8),INTENT(out)       :: lons(1:nx,1:ny),lats(1:nx,1:ny)

    REAL(r8)                   :: lat0,lon0,rlon1,rlat1,rlon2,rlat2,lon,lat
    REAL(r8)                   :: rlats(1:nx,1:ny),rlons(1:nx,1:ny)
    INTEGER                    :: ix,iy

    rlat1=FLOAT(concat_bytes1_sign(gds(11:13),3,.FALSE.))/1000.
    rlon1=FLOAT(concat_bytes1_sign(gds(14:16),3,.FALSE.))/1000.

    rlat2=FLOAT(concat_bytes1_sign(gds(18:20),3,.FALSE.))/1000.
    rlon2=FLOAT(concat_bytes1_sign(gds(21:23),3,.FALSE.))/1000.

    lat0=FLOAT(concat_bytes1_sign(gds(33:35),3,.FALSE.))/1000.
    lon0=FLOAT(concat_bytes1_sign(gds(36:38),3,.FALSE.))/1000.

    DO iy=1,ny
      rlats(:,iy)=rlat1+(iy-1)*(rlat2-rlat1)/FLOAT((ny-1))
    END DO
    DO ix=1,nx
      rlons(ix,:)=rlon1+(ix-1)*(rlon2-rlon1)/FLOAT((nx-1))
    END DO

    CALL rlatlon2latlon(rlats,rlons,lat0, lon0, lats, lons,nx,ny)
    RETURN

  END SUBROUTINE read_lonlat

  SUBROUTINE rlatlon2latlon(rlat,rlon,lat0, lon0, lat, lon, nx, ny)
    IMPLICIT NONE
    INTEGER,INTENT(in)    :: nx,ny
    REAL(r8),INTENT(in)   :: rlat(1:nx,1:ny),rlon(1:nx,1:ny),lon0,lat0
    REAL(r8),INTENT(out)  :: lon(1:nx,1:ny),lat(1:nx,1:ny)

    REAL(r8)              :: a1(1:nx,1:ny),a2(1:nx,1:ny),a3(1:nx,1:ny)
    REAL(r8)              :: lat00,lon00,rlon1(1:nx,1:ny)
    REAL(r8)              :: sinlat0,coslat0,sinlon0,coslon0
    REAL(r8)              :: sinlon(1:nx,1:ny),sinlat(1:nx,1:ny),coslat(1:nx,1:ny),coslon(1:nx,1:ny)

    WHERE (rlon>180.)
      rlon1=rlon-360.
    ELSEWHERE
      rlon1=rlon   
    END WHERE

    lat00 = lat0
    lon00 = lon0
    IF (lon00>180.) lon00 = lon00-360.

    sinlat0 = SIN(lat00*deg2rad)
    coslat0 = COS(lat00*deg2rad)
    sinlon0 = SIN(lon00*deg2rad)
    coslon0 = COS(lon00*deg2rad)
    sinlat = SIN(rlat*deg2rad)
    coslat = COS(rlat*deg2rad)
    sinlon = SIN(rlon1*deg2rad)
    coslon = COS(rlon1*deg2rad)
    a1 = sinlon0 * (-sinlat0*coslon*coslat + &
                     coslat0       *sinlat)- &
         coslon0 *           sinlon*coslat

    a2 = coslon0 * (-sinlat0*coslon*coslat + &
                     coslat0       *sinlat)+ &
         sinlon0 *           sinlon*coslat
    WHERE (abs(a2)<1.e-20)
      a2 = 1.e-20
    END WHERE

    lon = ATAN2(a1,a2)/deg2rad

    a3 = coslat0*coslat*coslon + sinlat0*sinlat
    lat = ASIN(a3)/deg2rad

    RETURN
  END SUBROUTINE rlatlon2latlon

  FUNCTION get_data_from_binary(filename,header,nx,ny) RESULT (data)

    REAL(r8)                      :: data(1:nx,1:ny)
    REAL                          :: data2(1:nx,1:ny)
    CHARACTER(len=256),INTENT(in) :: filename
    TYPE(grib_header_type),INTENT(in)  :: header

    INTEGER,INTENT(in)            :: nx,ny
    INTEGER(kind=1)               :: bin1(1:2),bin4(1:4)

    INTEGER                       :: idsf,ibsf,dval,irec,ix,iy,istat,ibyte,is,ie
    REAL(r8)                      :: dsf,bsf,ref_value
    INTEGER(kind=1),ALLOCATABLE   :: bytearr(:)

    ALLOCATE(bytearr(1:header%data_length+header%data_offset+8))

    CALL byte_to_word_signed(header%pds(27:28),idsf,2)
    dsf=FLOAT(idsf)

    OPEN(10,FILE=TRIM(filename),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)
    is=header%data_record
    ie=header%data_record+CEILING((header%data_length+header%data_offset)/4.)
    ibyte=1
    DO irec=is,ie
      READ(10,rec=irec,iostat=istat) bin4
      bytearr(ibyte:ibyte+3)=bin4
      ibyte=ibyte+4
    END DO
    CLOSE(10)

    bytearr(1:header%data_length)=bytearr(1+header%data_offset:header%data_length+header%data_offset)

    CALL byte_to_word_signed(bytearr(5:6),ibsf,2)
    bsf=FLOAT(ibsf)

    ref_value=to_float1(bytearr(7:10))

    ibyte=12
    DO iy=1,ny
      DO ix=1,nx
        dval=byte_to_word_data(bytearr(ibyte:ibyte+1))
        data(ix,iy)=(ref_value+FLOAT(dval)*(2.**bsf))/(10.**dsf)
        ibyte=ibyte+2
      END DO
    END DO

    RETURN
  END FUNCTION get_data_from_binary

  SUBROUTINE set_vertical_coords(filename,header,nsv,ivctype,sv,pp_index)

    TYPE(cosmo_non_state_data),INTENT(inout) :: nsv
    CHARACTER(len=256),INTENT(in)            :: filename
    TYPE(grib_header_type),INTENT(in)             :: header
    INTEGER,INTENT(in)                       :: ivctype
    REAL(r8),INTENT(in)                      :: sv(:)
    INTEGER,INTENT(in)                       :: pp_index(:)  

    REAL(r8),PARAMETER                       :: g = 9.80665
    REAL(r8),PARAMETER                       :: r = 287.05
    INTEGER                                  :: nhl,nfl,pos
    INTEGER(kind=1)                          :: w(4)
    REAL(r8),ALLOCATABLE                     :: ak(:),bk(:)
    REAL(r8),ALLOCATABLE                     :: vcoord(:)
    REAL(r8)                                 :: t0sl,p0sl,dt0lp,vcflat
    REAL(r8)                                 :: zgdrt,ztdbe,zbetf
    INTEGER                                  :: k,kflat,nx,ny
    REAL(r8),ALLOCATABLE                     :: hhl(:,:,:),p0hl(:,:,:)
    REAL(r8),ALLOCATABLE                     :: pp1(:,:),pp2(:,:)

    ! This routine calculates the vertical coordinate for the 3D grid

    IF (allocated(nsv%surface_orography).AND.(.NOT. nsv%vertical_coords_set)) THEN

      nx=SIZE(nsv%surface_orography,1)
      ny=SIZE(nsv%surface_orography,2)
      
      ! read constants and parameters from the binary data
      
      nhl=to_positive(header%gds(4))-4
      nfl=to_positive(header%gds(4))-5
      
      ALLOCATE(vcoord(1:nhl))
      ALLOCATE(ak(1:nhl))
      ALLOCATE(bk(1:nhl))
      ALLOCATE(hhl(1:nx,1:ny,1:nhl))
      ALLOCATE(p0hl(1:nx,1:ny,1:nhl))
      
      pos=to_positive(header%gds(5))
      w(1:4)=header%gds(pos:pos+3)
      p0sl=to_float1(w)
      pos=pos+4
      w(1:4)=header%gds(pos:pos+3)
      t0sl=to_float1(w)
      pos=pos+4
      w(1:4)=header%gds(pos:pos+3)
      dt0lp=to_float1(w)
      pos=pos+4
      w(1:4)=header%gds(pos:pos+3)
      vcflat=to_float1(w)
      pos=pos+4
      
      DO k=1,nhl
        w(1:4)=header%gds(pos:pos+3)
        vcoord(k)=to_float1(w)
        pos=pos+4
      END DO

      ! calculate needed parameters       
      zgdrt = g/r/t0sl
      IF (dt0lp /= 0.) THEN
        ztdbe = t0sl/dt0lp
      ELSE
        ztdbe = 0.0
      ENDIF
      zbetf = 2.0*dt0lp*zgdrt/t0sl
      
      IF (ivctype==1) THEN
        
        ! Calculate the inverse coordinate transformation, i.e. the ak's and bk's
        kflat = 0
        DO k = 1, nhl
          IF( vcoord(k) <= vcflat ) THEN
            ak(k) = vcoord(k)*p0sl
            bk(k) = 0.0
            kflat = k
          ELSE
            ak(k) = vcflat*p0sl*(1.0 - vcoord(k))/(1.0 - vcflat)
            bk(k) = (vcoord(k) - vcflat)/(1.0 - vcflat)
          ENDIF
        ENDDO
        
        ! Compute the surface reference pressure from surface topography
        hhl(:,:,nhl) = nsv%surface_orography(:,:)
        IF (dt0lp == 0.0) THEN
          p0hl (:,:,nhl) = p0sl*EXP ( - zgdrt*hhl(:,:,nhl) )
        ELSE
          p0hl (:,:,nhl) = p0sl*EXP ( - ztdbe &
           * (1.0 - SQRT(1.0 - zbetf*hhl(:,:,nhl))) )
        ENDIF
        ! Compute the reference pressure at half levels from surface topography
        ! and vertical coordinate parameters ak and bk as well as the
        ! height of half levels from the hydrostatic equation
        
        DO  k = 1, nhl-1
          p0hl(:,:,k) = ak(k) + bk(k)*p0hl(:,:,nhl)
          hhl (:,:,k) = (r/g)*LOG(p0sl/p0hl(:,:,k)) &
           *( t0sl - 0.5*dt0lp*LOG(p0sl/p0hl(:,:,k)) )
        ENDDO
      END IF
      
      IF (ivctype==2) THEN

        ! Calculate the inverse coordinate transformation, i.e. the ak's and bk's
        kflat = 0
        DO k = 1, nhl
          IF( vcoord(k) >= vcflat ) THEN
            ak(k) = vcoord(k)
            bk(k) = 0.0
            kflat = k
          ELSE
            ak(k) = vcoord(k)
            bk(k) = (vcflat - vcoord(k))/ vcflat
          ENDIF
        ENDDO
        
        ! Compute the height of the model half-levels
        hhl(:,:,nhl) = nsv%surface_orography(:,:)
        DO  k = 1, nhl-1
          hhl(:,:,k) = ak(k) + bk(k)*hhl(:,:,nhl)
        ENDDO
        
        ! Compute the reference pressure at half levels
        DO  k = 1, nhl
          IF (dt0lp == 0.0) THEN
            p0hl (:,:,k) = p0sl * EXP ( - zgdrt*hhl(:,:,k) )
          ELSE
            p0hl (:,:,k) = p0sl * EXP ( - ztdbe*(1.0         &
             - SQRT(1.0 - zbetf*hhl(:,:,k))) )
          ENDIF
        ENDDO

      END IF

      ! set the vertical coordinate information in the non-state variable

      nsv%nfl=nfl
      nsv%nhl=nhl
      nsv%p0sl=p0sl
      nsv%t0sl=t0sl
      nsv%dt0lp=dt0lp
      nsv%vcflat=vcflat

      ALLOCATE(nsv%vct_a(1:nhl))
      ALLOCATE(nsv%vct_b(1:nhl))
      ALLOCATE(nsv%hhl(1:nx,1:ny,1:nhl))
      ALLOCATE(nsv%hfl(1:nx,1:ny,1:nfl))
      ALLOCATE(nsv%p0hl(1:nx,1:ny,1:nhl))
      ALLOCATE(nsv%p0fl(1:nx,1:ny,1:nfl))

      nsv%hhl(1:nx,1:ny,1:nhl)=hhl
      nsv%p0hl(1:nx,1:ny,1:nhl)=p0hl
      
      DO  k = 1, nfl 
        nsv%hfl(:,:,k)=0.5*(nsv%hhl(:,:,k)+nsv%hhl(:,:,k+1))
        nsv%p0fl(:,:,k)=0.5*(nsv%p0hl(:,:,k)+nsv%p0hl(:,:,k+1))
      END DO

      if (pp_index(1)>-1) then
        ALLOCATE(nsv%phl(1:nx,1:ny,1:nhl))
        ALLOCATE(nsv%pfl(1:nx,1:ny,1:nfl))
        ALLOCATE(pp1(1:nx,1:ny))
        ALLOCATE(pp2(1:nx,1:ny))
        DO  k=1,nfl
          pp1=RESHAPE(sv(pp_index(k):pp_index(k)+nx*ny-1),(/ nx,ny /))
          nsv%pfl(1:nx,1:ny,k)=nsv%p0fl(1:nx,1:ny,k)+pp1(1:nx,1:ny)
          if (k==1) then
            nsv%phl(1:nx,1:ny,k)=nsv%p0hl(1:nx,1:ny,k)+pp1
          else
            nsv%phl(1:nx,1:ny,k)=nsv%p0hl(1:nx,1:ny,k)+0.5*( &
             pp2*(nsv%hfl(1:nx,1:ny,k-1)-nsv%hhl(1:nx,1:ny,k))/(nsv%hfl(1:nx,1:ny,k-1)-nsv%hfl(1:nx,1:ny,k))+&
             pp1*(nsv%hhl(1:nx,1:ny,k)-nsv%hfl(1:nx,1:ny,k))/(nsv%hfl(1:nx,1:ny,k-1)-nsv%hfl(1:nx,1:ny,k)))
          end if
          pp2=pp1
        END DO
        nsv%phl(1:nx,1:ny,nsv%nhl)=nsv%p0hl(1:nx,1:ny,nsv%nhl)+pp2
        DEALLOCATE(pp1)
        DEALLOCATE(pp2)
      end if

      nsv%vertical_coords_set=.TRUE.

    END IF
    
  END SUBROUTINE set_vertical_coords

  SUBROUTINE read_time(header,tf)
    
    TYPE(grib_header_type),INTENT(in)   :: header
    TYPE(time_type),INTENT(out)    :: tf

    TYPE(time_type)                :: ta
    INTEGER                        :: doy1(13),doy2(13),jd
    INTEGER                        :: ye,mo,da,ho,mi,dt,f1,f2,dc

    ye=header%pds(13)
    mo=header%pds(14)
    da=header%pds(15)
    ho=header%pds(16)
    mi=header%pds(17)

    doy1=(/ 0,31,59,90,120,151,181,212,243,273,304,334,365 /)
    doy2=(/ 0,31,60,91,121,152,182,213,244,274,305,335,366 /)

    SELECT CASE (header%pds(18))
      CASE (0)
        dt=60
      CASE (1)
        dt=60*60
      CASE (2)
        dt=60*60*24
      CASE DEFAULT
        dt=1
    END SELECT

    f1=header%pds(19)
    f2=header%pds(20)

    dc=(header%pds(25)-1)*100
    ye=ye+dc

    IF (f2==0) THEN
      tf=set_time(f1*dt)
    ELSE
      tf=set_time(f2*dt)
    END IF

    call set_calendar_type(GREGORIAN)
    jd=julian_day(ye,mo,da)

!    if is_leap_year(ye) then
!      ta=set_time(ho*60*60+mi*60,(ye-1)*365+doy2(mo)+mo)
!    else
!      ta=set_time(ho*60*60+mi*60,(ye-1)*365+doy1(mo)+mo)
!    end if

    tf=set_time(ho*60*60+mi*60,jd)
    
    RETURN
  END SUBROUTINE read_time
  
END MODULE cosmo_data_mod
