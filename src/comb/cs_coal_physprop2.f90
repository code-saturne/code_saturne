!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

!===============================================================================
! Function:
! --------
!> \file cs_coal_physprop2.f90
!>
!> \brief Calculation of the physical properties of the dispersed phase
!>        (classes of particules)
!>
!> Cell values
!> -----------
!> - Mass fraction of solid
!>   and eventual clippings
!> - Diameter
!> - Mass density
!>   and eventual clippings
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!______________________________________________________________________________!

subroutine cs_coal_physprop2 &
 ( ncelet , ncel )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use entsor
use cstnum
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use field

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel

! Local variables
character(len=80) :: name

integer          nbrint
parameter       (nbrint=8)
integer          iel    , icla
integer          n1     , n2     , n3     , n4     , n5    , n6
integer          n7     , n8
integer          inttmp(nbrint)

double precision xch    , dch    , xnp    , xck    , dck , d1s3
double precision xashcl , xuash
double precision x2min  , x2max  , dckmin , dckmax
double precision dchmin , dchmax , romin  , romax , coedmi
double precision ro2ini , roh2o

double precision, dimension(:), pointer :: nagcpi, agecpi
double precision, dimension(:), pointer :: cvar_xchcl, cvar_xckcl, cvar_xnpcl
double precision, dimension(:), pointer :: cvar_xwtcl
double precision, dimension(:), pointer :: cpro_x2, cpro_rom2, cpro_diam2

!===============================================================================

!===============================================================================
! 1. Initializations
!===============================================================================

d1s3 = 1.d0/3.d0

! --> Coefficient relating to the coke diameter
coedmi = 1.2d0
!
!===============================================================================
! 2. Calculation for each class
!    - of the solid mass fraction
!    - of the coke diameter
!    - of the coal mass density
!===============================================================================
!
do icla = 1, nclacp
  n1 = 0
  n2 = 0
  n3 = 0
  n4 = 0
  n5 = 0
  n6 = 0
  n7 = 0
  n8 = 0
  x2min  =  grand
  x2max  = -grand
  dchmin =  grand
  dchmax = -grand
  dckmin =  grand
  dckmax = -grand
  romin  =  grand
  romax  = -grand

  if (i_comb_drift.ge.1) then
    write(name,'(a,i2.2)') 'n_p_age_', icla
    call field_get_val_s_by_name(name,nagcpi)

    write(name,'(a,i2.2)') 'age_p_', icla
    call field_get_val_s_by_name(name,agecpi)
  endif

  call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xckcl)
  call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xchcl)
  call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnpcl)
  if ( ippmod(iccoal) .ge. 1 ) then
    call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwtcl)
  endif

  call field_get_val_s(ix2(icla),cpro_x2)
  call field_get_val_s(irom2(icla),cpro_rom2)
  call field_get_val_s(idiam2(icla),cpro_diam2)
  xashcl = xashch(ichcor(icla))

  do iel = 1, ncel
    xck    = cvar_xckcl(iel)
    xch    = cvar_xchcl(iel)
    xnp    = cvar_xnpcl(iel)
    xuash  = xnp*xmp0(icla)*(1.d0-xashcl)
    ! --- Calculation of the solid mass fraction
    cpro_x2(iel) = xch + xck + xnp*xmash(icla)
    !     Taking into account the humidity
    if ( ippmod(iccoal) .ge. 1 ) then
      cpro_x2(iel) = cpro_x2(iel)                     &
                          +cvar_xwtcl(iel)
    endif
    ! ---- Eventual clipping for the solid mass fraction
    if ( cpro_x2(iel) .gt. (1.d0+epsicp) ) then
      n1 = n1 + 1
      x2max = max(cpro_x2(iel),x2max)
      cpro_x2(iel) = 1.d0
    else if ( cpro_x2(iel) .lt. (zero-epsicp) ) then
      n2 = n2 + 1
      x2min = min(cpro_x2(iel),x2min)
      cpro_x2(iel) = zero
    endif


    ! --- Initialization

    cpro_rom2(iel) = rho20(icla)
    cpro_diam2(iel) = diam20(icla)

    if ( xuash.gt.epsicp ) then

      ! --- Calculation of the reactive coal diameter: Dch

      dch = diam20(icla)*(xch/xuash)**d1s3

      ! ---- Eventual clipping for the reactive coal diameter

      if ( dch .gt. (diam20(icla)+epsicp) ) then
        n3 = n3 + 1
        dchmax = max(dch,dchmax)
        dch = diam20(icla)
      else if ( dch .lt. (zero-epsicp) ) then
        n4 = n4 + 1
        dchmin = min(dch,dchmin)
        dch = zero
      endif

      ! --- Calculation of the coke diamter: Dck stores in cpro_diam2(iel)

      dck = ( (xch/rho20(icla)+xck/rhock(ichcor(icla)))/          &
              ((1.d0-xashcl)*pi/6.d0*xnp) )**d1s3

      ! ---- Eventual clipping for the coke diameter

      if ( dck .gt. coedmi*diam20(icla) ) then
        n5 = n5 + 1
        dckmax = max(dck,dckmax)
        dck = diam20(icla)*coedmi
      else if ( dck .lt. (zero-epsicp) ) then
        n6 = n6 + 1
        dckmin = min(dck,dckmin)
        dck = zero
      endif
      cpro_diam2(iel) = dck

      ! --- Density

      ro2ini = rho20(icla)
      !     Taking into account humidity
      if ( ippmod(iccoal) .eq. 1 ) then
      !     at the moment we asume that ROH2O is constant
        roh2o = 998.203
        ro2ini = rho20(icla)+ cvar_xwtcl(iel)                     &
                             *roh2o
      endif

      cpro_rom2(iel) =                                        &
        ( xashcl*diam20(icla)**3*rho20(icla) +                    &
          (1.d0-xashcl)*(dck**3-dch**3)*rhock(ichcor(icla)) +     &
          (1.d0-xashcl)*dch**3*ro2ini ) /                         &
        ( xashcl*diam20(icla)**3 +                                &
          (1.d0-xashcl)*dck**3 )

      ! ---- Clipping for density

      if ( cpro_rom2(iel) .gt. (ro2ini+epsicp) ) then
        n7 = n7 + 1
        romax = max(cpro_rom2(iel),romax)
        cpro_rom2(iel) = rho20(icla)
      endif
      if ( cpro_rom2(iel) .lt. (rhock(ichcor(icla))-epsicp) ) &
                              then
        n8 = n8 + 1
        romin = min(cpro_rom2(iel),romin)
        cpro_rom2(iel) = rhock(ichcor(icla))
      endif
    endif

    ! Particles' age of each particle class
    if(i_comb_drift.ge.1) then
      if (xnp.ge.epsicp) then
        agecpi(iel) = nagcpi(iel)/xnp
      else
        agecpi(iel) = 0.d0
      endif
    endif
  enddo

  if (irangp.ge.0) then

    inttmp(1) = n1
    inttmp(2) = n2
    inttmp(3) = n3
    inttmp(4) = n4
    inttmp(5) = n5
    inttmp(6) = n6
    inttmp(7) = n7
    inttmp(8) = n8
    call parism (nbrint,inttmp)
    !==========
    n1 = inttmp(1)
    n2 = inttmp(2)
    n3 = inttmp(3)
    n4 = inttmp(4)
    n5 = inttmp(5)
    n6 = inttmp(6)
    n7 = inttmp(7)
    n8 = inttmp(8)

    call parmax (x2max )
    !==========
    call parmax (dchmax)
    !==========
    call parmax (dckmax)
    !==========
    call parmax (romax )
    !==========

    call parmin (x2min )
    !==========
    call parmin (dchmin)
    !==========
    call parmin (dckmin)
    !==========
    call parmin (romin )
    !==========

    call synsca(cpro_x2)

  endif

  if ( n1 .gt. 0 ) then
    write(nfecra,1001) icla, n1, x2max
  endif
  if ( n2 .gt. 0 ) then
    write(nfecra,1002) icla, n2, x2min
  endif
  if ( n3 .gt. 0 ) then
    write(nfecra,1003) icla, n3, dchmax
  endif
  if ( n4 .gt. 0 ) then
    write(nfecra,1004) icla, n4, dchmin
  endif
  if ( n5 .gt. 0 ) then
    write(nfecra,1005) icla, n5, dckmax
  endif
  if ( n6 .gt. 0 ) then
    write(nfecra,1006) icla, n6, dckmin
  endif
  if ( n7 .gt. 0 ) then
    write(nfecra,1007) icla, n7, romax
  endif
  if ( n8 .gt. 0 ) then
    write(nfecra,1008) icla, n8, romin
  endif

enddo

!--------
! Formats
!--------

!===============================================================================
 1001 format(/,1X,' clipping in max for solid mass frac. for class',    &
        I3,/,10X,' Number of points : ',I8,                       &
           /,10X,' Max value        : ',G15.7)
 1002 format(/,1X,' clipping in min for solid mass frac. for class',    &
        I3,/,10X,' Number of points: ',I8,                       &
           /,10X,' Max value        : ',G15.7)
 1003 format(/,1X,' clipping in max of coal diameter for class    ',    &
        I3,/,10X,' Number of points: ',I8,                       &
           /,10X,' Max value        : ',G15.7)
 1004 format(/,1X,' clipping in max of coal diameter for class    ',    &
        I3,/,10X,' Number of points: ',I8,                       &
           /,10X,' Min value        : ',G15.7)
 1005 format(/,1X,' clipping in max for coke diameter for class   ',    &
        I3,/,10X,' Number of points: ',I8,                       &
           /,10X,' Max value        : ',G15.7)
 1006 format(/,1X,' clipping in min for coke diameter for class   ',    &
        I3,/,10X,' Number of points: ',I8,                       &
           /,10X,' Min value        : ',G15.7)
 1007 format(/,1X,' clipping in max of mass density for class       ',  &
        I3,/,10X,' Number of points: ',I8,                       &
           /,10X,' Max value        : ',G15.7)
 1008 format(/,1X,' clipping in min for mass density for class      ',  &
        I3,/,10X,' Number of points: ',I8,                       &
           /,10X,' Min value        : ',G15.7)
!===============================================================================

!----
! End
!----

return
end subroutine
