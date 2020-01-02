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
!> \file cs_coal_thfieldconv2.f90
!> \brief Calculation of the particles temperature
!>        Function with the solid enthalpy and concentrations
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!   mode          name          role
!______________________________________________________________________________!
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!______________________________________________________________________________!

subroutine cs_coal_thfieldconv2 &
 ( ncelet , ncel )

!==============================================================================
! Module files
!==============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use field

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel

! Local variables

integer          i      , icla   , icha   , iel
integer          ihflt2
double precision h2     , x2     , xch    , xck
double precision xash   , xnp    , xtes   , xwat

integer          iok1
double precision , dimension ( : )     , allocatable :: eh0,eh1

double precision, dimension(:), pointer :: cvar_xchcl, cvar_xckcl, cvar_xnpcl
double precision, dimension(:), pointer :: cvar_xwtcl
double precision, dimension(:), pointer :: cvar_h2cl
double precision, dimension(:), pointer :: cpro_temp2, cpro_temp1

!===============================================================================
! Conversion method
!
ihflt2 = 1
!
!===============================================================================
! 1. Preliminary calculations
!===============================================================================

!===============================================================================
! Deallocation dynamic arrays
!----
allocate(eh0(1:ncel),eh1(1:ncel),STAT=iok1)
!----
if ( iok1 > 0 ) THEN
  write(nfecra,*) ' Memory allocation error inside: '
  write(nfecra,*) '    cs_coal_thfieldconv2         '
  call csexit(1)
endif
!===============================================================================

! --- Tables initialization
eh0( : ) = zero
eh1( : ) = zero

! --- Initialization from T2 to T1

call field_get_val_s(itemp1,cpro_temp1)
do icla = 1, nclacp
  call field_get_val_s(itemp2(icla),cpro_temp2)
  do iel = 1, ncel
    cpro_temp2(iel) = cpro_temp1(iel)
  enddo
enddo

!===============================================================================
! 2. Particles temperature calculation
!===============================================================================

if ( ihflt2.eq.0 ) then

! --> H2 linear function of T2

  do icla = 1, nclacp
    call field_get_val_s(ivarfl(isca(ih2(icla))), cvar_h2cl)
    call field_get_val_s(itemp2(icla),cpro_temp2)
    icha = ichcor(icla)
    do iel = 1, ncel
      cpro_temp2(iel) =                                       &
            (cvar_h2cl(iel)-h02ch(icha))                         &!FIXME divide by x2
            / cp2ch(icha) + trefth
    enddo
  enddo

else

! --> H2 tabule

  do icla = 1, nclacp

    call field_get_val_s(ivarfl(isca(ixch(icla))), cvar_xchcl)
    call field_get_val_s(ivarfl(isca(ixck(icla))), cvar_xckcl)
    call field_get_val_s(ivarfl(isca(inp(icla))), cvar_xnpcl)
    if ( ippmod(iccoal) .eq. 1 ) then
      call field_get_val_s(ivarfl(isca(ixwt(icla))), cvar_xwtcl)
    endif
    call field_get_val_s(ivarfl(isca(ih2(icla))), cvar_h2cl)
    call field_get_val_s(itemp2(icla),cpro_temp2)

    i = npoc-1
    do iel = 1, ncel
      xch  = cvar_xchcl(iel)
      xck  = cvar_xckcl(iel)
      xnp  = cvar_xnpcl(iel)
      xash = xmash(icla)*xnp
      if ( ippmod(iccoal) .eq. 1 ) then
        xwat = cvar_xwtcl(iel)
      else
        xwat = 0.d0
      endif

      x2   = xch + xck + xash + xwat

      xtes = xmp0(icla)*xnp

      if ( xtes.gt.epsicp .and. x2.gt.epsicp*100.d0 ) then

        h2   = cvar_h2cl(iel)/x2

        eh1(iel) = xch /x2 * ehsoli(ich(ichcor(icla) ),i+1)    &
                  + xck /x2 * ehsoli(ick(ichcor(icla) ),i+1)    &
                  + xash/x2 * ehsoli(iash(ichcor(icla)),i+1)    &
                  + xwat/x2 * ehsoli(iwat(ichcor(icla)),i+1)
        if ( h2.ge.eh1(iel) ) then
          cpro_temp2(iel) = thc(i+1)

        endif

      endif

    enddo

    i = 1
    do iel = 1, ncel
      xch  = cvar_xchcl(iel)
      xck  = cvar_xckcl(iel)
      xnp  = cvar_xnpcl(iel)
      xash = xmash(icla)*xnp
      if ( ippmod(iccoal) .eq. 1 ) then
        xwat = cvar_xwtcl(iel)
      else
        xwat = 0.d0
      endif

      x2   = xch + xck + xash + xwat
      xtes = xmp0(icla)*xnp

      if ( xtes.gt.epsicp .and. x2.gt.epsicp*100.d0 ) then

        h2   = cvar_h2cl(iel)/x2


        eh0(iel) = xch /x2 * ehsoli(ich(ichcor(icla) ),i)      &
                  + xck /x2 * ehsoli(ick(ichcor(icla) ),i)      &
                  + xash/x2 * ehsoli(iash(ichcor(icla)),i)      &
                  + xwat/x2 * ehsoli(iwat(ichcor(icla)),i)
        if ( h2.le.eh0(iel) ) then
          cpro_temp2(iel) = thc(i)
        endif

      endif

    enddo

    do i = 1, npoc-1
      do iel = 1, ncel
        xch  = cvar_xchcl(iel)
        xck  = cvar_xckcl(iel)
        xnp  = cvar_xnpcl(iel)
        xash = xmash(icla)*xnp
        if ( ippmod(iccoal) .eq. 1 ) then
          xwat = cvar_xwtcl(iel)
        else
          xwat = 0.d0
        endif

        x2   = xch + xck + xash + xwat
        xtes = xmp0(icla)*xnp

        if ( xtes.gt.epsicp .and. x2.gt.epsicp*100.d0 ) then

          h2   = cvar_h2cl(iel)/x2

          eh0(iel) = xch /x2 * ehsoli(ich(ichcor(icla) ),i  )  &
                    + xck /x2 * ehsoli(ick(ichcor(icla) ),i  )  &
                    + xash/x2 * ehsoli(iash(ichcor(icla)),i  )  &
                    + xwat/x2 * ehsoli(iwat(ichcor(icla)),i  )

          eh1(iel) = xch /x2 * ehsoli(ich(ichcor(icla) ),i+1)  &
                    + xck /x2 * ehsoli(ick(ichcor(icla) ),i+1)  &
                    + xash/x2 * ehsoli(iash(ichcor(icla)),i+1)  &
                    + xwat/x2 * ehsoli(iwat(ichcor(icla)),i+1)

          if ( h2.ge.eh0(iel) .and. h2.le.eh1(iel) ) then
            cpro_temp2(iel) = thc(i) + (h2-eh0(iel)) *     &
                  (thc(i+1)-thc(i))/(eh1(iel)-eh0(iel))
          endif

        endif

      enddo
    enddo

  enddo

endif

!===============================================================================
! Deallocation dynamic arrays
!----
deallocate(eh0,eh1,STAT=iok1)
!----
if ( iok1 > 0 ) then
  write(nfecra,*) ' Memory deallocation error inside: '
  write(nfecra,*) '    cs_coal_thfieldconv2           '
  call csexit(1)
endif
!===============================================================================

!----
! End
!----

return
end subroutine
