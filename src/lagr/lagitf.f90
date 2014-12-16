!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine lagitf &
!================

 ( iprev, propce )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!     INTEGRATION DES EDS POUR LA TEMPERATURE FLUIDE
!      VU PAR LES PARTICULES.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! iprev            ! e  ! <-- ! time step indicator for fields                 !
!                  !    !     !   0: use fields at current time step           !
!                  !    !     !   1: use fields at previous time step          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstphy
use cstnum
use optcal
use entsor
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          iprev

double precision propce(ncelet,*)

! Local variables

logical          ltsvar
integer          npt   , iel   , mode

double precision ct    , aux1  , aux2   , ter1   , ter2
double precision energ , dissip

double precision, allocatable, dimension(:) :: tempf

double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_omg
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), allocatable :: auxl1

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! Allocate a temporary array
allocate(auxl1(nbpart))
allocate(tempf(ncelet))

! Initialize variables to avoid compiler warnings

dissip = 0.d0
energ = 0.d0

ct = 1.d0

mode = 1

if (associated(ptsvar)) then
  ltsvar = .true.
else
  ltsvar = .false.
endif

if (iprev.eq.0) then

  if (iscalt.gt.0) call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

  if (itytur.eq.2.or.iturb.eq.50) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iep), cvar_ep)
  elseif (itytur.eq.3) then
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
    call field_get_val_s(ivarfl(iep), cvar_ep)
  elseif (iturb.eq.60) then
    call field_get_val_s(ivarfl(ik), cvar_k)
    call field_get_val_s(ivarfl(iomg), cvar_omg)
  endif

else if (iprev.eq.1) then

  if (iscalt.gt.0) call field_get_val_prev_s(ivarfl(isca(iscalt)), cvar_scalt)

  if (itytur.eq.2.or.iturb.eq.50) then
    call field_get_val_prev_s(ivarfl(ik), cvar_k)
    call field_get_val_prev_s(ivarfl(iep), cvar_ep)
  elseif (itytur.eq.3) then
    call field_get_val_prev_s(ivarfl(ir11), cvar_r11)
    call field_get_val_prev_s(ivarfl(ir22), cvar_r22)
    call field_get_val_prev_s(ivarfl(ir33), cvar_r33)
    call field_get_val_prev_s(ivarfl(iep), cvar_ep)
  elseif (iturb.eq.60) then
    call field_get_val_prev_s(ivarfl(ik), cvar_k)
    call field_get_val_prev_s(ivarfl(iomg), cvar_omg)
  endif

endif

!===============================================================================
! 2. Temperature moyenne Fluide en degres Celsius
!===============================================================================

if ( ippmod(iccoal).ge.0 .or.                                     &
     ippmod(icpl3c).ge.0 .or.                                     &
     ippmod(icfuel).ge.0      ) then

  do iel = 1,ncel
    tempf(iel) = propce(iel,ipproc(itemp1)) - tkelvi
  enddo

else if ( ippmod(icod3p).ge.0 .or.                                &
          ippmod(icoebu).ge.0 .or.                                &
          ippmod(ielarc).ge.0 .or.                                &
          ippmod(ieljou).ge.0      ) then

  do iel = 1,ncel
    tempf(iel) = propce(iel,ipproc(itemp)) - tkelvi
  enddo

else if (itherm.eq.1 .and. itpscl.eq.2) then
  do iel = 1,ncel
    tempf(iel) = cvar_scalt(iel)
  enddo

else if (itherm.eq.1 .and. itpscl.eq.1) then
  do iel = 1,ncel
    tempf(iel) = cvar_scalt(iel) - tkelvi
  enddo

else if (itherm.eq.2) then
  do iel = 1,ncel
    call usthht (mode, cvar_scalt(iel), tempf(iel))
    !==========
  enddo
endif

!===============================================================================
! 3. INTEGRATION DE L'EDS SUR LES PARTICULES
!===============================================================================

do npt = 1,nbpart

  if ( ipepa(jisor,npt).gt.0 ) then

    iel = ipepa(jisor,npt)

    if (itytur.eq.2 .or. itytur.eq.3 .or.           &
        iturb.eq.50 .or. iturb.eq.60) then

      if ( itytur.eq.2 .or. iturb.eq.50 ) then

        energ  = cvar_k(iel)
        dissip = cvar_ep(iel)

      else if ( itytur.eq.3 ) then

        energ  = 0.5d0 * ( cvar_r11(iel)                   &
                         + cvar_r22(iel)                   &
                         + cvar_r33(iel) )
        dissip = cvar_ep(iel)

      else if (iturb.eq.60) then

        energ  = cvar_k(iel)
        dissip = cmu*cvar_k(iel)*cvar_omg(iel)

      endif

      auxl1(npt) = energ / ( ct*dissip )
      auxl1(npt) = max( auxl1(npt),epzero )

    else

      auxl1(npt) = epzero

    endif

  endif
enddo

if (nor.eq.1) then

  do npt = 1,nbpart

    if (ipepa(jisor,npt).gt.0) then

      iel = ipepa(jisor,npt)

      aux1 = -dtp/auxl1(npt)
      aux2 = exp(aux1)

      ter1 = eptpa(jtf,npt) * aux2
      ter2 = tempf(iel) * ( 1.d0-aux2 )

      eptp(jtf,npt) = ter1 + ter2

      ! Pour le cas NORDRE= 2, on calcule en plus TSVAR pour NOR= 2

      if (ltsvar) then
        ptsvar(jtf,npt) =    0.5d0 * ter1                                 &
                           + tempf(iel) * ( -aux2 +(aux2-1.d0) / aux1 )
      endif

    endif
  enddo

else if (nor.eq.2) then

  do npt = 1,nbpart

    if (ipepa(jisor,npt).gt.0 .and. ipepa(jord1,npt).eq.0) then

      iel = ipepa(jisor,npt)

      aux1 = -dtp/auxl1(npt)
      aux2 = exp(aux1)

      ter1 = 0.5d0 * eptpa(jtf,npt) * aux2
      ter2 = tempf(iel) * (1.d0 - (aux2-1.d0) / aux1)

      eptp(jtf,npt) = ptsvar(jtf,npt) + ter1 + ter2
    endif
  enddo

endif

! Free memory
deallocate(tempf)
deallocate(auxl1)

!===============================================================================

!=======
! FORMAT
!=======

!----
! FIN
!----

end subroutine
