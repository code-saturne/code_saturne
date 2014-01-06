!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine ecrlis &
!================

 ( nvar   , ncelet , ncel   ,                                     &
   rtp    , rtpa   , dt     , volume )

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE D'ECRITURE DES INFOS LISTING

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! e  ! <-- ! nombre de variables                            !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! irtp             ! e  ! <-- ! indice de rtp dans ra                          !
! rtp              ! tr ! <-- ! tableaux des variables au pdt courant          !
! (ncelet,nvar)    !    !     !                                                !
! rtpa             ! tr ! <-- ! tableaux des variables au pdt prec             !
! (ncelet,nvar)    !    !     !                                                !
! dt   (ncelet)    ! tr ! <-- ! valeur du pas de temps                         !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet)         !    !     !                                                !
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
use entsor
use optcal
use cstnum
use cstphy
use albase
use parall
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

integer          nvar, ncelet, ncel
double precision rtpa(ncelet,nvar), rtp(ncelet,nvar)
double precision dt(ncelet), volume(ncelet)

! Local variables

integer          ic, icel, ivar, ipp
character*200    chain, chainc

!==================================================================
! 1. DERIVE POUR LES VARIABLES TRANSPORTEES (sauf pression)
!==================================================================

do ivar = 1, nvar
  if (ivar.eq.ipr.and.ippmod(icompf).lt.0) cycle
  ipp = ipprtp(ivar)
  if (ilisvr(ipp).gt.0) then
    dervar(ipp) = 0
    do icel = 1, ncel
      dervar(ipp) = dervar(ipp)                                 &
              + (rtp(icel,ivar)-rtpa(icel,ivar))**2             &
              *  volume(icel)/dt(icel)
    enddo
    if (irangp.ge.0) call parsom (dervar(ipp))
    dervar(ipp) = dervar(ipp) / voltot
  endif
enddo

if (ippmod(icompf).lt.0) then
  ipp = ipprtp(ipr)
  if (dervar(ipp).lt.epzero) then
    dervar(ipp) = -1.d0
  endif
  dervar(ipp) = rnsmbr(ipp) / dervar(ipp)
endif

!===============================================================================
! 2. ECRITURE DES CRITERES DE CONVERGENCE
!===============================================================================

write(nfecra,1000)
write(nfecra,1010)
write(nfecra,1011)
write(nfecra,1010)

do ivar = 1, nvar
  ipp = ipprtp(ivar)
  if (ilisvr(ipp).gt.0) then
    chainc = 'c'
    chain = ' '
    ic = 4
    chain = nomvar(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+12
    chain = ' '
    write(chain,3000) rnsmbr(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+13
    chain = ' '
    write(chain,4000) nbivar(ipp)
    chainc(ic:ic+7) = chain(1:7)
    ic=ic+9
    chain = ' '
    write(chain,3000) resvar(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+14
    chain = ' '
    write(chain,3000) dervar(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+12
    write(nfecra,'(a)') chainc(1:ic)
  endif
enddo

write(nfecra,1010)
write(nfecra,*) ' '
write(nfecra,*) ' '

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format (/,3X,'** INFORMATIONS SUR LA CONVERGENCE',/,        &
          3X,'   -------------------------------')
 1011 format ('   Variable    Norm 2nd mb.',                      &
        '  Nbiter  Residu norme        derive')
 1010 format ('---------------------------',                      &
        '------------------------------------')

 3000 format (e12.5)
 4000 format (i7)

#else

 1000 format (/,3X,'** INFORMATION ON CONVERGENCE',/,             &
          3X,'   --------------------------')
 1011 format ('   Variable    Rhs norm    ',                      &
        '  N_iter  Norm. residual      derive')
 1010 format ('---------------------------',                      &
        '------------------------------------')

 3000 format (e12.5)
 4000 format (i7)

#endif
return
end subroutine
