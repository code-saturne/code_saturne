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

 ( nvar   , ndim   , ncelet , ncel   ,                            &
   irtp   ,                                                       &
   rtp    , rtpa   , dt     , volume , xyzcen ,                   &
   ra     )

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
! ndim             ! i  ! <-- ! spatial dimension                              !
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
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! ra(*)            ! ra ! --- ! main real work array                           !
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

integer          nvar, ndim, ncelet, ncel
integer          irtp
double precision rtpa(ncelet,nvar), rtp(ncelet,nvar)
double precision dt(ncelet), volume(ncelet)
double precision xyzcen(ndim,ncelet)
double precision, dimension(*), target :: ra

! Local variables

integer          ii, jj, ic, icel, ipp, ira, ivrtp, iok
integer          ipuvw
integer          icmin, icmax
integer          nbrval
integer          idivdt, ixmsdt, iel
double precision petit,xyzmin(3),xyzmax(3)
character*200    chain, chainc

double precision, dimension(:), allocatable, target :: momtmp
double precision, dimension(:), pointer :: varptr => null()

!===============================================================================
! 0. INITIALISATIONS LOCALES
!===============================================================================

petit  =-grand

!==================================================================
! 1. DERIVE POUR LES VARIABLES TRANSPORTEES (sauf pression)
!==================================================================

do ipp = 2, nvppmx
  iok = 1
  if(ipp.eq.ipprtp(ipr)) then
    iok = 0
  endif
  if(ilisvr(ipp).eq.1.and.itrsvr(ipp).ge.1) then
    if(iok.eq.1) then
      ira = abs(ipp2ra(ipp))
      ivrtp = (ira-irtp)/ncelet+1
      dervar(ipp) = 0.d0
      do icel = 1, ncel
        dervar(ipp) = dervar(ipp)                                 &
                + (rtp(icel,ivrtp)-rtpa(icel,ivrtp))**2           &
                *  volume(icel)/dt(icel)
      enddo
      if(irangp.ge.0) call parsom (dervar(ipp))
      !==========
      dervar(ipp) = dervar(ipp) / voltot
    endif
  endif
enddo

ipp = ipprtp(ipr)
if(dervar(ipp).lt.epzero) then
  dervar(ipp) = -1.d0
endif
dervar(ipp) = rnsmbr(ipp) / dervar(ipp)

!===============================================================================
! 2. ECRITURE DES CRITERES DE CONVERGENCE
!===============================================================================

write(nfecra,1000)
write(nfecra,1010)
write(nfecra,1011)
write(nfecra,1010)

do ipp = 2, nvppmx
  if(ilisvr(ipp).eq.1.and.itrsvr(ipp).ge.1) then

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

    write(nfecra,'(A)') chainc(1:ic)

  endif
enddo

write(nfecra,1010)
write(nfecra,*) ' '
write(nfecra,*) ' '

!===============================================================================
! 3.  SUIVI DES VARIABLES , MIN - MAX - CLIPPING
!===============================================================================

write(nfecra,1110)
write(nfecra,*) ' '
write(nfecra,*) ' '

write(nfecra,1150)
write(nfecra,1160)
write(nfecra,1161)
write(nfecra,1160)

do ipp = 2, nvppmx
  ipuvw = 0
  if (ipp.eq.ipprtp(ipr) .or.       &
      ipp.eq.ipprtp(iu ) .or.       &
      ipp.eq.ipprtp(iv ) .or.       &
      ipp.eq.ipprtp(iw )) then
    ipuvw = 1
  endif
  !   En v2f on ne clippe jamais f_barrre, on ne l'affiche donc pas
  if (iturb.eq.50) then
    if (ipp.eq.ipprtp(ifb)) ipuvw = 1
  endif
  !   En ALE on ne clippe pas la vitesse de maillage
  if (iale.eq.1) then
    if (ipp.eq.ipprtp(iuma) .or.                                &
        ipp.eq.ipprtp(ivma) .or.                                &
        ipp.eq.ipprtp(iwma)) ipuvw = 1
  endif
  !   Compressible
  if(ippmod(icompf).ge.0) then
    if(ipp.eq.ipprtp(isca(itempk))) then
      ipuvw = 1
    endif
  endif
  if(ilisvr(ipp).eq.1.and.itrsvr(ipp).gt.0.and.ipuvw.eq.0) then
    chainc = 'a'
    chain = ' '
    ic = 4
    chain = nomvar(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+12
    chain = ' '
    write(chain,3000) varmna(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+14
    chain = ' '
    write(chain,3000) varmxa(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+16

    chain = ' '
    write(chain,4000) iclpmn(ipp)
    chainc(ic:ic+7) = chain(1:7)
    ic=ic+11
    chain = ' '
    write(chain,4000) iclpmx(ipp)
    chainc(ic:ic+7) = chain(1:7)
    ic=ic+7

    write(nfecra,'(A)') chainc(1:ic)
    iclpmn(ipp) = 0
    iclpmx(ipp) = 0
  endif
enddo

write(nfecra,1160)
write(nfecra,*) ' '
write(nfecra,*) ' '

!===============================================================================
! 4.  PAS DE TEMPS LOCAL
!===============================================================================

if (idtvar.ge.0 .and. iptlro.eq.1) then

  write(nfecra,1200)
  write(nfecra,1210)

  chainc = ' '
  chainc = 'Clips a Dtrho'
  ic = 57
  chain = ' '
  write(chain,4000) nclptr
  chainc(ic:ic+7) = chain(1:7)
  ic = ic+7
  write(nfecra,'(a)') chainc(1:ic)
  write(nfecra,1210)

  write(nfecra,*) ' '
  write(nfecra,*) ' '

endif

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


 1110 format ('-----------------------------------------',        &
        '----------------------')

 1150 format (/,3X,'** INFORMATIONS SUR LES CLIPPINGS',/,         &
          3X,'   ------------------------------')
 1161 format ('   Variable    Min ss clips  Max ss clips',        &
        '   Clip min   Clip max')
 1160 format ('-----------------------------------------',        &
        '----------------------')


 1200 format (/,3X,'** INFORMATIONS SUR LE PAS DE TEMPS',/,       &
          3X,'   --------------------------------')
 1210 format ('---------------------------',                      &
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


 1110 format ('-----------------------------------------',        &
        '----------------------')

 1150 format (/,3X,'** INFORMATION ON CLIPPINGS',/,               &
          3X,'   ------------------------')
 1161 format ('   Variable    Min wo clips  Max wo clips',        &
        '   Min clip   Max clip')
 1160 format ('-----------------------------------------',        &
        '----------------------')


 1200 format (/,3X,'** INFORMATION ON THE TIME STEP',/,           &
         3X,'   -----------------------------')
 1210 format ('---------------------------',                      &
        '------------------------------------')

 3000 format (e12.5)
 4000 format (i7)

#endif
return
end subroutine
