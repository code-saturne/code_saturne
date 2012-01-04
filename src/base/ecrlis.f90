!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

!================================================
! 1. CALCUL DES MIN ET MAX DES VARIABLES
!================================================

do ipp = 2, nvppmx
  if(ilisvr(ipp).eq.1) then
    varmin(ipp) = grand
    varmax(ipp) = petit
    ira = abs(ipp2ra(ipp))

    ! For moments, we must divide by the cumulative time
    idivdt = ippmom(ipp)
    if (idivdt.eq.0) then
      varptr => ra(ira:ira+ncel)
    else
      allocate(momtmp(ncel))
      varptr => momtmp
      if (idivdt.gt.0) then
        do iel = 1, ncel
          momtmp(iel) = ra(ira+iel-1)/max(ra(idivdt+iel-1),epzero)
        enddo
      elseif (idivdt.lt.0) then
        do iel = 1, ncel
          momtmp(iel) = ra(ira+iel-1)/max(dtcmom(-idivdt),epzero)
        enddo
      endif
    endif

    do icel = 1, ncel
      if(varptr(icel).lt.varmin(ipp)) varmin(ipp) = varptr(icel)
      if(varptr(icel).gt.varmax(ipp)) varmax(ipp) = varptr(icel)
    enddo

    if (idivdt.ne.0) then
      deallocate(momtmp)
    endif

    if (irangp.ge.0) then
      call parmin (varmin(ipp))
      !==========
      call parmax (varmax(ipp))
      !==========
    endif

  endif
enddo

!==================================================================
! 2. DERIVE POUR LES VARIABLES TRANSPORTEES (sauf pression)
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


!==================================================================
! 3. MIN MAX DT + LOCALISATION
!==================================================================

if (idtvar.ge.0) then

!     -> min max dt
  ptploc(7,1) = grand
  ptploc(8,1) = petit
  icmin = 1
  icmax = 1
  do icel = 1, ncel
    if(dt(icel).lt.ptploc(7,1)) then
      ptploc(7,1) = dt(icel)
      icmin = icel
    endif
    if(dt(icel).gt.ptploc(8,1)) then
      ptploc(8,1) = dt(icel)
      icmax = icel
    endif
  enddo

  xyzmin(1) = xyzcen(1,icmin)
  xyzmin(2) = xyzcen(2,icmin)
  xyzmin(3) = xyzcen(3,icmin)
  xyzmax(1) = xyzcen(1,icmax)
  xyzmax(2) = xyzcen(2,icmax)
  xyzmax(3) = xyzcen(3,icmax)

  if (irangp.ge.0) then
    nbrval = 3
    call parmnl (nbrval, ptploc(7,1), xyzmin)
    !==========
    call parmxl (nbrval, ptploc(8,1), xyzmax)
    !==========
  endif

  ptploc(7,2) = xyzmin(1)
  ptploc(7,3) = xyzmin(2)
  ptploc(7,4) = xyzmin(3)
  ptploc(8,2) = xyzmax(1)
  ptploc(8,3) = xyzmax(2)
  ptploc(8,4) = xyzmax(3)

endif


!===============================================================================
! 4. ECRITURE DES CRITERES DE CONVERGENCE
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
! 5.  SUIVI DES VARIABLES , MIN - MAX - CLIPPING
!===============================================================================


write(nfecra,1100)
write(nfecra,1110)
write(nfecra,1111)
write(nfecra,1110)

do ipp = 2, nvppmx
  if(ilisvr(ipp).eq.1) then

    chainc = 'v'
    chain = ' '
    ic = 4
    chain = nomvar(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+12
    chain = ' '
    write(chain,3000) varmin(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+14
    chain = ' '
    write(chain,3000) varmax(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+16
    ipuvw = 0
    if(ipp.eq.ipprtp(ipr) .or.                           &
         ipp.eq.ipprtp(iu ) .or.                           &
         ipp.eq.ipprtp(iv ) .or.                           &
         ipp.eq.ipprtp(iw ) ) then
      ipuvw = 1
    endif
    !   En v2f phi-fbar on ne clippe jamais f_barre, on ne l'affiche donc pas
    if (iturb.eq.50) then
      if (ipp.eq.ipprtp(ifb)) ipuvw = 1
    endif
    !   En ALE on ne clippe pas la vitesse de maillage
    if (iale.eq.1) then
      if (ipp.eq.ipprtp(iuma) .or.                                &
          ipp.eq.ipprtp(ivma) .or.                                &
          ipp.eq.ipprtp(iwma)) ipuvw = 1
    endif
    if(ipuvw.eq.1) then
      chain = '     --         --'
      chainc(ic:ic+18) = chain(1:18)
      ic = ic+18
    else
      chain = ' '
      write(chain,4000) iclpmn(ipp)
      chainc(ic:ic+7) = chain(1:7)
      ic=ic+11
      chain = ' '
      write(chain,4000) iclpmx(ipp)
      chainc(ic:ic+7) = chain(1:7)
      ic=ic+7
    endif

    write(nfecra,'(A)') chainc(1:ic)
!MO          ICLPMN(IPP) = 0
!MO          ICLPMX(IPP) = 0
  endif
enddo

write(nfecra,1110)
write(nfecra,*) ' '
write(nfecra,*) ' '




write(nfecra,1150)
write(nfecra,1160)
write(nfecra,1161)
write(nfecra,1160)

do ipp = 2, nvppmx
  ipuvw = 0
  if(ipp.eq.ipprtp(ipr) .or.                             &
       ipp.eq.ipprtp(iu ) .or.                             &
       ipp.eq.ipprtp(iv ) .or.                             &
       ipp.eq.ipprtp(iw ) ) then
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
! 6.  PAS DE TEMPS LOCAL
!===============================================================================

if (idtvar.ge.0) then

  write(nfecra,1200)
  write(nfecra,1210)
  write(nfecra,1211)
  write(nfecra,1210)
  do ii = 1, 6
    chainc = ' '
    if(ii.eq.1) chainc = 'Courant min'
    if(ii.eq.2) chainc = 'Courant max'
    if(ii.eq.3) chainc = 'Fourier min'
    if(ii.eq.4) chainc = 'Fourier max'
    if(ii.eq.5) chainc = 'Cou_Fou min'
    if(ii.eq.6) chainc = 'Cou_Fou max'
!    Compressible
    if(ippmod(icompf).ge.0) then
      if(ii.eq.5) chainc = 'CFL/Mas min'
      if(ii.eq.6) chainc = 'CFL/Mas max'
    endif

    ic = 13
    do jj = 1, 4
      chain = ' '
      write(chain,3000) ptploc(ii,jj)
      chainc(ic:ic+13) = chain(1:13)
      ic = ic+13
    enddo
    write(nfecra,'(A)') chainc(1:ic)
  enddo
  write(nfecra,1210)

  if(idtvar.gt.0) then
    do ii = 7, 8
      chainc = ' '
      if(ii.eq.7) chainc = 'Dt min'
      if(ii.eq.8) chainc = 'Dt max'
      ic = 13
      do jj = 1, 4
        chain = ' '
        write(chain,3000) ptploc(ii,jj)
        chainc(ic:ic+13) = chain(1:13)
        ic = ic+13
      enddo
      write(nfecra,'(A)') chainc(1:ic)
    enddo
    write(nfecra,1210)
  endif

  if (iptlro.eq.1) then
    chainc = ' '
    chainc = 'Dt/Dtrho max'
    ic = 13
    do jj = 1, 4
      chain = ' '
      write(chain,3000) rpdtro(jj)
      chainc(ic:ic+13) = chain(1:13)
      ic = ic+13
    enddo
    write(nfecra,'(A)') chainc(1:ic)
    chainc = ' '
    chainc = 'Clips a Dtrho'
    ic = 57
    chain = ' '
    write(chain,4000) nclptr
    chainc(ic:ic+7) = chain(1:7)
    ic = ic+7
    write(nfecra,'(A)') chainc(1:ic)
    write(nfecra,1210)
  endif

  write(nfecra,*) ' '
  write(nfecra,*) ' '

endif

!===============================================================================
! 3.
!===============================================================================
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


 1100 format (/,3X,'** INFORMATIONS SUR LES VARIABLES',/,         &
          3X,'   ------------------------------')
 1111 format ('   Variable      Valeur min    Valeur max',        &
        '   Clip min   Clip max')
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
 1211 format ('Critere      Valeur     en         xc',            &
        '           yc           zc')

 3000 format (e12.5)
 4000 format (i7)

#else

 1000 format (/,3X,'** INFORMATION ON CONVERGENCE',/,             &
          3X,'   --------------------------')
 1011 format ('   Variable    Rhs norm    ',                      &
        '  N_iter  Norm. residual      derive')
 1010 format ('---------------------------',                      &
        '------------------------------------')


 1100 format (/,3X,'** INFORMATION ON VARIABLES',/,               &
          3X,'   ------------------------')
 1111 format ('   Variable      Min. value    Max. value',        &
        '   Min clip   Max clip')
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
 1211 format ('Criterion    Value      at         xc',            &
        '           yc           zc')

 3000 format (e12.5)
 4000 format (i7)

#endif
return
end subroutine
