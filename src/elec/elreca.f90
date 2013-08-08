!-------------------------------------------------------------------------------

!VERS

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

subroutine elreca &
!================

 ( nvar   , nscal  ,                                              &
   dt     , rtpa   , rtp    , propce , propfb ,                   &
   coefa  , coefb  )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE POUR LE MODULE ELECTRIQUE

!             CALCULS DU COEFFICIENT DE RECALAGE
!               POUR LES VARIABLES ELECTIQUES
!             RECALAGE DES VARIABLES ELECTRIQUES
!               EN FONCTION DE CE COEFFICIENT

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use elincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables

integer          iel    , ifac
integer          ipcefj , ipcdc1 , ipcdc2 , ipcdc3 , ipcsig
integer          ipdcrp , idimve, idir

double precision somje , coepoa , coefav , coepot
double precision emax  , aiex   , amex
double precision rayo  , econs  , z1     , z2   , posi
double precision dtj   , dtjm   , delhsh , cdtj , cpmx
double precision xelec , yelec  , zelec

double precision, allocatable, dimension(:) :: w1

logical          ok

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================



!===============================================================================
! 2.  ARC ELECTRIQUE
!===============================================================================



if ( ippmod(ielarc).ge.1 ) then

! 2.1 :  cas general
! ===============================

  if ( modrec .eq. 1) then

!       CALCUL DU COEFFICIENT DE RECALAGE
!       -------------------------------

!  Calcul de l'integrale sur le Volume de J.E
!     (c'est forcement positif ou nul)

    ipcefj = ipproc(iefjou)
    somje = 0.d0
    do iel = 1, ncel
      somje = somje+propce(iel,ipcefj)*volume(iel)
    enddo

    if(irangp.ge.0) then
      call parsom (somje)
    endif

    coepot = couimp*dpot/max(somje,epzero)

    coepoa = coepot

!  On impose COEPOT >= 0.75 et COEPOT <= 1.5

    if ( coepot .gt. 1.50d0 ) coepot = 1.50d0
    if ( coepot .lt. 0.75d0 ) coepot = 0.75d0

    write(nfecra,1000)coepoa,coepot
 1000     format(/,                                               &
 ' Courant impose/Courant= ',E14.5,', Coeff. recalage= ',E14.5)

!       RECALAGE DES VARIABLES ELECTRIQUES
!       ---------------------------------------

!        Valeur de DPOT
!        --------------

    dpot = dpot*coepot

!        Potentiel Electrique (on pourrait eviter ; c'est pour le post)
!        --------------------

    do iel = 1, ncel
      rtp(iel,isca(ipotr)) = rtp(iel,isca(ipotr))*coepot
    enddo


!        Densite de courant (sert pour A et pour jXB)
!        ------------------

    if(ippmod(ielarc).ge.1 ) then
      do idimve = 1, ndimve
        ipdcrp = ipproc(idjr(idimve))
        do iel = 1, ncel
          propce(iel,ipdcrp) = propce(iel,ipdcrp) * coepot
        enddo
      enddo
    endif

!        Effet Joule (sert pour H au pas de temps suivant)
!        -----------

    ipcefj = ipproc(iefjou)
    do iel = 1, ncel
      propce(iel,ipcefj) = propce(iel,ipcefj)*coepot**2
    enddo

  else if ( modrec .eq. 2) then

! 2.2 : 2eme exemple : Autre methode de recalage
! ==============================================
!    Ceci est un cas particulier et doit etre adapte en fonction
!    du cas et du maillage (intervenir aussi dans uselcl)

!        Calcul de l'integrale sur le Volume de J.E
!        -----------------------------------
!        (c'est forcement positif ou nul)

    call uielrc(ncelet, izreca, crit_reca)

    ipcefj = ipproc(iefjou)
    somje = 0.d0
    do iel = 1, ncel
      somje = somje+propce(iel,ipcefj)*volume(iel)
    enddo

    if(irangp.ge.0) then
      call parsom (somje)
    endif

    if (somje .ne. 0) then
      coepot = couimp*dpot/max(somje,epzero)
    endif
    write(nfecra,1001) couimp,dpot,somje

!        Calcul de l'intensite du courant d'arc
!        --------------------------------------
!          Calcul de l'integrale de J sur une surface plane
!          perpendiculaire a l'axe de l'arc

!       ATTENTION : changer la valeur des tests sur CDGFAC(3,IFAC)
!                   en fonction du maillage

    ipcdc3 = ipproc(idjr(idreca))
    elcou  = 0.d0
    do ifac = 1, nfac
      if (izreca(ifac) .gt. 0) then
        ok = .true.
        do idir = 1, 3
          if (abs(surfac(idir, ifac)) .gt. 0.d0 .and. idir.ne.idreca) then
            ok = .false.
          endif
        enddo
        if (ok .eqv. .true.) then
          iel = ifacel(1,ifac)
          elcou = elcou + propce(iel,ipcdc3) * surfac(idreca,ifac)
        endif
      endif
    enddo

    if(irangp.ge.0) then
      call parsom (elcou)
    endif

    if ( abs(elcou).ge.1.d-06 ) then
      elcou=abs(elcou)
    else
      elcou=0.d0
    endif
    if(elcou.ne.0.d0) coepoa = couimp/elcou
    coepot = coepoa

    WRITE(NFECRA,*) ' ELCOU = ',ELCOU

    dtj = 1.d15
    dtjm =dtj
    delhsh = 0.d0
    cdtj= 20.d0

    do iel = 1, ncel
      if(propce(iel,ipproc(irom)).ne.0.d0)                     &
           delhsh =  propce(iel,ipcefj) * dt(iel)                 &
           /propce(iel,ipproc(irom))

      if(delhsh.ne.0.d0) then
        dtjm= rtp(iel,isca(iscalt))/delhsh
      else
        dtjm= dtj
      endif
      dtjm=abs(dtjm)
      dtj =min(dtj,dtjm)
    enddo
    if(irangp.ge.0) then
      call parmin (dtj)
    endif
    WRITE(NFECRA,*) ' DTJ = ',DTJ

    cpmx= sqrt(cdtj*dtj)
    coepot=cpmx
    if(ntcabs.gt.2) then
      if(coepoa.ge.1.05d0) then
        coepot=cpmx
      else
        coepot=coepoa
      endif
    endif

    write(nfecra,1008)cpmx,coepoa,coepot
    write(nfecra,1009)elcou,dpot*coepot

!        RECALAGE DES VARIABLES ELECTRIQUES
!        ----------------------------------

!         Valeur de DPOT
!         --------------

    dpot = dpot*coepot

!         Potentiel Electrique (on pourrait eviter ; c'est pour le post)
!         --------------------

    do iel = 1, ncel
      rtp(iel,isca(ipotr)) = rtp(iel,isca(ipotr))*coepot
    enddo


!      Densite de courant (sert pour A et pour jXB)
!      ------------------

    if(ippmod(ielarc).ge.1 ) then
      do idimve = 1, ndimve
        do iel = 1, ncel
          ipdcrp = ipproc(idjr(idimve))
          propce(iel,ipdcrp) = propce(iel,ipdcrp) * coepot
        enddo
      enddo
    endif

!      Effet Joule (sert pour H au pas de temps suivant)
!      -----------

    ipcefj = ipproc(iefjou)
    do iel = 1, ncel
      propce(iel,ipcefj) = propce(iel,ipcefj)*coepot**2
    enddo
  endif
endif

!===============================================================================
! 3.  EFFET JOULE
!===============================================================================

if ( ippmod(ieljou).ge.1 ) then

! 3.1  CALCUL DU COEFFICIENT DE RECALAGE
! --------------------------------------

!  Calcul de l'integrale sur le Volume de J.E
!     (c'est forcement positif ou nul)

  ipcefj = ipproc(iefjou)
  somje = 0.d0
  do iel = 1, ncel
    somje = somje+propce(iel,ipcefj)*volume(iel)
  enddo

  if(irangp.ge.0) then
    call parsom (somje)
  endif

  coepot = sqrt(puisim/max(somje,epzero))

  coefav = coepot

!  On impose COEF >= 0.75 et COEF <= 1.5

  if ( coepot .gt. 1.50d0 ) coepot = 1.50d0
  if ( coepot .lt. 0.75d0 ) coepot = 0.75d0

  write(nfecra,2000)coefav,coejou
 2000   format(/,                                                 &
 ' Puissance impose/Somme jE= ',E14.5,', Coeff. recalage= ',E14.5)


! 3.2  RECALAGE DES VARIABLES JOULE
! ---------------------------------

!       Valeur de DPOT (au cas ou utile)
!       --------------

  dpot = dpot*coepot

!       Coefficient correcteur COEJOU cumule
!       ------------------------------------

  coejou = coejou*coepot

!       Potentiel Electrique (on pourrait eviter ; c'est pour le post)
!       --------------------

  if ( ippmod(ieljou).ne.3 .and. ippmod(ieljou).ne.4 ) then
    do iel = 1, ncel
      rtp(iel,isca(ipotr)) = rtp(iel,isca(ipotr))*coepot
    enddo
  endif

!      Potentiel complexe (on pourrait eviter ; c'est pour le post)
!      -----------------

  if ( ippmod(ieljou).eq.2 ) then
    do iel = 1, ncel
      rtp(iel,isca(ipoti)) = rtp(iel,isca(ipoti))*coepot
    enddo
  endif

!      Effet Joule (sert pour H au pas de temps suivant)
!      -----------

  ipcefj = ipproc(iefjou)
  do iel = 1, ncel
    propce(iel,ipcefj) = propce(iel,ipcefj)*coepot**2
  enddo

endif

!--------
! FORMATS
!--------

 1001  format(/, ' Courant impose= ',E14.5, /,                    &
              ' Dpot= ',E14.5,/,                            &
              ' Somje= ',E14.5)

 1008  format(/,' Cpmx   = ',E14.5,/,                             &
          ' COEPOA = ',E14.5,/,                             &
          ' COEPOT = ',E14.5)

 1009  format(/,' Courant calcule = ',E14.5,/,                    &
          ' Dpot recale     = ',E14.5)

!----
! FIN
!----

return
end subroutine
