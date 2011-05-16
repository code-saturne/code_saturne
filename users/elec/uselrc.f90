!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine uselrc &
!================

 ( nvar   , nscal  ,                                              &
   ia     ,                                                       &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , viscf  , viscb  ,                            &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )

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
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! viscf(nfac)      ! tr ! --- ! tableau de travail    faces internes           !
! viscb(nfabor     ! tr ! --- ! tableau de travail    faces de bord            !
! w1..9(ncelet     ! tr ! --- ! tableau de travail    cellules                 !
! ra(*)            ! ra ! --- ! main real work array                           !
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

integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision viscf(nfac), viscb(nfabor)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision ra(*)

! Local variables

integer          iel    , ifac   , iutile
integer          ipcefj , ipcdc1 , ipcdc2 , ipcdc3 , ipcsig
integer          ipdcrp , idimve

double precision somje , coepoa , coefav , coepot
double precision emax  , aiex   , amex
double precision rayo  , econs  , z1     , z2   , posi
double precision dtj   , dtjm   , delhsh , cdtj , cpmx
double precision xelec , yelec  , zelec

!===============================================================================
!===============================================================================
! 1. INITIALISATION
!===============================================================================



!===============================================================================
! 2.  ARC ELECTRIQUE
!===============================================================================



if ( ippmod(ielarc).ge.1 ) then

!       3 Exemples : pour activer l'un des 3 mettre IUTILE= 1, 2 ou 3
!       Par defaut IUTILE = 1

  iutile = 1

! 2.1 : 1er exemple : cas general
! ===============================

  if ( iutile .eq. 1) then

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

! Fin 1er exemple

  else if ( iutile .eq. 2) then

! 2.2 : 2eme exemple : Autre methode de recalage
! ==============================================
!    Ceci est un cas particulier et doit etre adapte en fonction
!    du cas et du maillage (intervenir aussi dans uselcl)

!        Calcul de l'integrale sur le Volume de J.E
!        -----------------------------------
!        (c'est forcement positif ou nul)

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

    ipcdc3 = ipproc(idjr(3))
    elcou  = 0.d0
    do ifac = 1, nfac
      if( surfac(1,ifac).eq.0.d0 .and. surfac(2,ifac).eq.0.d0     &
           .and. cdgfac(3,ifac) .lt. 0.7d-2                       &
           .and. cdgfac(3,ifac) .gt. 0.65d-2 ) then
        iel = ifacel(1,ifac)
        elcou = elcou + propce(iel,ipcdc3) * surfac(3,ifac)
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
    cdtj= 2.0d2

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
    if(ntcabs.gt.5) then
      if(coepoa.ge.1.05d0 .and. coepot.le.cpmx) then
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

  else if ( iutile .eq. 3) then

! 2.3 :  3eme exemple : cas avec claquage
! =======================================
!    Ceci est un cas particulier et doit etre adapte en fonction
!    du cas et du maillage (intervenir aussi dans uselcl)


!        Utilisation d'une rampe d'intensite
!        -----------------------------------

    if ( ntcabs.le.200 ) then
      couimp = 200.d0
    endif

    if ( ntcabs.gt.200.and.ntcabs.le.400 ) then
      couimp = 200.d0 + 2 * (ntcabs-200)
    endif

    if ( ntcabs.gt.400 ) then
      couimp = 600.d0
    endif

!        UTILISANT D'UN CLAQUAGE AUTO
!        ----------------------------

    if(ntcabs.le.400.or.ntcabs.eq.ntpabs+1) iclaq = 0

    econs = 1.5d5

!        ON REPERE SI IL Y A CLAQUAGE ET SI OUI OU
!        -----------------------------------------

    if(ntcabs.ge.400 .and. iclaq .eq. 0 ) then
      amex = 1.d30
      aiex = -1.d30
      emax = 0.d0

!     les composantes du champ electrique : J/SIGMA

      ipcdc1 = ipproc(idjr(1))
      ipcdc2 = ipproc(idjr(2))
      ipcdc3 = ipproc(idjr(3))
      ipcsig = ipproc(ivisls(ipotr))

      do iel = 1, ncel

        xelec = propce(iel,ipcdc1)/propce(iel,ipcsig)
        yelec = propce(iel,ipcdc2)/propce(iel,ipcsig)
        zelec = propce(iel,ipcdc3)/propce(iel,ipcsig)

        w1(iel) = sqrt ( xelec**2 + yelec**2 + zelec**2 )

        amex =  min(amex,w1(iel))
        aiex =  max(aiex,w1(iel))
        if( w1(iel) .ge. econs) then
          WRITE(NFECRA,*) 'claquage ', NTCABS, W1(IEL)
          iclaq = 1
          ntdcla = ntcabs
          if(w1(iel).gt.emax) then
            xclaq =  xyzcen(1,iel)
            yclaq =  xyzcen(2,iel)
            zclaq =  xyzcen(3,iel)
            emax = w1(iel)
          endif
        endif
      enddo

      write(nfecra,*)
      WRITE(NFECRA,*) ' NT min et max de E = ',NTCABS,AMEX,AiEX
      write(nfecra,*)
      write(nfecra,*) xclaq,yclaq,zclaq,ntdcla

    endif

!        SI IL Y A CLAQUAGE : ON IMPOSE COLONNE CHAUDE DU CENTRE VERS
!        LE POINT DE CLAQUAGE
!        =============================================================

    if(iclaq .eq. 1) then
      if(ntcabs.le.ntdcla+30) then
        z1 = zclaq - 3.d-4
        if(z1.le.0.d0) z1 = 0.d0
        z2 = zclaq + 3.d-4
        if(z2.ge.2.d-2) z2 = 2.d-2

        do iel = 1, ncel

          if( xyzcen(3,iel).ge.z1 .and. xyzcen(3,iel).le.z2) then
            rayo = sqrt((xclaq*xyzcen(1,iel)-yclaq*xyzcen(2,iel)  &
                 /sqrt(xclaq**2+yclaq**2))**2+(xyzcen(3,iel)      &
                 -zclaq)**2)
            posi=xclaq*xyzcen(1,iel)
            if( rayo.le.5d-4 .and. posi.ge.0d0 ) then
!                    RTP(IEL,ISCA(IHM)) = 16.D6
              rtp(iel,isca(ihm)) = 8.d7
            endif
          endif
        enddo
      else
        iclaq = 0
      endif
    endif

!        Calcul de l'integrale sur le Volume de J.E
!        -----------------------------------
!        (c'est forcement positif ou nul)

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

    ipcdc3 = ipproc(idjr(3))
    elcou = 0.d0
    do ifac = 1, nfac
      if( surfac(1,ifac).eq.0.d0 .and. surfac(2,ifac).eq.0.d0     &
           .and. cdgfac(3,ifac) .gt. 0.05d-2                      &
           .and. cdgfac(3,ifac) .lt. 0.08d-2 ) then
        iel = ifacel(1,ifac)
        elcou = elcou + propce(iel,ipcdc3) * surfac(3,ifac)
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
    cdtj= 2.0d2

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

    cpmx= sqrt(cdtj*dtj)
    coepot=cpmx
    if(ntcabs.gt.5) then
      if(coepoa.ge.1.05d0 .and. coepot.le.cpmx) then
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

  else
    write(nfecra,5000) iutile
    call csexit(1)
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

 5000  format(/,' ERREUR DANS USELRC :',/,                        &
          ' VALEUR NON PERMISE DE IUTILE ',/,               &
          ' VERIFIER VOS DONNEES ')

!----
! FIN
!----

return
end subroutine
