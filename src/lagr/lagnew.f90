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

subroutine lagnew &
!================

 ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   npt    , nptnew , new    ,                                     &
   izone  ,                                                       &
   ifrlag , isorti , iworkp ,                                     &
   ettp   )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!   TIRAGE ALEATOIRE DE LA POSITION DES PARTICULES A INJECTER
!   DANS LE DOMAINE DE CALCUL AU NIVEAU DES FACES
!   D'UNE ZONE DE COULEUR : POSITIONS + REPERAGE DE LA CELLULE
!   ATTENTION : CE TIRAGE ALEATOIRE NE FONCTIONNE QU'AVEC DES FACES
!   TRIANGULAIRE OU QUADRANGULAIRES.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! npt              ! e  ! --> ! nombre courant de particules                   !
! nptnew           ! e  ! <-- ! nombre total de nouvelles particules           !
!                  !    !     ! pour toutes les zones d'injection              !
! new              ! e  ! <-- ! nombre de nouvelles part a injecter            !
!                  !    !     ! pour la zone d'injection courante              !
! izone            ! e  ! <-- ! numero  de la zone d'injection                 !
! ifrlag           ! te ! <-- ! numero de zone de la face de bord              !
!   (nfabor)       !    !     !  pour le module lagrangien                     !
! isorti           ! te ! --> ! pour chaque particule :                        !
!   (nbpmax)       !    !     !    * numero de sa cellule                      !
!                  !    !     !    * 0 si sortie du domaine                    !
! iworkp(nbpmax    ! te ! --> ! numero de la face d'injection                  !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!==============================================================================
! Module files
!==============================================================================

use paramx
use numvar
use optcal
use entsor
use lagpar
use lagran
use mesh

!==============================================================================

implicit none

! Arguments

integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          npt    , nptnew , new
integer          izone

integer          ifrlag(nfabor) , isorti(nbpmax) , iworkp(nbpmax)

double precision ettp(nbpmax,nvp)

! Local variables

integer          nnn , nnnn , ifac , np , ii
integer          ifrbr , minfac , maxfac
integer          iconfo(100)

double precision surfm , rda(1) , rd , eps , pm1 , pm6
double precision ctr(6,3) , vec(3) , are(2,3) , surftr(2)

!===============================================================================

eps = 1.d-3

!     CALCUL DE LA SURFACE MAX DE L'ENTREE :

surfm = -10.d0
minfac = nfabor+1
maxfac = 0

do ifac = 1,nfabor
  ifrbr = ifrlag(ifac)
  if (ifrbr.eq.izone) then
    surfm = max( surfm,surfbn(ifac) )
    minfac = min (ifac,minfac)
    maxfac = max (ifac,maxfac)
  endif
enddo

!     STOP PAR SECURITE

if (maxfac.eq.0.or.minfac.eq.nfabor+1) then
  write(nfecra,9000) izone
  call csexit (1)
  !==========
endif

!     BOUCLE SUR LES NOUVELLES PARTICULES :

do np = 1,new

!       incrementation du pointeur sur les particules

  npt = npt + 1

!       tirage aleatoire d'une face :

 100    continue

    nnn = 1
    call zufall(nnn,rda)
    rd=rda(1)

    rd = rd*(dble(maxfac-minfac+1)-eps)
    ifac = minfac + int(rd)

    if (ifac .lt. minfac .or. ifac .gt. maxfac)  goto 100

!         numero de la face :

    ifrbr = ifrlag(ifac)

    if (ifrbr.ne.izone) goto 100

!         tirage aleatoire pour determiner si cette face convient
!         plus la fa7 est grande plus elle a des chance d'etre choisie

    nnn = 1
    call zufall(nnn,rda)
    rd=rda(1)

    if (rd.gt.(surfbn(ifac)/surfm)) goto 100

!    ATTENTION :

!         type de face : 3 ou 4 points supports
!         pour l'instant je ne sais pas traiter les autres
!         avec plus de points supports...

    ii = ipnfbr(ifac+1)-ipnfbr(ifac)
    if (ii.gt.4) goto 100

!        si face a 4 points, on choisit l'un des deux triangles

  if (ii.eq.4) then

    ii = 0
    do nnnn = ipnfbr(ifac),ipnfbr(ifac+1)-1
      ii = ii+1
      iconfo(ii) = nodfbr(nnnn)
    enddo

!        longueur des arretes 1 et 2 du premier triangle :

    are(1,1) = xyznod(1,iconfo(2)) - xyznod(1,iconfo(1))
    are(1,2) = xyznod(2,iconfo(2)) - xyznod(2,iconfo(1))
    are(1,3) = xyznod(3,iconfo(2)) - xyznod(3,iconfo(1))
    are(2,1) = xyznod(1,iconfo(3)) - xyznod(1,iconfo(1))
    are(2,2) = xyznod(2,iconfo(3)) - xyznod(2,iconfo(1))
    are(2,3) = xyznod(3,iconfo(3)) - xyznod(3,iconfo(1))

!        surface du premier triangle

    vec(1) = are(1,2)*are(2,3)-are(1,3)*are(2,2)
    vec(2) = are(1,3)*are(2,1)-are(1,1)*are(2,3)
    vec(3) = are(1,1)*are(2,2)-are(1,2)*are(2,1)
    surftr(1) = sqrt(vec(1)**2+vec(2)**2+vec(3)**2)

!        longueur des arretes 1 et 2 du deuxieme triangle :

    are(1,1) = xyznod(1,iconfo(3)) - xyznod(1,iconfo(1))
    are(1,2) = xyznod(2,iconfo(3)) - xyznod(2,iconfo(1))
    are(1,3) = xyznod(3,iconfo(3)) - xyznod(3,iconfo(1))
    are(2,1) = xyznod(1,iconfo(4)) - xyznod(1,iconfo(1))
    are(2,2) = xyznod(2,iconfo(4)) - xyznod(2,iconfo(1))
    are(2,3) = xyznod(3,iconfo(4)) - xyznod(3,iconfo(1))

!        surface du deuxieme triangle

    vec(1) = are(1,2)*are(2,3) - are(1,3)*are(2,2)
    vec(2) = are(1,3)*are(2,1) - are(1,1)*are(2,3)
    vec(3) = are(1,1)*are(2,2) - are(1,2)*are(2,1)
    surftr(2) = sqrt(vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3))

!        tirage d'un nombre aleatoire entre 0 et 1
!        pour determiner quel triangle choisir

    nnn = 1
    call zufall(nnn,rda)
    rd=rda(1)

!        si le deuxieme triangle est choisit, on reorganise
!        les points : 4 <--> 2

    if (rd.le.( surftr(2) / (surftr(1)+surftr(2)) )) then
      nnnn = iconfo(4)
      iconfo(4) = iconfo(2)
      iconfo(2) = nnnn
    endif

! dans le cas ou la face est un triangle...

  else if (ii.eq.3) then

    ii = 0
    do nnnn = ipnfbr(ifac),ipnfbr(ifac+1)-1
      ii = ii+1
      iconfo(ii) = nodfbr(nnnn)
    enddo

  endif

!        constitution des coordonnees du triangle

  do nnnn = 1,3
    ctr(nnnn,1) = xyznod(1,iconfo(nnnn))
    ctr(nnnn,2) = xyznod(2,iconfo(nnnn))
    ctr(nnnn,3) = xyznod(3,iconfo(nnnn))
  enddo

!        tirage aleatoire d'un point dans le triangle
!        constitue des points 1,2,3

 200    continue

!        1) tirage du point 4 sur l'arete 12

 300      continue
    nnn = 1
    call zufall(nnn,rda)
    rd=rda(1)
    if (rd.eq.0.d0 .or. rd.eq.1.d0) goto 300

    do nnnn = 1,3
      ctr(4,nnnn) = rd*ctr(1,nnnn) + (1.d0-rd)*ctr(2,nnnn)
    enddo

!        2) tirage du point 5 sur l'arete 13

 400      continue
    nnn = 1
    call zufall(nnn,rda)
    rd=rda(1)
    if (rd.eq.0.d0 .or. rd.eq.1.d0) goto 400

    do nnnn = 1,3
      ctr(5,nnnn) = rd*ctr(1,nnnn) +( 1.d0-rd)*ctr(3,nnnn)
    enddo

!        3) le point 6 est le sommet du parallelogramme 1465

    do nnnn = 1,3
      ctr(6,nnnn) = ctr(4,nnnn) + ctr(5,nnnn) - ctr(1,nnnn)
    enddo

!        4) reste a verifier que le point 6 appartient au triangle 123

!        4.1) vecteur normal au triangle : 12^13

    vec(1) = (ctr(2,2)-ctr(1,2))*(ctr(3,3)-ctr(1,3))              &
           - (ctr(2,3)-ctr(1,3))*(ctr(3,2)-ctr(1,2))
    vec(2) = (ctr(2,3)-ctr(1,3))*(ctr(3,1)-ctr(1,1))              &
           - (ctr(2,1)-ctr(1,1))*(ctr(3,3)-ctr(1,3))
    vec(3) = (ctr(2,1)-ctr(1,1))*(ctr(3,2)-ctr(1,2))              &
           - (ctr(2,2)-ctr(1,2))*(ctr(3,1)-ctr(1,1))

!        4.2) produit mixte pour le point 1 :

    pm1 = 0.d0
    pm1 = pm1 + vec(1) *                                          &
        ( (ctr(2,2)-ctr(1,2))*(ctr(3,3)-ctr(2,3))                 &
         -(ctr(2,3)-ctr(1,3))*(ctr(3,2)-ctr(2,2)) )
    pm1 = pm1 + vec(2) *                                          &
        ( (ctr(2,3)-ctr(1,3))*(ctr(3,1)-ctr(2,1))                 &
         -(ctr(2,1)-ctr(1,1))*(ctr(3,3)-ctr(2,3)) )
    pm1 = pm1 + vec(3) *                                          &
        ( (ctr(2,1)-ctr(1,1))*(ctr(3,2)-ctr(2,2))                 &
         -(ctr(2,2)-ctr(1,2))*(ctr(3,1)-ctr(2,1)) )

!        4.3) produit mixte pour le point 6 :

    pm6 = 0.d0
    pm6 = pm6 + vec(1) *                                          &
        ( (ctr(2,2)-ctr(6,2))*(ctr(3,3)-ctr(2,3))                 &
         -(ctr(2,3)-ctr(6,3))*(ctr(3,2)-ctr(2,2)) )
    pm6 = pm6 + vec(2) *                                          &
        ( (ctr(2,3)-ctr(6,3))*(ctr(3,1)-ctr(2,1))                 &
         -(ctr(2,1)-ctr(6,1))*(ctr(3,3)-ctr(2,3)) )
    pm6 = pm6 + vec(3) *                                          &
        ( (ctr(2,1)-ctr(6,1))*(ctr(3,2)-ctr(2,2))                 &
         -(ctr(2,2)-ctr(6,2))*(ctr(3,1)-ctr(2,1)) )

!        4.4) 6 est dans le triangle si PM1*PM6>=0

    if (pm1*pm6.lt.0.d0) goto 200

!        5) POUR PLUS DE SECURITE, ON DEPLACE LE POINT
!           D'UN EPSILON EN DIRECTION DU CENTRE CELLULE

! ATTENTION : CE DECALAGE PEUT ETRE DANGEREUX DANS LE CAS
!             DE CELULES CONCAVES

  ctr(6,1) = ctr(6,1) + (xyzcen(1,ifabor(ifac))-ctr(6,1))*eps
  ctr(6,2) = ctr(6,2) + (xyzcen(2,ifabor(ifac))-ctr(6,2))*eps
  ctr(6,3) = ctr(6,3) + (xyzcen(3,ifabor(ifac))-ctr(6,3))*eps


!        LE TRAITEMENT EST TERMINE POUR LE POINT NPT,
!        ON REMPLIT LES TABLEAUX POUR LE LAGRANGIEN :

  ettp(npt,jxp) = ctr(6,1)
  ettp(npt,jyp) = ctr(6,2)
  ettp(npt,jzp) = ctr(6,3)

  isorti(npt) = ifabor(ifac)
  iworkp(npt) = ifac

enddo

!===============================================================================

!-------
! FORMAT
!-------

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ',/,&
'@    =========   (LAGNEW).                                   ',/,&
'@                                                            ',/,&
'@    PROBLEME DANS LA GESTION DE NOUVELLES PARTICULES        ',/,&
'@                                                            ',/,&
'@  Le nombre de faces de la zone ',I10                        ,/,&
'@    est egal a zero.                                        ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Contacter l''equipe de developpement.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return

end subroutine
