!-------------------------------------------------------------------------------

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

subroutine lagnew &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndnod , lndfac , lndfbr , ncelbr ,                   &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   npt    , nptnew , new    ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   izone  ,                                                       &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   ifrlag , isorti , iworkp ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   surfbn , ettp   ,                                              &
   rdevel , rtuser , ra     )


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
!                  !    !     !                                                !
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndnod           ! e  ! <-- ! dim. connectivite cellules->faces              !
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
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
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! izone            ! e  ! <-- ! numero  de la zone d'injection                 !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
!  (nfml,nprfml    !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac                      !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !                                                !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr                      !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !                                                !
! ifrlag           ! te ! <-- ! numero de zone de la face de bord              !
!   (nfabor)       !    !     !  pour le module lagrangien                     !
! isorti           ! te ! --> ! pour chaque particule :                        !
!   (nbpmax)       !    !     !    * numero de sa cellule                      !
!                  !    !     !    * 0 si sortie du domaine                    !
! iworkp(nbpmax    ! te ! --> ! numero de la face d'injection                  !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! volume(ncelet    ! tr ! <-- ! volume d'un des ncelet elements                !
! surfbn(nfabor    ! tr ! <-- ! surface des faces de bord                      !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

implicit none

!==============================================================================
! Common blocks
!==============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"

!==============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndnod , lndfac , lndfbr , ncelbr
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          npt    , nptnew , new
integer          nideve , nrdeve , nituse , nrtuse
integer          izone
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          ifrlag(nfabor) , isorti(nbpmax) , iworkp(nbpmax)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision surfbn(nfabor)
double precision ettp(nbpmax,nvp)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

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
