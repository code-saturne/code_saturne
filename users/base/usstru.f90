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

! FONCTION :
! ----------

! GESTION DES STRUCTURES MOBILES EN ALE AVEC COUPLAGE INTERNE

! ON TROUVERA ICI DEUX ROUTINES DIFFERENTES :

! - USSTR1 : APPELEE A L'INITIALISATION, POUR DEFINIR LES STRUCTURES
!               ET LEURS PARAMETRES INITIAUX (VITESSE, DEPLACEMENT)

! - USSTR2 : APPELE A CHAQUE PAS DE TEMPS POUR DEFINIR LES
!               CARACTERISTIQUES (POTENTIELLEMENT VARIABLES) DES
!               STRUCTURES


! Boundary faces identification
! =============================

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


!-------------------------------------------------------------------------------
subroutine usstr1 &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr , idfstr ,                   &
   idevel , ituser , ia     ,                                     &
   aexxst , bexxst , cfopre ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   xstr0  , vstr0  , xstreq ,                                     &
   rdevel , rtuser , ra     )



!===============================================================================
! FONCTION :
! ----------

! GESTION DES STRUCTURES MOBILES EN ALE AVEC COUPLAGE INTERNE

!   DEFINITION DES STRUCTURES
!   DEPLACEMENT ET VITESSE INITIAUX

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! maxelt           ! i  ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfac+1)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (lndfac)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfabor+1)     !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (lndfbr)       !    !     !  (optionnel)                                   !
! idfstr(nfabor    ! te ! <-- ! definition des structures                      !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! aexxst,bexxst    ! r  ! <-- ! coefficients de prediction du deplact          !
! cfopre           ! r  ! <-- ! coeff de prediction des efforts                !
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
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! xstr0(ndim,      ! tr ! <-- ! deplacement initial des structures             !
!       nbstru)    !    !     !                                                !
! vstr0(ndim,      ! tr ! <-- ! vitesse initiale des structures                !
!       nbstru)    !    !     !                                                !
! xstreq(ndim,     ! tr ! <-- ! deplacement du maillage initial par            !
!       nbstru)    !    !     ! rapport a l'equilibre                          !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "cstnum.h"
include "optcal.h"
include "entsor.h"
include "pointe.h"
include "albase.h"
include "period.h"
include "parall.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nbstru
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idfstr(nfabor)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision aexxst, bexxst, cfopre
double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision xstr0(3,nstrmx), xstreq(3,nstrmx)
double precision vstr0(3,nstrmx)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac
integer          ilelt, nlelt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START


if(1.eq.1) then
  nbstru = 0
  return
endif

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  INITIALISATIONS

!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2.  DEFINITION DES STRUCTURES (appel a l'initialisation)
!===============================================================================

!     On remplit le tableau IDFSTR(NFABOR)
!     IDFSTR(IFAC) est le numero de la structure a laquelle appartient
!       la face de bord IFAC (0 si elle n'appartient a aucune structure)
!     Le nombre de structure est automatiquement determine a partir du
!       plus grand element de IDFSTR (les numeros des structures doivent
!       donc etre affectes de maniere sequentielle sans trou en commencant
!       par 1).

! Dans l'exemple ci-dessous la structure 1 est bordee par les faces de
!   couleur 4, la structure 2 par les faces de couleur 6


CALL GETFBR('4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = 1

enddo

CALL GETFBR('6',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  idfstr(ifac) = 2

enddo


! Pour chaque structure, on definit si necessaire
! - une vitesse initiale VSTR0
! - un deplacement initial XSTR0 (i.e. la valeur du deplacement XSTR
!     au temps t=0 par rapport au maillage initial)
! - un deplacement par rapport a l'equilibre XSTREQ (i.e. le deplacement
!     du maillage initial par rapport a la position d'equilibre de la
!     structure ; la force de retour exercee par le ressort a un temps
!     donne pour un deplacement XSTR sera donc -k*(XSTR+XSTREQ) ).
! Toutes les composantes de XSTR0, XSTREQ et VSTR0 sont initialisees a 0

!     Exemple : deplacement initial y=2 pour la structure 1
!               deplacement par rapport a l'equilibre yeq=1 pour la
!                  structure 1
!               vitesse initiale uz=0.5 pour la structure 2

! Dans le cas d'un calcul initial, ou d'une suite d'un calcul sans ALE,
!   une iteration 0 est automatiquement realisee pour gerer un eventuel
!   deplacement initial des structures. Si necessaire, positionner
!   ITALIN a 1 dans usalin pour activer une iteration 0 dans les autres
!   cas.

xstr0(2,1)  = 2.d0
xstreq(2,1) = 1.d0
vstr0(3,2)  =-0.5d0

! --- Si necessaire on definit les coefficients d'extrapolation utiles
!       en couplage explicite :
!       deplacement predit = X(n) + AEXXST.DT.X'(n)
!                                 + BEXXST.DT.( X'(n)-X'(n-1) )
!       force envoyee a la structure = CFOPRE.F(n) + (1.D0-CFOPRE).F(n-1)

aexxst =  0.5d0
bexxst =  0.0d0
cfopre =  2.d0

! --- Ecriture des fichiers historiques des structures mobiles
!     (deplacement, vitesse, acceleration, force)
!     La periodicite de sortie est la meme que pour les historiques
!     standards (NTHIST)
ihistr = 1

return

end subroutine


!===============================================================================


subroutine usstr2 &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nbstru ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , idfstr ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dtcel  ,                                                       &
   xmstru , xcstru , xkstru , xstreq , xstr   , vstr   , forstr , &
   dtstr  ,                                                       &
   rdevel , rtuser ,                                              &
   ra     )



!===============================================================================
! FONCTION :
! ----------

! GESTION DES STRUCTURES MOBILES EN ALE AVEC COUPLAGE INTERNE

! SPECIFICATION DES CARACTERISTIQUES DES STRUCTURES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nbstru           ! e  ! <-- ! nombre de structures definies                  !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfac+1)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (lndfac)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfabor+1)     !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (lndfbr)       !    !     !  (optionnel)                                   !
! idfstr(nfabor    ! te ! <-- ! definition des structures                      !
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
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dtcel(ncelet)    ! tr ! <-- ! pas de temps dans les cellules                 !
! xmstru(ndim,     ! tr ! --> ! matrice de masse des structures                !
!  ndim,nbstru)    !    !     !                                                !
! xcstru(ndim,     ! tr ! --> ! matrice de friction des structures             !
!  ndim,nbstru)    !    !     !                                                !
! xkstru(ndim,     ! tr ! --> ! matrice de raideur des structures              !
!  ndim,nbstru)    !    !     !                                                !
! xstreq(ndim,     ! tr ! <-- ! deplacement du maillage initial par            !
!       nbstru)    !    !     ! rapport a l'equilibre                          !
! xstr(ndim,       ! tr ! <-- ! deplacement des structures                     !
!       nbstru)    !    !     !                                                !
! vstr(ndim,       ! tr ! <-- ! vitesse  des structures                        !
!       nbstru)    !    !     !                                                !
! forstr(ndim      ! tr ! <-- ! effort sur les structures (contient            !
!       nbstru)    !    !     !            les efforts dus au fluide)          !
! dtstr(nbstru)    ! tr ! --> ! pas de temps des structures                    !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "cstnum.h"
include "optcal.h"
include "pointe.h"
include "albase.h"
include "period.h"
include "parall.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nbstru
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idfstr(nfabor)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dtcel(ncelet)
double precision xmstru(3,3,nstrmx)
double precision xcstru(3,3,nstrmx)
double precision xkstru(3,3,nstrmx)
double precision xstreq(3,nstrmx)
double precision xstr(3,nstrmx)
double precision vstr(3,nstrmx)
double precision forstr(3,nstrmx)
double precision dtstr(nstrmx)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ii, jj, istr
double precision theta, sint, cost, xm, xc, xk, fx, fy

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START


if(1.eq.1) return

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  INITIALISATIONS

!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2.  COEFFICIENTS DES STRUCTURES (appel a chaque pas de temps)
!===============================================================================


!     On remplit ici les coefficients de definissant la structure.
!      - sa masse M                     (XMSTRU)
!      - son coefficient de friction C  (XCSTRU)
!      - son coefficient de raideur K   (XKSTRU)

!     Le tableau FORSTR contient les efforts exerces par le fluide
!       sur chacune des structures. Il est possible d'y ajouter une
!       composante de force exterieure (gravite par exemple)

!     Le tableau XSTR contient le deplacement des structures par rapport
!       au maillage initial
!     Le tableau XSTR0 contient le deplacement des structures dans le
!       maillage initial par rapport a leur position d'equilibre
!     Le tableau XPSTR contient la vitesse des structures

!     XSTR, XSTR0 et VSTR sont des DONNEES pouvant servir eventuellement
!     a calculer M, C et K. ILS NE DOIVENT PAS ETRE MODIFIES.

!     L'evolution du systeme est resolue de maniere tridimensionnelle,
!       ces trois coefficients sont donc en fait des matrices 3x3.

!     L'equation resolue est

!       M.X'' + C.X' + K.(X+X0) = F
!       = -     = -    =  - --    -

!       (X est le deplacement par rapport a la position initiale du maillage,
!        X0 est le deplacement de la position dans le maillage initial par
!        rapport a l aposition d'equilibre)

!     La resolution est effectuee par la methode de Newmark HHT.
!     Le pas de temps utilise peut etre different du pas de temps
!       fluide (a definir dans le tableau DTSTR, initialise par defaut
!       au pas de temps fluide).


!     On met a zero tous les coefficients
do istr = 1, nbstru

  do ii = 1, 3
    do jj = 1, 3
      xmstru(ii,jj,istr) = 0.d0
      xcstru(ii,jj,istr) = 0.d0
      xkstru(ii,jj,istr) = 0.d0
    enddo
  enddo

enddo

! Dans l'exemple ci-dessous, la structure 1, de masse 5 kg est retenue par
!   ressort isotrope de raideur 2 N/m et de coefficient de friction 3 kg.s

do ii = 1, 3
  xmstru(ii,ii,1) = 5.d0
  xcstru(ii,ii,1) = 2.d0
  xkstru(ii,ii,1) = 3.d0
enddo


! Dans l'exemple ci-dessous, la structure 2 est contrainte a un mouvement
!   decompose en deux :
!  - dans le plan xOy, le mouvement est force dans une direction X, avec
!    une masse xm, une friction xc et une raideur xk (et la composante
!    normale Y est donc forcee a 0). L'axe X est incline d'un angle THETA
!    par rapport a l'axe x du repere global.
!  - dans la direction z, le mouvement est un mouvement d'oscillation
!    harmonique de masse 1 et de raideur 1 (et de friction nulle) avec un
!    forcage externe en 3.cos(4.t) (en plus des efforts fluides).

theta = pi/6.d0
cost = cos(theta)
sint = sin(theta)


! Dans le repere local on a donc
!      xm.X'' + xc.X' + xk.X = FX
!                          Y = 0
!         Z''         +    Z = FZ + 3.cos(4.t)

!   FX, FY et FZ sont les composantes des efforts fluides dans le repere
!     local. Soit, a partir des composantes dans le repere global :
!     FX =  COST*Fx + SINT*Fy
!     FY = -SINT*Fx + COST*Fy
!     FZ = Fz

! Apres changement de repere, on obtient donc :

xm = 1.d0
xc = 3.d-1
xk = 2.d0
fx = forstr(1,2)
fy = forstr(2,2)

xmstru(1,1,2) = xm*cost**2
xmstru(1,2,2) = xm*cost*sint
xmstru(2,1,2) = xm*cost*sint
xmstru(2,2,2) = xm*sint**2
xmstru(3,3,2) = 1.d0

xcstru(1,1,2) = xc*cost**2
xcstru(1,2,2) = xc*cost*sint
xcstru(2,1,2) = xc*cost*sint
xcstru(2,2,2) = xc*sint**2

xkstru(1,1,2) = (xk-1.d0)*cost**2 + 1.d0
xkstru(1,2,2) = (xk-1.d0)*cost*sint
xkstru(2,1,2) = (xk-1.d0)*cost*sint
xkstru(2,2,2) = (xk-1.d0)*sint**2 + 1.d0
xkstru(3,3,2) = 1.d0

forstr(1,2) = fx*cost**2   + fy*sint*cost
forstr(2,2) = fx*sint*cost + fy*sint**2
forstr(3,2) = forstr(3,2) + 3.d0*cos(4.d0*ttcabs)

do istr = 1, nbstru
  dtstr(istr) = dtcel(1)
enddo


return

end subroutine
