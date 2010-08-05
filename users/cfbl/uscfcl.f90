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

subroutine uscfcl &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icodcl , itrifb , itypfb , izfppp ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! --------

!    ROUTINE UTILISATEUR
!    REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
!    (ICODCL,RCODCL) POUR LES VARIABLES INCONNUES



!    CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE
!    =============================================


! Introduction
! ============

! Here one defines boundary conditions on a per-face basis.

! Boundary faces may be identified using the 'getfbr' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


! Boundary condition types
! ========================

! Boundary conditions setup for standard variables (pressure, velocity,
! turbulence, scalars) is described precisely in the 'usclim' subroutine.

! Detailed explanation will be found in the theory guide.


! TYPE DE CONDITIONS AUX LIMITES
! ==============================

! En compressible, on ne peut affecter que les conditions aux
!  limites predefinies

!    IPAROI, ISYMET, IESICF, ISSPCF, ISOPCF, IERUCF, IEQHCF

!    IPAROI : paroi standard
!    ISYMET : symetrie standard

!    IESICF, ISSPCF, ISOPCF, IERUCF, IEQHCF : entree/sortie

! Pour les entrees/sorties, on peut
!  imposer une valeur pour la turbulence et les scalaires
!  passifs dans RCODCL(.,.,1) pour le cas ou le flux de masse
!  serait entrant. Si on ne le fait pas, une condition de
!  nul est appliquee.

! IESICF : entree sortie imposee (par exemple entree supersonique)
!         l'utilisateur impose la vitesse et toutes les
!           variables thermodynamiques
! ISSPCF : sortie supersonique
!         l'utilisateur n'impose rien
! ISOPCF : sortie subsonique a pression imposee
!         l'utilisateur impose la pression
! IERUCF : entree subsonique a vitesse et rho imposes
!         l'utilisateur impose la vitesse et la masse volumique
! IEQHCF : entree subsonique a debit et debit enthalpique imposes
!         a implementer

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
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! icodcl           ! ia ! --> ! boundary condition code                        !
!  (nfabor, nvar)  !    !     ! = 1  -> Dirichlet                              !
!                  !    !     ! = 2  -> flux density                           !
!                  !    !     ! = 4  -> sliding wall and u.n=0 (velocity)      !
!                  !    !     ! = 5  -> friction and u.n=0 (velocity)          !
!                  !    !     ! = 6  -> roughness and u.n=0 (velocity)         !
!                  !    !     ! = 9  -> free inlet/outlet (velocity)           !
!                  !    !     !         inflowing possibly blocked             !
! itrifb           ! ia ! <-- ! indirection for boundary faces ordering        !
!  (nfabor, nphas) !    !     !                                                !
! itypfb           ! ia ! --> ! boundary face types                            !
!  (nfabor, nphas) !    !     !                                                !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
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
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! rcodcl           ! ra ! --> ! boundary condition values                      !
!  (nfabor,nvar,3) !    !     ! rcodcl(1) = Dirichlet value                    !
!                  !    !     ! rcodcl(2) = exterior exchange coefficient      !
!                  !    !     !  (infinite if no exchange)                     !
!                  !    !     ! rcodcl(3) = flux density value                 !
!                  !    !     !  (negative for gain) in w/m2 or                !
!                  !    !     !  roughness height (m) if icodcl=6              !
!                  !    !     ! for velocities           ( vistl+visct)*gradu  !
!                  !    !     ! for pressure                         dt*gradp  !
!                  !    !     ! for scalars    cp*(viscls+visct/sigmas)*gradt  !
! w1,2,3,4,5,6     ! ra ! --- ! work arrays                                    !
!  (ncelet)        !    !     !  (computation of pressure gradient)            !
! coefu            ! ra ! --- ! work array                                     !
!  (nfabor, 3)     !    !     !  (computation of pressure gradient)            !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
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
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          izfppp(nfabor)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          ifac  , iel   , ii    , iphas
integer          izone , iutile
integer          ilelt, nlelt

double precision uref2 , dhyd  , rhomoy
double precision ustar2, xkent , xeent , d2s3

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if(1.eq.1) then
  write(nfecra,9001)
  call csexit (1)
  !==========
endif

 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@     COMPRESSIBLE                                           ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uscfcl DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1.  INITIALISATIONS

!===============================================================================

idebia = idbia0
idebra = idbra0

d2s3 = 2.d0/3.d0

!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR LES FACES DE BORD
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LA CONDITION LIMITE

!          IMPOSER ICI LES CONDITIONS LIMITES SUR LES FACES DE BORD

!===============================================================================

iphas = 1


! --- Exemple d'entree/sortie pour laquelle tout est connu

!       sans presumer du caractere subsonique ou supersonique,
!         l'utilisateur souhaite imposer toutes les caracteristiques
!         de l'ecoulement
!       une entree supersonique est un cas particulier

!       La turbulence et les scalaires utilisateur prennent un flux nul
!         si la vitesse est sortante.

CALL GETFBR('1 and X <= 1.0 ',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!      ELEMENT ADJACENT A LA FACE DE BORD

  iel = ifabor(ifac)

!       On numerote les zones de 1 a n...
  izone = 1
  izfppp(ifac) = izone

  itypfb(ifac,iphas) = iesicf

!   - Vitesse
  rcodcl(ifac,iu(iphas),1) = 5.0d0
  rcodcl(ifac,iv(iphas),1) = 0.0d0
  rcodcl(ifac,iw(iphas),1) = 0.0d0

!   - Pression, Masse Volumique, Temperature, Energie Totale Specifique

!       Seules 2 variables sur les 4 sont independantes
!       On peut donc fixer le couple de variables que l'on veut
!       (sauf Temperature-Energie) et les 2 autres variables seront
!       calculees automatiquement

!    ** Choisir les 2 variables a imposer et effacer les autres
!       (elles sont calculées par les lois thermodynamiques dans uscfth)

!     Pression (en Pa)
  rcodcl(ifac,ipr(iphas),1) = 5.d5

!     Masse Volumique (en kg/m3)
!        RCODCL(IFAC,ISCA(IRHO  (IPHAS)),1) = 1.D0

!     Temperature (en K)
  rcodcl(ifac,isca(itempk(iphas)),1) = 300.d0

!     Energie Totale Specifique (en J/kg)
!        RCODCL(IFAC,ISCA(IENERG(IPHAS)),1) = 355.D3


!   - Turbulence

  uref2 = rcodcl(ifac,iu(iphas),1)**2                             &
         +rcodcl(ifac,iv(iphas),1)**2                             &
         +rcodcl(ifac,iw(iphas),1)**2
  uref2 = max(uref2,1.d-12)


!       Exemple de turbulence calculee a partir
!         de formules valables pour une conduite

!       On veillera a specifier le diametre hydraulique
!         adapte a l'entree courante.

!       On s'attachera egalement a utiliser si besoin une formule
!         plus precise pour la viscosite dynamique utilisee dans le
!         calcul du nombre de Reynolds (en particulier, lorsqu'elle
!         est variable, il peut etre utile de reprendre ici la loi
!         imposee dans USCFPV. On utilise ici par defaut la valeur
!         VISCL0 donnee dans USINI1
!       En ce qui concerne la masse volumique, on peut utiliser directement
!         sa valeur aux faces de bord si elle est connue (et imposée
!         ci-dessus). Dans le cas général, on propose d'utiliser la valeur
!         de la cellule adjacente.

!         Diametre hydraulique
  dhyd   = 0.075d0

!         Vitesse de frottement (au carre = USTAR2)
  rhomoy = rtp(iel,isca(irho(iphas)))
  ustar2 = 0.d0
  xkent  = epzero
  xeent  = epzero

  call keendb                                                     &
  !==========
    ( uref2, dhyd, rhomoy, viscl0(iphas), cmu, xkappa,            &
      ustar2, xkent, xeent )

  if    (itytur(iphas).eq.2) then

    rcodcl(ifac,ik(iphas),1)  = xkent
    rcodcl(ifac,iep(iphas),1) = xeent

  elseif(itytur(iphas).eq.3) then

    rcodcl(ifac,ir11(iphas),1) = 2.d0/3.d0*xkent
    rcodcl(ifac,ir22(iphas),1) = 2.d0/3.d0*xkent
    rcodcl(ifac,ir33(iphas),1) = 2.d0/3.d0*xkent
    rcodcl(ifac,ir12(iphas),1) = 0.d0
    rcodcl(ifac,ir13(iphas),1) = 0.d0
    rcodcl(ifac,ir23(iphas),1) = 0.d0
    rcodcl(ifac,iep(iphas),1)  = xeent

  elseif(iturb(iphas).eq.50) then

    rcodcl(ifac,ik(iphas),1)   = xkent
    rcodcl(ifac,iep(iphas),1)  = xeent
    rcodcl(ifac,iphi(iphas),1) = d2s3
    rcodcl(ifac,ifb(iphas),1)  = 0.d0

  elseif(iturb(iphas).eq.60) then

    rcodcl(ifac,ik(iphas),1)   = xkent
    rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

  endif

!   - On traite les scalaires rattaches a la phase courante
!       (ne pas boucler sur NSCAL sous peine de modifier rho et energie)
  if(nscaus.gt.0) then
    do ii = 1, nscaus
      if(iphsca(ii).eq.iphas) then
        rcodcl(ifac,isca(ii),1) = 1.d0
      endif
    enddo
  endif

enddo

! --- Exemple de sortie supersonique

!       toutes les caractéristiques sortent
!       on ne doit rien imposer (ce sont les valeurs internes qui sont
!         utilisees pour le calcul des flux de bord)

!       pour la turbulence et les scalaires, si on fournit ici des
!         valeurs de RCODCL, on les impose en Dirichlet si le flux
!         de masse est entrant ; sinon, on impose un flux nul (flux de
!         masse sortant ou RCODCL renseigné ici).
!         Noter que pour la turbulence, il faut renseigner RCODCL pour
!         toutes les variables turbulentes (sinon, on applique un flux
!         nul).


CALL GETFBR('2',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!       On numerote les zones de 1 a n...
  izone = 2
  izfppp(ifac) = izone

  itypfb(ifac,iphas) = isspcf

enddo

! --- Exemple d'entree subsonique (debit, debit enthalpique)

!       2 caracteristiques sur 3 entrent : il faut donner 2 informations
!         la troisieme est deduite par un scenario de 2-contact et
!         3-detente dans le domaine
!       ici on choisit de donner (rho*(U.n), rho*(U.n)*H)
!         avec H = 1/2 U*U + P/rho + e
!              n la normale unitaire entrante

!       ATTENTION, on donne des DENSITES de debit (par unite de surface)

CALL GETFBR('3',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!       On numerote les zones de 1 a n...
  izone = 3
  izfppp(ifac) = izone

  itypfb(ifac,iphas) = ieqhcf

!   - Densite de debit massique (en kg/(m2 s))
  rcodcl(ifac,irun(iphas),1) = 5.d5

!   - Densite de debit enthalpique (en J/(m2 s))
  rcodcl(ifac,irunh(iphas),1) = 5.d5


!     Condition non disponible dans la version presente
  call csexit (1)
  !==========

enddo

! --- Exemple d'entree subsonique (masse volumique, vitesse)

!       2 caracteristiques sur 3 entrent : il faut donner 2 informations
!         la troisieme est deduite par un scenario de 2-contact et
!         3-detente dans le domaine
!       ici on choisit de donner (rho, U)

CALL GETFBR('4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!       On numerote les zones de 1 a n...
  izone = 4
  izfppp(ifac) = izone

  itypfb(ifac,iphas) = ierucf

!   - Vitesse d'entree
  rcodcl(ifac,iu(iphas),1) = 5.0d0
  rcodcl(ifac,iv(iphas),1) = 0.0d0
  rcodcl(ifac,iw(iphas),1) = 0.0d0

!   - Masse Volumique (en kg/m3)
  rcodcl(ifac,isca(irho  (iphas)),1) = 1.d0


!   - Turbulence

  uref2 = rcodcl(ifac,iu(iphas),1)**2                             &
         +rcodcl(ifac,iv(iphas),1)**2                             &
         +rcodcl(ifac,iw(iphas),1)**2
  uref2 = max(uref2,1.d-12)


!       Exemple de turbulence calculee a partir
!         de formules valables pour une conduite

!       On veillera a specifier le diametre hydraulique
!         adapte a l'entree courante.

!       On s'attachera egalement a utiliser si besoin une formule
!         plus precise pour la viscosite dynamique utilisee dans le
!         calcul du nombre de Reynolds (en particulier, lorsqu'elle
!         est variable, il peut etre utile de reprendre ici la loi
!         imposee dans USCFPV. On utilise ici par defaut la valeur
!         VISCL0 donnee dans USINI1
!       En ce qui concerne la masse volumique, on peut utiliser directement
!         sa valeur aux faces de bord si elle est connue (et imposée
!         ci-dessus). Dans le cas général, on propose d'utiliser la valeur
!         de la cellule adjacente.

!         Diametre hydraulique
  dhyd   = 0.075d0

!         Calcul de la vitesse de frottement au carre (USTAR2)
!           et de k et epsilon en entree (XKENT et XEENT) a partir
!           de lois standards en conduite circulaire
!           (leur initialisation est inutile mais plus propre)
  rhomoy = propfb(ifac,ipprob(irom(iphas)))
  ustar2 = 0.d0
  xkent  = epzero
  xeent  = epzero

  call keendb                                                     &
  !==========
    ( uref2, dhyd, rhomoy, viscl0(iphas), cmu, xkappa,            &
      ustar2, xkent, xeent )

  if    (itytur(iphas).eq.2) then

    rcodcl(ifac,ik(iphas),1)  = xkent
    rcodcl(ifac,iep(iphas),1) = xeent

  elseif(itytur(iphas).eq.3) then

    rcodcl(ifac,ir11(iphas),1) = 2.d0/3.d0*xkent
    rcodcl(ifac,ir22(iphas),1) = 2.d0/3.d0*xkent
    rcodcl(ifac,ir33(iphas),1) = 2.d0/3.d0*xkent
    rcodcl(ifac,ir12(iphas),1) = 0.d0
    rcodcl(ifac,ir13(iphas),1) = 0.d0
    rcodcl(ifac,ir23(iphas),1) = 0.d0
    rcodcl(ifac,iep(iphas),1)  = xeent

  elseif(iturb(iphas).eq.50) then

    rcodcl(ifac,ik(iphas),1)   = xkent
    rcodcl(ifac,iep(iphas),1)  = xeent
    rcodcl(ifac,iphi(iphas),1) = d2s3
    rcodcl(ifac,ifb(iphas),1)  = 0.d0

  elseif(iturb(iphas).eq.60) then

    rcodcl(ifac,ik(iphas),1)   = xkent
    rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

  endif

!   - On traite les scalaires rattaches a la phase courante
!       (ne pas boucler sur NSCAL sous peine de modifier rho et energie)
  if(nscaus.gt.0) then
    do ii = 1, nscaus
      if(iphsca(ii).eq.iphas) then
        rcodcl(ifac,isca(ii),1) = 1.d0
      endif
    enddo
  endif

enddo

! --- Exemple de sortie subsonique

!       1 caracteristique sur 3 sort : il faut donner 1 information
!         les deux autres sont deduites par un scenario de 2-contact et
!         3-detente dans le domaine
!       ici on choisit de donner P

!       La turbulence et les scalaires utilisateur prennent un flux nul.

CALL GETFBR('5',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!       On numerote les zones de 1 a n...
  izone = 5
  izfppp(ifac) = izone

  itypfb(ifac,iphas) = isopcf

!     Pression (en Pa)
  rcodcl(ifac,ipr(iphas),1) = 5.d5

enddo

! --- Exemple de paroi

CALL GETFBR('7',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!       On numerote les zones de 1 a n...
  izone = 7
  izfppp(ifac) = izone

  itypfb(ifac,iphas) = iparoi


!     Paroi défilante
!       Par défaut, la paroi n'est pas défilante
!       Si la paroi est défilante, donner les composantes non nulles
!         de la vitesse. La vitesse sera projetée dans le plan
!         tangeant à la paroi. Dans l'exemple suivant, on impose
!         Ux = 1. (l'exemple est activé si IUTILE=1)

  iutile = 0
  if(iutile.eq.1) then
    rcodcl(ifac,iu(iphas),1) = 1.d0
  endif

!     Température imposée
!       Par défaut, la paroi est adiabatique
!       Si la paroi est à température imposée, l'indiquer par
!         ICODCL = 5 et donner la valeur en Kelvin dans RCODCL(.,.,1)
!         Dans l'exemple suivant, on impose T = 293.15 K (l'exemple
!         est activé si IUTILE=1)

  iutile = 0
  if(iutile.eq.1) then
    icodcl(ifac,isca(itempk(iphas)))   = 5
    rcodcl(ifac,isca(itempk(iphas)),1) = 20.d0 + 273.15d0
  endif

!     Flux imposé
!       Par défaut, la paroi est adiabatique
!       Si la paroi est à flux imposé, l'indiquer par
!         ICODCL = 3 et donner la valeur en Watt/m2 dans RCODCL(.,.,3)
!         Dans l'exemple suivant, on impose un flux de 1000 W/m2
!         - la plage en été - (l'exemple est activé si IUTILE=1)

  iutile = 0
  if(iutile.eq.1) then
    icodcl(ifac,isca(itempk(iphas)))   = 3
    rcodcl(ifac,isca(itempk(iphas)),3) = 1000.d0
  endif

enddo

! --- Exemple de symetrie

CALL GETFBR('8',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!       On numerote les zones de 1 a n...
  izone = 8
  izfppp(ifac) = izone

  itypfb(ifac,iphas) = isymet

enddo

!     Il est deconseille d'utiliser d'autres types de conditions
!     aux limites que ceux proposes ci-dessus.

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
