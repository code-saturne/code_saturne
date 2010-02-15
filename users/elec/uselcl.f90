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

subroutine uselcl &
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

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE
!                MODULE ELECTRIQUE
!   (Effet Joule, Arc Electrique, Conduction ionique)
!    REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
!    (ICODCL,RCODCL) POUR LES VARIABLES INCONNUES
!    PENDANT DE USCLIM.F



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

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "elincl.h"

!===============================================================================

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
integer          ifac, ii, iphas , iel
integer          idim
integer          izone,iesp
integer          ilelt, nlelt

double precision uref2, d2s3
double precision rhomoy, dhy, ustar2
double precision xkent, xeent
double precision z1   , z2

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
endif

 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@                      MODULE ELECTRIQUE                     ',/,&
'@                                                            ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uselcl DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@     Ce sous-programme utilisateur permet de definir les    ',/,&
'@       conditions aux limites. Il est indispensable.        ',/,&
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


! --- On impose en couleur 1 une entree ; exemple de Cathode
!     ======================================================

CALL GETFBR('1',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac,iphas) = ientre

!      - Numero de zone (on numerote de 1 a n)
  izone = 1

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

  rcodcl(ifac,iu(iphas),1) = 0.d0
  rcodcl(ifac,iv(iphas),1) = 0.d0
  rcodcl(ifac,iw(iphas),1) = 0.d0

!         Turbulence

!     (ITYTUR est un indicateur qui vaut ITURB/10)
  if (itytur(iphas).eq.2 .or. itytur(iphas).eq.3                  &
       .or. iturb(iphas).eq.50 .or. iturb(iphas).eq.60) then

    uref2 = rcodcl(ifac,iu(iphas),1)**2                           &
           +rcodcl(ifac,iv(iphas),1)**2                           &
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
!         imposee dans USPHYV. On utilise ici par defaut la valeur
!         VISCL0 donnee dans USINI1
!       En ce qui concerne la masse volumique, on dispose directement
!         de sa valeur aux faces de bord (ROMB) et c'est celle que
!         utilise donc ici (elle est en particulier coherente avec
!         le traitement implante dans USPHYV, en cas de masse
!         volumique variable)

!         Diametre hydraulique
    dhy     = 0.075d0

!         Calcul de la vitesse de frottement au carre (USTAR2)
!           et de k et epsilon en entree (XKENT et XEENT) a partir
!           de lois standards en conduite circulaire
!           (leur initialisation est inutile mais plus propre)
    rhomoy = propfb(ifac,ipprob(irom(iphas)))
    ustar2 = 0.d0
    xkent  = epzero
    xeent  = epzero

    call keendb                                                   &
    !==========
     ( uref2, dhy, rhomoy, viscl0(iphas), cmu, xkappa,            &
       ustar2, xkent, xeent )

    if (itytur(iphas).eq.2) then

      rcodcl(ifac,ik(iphas),1)  = xkent
      rcodcl(ifac,iep(iphas),1) = xeent

    elseif(itytur(iphas).eq.3) then

      rcodcl(ifac,ir11(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir22(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir33(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir12(iphas),1) = 0.d0
      rcodcl(ifac,ir13(iphas),1) = 0.d0
      rcodcl(ifac,ir23(iphas),1) = 0.d0
      rcodcl(ifac,iep(iphas),1)  = xeent

    elseif (iturb(iphas).eq.50) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iep(iphas),1)  = xeent
      rcodcl(ifac,iphi(iphas),1) = d2s3
      rcodcl(ifac,ifb(iphas),1)  = 0.d0

    elseif (iturb(iphas).eq.60) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

    endif

  endif

! --- On traite les scalaires

!      Enthalpie en J/kg

  ii = ihm
  icodcl(ifac,isca(ii))   = 1
  rcodcl(ifac,isca(ii),1) = 1.d6

!  Potentiel electrique reel impose a 0. volts (exemple de Cathode en arc)

  ii = ipotr
  icodcl(ifac,isca(ii))   = 1
  rcodcl(ifac,isca(ii),1) = 0.d0

!  Fraction massique des (N-1) constituants

  if ( ngazg .gt. 1 ) then
    do iesp=1,ngazg-1
      ii = iycoel(iesp)
      icodcl(ifac,isca(ii))   = 1
      rcodcl(ifac,isca(ii),1) = 0.d0
    enddo
  endif

!  Specifique Version Effet Joule :

!       Potentiel Imaginaire impose a 0

  if ( ippmod(ieljou).ge. 2 ) then
    ii = ipoti
    icodcl(ifac,isca(ii))   = 1
    rcodcl(ifac,isca(ii),1) = 0.d0
  endif

!  Specifique Version Arc Electrique :

!       Potentiel vecteur : Flux nul

  if ( ippmod(ielarc).ge.2 ) then
    do idim= 1,ndimve
      ii = ipotva(idim)
      icodcl(ifac,isca(ii))   = 3
      rcodcl(ifac,isca(ii),3) = 0.d0
    enddo
  endif

enddo

! --- On impose en couleur 5 une entree/sortie ;
!     ====================================== exemple d'Electrode en Joule
!                                            ============================

CALL GETFBR('5',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SORTIE : FLUX NUL VITESSE ET TEMPERATURE, PRESSION IMPOSEE
!            Noter que la pression sera recalee a P0
!                sur la premiere face de sortie libre (ISOLIB)

  itypfb(ifac,iphas)   = isolib

!      - Numero de zone (on numerote de 1 a n)
  izone = 2

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

! --- On traite les scalaires rattaches a la phase courante

!  Enthalpie en J/kg  (Par defaut flux nul avec ISOLIB)
!     Rien a faire

!  Fraction massique des (N-1) constituants (Par defaut flux nul avec ISOLIB)
!     Rien a faire

!  Specifique Version Effet Joule :

!     En effet Joule,
!       si l'on souhaite faire un calcul en recalant les conditions
!         aux limites (utiliser IELCOR=1 dans useli1)
!       pour atteindre la valeur de la puissance PUISIM
!         (a imposer dans useli1 en Ampere.Volt)
!       on multiplie la condition limite initiale sur le potentiel
!          reel (et sur le potentiel imaginaire s'il est pris en
!          compte) par le coefficient COEJOU.
!       COEJOU est determine automatiquement pour que la puissance
!          dissipee par effet Joule (partie reelle et partie
!          imaginaire si besoin) soit PUISIM
!       au debut du calcul, COEJOU vaut 1 ; COEJOU est transmis dans
!          les fichiers suites.

!     Si on ne souhaite pas faire un calcul avec recalage, on impose
!       directement une valeur adaptee.

  if ( ippmod(ieljou).ge. 1 ) then
    ii = ipotr
    icodcl(ifac,isca(ii))   = 1
    if(ielcor.eq.1) then
      rcodcl(ifac,isca(ii),1) = 500.d0*coejou
    else
      rcodcl(ifac,isca(ii),1) = 500.d0
    endif
  endif

  if ( ippmod(ieljou).ge. 2 ) then
    ii = ipoti
    icodcl(ifac,isca(ii))   = 1
    if(ielcor.eq.1) then
      rcodcl(ifac,isca(ii),1) = sqrt(3.d0)*500.d0*coejou
    else
      rcodcl(ifac,isca(ii),1) = sqrt(3.d0)*500.d0
    endif
  endif

enddo

! --- On impose en couleur 2 une entree/sortie ;
!     ============================== exemple d'Anode en arc electrique
!                                    =================================

CALL GETFBR('2',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SORTIE : FLUX NUL VITESSE ET TEMPERATURE, PRESSION IMPOSEE
!            Noter que la pression sera recalee a P0
!                sur la premiere face de sortie libre (ISOLIB)

  itypfb(ifac,iphas)   = isolib

!      - Numero de zone (on numerote de 1 a n)
  izone = 3

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

! --- On traite les scalaires rattaches a la phase courante

!  Enthalpie en J/kg  (Par defaut flux nul avec ISOLIB)
!     Rien a faire

!  Potentiel electrique reel

!     En arc electrique,
!       si l'on souhaite faire un calcul en recalant le potentiel
!          de l'anode (utiliser IELCOR=1 dans useli1)
!       pour atteindre la valeur du courant COUIMP
!         (a imposer dans useli1 en Amperes)
!       on utilise alors la valeur DPOT comme condition limite
!       DPOT est en effet automatiquement adaptee par le calcul
!          pour que (j.E Volume/DPOT) = COUIMP
!          (initialiser DPOT dans useli1 en Volts avec une valeur
!           representative de la difference de potentiel imposee)

!     Si on ne souhaite pas faire un calcul avec recalage,  on impose
!       directement une valeur adaptee au cas
!       (par exemple, ici 1000 Volts ).

  ii = ipotr
  icodcl(ifac,isca(ii))   = 1

  if ( ippmod(ielarc).ge.1 .and. ielcor .eq.1) then
    rcodcl(ifac,isca(ii),1) = dpot
  else
    rcodcl(ifac,isca(ii),1) = 1000.d0
  endif


!  Fraction massique des (N-1) constituants (Par defaut flux nul avec ISOLIB)

!  Specifique Version Arc Electrique :
!      Potentiel vecteur : flux nul (par defaut)


enddo

! --- On impose en couleur 3 une paroi
!     ================================

CALL GETFBR('3',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          PAROI : DEBIT NUL (FLUX NUL POUR LA PRESSION)
!                  FROTTEMENT POUR LES VITESSES (+GRANDEURS TURB)
!                  FLUX NUL SUR LES SCALAIRES

!          Pour un calcul arc electrique 3D, on cale le potentiel
!            vecteur avec une condition de Dirichlet issue des valeurs
!            du potentiel vecteur au pas de temps precedent
!            dans une zone de paroi choisie
!            Par defaut, ailleurs, un flux nul s'applique (paroi isolee).

  itypfb(ifac,iphas)   = iparoi

!      - Numero de zone (on numerote de 1 a n)
  izone = 4

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

  if ( ippmod(ielarc).ge.2 ) then
    if ( cdgfbo(1,ifac) .le.  2.249d-2  .or.                      &
         cdgfbo(1,ifac) .ge.  2.249d-2  .or.                      &
         cdgfbo(3,ifac) .le. -2.249d-2  .or.                      &
         cdgfbo(3,ifac) .ge.  2.249d-2       ) then
      iel = ifabor(ifac)
      do idim = 1, ndimve
        ii = ipotva(idim)
        icodcl(ifac,isca(ii))   = 1
        rcodcl(ifac,isca(ii),1) = rtpa(iel,isca(ii))
      enddo
    endif
  endif

enddo

! --- On impose en couleur 51 : anode avec claquage
!     =============================================

CALL GETFBR('51',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

  itypfb(ifac,iphas)   = iparoi

!      - Numero de zone (on numerote de 1 a n)
  izone = 5

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

! ---- Enthalpie (J/kg ) : coef echange impose

  ii=ihm
  icodcl(ifac,isca(ii))   = 1
  rcodcl(ifac,isca(ii),1) = 2.d4
  rcodcl(ifac,isca(ii),2) = 1.d5

!  Potentiel electrique reel

  ii = ipotr
  icodcl(ifac,isca(ii))   = 1

  if ( ippmod(ielarc).ge.1 .and. ielcor .eq.1) then
    rcodcl(ifac,isca(ii),1) = dpot
  else
    rcodcl(ifac,isca(ii),1) = 100.d0
  endif

!       Si CLAQUAGE : a adapter en fonction du cas et du
!                     sous-programme USELRC

  if ( ippmod(ielarc).ge.1 .and. ielcor .eq.1) then
    if(iclaq.eq.1 .and. ntcabs.le.ntdcla+30) then

      z1 = zclaq - 2.d-4
      if(z1.le.0.d0) z1 = 0.d0
      z2 = zclaq + 2.d-4
      if(z2.ge.2.d-2) z2 = 2.d-2

      if( cdgfbo(3,ifac).ge.z1 .and.                              &
           cdgfbo(3,ifac).le.z2       ) then
        icodcl(ifac,isca(ii))   = 1
        rcodcl(ifac,isca(ii),1) = dpot
      else
        icodcl(ifac,isca(ii))   = 3
        rcodcl(ifac,isca(ii),3) = 0.d0
      endif
    endif
  endif

!       Potentiel vecteur : Flux nul

  if ( ippmod(ielarc).ge.2 ) then
    do idim= 1,ndimve
      ii = ipotva(idim)
      icodcl(ifac,isca(ii))   = 3
      rcodcl(ifac,isca(ii),3) = 0.d0
    enddo
  endif

enddo

! --- On impose en couleur 4 une symetrie
!     ===================================

CALL GETFBR('4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SYMETRIES

  itypfb(ifac,iphas)   = isymet

!      - Numero de zone (on numerote de 1 a n)
  izone = 6

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

!     Par defaut tous les scalaires (potentiels en particulier)
!       recoivent une condition de flux nul.
!     En effet Joule, on peut souhaiter imposer une condition
!       d'antisymetrie sur le potentiel imaginaire selon la
!       configuration des electrodes :
  if ( ippmod(ieljou).ge. 2 ) then
    ii = ipoti
    icodcl(ifac,isca(ii))   = 1
    rcodcl(ifac,isca(ii),1) = 0.d0
  endif

enddo

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
