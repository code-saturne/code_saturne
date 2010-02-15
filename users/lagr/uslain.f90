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

subroutine uslain &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   ntersl , nvlsta , nvisbr ,                                     &
   nptnew ,                                                       &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   itypfb , itrifb , itepa  , ifrlag , injfac ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   ettp   , tepa   , vagaus , w1     , w2     , w3     ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    SOUS-PROGRAMME UTILISATEUR (INTERVENTION NON OBLIGATOIRE)

!    ROUTINE UTILISATEUR POUR LES CONDITIONS AUX LIMITES RELATIVES
!      AUX PARTICULES (ENTREE ET TRAITEMENT AUX AUTRES BORDS)

!    CE SOUS-PROGRAMME EST APPELE
!      APRES INITIALISATION DES TABLEAUX ETTP TEPA ET ITEPA
!      POUR LES NOUVELLES PARTICULES AFIN DE LES MODIFIER POUR
!      INJECTER DES PROFILS DE PARTICULES.


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
! lndfac           ! e  ! <-- ! longueur du tableau nodfac                     !
! lndfbr           ! e  ! <-- ! longueur du tableau nodfbr                     !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! ntersl           ! e  ! <-- ! nbr termes sources de couplage retour          !
! nvlsta           ! e  ! <-- ! nombre de var statistiques lagrangien          !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! nptnew           ! e  ! <-- ! nombre total de nouvelles particules           !
!                  !    !     ! pour toutes les zones d'injection              !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
!                  !    !     ! le module lagrangien                           !
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
! itrifb           ! ia ! <-- ! indirection for boundary faces ordering        !
!  (nfabor, nphas) !    !     !                                                !
! itypfb(nfabor    ! te ! <-- ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
! ifrlag(nfabor    ! te ! --> ! type des faces de bord lagrangien              !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! injfac(nbpmax    ! te ! <-- ! numero de la face de bord d'injection          !
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
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec)                     !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules etape courante                !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
! vagaus           ! tr ! --> ! variables aleatoires gaussiennes               !
!(nbpmax,nvgaus    !    !     !                                                !
! w1..w3(ncelet    ! tr ! --- ! tableaux de travail                            !
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
include "cstnum.h"
include "cstphy.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"
include "ppppar.h"
include "ppthch.h"
include "cpincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          ntersl , nvlsta , nvisbr
integer          nptnew
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1) , nodfac(lndfac)
integer          ipnfbr(nfabor+1) , nodfbr(lndfbr)
integer          itypfb(nfabor,nphas) , itrifb(nfabor,nphas)
integer          itepa(nbpmax,nivep) , ifrlag(nfabor)
integer          injfac(nbpmax)
integer          idevel(nideve) , ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac) , surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac) , cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod) , volume(ncelet)
double precision dt(ncelet) , rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*) , propfb(nfabor,*)
double precision coefa(nfabor,*) , coefb(nfabor,*)
double precision ettp(nbpmax,nvp) , tepa(nbpmax,nvep)
double precision vagaus(nbpmax,*)
double precision w1(ncelet) ,  w2(ncelet) ,  w3(ncelet)
double precision rdevel(nrdeve) , rtuser(nrtuse)
double precision ra(*)

! Local variables

integer          idebia , idebra
integer          iclas , izone , ifac
integer          ii , ip , npt , npar1 , npar2, ipnorm

! Local variables UTILISATEUR
!     (VGAUSS est dimensionne a 3, mais 2 suffirait ici)

double precision vgauss(3)

!===============================================================================


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

!     Par defaut on ne les modifie pas.

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

if (nbpnew.eq.0) return

!===============================================================================
! 1. GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 2. INITIALISATIONS
!===============================================================================



!===============================================================================
! 3. MODIFICATION DES DONNEES PARTICULAIRES DES NOUVELLES PARTICULES
!    (PROFILS D'INJECTION, REPLACEMENT DU POINT D'INJECTION,
!     MODIFICATION DES POIDS STATISTIQUES, CORRECTION DU DIAMETRE
!     SI OPTION ECART-TYPE ACTIVEE...)
!===============================================================================

!       LES MODIFICATIONS DES DONNEES PARTICULAIRES
!       INTERVIENNENT APRES TOUTES LES INITIALISATIONS LIEES
!       A L'INJECTION DES NOUVELLES PARTICULES, MAIS AVANT LE
!       TRAITEMENT DE L'INJECTION CONTINUE (IL EST DONC POSSIBLE
!       D'IMPOSER UN PROFIL D'INJECTION AVEC OPTION D'INJECTION
!       CONTINUE).

!   reinitialisation du compteur de nouvelles particules
npt = nbpart

!     pour chaque zone de bord:
do ii = 1,nfrlag
  izone = ilflag(ii)

!       pour chaque classe :
  do iclas = 1, iusncl(izone)

!         si de nouvelles particules doivent entrer :
    if (mod(ntcabs,iuslag(iclas,izone,ijfre)).eq.0) then

      do ip = npt+1 , npt+iuslag(iclas,izone,ijnbp)

!         NUMERO DE LA FACE DE BORD D'INJECTION D'ORIGINE

      ifac = injfac(ip)

!         EXEMPLE DE MODIFICATION DES VITESSES D'INJECTION
!           EN FONCTION DE LA POSITION D'INJECTION

!     Appeler par exemple votre propre sous-programme
!       qui fournirait les trois composantes de la vitesse instantanee
!       ETTP(IP,JUP),ETTP(IP,JVP),ETTP(IP,JWP)
!       en fonction de ETTP(IP,JZP) (interpolation par exemple)
!     Plus simplement on peut aussi imaginer de fournir les trois
!       composantes de la vitesse instantanee, sous la forme d'une
!       valeur moyenne (prise arbitrairement egale a (2,0,0) m/s ici)
!       et d'une valeur fluctuante calculee a partir d'une fluctuation
!       (prise arbitrairement egale a 0,2 m/s ici pour les composantes
!       1 et 3) :
        ipnorm = 2
        call normalen(ipnorm,vgauss)
        ettp(ip,jup) = 2.d0 + vgauss(1) * 0.2d0
        ettp(ip,jvp) = 0.d0
        ettp(ip,jwp) = 0.d0 + vgauss(2) * 0.2d0

      enddo

      npt = npt + iuslag(iclas,izone,ijnbp)

    endif

  enddo
enddo

!===============================================================================
! 4. SIMULATION DES VITESSES TURBULENTES FLUIDES INSTANTANNEES VUES
!    PAR LES PARTICULES SOLIDES LE LONG DE LEUR TRAJECTOIRE.
!===============================================================================

!    En entrant dans ce sous-programme, les tableaux :
!         ETTP(IP,JUF)
!         ETTP(IP,JVF)
!         ETTP(IP,JWF)
!      contiennent les composantes de la vitesse instantanee
!      (partie fluctuante + partie moyenne) du fluide vu
!      par les particules.

!    Lorsque la vitesse du fluide vu est modifiee ci-dessus,
!      le plus souvent l'utilisateur n'en connait que
!      la partie moyenne. Dans certaines configurations d'ecoulement
!      et d'injection des particules, il peut parfois s'averer
!      necessaire d'en reconstruire la partie turbulente.
!      C'est l'objet de l'appel au sous-programme suivant.

!    Attention : il ne faut reconstruire cette composante turbulente
!      que sur les vitesses du fluide vu modifiees.

!    La reconstruction est desactivee ici et doit etre adaptee au
!      cas traite.

if ( 1.eq.0 ) then

  npar1 = nbpart+1
  npar2 = nbpart+nbpnew

  call lagipn                                                     &
  !==========
  ( idebia , idebra ,                                             &
    ncelet , ncel   ,                                             &
    nbpmax , nvp    , nvp1   , nvep   , nivep  ,                  &
    npar1  , npar2  ,                                             &
    nideve , nrdeve , nituse , nrtuse ,                           &
    itepa  ,                                                      &
    idevel , ituser , ia     ,                                    &
    rtpa   ,                                                      &
    ettp   , tepa   , vagaus ,                                    &
    w1     , w2     , w3     ,                                    &
    rdevel , rtuser , ra     )

endif

!===============================================================================

!--------
! FORMATS
!--------

!----
! FIN
!----

return

end subroutine
