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

subroutine mtkpdc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ncepdp , iphas  , iappel ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc ,                   &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc ,                                     &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DES PERTES DE CHARGE POUR MATISSE (COPIE DE USKPDC)

!   TRAITEMENT DES REGISTRES UNIQUEMENT (ENTREE ET SORTIE)

!   LES AUTRES PERTES DE CHARGES (ISOTROPES)
!                                 SONT TRAITEES DANS MTTSNS



! IAPPEL = 1 :
!             CALCUL DU NOMBRE DE CELLULES OU L'ON IMPOSE UNE PDC
! IAPPEL = 2 :
!             REPERAGE DES CELLULES OU L'ON IMPOSE UNE PDC
! IAPPEL = 3 :
!             CALCUL DES VALEURS DES COEFS DE PDC


! CKUPDC EST LE COEFF DE PDC CALCULE.

!  IL INTERVIENT DANS LA QDM COMME SUIT :
!    RHO DU/DT = - GRAD P + TSPDC        (+ AUTRES TERMES)
!                      AVEC TSPDC = - RHO CKUPDC U ( en kg/(m2 s))


!  POUR UNE PDC REPARTIE,

!    SOIT KSIL = DHL/(0.5 RHO U**2) DONNE DANS LA LITTERATURE
!    (DHL EST LA PERTE DE CHARGE PAR UNITE DE LONGUEUR)

!    LE TERME SOURCE TSPDC VAUT DHL = - KSIL *(0.5 RHO U**2)

!    ON A CKUPDC = 0.5 KSIL ABS(U)


!  POUR UNE PDC SINGULIERE,

!    SOIT KSIS = DHS/(0.5 RHO U**2) DONNE DANS LA LITTERATURE
!    (DHS EST LA PERTE DE CHARGE SINGULIERE)

!    LE TERME SOURCE TSPDC VAUT DHS/L = - KSIS/L *(0.5 RHO U**2)

!    ON A CKUPDC = 0.5 KSIS/L ABS(U)

!    OU L DESIGNE LA LONGUEUR SUR LAQUELLE
!           ON A CHOISI DE REPRESENTER LA ZONE DE PDC SINGULIERE


!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
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
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iphas            ! i  ! <-- ! phase number                                   !
! iappel           ! e  ! <-- ! indique les donnes a renvoyer                  !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
! icepdc(ncepdp    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
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
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
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

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.h"
include "pointe.h"
include "numvar.h"
include "optcal.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"
include "period.h"
include "matiss.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas  , iappel

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)


! Local variables

integer          idebia, idebra
integer          iel   , ielpdc, ikpdc
integer          ifml  , icoul
double precision alpha , cosalp, sinalp
double precision vit2  , vit3  , ck2   , ck3

!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0


if(iappel.eq.1.or.iappel.eq.2) then

!===============================================================================

! 1. POUR CHAQUE PHASE : UN OU DEUX APPELS

!      PREMIER APPEL :

!        IAPPEL = 1 : NCEPDP : CALCUL DU NOMBRE DE CELLULES
!                                AVEC PERTES DE CHARGE


!      DEUXIEME APPEL (POUR LES PHASES AVEC NCEPDP > 0) :

!        IAPPEL = 2 : ICEPDC : REPERAGE DU NUMERO DES CELLULES
!                                AVEC PERTES DE CHARGE

! REMARQUES :

!        Ne pas utiliser CKUPDC dans cette section
!          (il est rempli au troisieme appel, IAPPEL = 3)

!        Ne pas utiliser ICEPDC dans cette section
!           au premier appel (IAPPEL = 1)

!        On passe ici a chaque pas de temps
!           (ATTENTION au cout calcul de vos developpements)

!===============================================================================


!  1.1 A completer par l'utilisateur : selection des cellules
!     Pour Matisse, c'est fait par defaut selon les couleurs
!  -----------------------------------------------------------

! --- Aucune pdc (initialisation)

  ielpdc = 0

! --- Pdc definies selon les couleurs du maillage
!     registres d'entree ICMTRI et de sortie ICMTRO
  do iel = 1, ncel
    ifml  = ifmcel(iel   )
    icoul = iprfml(ifml,1)
    if(icoul.eq.icmtri) then
      ielpdc = ielpdc + 1
      if (iappel.eq.2) icepdc(ielpdc) = iel
    elseif(icoul.eq.icmtro) then
      ielpdc = ielpdc + 1
      if (iappel.eq.2) icepdc(ielpdc) = iel
    endif
  enddo


!  1.2 Sous section generique a ne pas modifier
!  ---------------------------------------------

! --- Pour IAPPEL = 1,
!      Renseigner NCEPDP, nombre de cellules avec pdc
!      Le bloc ci dessous est valable pourles 2 exemples ci dessus

  if (iappel.eq.1) then
    ncepdp = ielpdc
  endif

!-------------------------------------------------------------------------------

elseif(iappel.eq.3) then

!===============================================================================

! 2. POUR CHAQUE PHASE AVEC NCEPDP > 0 , TROISIEME APPEL

!      TROISIEME APPEL (POUR LES PHASES AVEC NCEPDP > 0) :

!       IAPPEL = 3 : CKUPDC : CALCUL DES COEFFICIENTS DE PERTE DE CHARGE
!                             DANS LE REPERE DE CALCUL
!                             STOCKES DANS L'ORDRE
!                             K11, K22, K33, K12, K13, K23


!    REMARQUE :

!        Veillez a ce que les coefs diagonaux soient positifs.

!        Vous risquez un PLANTAGE si ce n'est pas le cas.

!        AUCUN controle ulterieur ne sera effectue.

!      ===========================================================


!  2.1 A completer par l'utilisateur : valeur des coefs
!     Pour Matisse, c'est fait par defaut selon les couleurs
! -----------------------------------------------------

! --- Attention
!   Il est important que les CKUPDC soient completes (par des valeurs
!     nulles eventuellement) dans la mesure ou ils seront utilises pour
!     calculer un terme source dans les cellules identifiees precedemment.
!   On les initialise tous par des valeurs nulles.

  do ikpdc = 1, 6
    do ielpdc = 1, ncepdp
      ckupdc(ielpdc,ikpdc) = 0.d0
    enddo
  enddo

! --- Tenseur diagonal : pas pour Matisse
!     On elimine la section

! --- Tenseur 3x3


  do ielpdc = 1, ncepdp

!     Identification de la zone traitee
    iel   = icepdc(ielpdc)
    ifml  = ifmcel(iel   )
    icoul = iprfml(ifml,1)

!     Donnees relatives aux registres "amont"
    if(icoul.eq.icmtri)then
      alpha = argamt*pi/180.d0
      ck2 = 0.5d0*pdcalg/epregi
      ck3 = 0.5d0*pdcatv/epregi
!     Donnees relatives aux registres "aval"
    elseif(icoul.eq.icmtro)then
      alpha = argavl*pi/180.d0
      ck2 = 0.5d0*pdcslg/epregi
      ck3 = 0.5d0*pdcstv/epregi
!     Par securite (la selection precedente des elements
!       a ete faite sur les couleurs ICMTRI et ICMTRO)
    else
      alpha = 0.d0
      ck2 = 0.d0
      ck3 = 0.d0
    endif

!     Calcul des pertes de charge
    cosalp = cos(alpha)
    sinalp = sin(alpha)
    vit2 = cosalp*rtpa(iel,iv(iphas))                             &
         - sinalp*rtpa(iel,iw(iphas))
    vit3 = sinalp*rtpa(iel,iv(iphas))                             &
         + cosalp*rtpa(iel,iw(iphas))

    ckupdc(ielpdc,1) = 0.d0
    ckupdc(ielpdc,2) =                                            &
         cosalp**2*ck2*abs(vit2)+sinalp**2*ck3*abs(vit3)
    ckupdc(ielpdc,3) =                                            &
         sinalp**2*ck2*abs(vit2)+cosalp**2*ck3*abs(vit3)
    ckupdc(ielpdc,4) = 0.d0
    ckupdc(ielpdc,5) = 0.d0
    ckupdc(ielpdc,6) =                                            &
         cosalp*sinalp*(-ck2*abs(vit2)+ck3*abs(vit3))

  enddo

!-------------------------------------------------------------------------------

endif

return

end subroutine
