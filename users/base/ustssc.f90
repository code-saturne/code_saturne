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

subroutine ustssc &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iscal  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , lstelt , &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   crvexp , crvimp ,                                              &
   viscf  , viscb  , xam    ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     , w10    , w11    ,          &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE UTILISATEUR
!   ON PRECISE LES TERMES SOURCES UTILISATEURS
!   POUR UN SCALAIRE SUR UN PAS DE TEMPS

! ON RESOUT RHO*VOLUME*D(VAR)/DT = CRVIMP*VAR + CRVEXP

! ON FOURNIT ICI CRVIMP ET CRVEXP (ILS CONTIENNENT RHO*VOLUME)
!    CRVEXP en kg variable/s :
!     ex : pour les temperatures      kg degres/s
!          pour les enthalpies        Joules/s
!    CRVIMP en kg /s :

! VEILLER A UTILISER UN CRVIMP NEGATIF
! (ON IMPLICITERA CRVIMP
!  IE SUR LA DIAGONALE DE LA MATRICE, LE CODE AJOUTERA :
!   MAX(-CRVIMP,0) EN SCHEMA STANDARD EN TEMPS
!       -CRVIMP    SI LES TERMES SOURCES SONT A L'ORDRE 2

! CES TABLEAUX SONT INITIALISES A ZERO AVANT APPEL A CE SOUS
!   PROGRAMME ET AJOUTES ENSUITE AUX TABLEAUX PRIS EN COMPTE
!   POUR LA RESOLUTION

! EN CAS D'ORDRE 2 DEMANDE SUR LES TERMES SOURCES, ON DOIT
!   FOURNIR CRVEXP A L'INSTANT N     (IL SERA EXTRAPOLE) ET
!           CRVIMP A L'INSTANT N+1/2 (IL EST  DANS LA MATRICE,
!                                     ON LE SUPPOSE NEGATIF)

! L'identification des faces de bord concernees se fait grace
! a la commande GETFBR.


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.


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
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iscal            ! i  ! <-- ! scalar number                                  !
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
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
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
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! crvexp(ncelet    ! tr ! --> ! tableau de travail pour part explicit          !
! crvimp(ncelet    ! tr ! --> ! tableau de travail pour terme instat           !
! viscf(nfac)      ! tr ! --- ! tableau de travail    faces internes           !
! viscb(nfabor     ! tr ! --- ! tableau de travail    faces de bord            !
! xam(nfac,2)      ! tr ! --- ! tableau de travail    faces de bord            !
! w1..11(ncelet    ! tr ! --- ! tableau de travail    cellules                 !
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
include "pointe.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "parall.h"
include "period.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iscal

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          maxelt, lstelt(maxelt)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision crvexp(ncelet), crvimp(ncelet)
double precision viscf(nfac), viscb(nfabor)
double precision xam(nfac,2)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet), w11(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

character*80     chaine
integer          idebia, idebra
integer          ivar, iiscvr, ipcrom, iel, iphas, iutile
integer          ilelt, nlelt

double precision tauf, prodf, volf, pwatt

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! --- Numero du scalaire a traiter : ISCAL

!     ISCAL est une donnee d'entree de ce sous programme
!       qui indique quel scalaire on est en train de traiter
!       en effet ce sous programme est appele successivement pour
!       tous les scalaires du calcul.
!     ISCAL ne doit pas etre modifie par l'utilisateur

!     Si le calcul comporte plusieurs scalaires (de 1 a n par exemple)
!       et que l'on souhaite imposer des termes sources pour un
!       scalaire p donne, on doit utiliser ISCAL ici :
!       CRVIMP et CRVEXP ne seront renseignes que "IF(ISCAL.EQ.p)"


! --- Numero de la variable associee au scalaire a traiter ISCAL
ivar = isca(iscal)

! --- Nom de la variable associee au scalaire a traiter ISCAL
chaine = nomvar(ipprtp(ivar))

! --- Indicateur de variance
!         Si ISCAVR = 0 :
!           le scalaire ISCAL n'est pas une variance
!         Si ISCAVR > 0 et ISCAVR < NSCAL + 1 :
!           le scalaire ISCAL est une variance associee
!           au scalaire ISCAVR
iiscvr = iscavr(iscal)

! --- Numero de phase associee au scalaire ISCAL
iphas = iphsca(iscal)

! --- Numero des grandeurs physiques (voir usclim)
ipcrom = ipproc(irom(iphas))

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif


!===============================================================================
! 2. EXEMPLE DE TERME SOURCE ARBITRAIRE, POUR LA VARIABLE F :

!                             S = A * F + B

!            APPARAISSANT DANS LES EQUATIONS SOUS LA FORME :

!                       RHO VOLUME D(F)/Dt = VOLUME*S


!   CE TERME A UNE PARTIE QU'ON VEUT IMPLICITER         : A
!           ET UNE PARTIE QU'ON VA TRAITER EN EXPLICITE : B


!   ICI PAR EXEMPLE ARBITRAIREMENT :

!     A = - RHO / TAUF
!     B =   RHO * PRODF
!        AVEC
!     TAUF   = 10.D0  [secondes  ] (TEMPS DE DISSIPATION DE F)
!     PRODF  = 100.D0 [variable/s] (PRODUCTION DE F PAR UNITE DE TEMPS)

!   ON A ALORS
!     CRVIMP(IEL) = VOLUME(IEL)* A = - VOLUME(IEL) (RHO / TAUF )
!     CRVEXP(IEL) = VOLUME(IEL)* B =   VOLUME(IEL) (RHO * PRODF)

!   ON REMPLIT CI-DESSOUS CRVIMP ET CRVEXP CORRESPONDANTS.

!===============================================================================



!     ATTENTION, L'EXEMPLE EST COMPLETEMENT ARBITRAIRE
!     =========
!         ET DOIT ETRE REMPLACE PAR LES TERMES UTILISATEURS ADEQUATS


! ----------------------------------------------

! Il est assez courant que l'on oublie d'eliminer cet exemple
!   de ce sous-programme.
! On a donc prevu le test suivant pour eviter les mauvaises surprises

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------

tauf  = 10.d0
prodf = 100.d0

do iel = 1, ncel
  crvimp(iel) = - volume(iel)*propce(iel,ipcrom)/tauf
enddo

do iel = 1, ncel
  crvexp(iel) =   volume(iel)*propce(iel,ipcrom)*prodf
enddo

!===============================================================================
! 3. EXEMPLE DE TERME SOURCE ARBITRAIRE, POUR LA VARIABLE F = ENTHALPIE :

!     IL S'AGIT D'UNE PARTICULARISATION DE L'EXEMPLE PRECEDENT AVEC

!                             S = B

!            APPARAISSANT DANS LES EQUATIONS SOUS LA FORME :

!                       RHO VOLUME D(F)/Dt = VOLUME*S


!   IL N'Y A RIEN A IMPLICITER, ON PEUT DONC IMPOSER

!     CRVIMP(IEL) = 0.D0


!   CE TERME A UNE PARTIE QU'ON VA TRAITER EN EXPLICITE : B

!     SI F EST UNE ENTHALPIE (J/kg), ON AURA ALORS B EN Watt/m3

!     SI ON CONNAIT LA PUISSANCE PWATT A INTRODUIRE UNIFORMEMENT DANS UN
!       VOLUME DONNE VOLF TEL QUE 0<X<1.2, 3.1<Y<4, ON POURRA ALORS
!       ECRIRE, POUR LES CELLULES DU MAILLAGE CONCERNEES :

!       CRVEXP(IEL) = VOLUME(IEL)* B =   VOLUME(IEL)*(PWATT/VOLF)


!===============================================================================



!     ATTENTION, L'EXEMPLE EST COMPLETEMENT ARBITRAIRE
!     =========
!         ET DOIT ETRE REMPLACE PAR LES TERMES UTILISATEURS ADEQUATS


! ----------------------------------------------

! Il est assez courant que l'on oublie d'eliminer cet exemple
!   de ce sous-programme.
! On a donc prevu le test suivant pour eviter les mauvaises surprises

iutile = 0

if(iutile.eq.0) return

! ----------------------------------------------

! PUISSANCE A INTRODUIRE DANS VOLF
!   ATTENTION on suppose qu'on travaille en enthalpie,
!             si on travaille en temperature, il ne faut pas oublier de
!               diviser PWATT par Cp.

pwatt = 100.d0

! CALCUL DE VOLF

volf  = 0.d0
CALL GETCEL('X > 0.0 and X < 1.2 and Y > 3.1 and'//               &
            'Y < 4.0',NLELT,LSTELT)

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
  volf = volf + volume(iel)
enddo

do ilelt = 1, nlelt
  iel = lstelt(ilelt)
! PAS DE TERME IMPLICITE
  crvimp(iel) = 0.d0
! PUISSANCE EN EXPLICITE
  crvexp(iel) = volume(iel)*pwatt/volf
enddo

!--------
! FORMATS
!--------

 1000 format(' TERMES SOURCES UTILISATEURS POUR LA VARIABLE ',A8,/)

!----
! FIN
!----

return

end subroutine
