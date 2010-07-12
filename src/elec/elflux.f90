!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2010 EDF S.A., France

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

subroutine elflux &
!================

 ( idbia0 , idbra0 , iappel ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , viscf  , viscb  ,                            &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE POUR LE MODULE ELECTRIQUE

!              CALCULS DES VARIABLES PROPCE
!        ELLES SONT UTILISEES POUR LE CALCUL DES TS

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! iappel           ! e  ! <-- ! numero d'appel                                 !
!                  !    !     ! 1 : j, e, j.e                                  !
!                  !    !     ! 2 : b,    jxb                                  !
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
! iprfml           ! te ! <-- ! proprietes d'une famille                       !
! ipnfac(nfac+1)   ! ia ! <-- ! interior faces -> vertices index (optional)    !
! nodfac(lndfac)   ! ia ! <-- ! interior faces -> vertices list (optional)     !
! ipnfbr(nfabor+1) ! ia ! <-- ! boundary faces -> vertices index (optional)    !
! nodfbr(lndfbr)   ! ia ! <-- ! boundary faces -> vertices list (optional)     !
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
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! viscf(nfac)      ! tr ! --- ! tableau de travail    faces internes           !
! viscb(nfabor     ! tr ! --- ! tableau de travail    faces de bord            !
! w1..9(ncelet     ! tr ! --- ! tableau de travail    cellules                 !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!-------------------------------------------------------------------------------
!===============================================================================

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "dimfbr.h"
include "paramx.h"
include "numvar.h"
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "cstnum.h"
include "parall.h"
include "period.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "elincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0 , iappel
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve)
integer          ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision viscf(nfac), viscb(nfabor)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra
integer          iel
integer          ipcefj, ipcsig, ipcsii
integer          ipcla1, ipcla2, ipcla3
integer          ipcdc1, ipcdc2, ipcdc3
integer          ipcdi1, ipcdi2, ipcdi3
integer          inc   , iccocg, nswrgp, imligp, iwarnp
integer          ivar0 , iclimv
integer          iphydp, ivar  , modntl

double precision epsrgp, climgp, extrap, vrmin, vrmax, var

!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

! --- Numero des grandeurs physiques
ipcsig = ipproc(ivisls(ipotr))

if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
  ipcsii = ipproc(ivisls(ipoti))
endif

ipcefj = ipproc(iefjou)
ipcdc1 = ipproc(idjr(1))
ipcdc2 = ipproc(idjr(2))
ipcdc3 = ipproc(idjr(3))

if ( ippmod(ieljou).eq.4 ) then
  ipcdi1 = ipproc(idji(1))
  ipcdi2 = ipproc(idji(2))
  ipcdi3 = ipproc(idji(3))
endif

! --- Necessite d'une impression (impression si MODNTL=0)
if(ntlist.gt.0) then
  modntl = mod(ntcabs,ntlist)
elseif(ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif

!===============================================================================
! 1. PREMIER APPEL : J, E => J.E
!===============================================================================

if(iappel.eq.1) then

!===============================================================================
! 1.1 PRISE EN COMPTE DES TERMES SOURCES ET VARIABLES STANDARD ET
!   COMMUNES A TOUTES LES VERSIONS ELECTRIQUES : EFFET JOULE, E, J
!         CAS IELJOU >= 1 ou IELARC >= 1 ou IELION >= 1
!===============================================================================

!    Pour toutes les versions electriques :
!      on doit calculer le terme joule (partie reelle) = j . E
!      produit de la densite de courant par le champ electrique


!   2.1 Calcul du grad (potR) (W4, W5, W6)
!  ---------------------------

  ivar = isca(ipotr)
  iclimv = iclrtp(ivar,icoef)

  inc = 1
  iccocg = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  ! En periodique et parallele, echange avant calcul du gradient
  !     C'est indispensable car on vient de calculer IVAR
  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(rtp(1,ivar))
    !==========
  endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
  ivar0 = 0

!    Sans prise en compte de la pression hydrostatique

  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra     , ra     , ra     ,                                     &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
!       POTR
   w4     , w5     , w6     ,                                     &
!       d POTR /dx   d POTR /dy   d POTR /dz
   w7     , w8     , w9    ,                                      &
   rdevel , rtuser , ra     )


!   2.2 Calcul du champ electrique E = - grad (potR) : (-W4, -W5, -W6)
!  -------------------------------------------------

!   2.3 Calcul de la densite de courant j = sig E  :
!  -------------------------------------------------
!                                   PROPCE(IEL,IPCSIG) (-W4, -W5, -W6)


  if( ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1 ) then
    do iel = 1, ncel
      propce(iel,ipcdc1)= - propce(iel,ipcsig) * w4(iel)
      propce(iel,ipcdc2)= - propce(iel,ipcsig) * w5(iel)
      propce(iel,ipcdc3)= - propce(iel,ipcsig) * w6(iel)
    enddo
  endif

!   2.4 Calcul de l'Effet Joule j . E : PROPCE( .,IPPROC(IEFJOU) )
!  -----------------------------------------------------------------
!                                                           sig E.E
  do iel = 1, ncel
    propce(iel,ipcefj)=                                           &
         propce(iel,ipcsig)*(w4(iel)**2+w5(iel)**2+w6(iel)**2)
  enddo


!   2.5 On imprime les extrema de E et j
!  -------------------------------------

  if(modntl.eq.0) then

    write(nfecra,1000)

!     Grad PotR = -E
    var    = w4(1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var    = w4(iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    WRITE(NFECRA,1010)'Gr_PotRX',VRMIN,VRMAX

    var    = w5(1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var    = w5(iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    WRITE(NFECRA,1010)'Gr_PotRY',VRMIN,VRMAX

    var    = w6(1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var    = w6(iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    WRITE(NFECRA,1010)'Gr_PotRZ',VRMIN,VRMAX

    var    = -propce(1,ipcsig) * w4(1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var = -propce(iel,ipcsig) * w4(iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    WRITE(NFECRA,1010)'Cour_ReX',VRMIN,VRMAX

    var    = -propce(1,ipcsig) * w5(1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var = -propce(iel,ipcsig) * w5(iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    WRITE(NFECRA,1010)'Cour_ReY',VRMIN,VRMAX

    var    = -propce(1,ipcsig) * w6(1)
    vrmin = var
    vrmax = var
    do iel = 1, ncel
      var = -propce(iel,ipcsig) * w6(iel)
      vrmin = min(vrmin,var)
      vrmax = max(vrmax,var)
    enddo
    if (irangp.ge.0) then
      call parmin (vrmin)
      call parmax (vrmax)
    endif
    WRITE(NFECRA,1010)'Cour_ReZ',VRMIN,VRMAX

  endif



!===============================================================================
! 1.2 PRISE EN COMPTE ET AJOUT DES TERMES SOURCES ET VARIABLES
!    RELATIVES A LA PRISE EN COMPTE DU POTENTIEL COMPLEXE
!                           CAS IELJOU = 2
!===============================================================================

  if(ippmod(ieljou).ge.2 .or. ippmod(ieljou).eq.4) then


!   3.1 Calcul du grad (potI) :  (W4, W5, W6)
!  ----------------------------

    ivar = isca(ipoti)
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    ! En periodique et parallele, echange avant calcul du gradient
    !     C'est indispensable car on vient de calculer IVAR
    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ivar))
      !==========
    endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
    ivar0 = 0


!    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE

    iphydp = 0

    call grdcel                                                   &
    !==========
  ( idbia0 , idbra0 ,                                             &
    ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,&
    nnod   , lndfac , lndfbr , ncelbr , nphas  ,                  &
    nideve , nrdeve , nituse , nrtuse ,                           &
    ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp ,&
    iwarnp , nfecra , epsrgp , climgp , extrap ,                  &
    ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                  &
    ipnfac , nodfac , ipnfbr , nodfbr ,                           &
    idevel , ituser , ia     ,                                    &
    xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume ,&
    ra     , ra     , ra     ,                                    &
    rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv) ,              &
    w4     , w5     , w6     ,                                    &
!       d POTI /dx   d POTI /dy   d POTI /dz
    w7     , w8     , w9     ,                                    &
    rdevel , rtuser , ra     )


!   3.2 Calcul du champ electrique Ei = - grad (potI) : (-W4, -W5, -W6)
!  -------------------------------------------------

!   3.3 Partie imaginaire de la densite de courant :
!  ------------------------------------------------
!                                   PROPCE(IEL,IPCSIG) (-W4, -W5, -W6)

  if ( ippmod(ieljou).eq.4 ) then
    do iel = 1, ncel
      propce(iel,ipcdi1)= -propce(iel,ipcsig)*w4(iel)
      propce(iel,ipcdi2)= -propce(iel,ipcsig)*w5(iel)
      propce(iel,ipcdi3)= -propce(iel,ipcsig)*w6(iel)
    enddo
  endif


!   3.4 Effet Joule total : PROPCE( .,IPPROC(IEFJOU) )
!  ----------------------------------------------------

    do iel = 1, ncel

!             ajout de la partie imaginaire et ...
      propce(iel,ipcefj) = propce(iel,ipcefj)                     &
         + propce(iel,ipcsii)*(w4(iel)**2+w5(iel)**2+w6(iel)**2)
!             .    ..division par 2
      propce(iel,ipcefj) = 0.5d0*propce(iel,ipcefj)

    enddo



!   3.5 On imprime les extrema de E et j
!  -------------------------------------

    if(modntl.eq.0) then

!     Grad PotI = -Ei
      var    = w4(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = w4(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      WRITE(NFECRA,1010)'Gr_PotIX',VRMIN,VRMAX

      var    = w5(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = w5(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      WRITE(NFECRA,1010)'Gr_PotIY',VRMIN,VRMAX

      var    = w6(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = w6(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      WRITE(NFECRA,1010)'Gr_PotIZ',VRMIN,VRMAX

!     j=sigma E
      var    = -propce(1,ipcsii) * w4(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var = -propce(iel,ipcsii) * w4(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      WRITE(NFECRA,1010)'Cour_ImX',VRMIN,VRMAX

      var    = -propce(1,ipcsii) * w5(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var = -propce(iel,ipcsii) * w5(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      WRITE(NFECRA,1010)'Cour_ImY',VRMIN,VRMAX

      var    = -propce(1,ipcsii) * w6(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var = -propce(iel,ipcsii) * w6(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      WRITE(NFECRA,1010)'Cour_ImZ',VRMIN,VRMAX

      write(nfecra,1001)

    endif

  endif

endif
!     Fin du test IAPPEL = 1

!===============================================================================
! 2. DEUXIEME APPEL : A, B, JXB
!===============================================================================

if (iappel.eq.2) then

!===============================================================================
! 2.1 PRISE EN COMPTE ET AJOUT DES TERMES SOURCES ET VARIABLES RELATIVES
!    A LA PRISE EN COMPTE DE L'ARC ELECTRIQUE
!                           CAS IELARC = 1 OU 2
!===============================================================================


  if(ippmod(ielarc).ge.1) then
    ipcla1 = ipproc(ilapla(1))
    ipcla2 = ipproc(ilapla(2))
    ipcla3 = ipproc(ilapla(3))
  endif

!   4.1 ARC ELECTRIQUE DIT 3D : IELARC = 2
!  ----------------------------------------
!         ON PASSE PAR LA RESOLUTION D'EQUATIONS DE POISSON SUR LES
!        -----------------------------------------------------------
!           COMPOSANTES DU POTENTIEL VECTEUR
!          ----------------------------------

  if( ippmod(ielarc).ge.2 ) then

! --> Calcul des composantes du champ magnetique B = (W1, W2, W3)
!===================================================


!    Sur Ax

    ivar = isca(ipotva(1))
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    ! En periodique et parallele, echange avant calcul du gradient
    !     C'est indispensable car on vient de calculer IVAR
    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ivar))
      !==========
    endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)

    ivar0 = 0

!    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE

    iphydp = 0

    call grdcel                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra     , ra     , ra     ,                                     &
   rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv)  ,              &
   w4     , w5     , w6     ,                                     &
!       d Ax /dx   d Ax /dy   d Ax /dz
   w7     , w8     , w9    ,                                      &
   rdevel , rtuser , ra     )

!       B = rot A

    do iel = 1, ncel
      w1(iel)=  zero
      w2(iel)=  w6(iel)
      w3(iel)= -w5(iel)
    enddo

!    Sur Ay

    ivar = isca(ipotva(2))
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    ! En periodique et parallele, echange avant calcul du gradient
    !     C'est indispensable car on vient de calculer IVAR
    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ivar))
      !==========
    endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)

    ivar0 = 0

!    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE

    iphydp = 0

    call grdcel                                                   &
    !==========
  ( idbia0 , idbra0 ,                                             &
    ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,&
    nnod   , lndfac , lndfbr , ncelbr , nphas  ,                  &
    nideve , nrdeve , nituse , nrtuse ,                           &
    ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp ,&
    iwarnp , nfecra , epsrgp , climgp , extrap ,                  &
    ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                  &
    ipnfac , nodfac , ipnfbr , nodfbr ,                           &
    idevel , ituser , ia     ,                                    &
    xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume ,&
    ra     , ra     , ra     ,                                    &
    rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv) ,              &
    w4     , w5     , w6     ,                                    &
!       d Ay /dx   d Ay /dy   d Ay /dz
    w7     , w8     , w9     ,                                    &
    rdevel , rtuser , ra     )

!       B = rot A

    do iel = 1, ncel
      w1(iel)= w1(iel) - w6(iel)
      w2(iel)= w2(iel) + zero
      w3(iel)= w3(iel) + w4(iel)
    enddo

!    Sur Az

    ivar = isca(ipotva(3))
    iclimv = iclrtp(ivar,icoef)

    inc = 1
    iccocg = 1
    nswrgp = nswrgr(ivar)
    imligp = imligr(ivar)
    iwarnp = iwarni(ivar)
    epsrgp = epsrgr(ivar)
    climgp = climgr(ivar)
    extrap = extrag(ivar)

    ! En periodique et parallele, echange avant calcul du gradient
    !     C'est indispensable car on vient de calculer IVAR
    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(rtp(1,ivar))
      !==========
    endif

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)

    ivar0 = 0

!    SANS PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE

    iphydp = 0

    call grdcel                                                   &
    !==========
  ( idbia0 , idbra0 ,                                             &
    ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml ,&
    nnod   , lndfac , lndfbr , ncelbr , nphas  ,                  &
    nideve , nrdeve , nituse , nrtuse ,                           &
    ivar0  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp ,&
    iwarnp , nfecra , epsrgp , climgp , extrap ,                  &
    ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                  &
    ipnfac , nodfac , ipnfbr , nodfbr ,                           &
    idevel , ituser , ia     ,                                    &
    xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume ,&
    ra     , ra     , ra     ,                                    &
    rtp(1,ivar), coefa(1,iclimv) , coefb(1,iclimv) ,              &
    w4     , w5     , w6     ,                                    &
!       d Az /dx   d Az /dy   d Az /dz
    w7     , w8     , w9     ,                                    &
    rdevel , rtuser , ra     )

!       B = rot A

    do iel = 1, ncel
      w1(iel)= w1(iel) + w5(iel)
      w2(iel)= w2(iel) - w4(iel)
      w3(iel)= w3(iel) + zero
    enddo

!   4.2 CAS D'ARC AXISYMETRIQUE : IELARC = 1
!  ------------------------------------------
!         ON PEUT UTILISER LE THEOREME D'AMPERE
!        ---------------------------------------

  else if(ippmod(ielarc).eq.1) then

!      Calcul du Champ magnetique calcule :
!===========================================

!       B = / j d S

!       Cette version n'a pas ete developpe pour l'instant :
!         - elle ne fonctionne que dans le cas d'un arc axisymetrique
!         - elle demande donc un reperage des plans perpendiculaires
!           a l'arc
!         - elle necessite la connaissance des coordonnees de chaque
!           point du maillage

    write(nfecra,2000)
    call csexit (1)

  endif

!   4.3 ARC ELECTRIQUE : CALCUL DES FORCES DE LAPLACE j X B
!  ---------------------------------------------------------
!         POUR TOUS LES CAS D'ARCS ELECTRIQUES
!        --------------------------------------

  if( ippmod(ielarc) .ge. 1 ) then
    do iel = 1, ncel
      propce(iel,ipcla1)= propce(iel,ipcdc2) * w3(iel)            &
                         -propce(iel,ipcdc3) * w2(iel)
      propce(iel,ipcla2)= propce(iel,ipcdc3) * w1(iel)            &
                         -propce(iel,ipcdc1) * w3(iel)
      propce(iel,ipcla3)= propce(iel,ipcdc1) * w2(iel)            &
                         -propce(iel,ipcdc2) * w1(iel)
    enddo
  endif


!   4.4 Impression de B en arc electrique
!  --------------------------------------

  if( ippmod(ielarc) .ge. 2 ) then

    if(modntl.eq.0) then

      write(nfecra,1000)

!     B=rot A
      var    = w1(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = w1(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      WRITE(NFECRA,1010)'Ch_MagX ',VRMIN,VRMAX

      var    = w2(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = w2(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      WRITE(NFECRA,1010)'Ch_MagY ',VRMIN,VRMAX

      var    = w3(1)
      vrmin = var
      vrmax = var
      do iel = 1, ncel
        var    = w3(iel)
        vrmin = min(vrmin,var)
        vrmax = max(vrmax,var)
      enddo
      if (irangp.ge.0) then
        call parmin (vrmin)
        call parmax (vrmax)
      endif
      WRITE(NFECRA,1010)'Ch_MagZ ',VRMIN,VRMAX

      write(nfecra,1001)

    endif

  endif

endif
!     Fin du test IAPPEL = 2

!--------
! FORMATS
!--------


 1000 format(/,                                                   &
'-----------------------------------------                    ',/,&
'   Variable         Minimum       Maximum                    ',/,&
'-----------------------------------------                    '  )
 1010 format(                                                           &
'v  ',A8,'    ',E12.5,'  ',E12.5                                 )
 1001 format(                                                           &
'-----------------------------------------                    '  )

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DANS ELTSSC            ',/,&
'@    =========                                               ',/,&
'@     LA VERSION ARC ELECTRIQUE AVEC THEOREME D''AMPERE      ',/,&
'@       (IELARC = 1)    N''EST PAS DISPONIBLE                ',/,&
'@                                                            ',/,&
'@            VEUILLEZ UTILISER LA VERSION DITE 3D            ',/,&
'@                      IELARC = 2                            ',/,&
'@   AVEC EQUATION DE TRANSPORT SUR LE POTENTIEL VECTEUR A    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


!----
! FIN
!----

return
end subroutine
