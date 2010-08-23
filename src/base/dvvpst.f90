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

subroutine dvvpst &
!================

 ( idbia0 , idbra0 , nummai , numtyp ,                            &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta , nvisbr ,                   &
   ncelps , nfacps , nfbrps ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   itypps , ifacel , ifabor , ifmfbr , ifmcel , iprfml ,          &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statce , stativ , statfb ,                   &
   tracel , trafac , trafbr , rdevel , rtuser , ra     )
!===============================================================================
! FONCTION :
! --------

! ROUTINE DEVELOPPEUR POUR LA SORTIE STANDARD
! DES VALEURS SUR LES MAILLAGES DE POST TRAITEMENT
!   (appelee apres usproj)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nummai           ! ec ! <-- ! numero du maillage post                        !
! numtyp           ! ec ! <-- ! numero de type de post-traitement              !
!                  !    !     ! (-1: volume, -2: bord, nummai par defaut)      !
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
! nvlsta           ! e  ! <-- ! nombre de variables stat. lagrangien           !
! nvisbr           ! e  ! <-- ! nombre de statistiques aux frontieres          !
! ncelps           ! e  ! <-- ! nombre de cellules du maillage post            !
! nfacps           ! e  ! <-- ! nombre de faces interieur post                 !
! nfbrps           ! e  ! <-- ! nombre de faces de bord post                   !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! itypps(3)        ! te ! <-- ! indicateur de presence (0 ou 1) de             !
!                  !    !     ! cellules (1), faces (2), ou faces de           !
!                  !    !     ! de bord (3) dans le maillage post              !
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
! lstcel(ncelps    ! te ! <-- ! liste des cellules du maillage post            !
! lstfac(nfacps    ! te ! <-- ! liste des faces interieures post               !
! lstfbr(nfbrps    ! te ! <-- ! liste des faces de bord post                   !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! te ! --- ! macro tableau entier                           !
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
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet)         !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! statce           ! tr ! <-- ! statistiques cellules (lagrangien)             !
!(ncelet,nvlsta    !    !     !                                                !
! stativ           ! tr ! <-- ! statistiques variance cellules                 !
!(ncelet,nvlsta    !    !     !                          (lagrangien)          !
! statfb           ! tr ! <-- ! statistiques faces bord (lagrangien)           !
!(nfabor,nvisbr    !    !     !                                                !
! tracel(*)        ! tr ! <-- ! tab reel valeurs cellules post                 !
! trafac(*)        ! tr ! <-- ! tab reel valeurs faces int. post               !
! trafbr(*)        ! tr ! <-- ! tab reel valeurs faces bord post               !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use optcal
use numvar
use parall
use period
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use radiat
use cplsat

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nummai , numtyp
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , nvlsta , nvisbr
integer          ncelps , nfacps , nfbrps
integer          nideve , nrdeve , nituse , nrtuse

integer          itypps(3)
integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          lstcel(ncelps), lstfac(nfacps), lstfbr(nfbrps)
integer          idevel(nideve), ituser(nituse), ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtpa(ncelet,*), rtp(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision statce(ncelet,nvlsta), statfb(nfabor,nvisbr)
double precision stativ(ncelet,nvlsta)
double precision tracel(ncelps*3)
double precision trafac(nfacps*3), trafbr(nfbrps*3)
double precision rdevel(nrdeve), rtuser(nrtuse)
double precision ra(*)

! Local variables

character*32     namevr, namev1, namev2
character*80     name80

integer          idebia, idebra, ifinia, ifinra
integer          igradx, igrady, igradz
integer          itravx, itravy, itravz, itreco
integer          iw1   , iw2
integer          inc   , iccocg, nswrgp, imligp, iwarnp, iphydp
integer          isorva, isaut
integer          ifac  , iloc  , iphas , ivar , iclvar
integer          ira   , idivdt, ineeyp
integer          ipp   , idimt , ii    , kk   , iel
integer          ivarl , iip   , iph
integer          iii, ivarl1 , ivarlm , iflu   , ilpd1  , icla
integer          iscal , ipcvsl, ipcvst, iflmab
integer          ientla, ivarpr
integer          iyplbp
integer          ipccp , ipcrom

double precision cp0iph, xcp   , xvsl  , surfbn, distbr
double precision visct , flumab, diipbx, diipby, diipbz
double precision epsrgp, climgp, extrap
double precision omgnrm, daxis2
double precision rbid(1)


!===============================================================================

! Initialize variables to avoid compiler warnings

ipp = 0
cp0iph = 0

! Memoire

idebia = idbia0
idebra = idbra0

!===============================================================================
!     1.1. TRAITEMENT POUR LE MAILLAGE FLUIDE
!===============================================================================

if (numtyp .eq. -1) then


!       1.1.1 TRAITEMENT DES VARIABLES POST TRAITABLES
!       ----------------------------------------------

  do ipp = 2, nvppmx

!         -> si chrono demande sur la variable
    if(ichrvr(ipp).eq.1) then

!           -> pointeur de la variable dans RA (en absolu)
!              (si negatif, c'est un vecteur)
      ira = ipp2ra(ipp)
!           -> si c'est un moment cumule, il faut diviser par le temps
!             (si   0 ce n'est pas un moment,
!              si > 0 c'est le pointeur dans RA,
!              si < 0 c'est le rang dans DTCMOM)
      idivdt = ippmom(ipp)
!           -> dimension de la variable a ecrire
      idimt = 1
      if(ira.lt.0) then
        idimt = 3
        ira = -ira
      endif

!           -> nom de la variable
      name80 = nomvar(ipp)
      namevr = name80(1:32)

!           -> si l'on a une variable vectorielle, on supprime
!              la partie X, Y, ou Z en derniere ou avant-derniere
!              position dans le nom

      if (idimt.eq.3) then
        name80 = nomvar(ipp+1)
        namev1 = name80(1:32)
        name80 = nomvar(ipp+2)
        namev2 = name80(1:32)
        call pstsnv ( namevr , namev1 , namev2 )
        !==========
      endif

!           -> si c'est un moment cumule, il faut diviser par le temps
      if(idivdt.gt.0) then
        do iloc = 1, ncelps
          iel = lstcel(iloc)
          tracel(iloc) = ra(ira+iel-1)/                           &
               max(ra(idivdt+iel-1),epzero)
        enddo
      elseif(idivdt.lt.0) then
        do iloc = 1, ncelps
          iel = lstcel(iloc)
          tracel(iloc) = ra(ira+iel-1)/                           &
               max(dtcmom(-idivdt),epzero)
        enddo
      endif

!           Ecriture effective des valeurs calculees

!           Valeurs non entrelacées
      ientla = 0

      if(idivdt.eq.0) then
        ivarpr = 1
        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
        !==========
                    ntcabs, ttcabs, ra(ira), rbid, rbid)

      else
        ivarpr = 0
        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
        !==========
                    ntcabs, ttcabs, tracel, rbid, rbid)

      endif

    endif
!         Fin du traitement en cas de sortie de la variable

  enddo
!       Fin de la boucle sur les variables

!       1.1.2 VARIABLES SUPPLEMENTAIRES AUTOMATIQUES
!       --------------------------------------------

!       Distance a la paroi (si LES+VanDriest ou Rij+Echo ou K-w SST)

  if (ineedy.eq.1 .and. abs(icdpar).eq.1) then

    NAMEVR = 'DistParoi'
    idimt = 1
    ientla = 0
    ivarpr = 1

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
    !==========
                ntcabs, ttcabs, ra(idipar), rbid, rbid)

  endif

!       Yplus (si LES+VanDriest)

  if (ineedy.eq.1 .and. abs(icdpar).eq.1) then

    ineeyp = 0
    do iphas = 1, nphas
      if(itytur(iphas).eq.4.and.idries(iphas).eq.1) then
        ineeyp = 1
      endif
    enddo

    if (ineeyp.eq.1) then

      NAMEVR = 'Yplus'
      idimt = 1
      ientla = 0
      ivarpr = 1

      call psteva(nummai, namevr, idimt, ientla, ivarpr,          &
      !==========
                  ntcabs, ttcabs, ra(iyppar), rbid, rbid)

    endif

  endif

  ! Vitesse et pression absolues en cas de calcul en repère relatif

  if (icorio.eq.1) then

    iphas = 1
    ipcrom = ipproc(irom(iphas))
    omgnrm = sqrt(omegax**2 + omegay**2 + omegaz**2)

    NAMEVR = 'Pressure'
    idimt = 1
    ientla = 0
    ivarpr = 0

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      daxis2 =   (omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))**2 &
               + (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))**2 &
               + (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))**2

      daxis2 = daxis2 / omgnrm**2

      tracel(iloc) = rtp(iel,ipr(iphas))                          &
          + 0.5d0*propce(iel,ipcrom)*omgnrm**2*daxis2

    enddo

    call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)


    NAMEVR = 'Velocity'
    idimt = 3
    ientla = 1
    ivarpr = 0

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      tracel(1 + (iloc-1)*idimt) = rtp(iel,iu(iphas)) &
          + (omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))

      tracel(2 + (iloc-1)*idimt) = rtp(iel,iv(iphas)) &
          + (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))

      tracel(3 + (iloc-1)*idimt) = rtp(iel,iw(iphas)) &
          + (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))

    enddo

    call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)
  endif

  ! Vitesse et pression relatives en cas de calcul en repère mobile

  if (imobil.eq.1) then

    iphas = 1
    ipcrom = ipproc(irom(iphas))
    omgnrm = sqrt(omegax**2 + omegay**2 + omegaz**2)

    NAMEVR = 'Rel Pressure'
    idimt = 1
    ientla = 0
    ivarpr = 0

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      daxis2 =   (omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))**2 &
               + (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))**2 &
               + (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))**2

      daxis2 = daxis2 / omgnrm**2

      tracel(iloc) = rtp(iel,ipr(iphas))                          &
          - 0.5d0*propce(iel,ipcrom)*omgnrm**2*daxis2

    enddo

    call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)


    NAMEVR = 'Rel Velocity'
    idimt = 3
    ientla = 1
    ivarpr = 0

    do iloc = 1, ncelps

      iel = lstcel(iloc)

      tracel(1 + (iloc-1)*idimt) = rtp(iel,iu(iphas)) &
          - (omegay*xyzcen(3,iel) - omegaz*xyzcen(2,iel))

      tracel(2 + (iloc-1)*idimt) = rtp(iel,iv(iphas)) &
          - (omegaz*xyzcen(1,iel) - omegax*xyzcen(3,iel))

      tracel(3 + (iloc-1)*idimt) = rtp(iel,iw(iphas)) &
          - (omegax*xyzcen(2,iel) - omegay*xyzcen(1,iel))

    enddo

    call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
    !==========
                ntcabs, ttcabs, tracel, rbid, rbid)

  endif


!===============================================================================
!     1.2. TRAITEMENT POUR LE MAILLAGE DE BORD
!===============================================================================

else if  (numtyp .eq. -2) then


! --    1.2.1 TRAITEMENT DE YPLUS AU BORD
!       ----------------------------------

  if(mod(ipstdv,ipstyp).eq.0) then

!       Phase
    do iphas = 1, nphas

!       Initialisation
      do ii = 1, 32
        NAMEVR (II:II) = ' '
      enddo

!       Nom de la variable
      if(nphas.gt.1) then
        write(namevr,1000)iphas
 1000         format('YplusPhase',I2.2)
      else
        NAMEVR = 'Yplus'
      endif

!       Dimension de la variable (3 = vecteur, 1=scalaire)
      idimt = 1

!       Calcul des valeurs de la variable sur les faces de bord

      iyplbp = iyplbr+(iphas-1)*nfabor
      do iloc = 1, nfbrps
        ifac = lstfbr(iloc)
        trafbr(1 + (iloc-1)*idimt) = ra(iyplbp+ifac-1)
      enddo

!           Valeurs non entrelacées, définies sur tableau de travail
      ientla = 0
      ivarpr = 0

      call psteva(nummai, namevr, idimt, ientla, ivarpr,          &
      !==========
                  ntcabs, ttcabs, rbid, rbid, trafbr)

    enddo
!     fin du test sur les phases

  endif
! fin du test sur sortie de yplus


! --    1.2.2 TRAITEMENT DES VARIABLES AU BORD SANS RECONSTRUCTION
!       -----------------------------------------------------------

  if(mod(ipstdv,ipstcl).eq.0) then

!       Le codage ci-dessous est relativement avance :
!         il accede aux variables directement au travers du macro tableau RA
!         il comprend un artifice "ISAUT" permettant de reperer les vecteurs
!           a posttraiter


!       Boucle sur les variables usuelles
!       ---------------------------------

    isaut = 0

    do ivar = 1, nvar


!         Variables post-traitables
!         (ISAUT est utilise pour ne pas postraiter plusieurs fois les
!          composantes 2 et 3 d'un vecteur, initialise a 0 avant la boucle
!          sur IVAR)

      ipp = ipprtp(ivar)

      isorva = 0
      if (isaut .gt. 0) then
        isaut = isaut - 1
        isorva = 0
      else if (ichrvr(ipp).eq.1) then
        isorva = 1
        name80 = nomvar(ipp)
        namevr = name80(1:32)
      endif

!         Traitement des variables definies aux centres cellules à sortir
!         ---------------------------------------------------------------
      if (isorva .eq. 1) then


!         -> on verifie le signe du pointeur de la variable dans RA
!         (si negatif, c'est un vecteur)

!         -> dimension de la variable a ecrire
        idimt = 1
        if(ipp2ra(ipp).lt.0) then
          idimt = 3
          isaut = 2
        endif

!           -> si l'on a une variable vectorielle, on supprime
!              la partie X, Y, ou Z en derniere ou avant-derniere
!              position dans le nom

      if (idimt.eq.3) then
        name80 = nomvar(ipp+1)
        namev1 = name80(1:32)
        name80 = nomvar(ipp+2)
        namev2 = name80(1:32)
        call pstsnv ( namevr , namev1 , namev2 )
        !==========
      endif

!         Calcul des valeurs (non-reconstruites) de la variable
!         sur les faces de bord

        do kk = 0, idimt-1

          iclvar = iclrtp(ivar+kk,icoef)
          do iloc = 1, nfbrps

            ifac = lstfbr(iloc)
            iel = ifabor(ifac)

            trafbr(kk + (iloc-1)*idimt + 1)                       &
                 =   coefa(ifac,iclvar)                           &
                   + coefb(ifac,iclvar)*rtp(iel,ivar+kk)

          enddo

        enddo

!             Valeurs entrelacées, définies sur tableau de travail
        ientla = 1
        ivarpr = 0

        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
        !==========
                    ntcabs, ttcabs, rbid, rbid, trafbr)

      endif
!         Fin du traitement en cas de sortie de la variable

    enddo
!       Fin de la boucle sur les variables

  endif
!     Fin du test sur sortie des variables



!       1.2.3 TRAITEMENT FLUX THERMIQUE AU BORD
!       ----------------------------------------
!           Si on travaille en enthalpie, on calcule un flux d'enthalpie

  if(mod(ipstdv,ipstft).eq.0) then

!       Phase
    do iphas = 1, nphas

      if(iscalt(iphas).gt.0 .and. nscal.gt.0 .and.                &
           iscalt(iphas).le.nscal) then

!       Initialisation
        do ii = 1, 32
          NAMEVR (II:II) = ' '
        enddo

!       Nom de la variable
        if(nphas.gt.1) then
          write(namevr,2000)iphas
 2000           format('Flux th. entrant W.m-2 Phase',I2.2)
        else
          NAMEVR = 'Flux thermique entrant W.m-2'
        endif

!       Dimension de la variable (3 = vecteur, 1=scalaire)
        idimt = 1

!       Numero de la variable

        iscal  = iscalt(iphas)
        ivar   = isca(iscal)
        iclvar = iclrtp(ivar,icoef)

!       Calcul des valeurs de la variable sur les faces de bord

!          Reservation de la memoire pour reconstruction

        ifinia = idebia

        igradx = idebra
        igrady = igradx+ncelet
        igradz = igrady+ncelet
        itravx = igradz+ncelet
        itravy = itravx+ncelet
        itravz = itravy+ncelet
        itreco = itravz+ncelet
        ifinra = itreco+nfabor

!          Verification de la disponibilite de la memoire

        call iasize('dvvpst',ifinia)
        call rasize('dvvpst',ifinra)


!          Calcul du gradient de la temperature / enthalpie


!      Pour calculer le gradient de Temperature
!        - dans les calculs paralleles, il est necessaire que
!          les cellules situees sur un bord de sous-domaine connaissent
!          la valeur de temperature dans les cellules situees en
!          vis-a-vis sur le sous-domaine voisin.
!        - dans les calculs periodiques, il est necessaire que
!          les cellules periodiques aient acces a la valeur de la
!          temperature des cellules periodiques correspondantes

!      Pour cela, il est necessaire d'appeler les routines de
!        de synchronisation des halos pour echanger les valeurs de temperature
!        avant de calculer le gradient.
!      En effet, on se situe ici a la fin du pas de temps n. Or,
!        les variables RTP ne seront echangees qu'en debut du pas de
!        temps n+1. Ici, seules les variables RTPA (obtenues a la fin
!        du pas de temps n-1) ont deja ete echangees.

!      Si le calcul n'est ni periodique, ni parallele, on peut conserver
!        appels (les tests sur IPERIO et IRANGP assurent la generalite)


!          Echange pour le parallelisme et la periodicite

        if (irangp.ge.0.or.iperio.eq.1) then
          call synsca(rtp(1,ivar))
          !==========
        endif


!          Calcul du gradient

        inc = 1
        iccocg = 1
        nswrgp = nswrgr(ivar)
        imligp = imligr(ivar)
        iwarnp = iwarni(ivar)
        epsrgp = epsrgr(ivar)
        climgp = climgr(ivar)
        extrap = extrag(ivar)
        iphydp = 0

        call grdcel                                               &
        !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   ra(itravx) , ra(itravx) , ra(itravx) ,                         &
   rtp(1,ivar) , coefa(1,iclvar) , coefb(1,iclvar) ,              &
   ra(igradx) , ra(igrady) , ra(igradz) ,                         &
!        ----------   ----------   ----------
   ra(itravx) , ra(itravy) , ra(itravz) ,                         &
   rdevel , rtuser , ra     )


!          Calcul de la valeur reconstruite dans les cellules de bord

        do ifac = 1, nfabor
          iel = ifabor(ifac)
          iii = idiipb-1+3*(ifac-1)
          diipbx = ra(iii+1)
          diipby = ra(iii+2)
          diipbz = ra(iii+3)
          ra(itreco+ifac-1) = rtp(iel,ivar)                       &
               + diipbx*ra(igradx+iel-1)                          &
               + diipby*ra(igrady+iel-1)                          &
               + diipbz*ra(igradz+iel-1)
        enddo

!          Calcul du flux (ouf          !) convectif et diffusif

        if(ivisls(iscal).gt.0) then
          ipcvsl = ipproc(ivisls(iscal))
        else
          ipcvsl = 0
        endif
        ipcvst = ipproc(ivisct(iphas))
        iflmab = ipprob(ifluma(ivar))

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          iel = ifabor(ifac)

          if(ipcvsl.gt.0) then
            xvsl = propce(iel,ipcvsl)
          else
            xvsl = visls0(iscal)
          endif
          surfbn = ra(isrfbn-1+ifac)
          distbr = ra(idistb-1+ifac)
          visct  = propce(iel,ipcvst)
          flumab = propfb(ifac,iflmab)

          trafbr(1 + (iloc-1)*idimt) =                            &
               (xvsl+visct/sigmas(iscal))/max(distbr,epzero)*     &
               (coefa(ifac,iclvar)+(coefb(ifac,iclvar)-1.d0)*     &
               rtp(iel,ivar))                                     &
               - flumab/max(surfbn,epzero**2)*                    &
               (coefa(ifac,iclvar)+ coefb(ifac,iclvar)*           &
               rtp(iel,ivar))

        enddo

!          Pour la temperature, on multiplie par CP
        if(abs(iscsth(iscal)).eq.1) then
          if(icp(iphas).gt.0) then
            ipccp  = ipproc(icp   (iphas))
          else
            ipccp  = 0
            cp0iph = cp0(iphas)
          endif
          do iloc = 1, nfbrps
            ifac = lstfbr(iloc)
            iel = ifabor(ifac)
            if(ipccp.gt.0) then
              xcp = propce(iel,ipccp)
            else
              xcp    = cp0iph
            endif
            trafbr(1 + (iloc-1)*idimt)                            &
                 = xcp*trafbr(1 + (iloc-1)*idimt)
          enddo
        endif

!             Valeurs entrelacées, définies sur tableau de travail
        ientla = 1
        ivarpr = 0

        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
        !==========
                    ntcabs, ttcabs, rbid, rbid, trafbr)

      endif
!         Fin du test sur variable thermique

    enddo
!       Fin de boucle sur les phases

  endif
!      Fin du test sur sortie des flux thermiques

! --    1.2.4 TRAITEMENT DES EFFORTS AUX BORDS
!       --------------------------------------

  if(mod(ipstdv,ipstfo).eq.0) then

!       Initialisation
    do ii = 1, 32
      NAMEVR (II:II) = ' '
    enddo

!       Nom de la variable
    NAMEVR = 'Efforts'

!       Dimension de la variable (3 = vecteur, 1=scalaire)
    idimt = 3

!       Calcul des valeurs de la variable sur les faces de bord

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      surfbn = ra(isrfbn-1+ifac)
      trafbr(1 + (iloc-1)*idimt ) =                               &
             ra(iforbr+(ifac-1)*idimt  )/surfbn
      trafbr(2 + (iloc-1)*idimt ) =                               &
             ra(iforbr+(ifac-1)*idimt+1)/surfbn
      trafbr(3 + (iloc-1)*idimt ) =                               &
             ra(iforbr+(ifac-1)*idimt+2)/surfbn
    enddo

!           Valeurs entrelacées, définies sur tableau de travail
      ientla = 1
      ivarpr = 0

      call psteva(nummai, namevr, idimt, ientla, ivarpr,          &
      !==========
                  ntcabs, ttcabs, rbid, rbid, trafbr)

  endif
! fin du test sur sortie des efforts

endif
!     Fin du test sur le numero de maillage post.

!===============================================================================
!     2.1. VARIABLES LAGRANGIENNES
!===============================================================================

if (nummai .eq. -1) then

  if (iilagr.gt.0 .and. istala.ge.1) then

!         Toutes les statistiques standard sont de dimension 1,
!         et sont definies ou calculees sur tableau de travail
!         de maniere non entrelacee (sans importance de toutes
!         manieres pour une variable scalaire)

    idimt  = 1
    ientla = 0

    iii = nvlsta-nvlsts

    do icla  = 0, nbclst

!     -> si IPAS = 0 : statistiques globales
!        si 0 < IPAS =< NBCLST : statistiques par groupe

      do ivarl = 1, nvlsta

        ivarl1 = icla*nvlsta +ivarl
        ivarlm = ivarl1
        ilpd1  = icla*nvlsta +ilpd
        iflu   = 0

        if (ivarl.le.iii) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlag(ivarl)
          else
            WRITE(NAME80,'(A8,A4,I3)') NOMLAG(IVARL),'_grp',ICLA
          endif
        else if (nvlsts.gt.0) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlag(ilvu(ivarl-iii))
          else
            WRITE(NAME80,'(A8,A4,I3)')                            &
                 NOMLAG(ILVU(IVARL-III)),'_grp',ICLA
          endif
        endif

        namevr = name80(1:32)

        call uslaen                                               &
        !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivarl  , ivarl1 , ivarlm , iflu   , ilpd1  , icla   ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statce , stativ , tracel ,                   &
   rdevel , rtuser , ra     )

!           La variable est deja definie sur le maillage volumique
!           global ; on utilise donc l'indirection  (donc IVARPR = 1)
        ivarpr = 1

        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
                    ntcabs, ttcabs, tracel, rbid, rbid)
      enddo

      do ivarl = 1, nvlsta-1

        ivarl1 = icla*(nvlsta-1)+ivarl
        ivarlm = icla*nvlsta+ivarl
        ilpd1  = icla*nvlsta +ilpd
        iflu   = 1

        if (ivarl.le.iii) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlav(ivarl)
          else
            WRITE(NAME80,'(A8,A4,I3)') NOMLAV(IVARL),'_grp',ICLA
          endif
        else if (nvlsts.gt.0) then
          if (ivarl.eq.ivarl1) then
            name80 = nomlav(ilvu(ivarl-iii))
          else
            WRITE(NAME80,'(A8,A4,I3)')                            &
                 NOMLAV(ILVU(IVARL-III)),'_grp',ICLA
          endif
        endif

        namevr = name80(1:32)

        call uslaen                                               &
        !==========
 ( ifinia , ifinra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , nvlsta ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivarl  , ivarl1 , ivarlm , iflu   , ilpd1  , icla   ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  , statce , stativ , tracel ,                   &
   rdevel , rtuser , ra     )

!           La variable est deja definie sur le maillage volumique
!           global ; on utilise donc l'indirection  (donc IVARPR = 1)
        ivarpr = 1

        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
                    ntcabs, ttcabs, tracel, rbid, rbid)
      enddo

    enddo

  endif

endif

if (nummai.eq.-2) then

  if (iilagr.gt.0 .and. iensi3.eq.1) then

    iii = nvisbr-nusbor

    do ivarl = 1,nvisbr

      if (ivarl.le.iii) then
        name80 = nombrd(ivarl)
      else if (nusbor.gt.0) then
        name80 = nombrd(iusb(ivarl-iii))
      endif
      namevr = name80(1:32)

      if (imoybr(ivarl).eq.2) then

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,inbr).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl)/statfb(ifac,inbr)
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      else if (imoybr(ivarl).eq.1) then

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,inbr).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl) / tstatp
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      else

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          if (statfb(ifac,inbr).gt.seuilf) then
            trafbr(iloc) = statfb(ifac,ivarl)
          else
            trafbr(iloc) = 0.d0
          endif
        enddo

      endif

      idimt  = 1
      ientla = 0
      ivarpr = 0

      call psteva(nummai, namevr, idimt, ientla, ivarpr,          &
                  ntcabs, ttcabs, rbid, rbid, trafbr)

    enddo

    NAME80 = 'lagrangian_boundary_zones'
    namevr = name80(1:32)

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      trafbr(iloc) = ia(iifrla+ifac-1)
    enddo

    idimt  = 1
    ientla = 0
    ivarpr = 0

    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
         ntcabs, ttcabs, rbid, rbid, trafbr)

  endif
endif
!     Fin du test sur le numero de maillage post.

!===============================================================================
!     2.2. VARIABLES RADIATIVES AUX FRONTIERES
!===============================================================================


if (nummai.eq.-2) then

  if (iirayo.gt.0) then

    do ivarl = 1,nbrayf

      if (irayvf(ivarl).eq.1) then

        name80 = nbrvaf(ivarl)
        namevr = name80(1:32)

        if (ivarl .eq. itparp)      then
          ipp =  ipprob(itparo)
        else if (ivarl .eq. iqincp) then
          ipp = ipprob(iqinci)
        else if (ivarl .eq. ixlamp)  then
          ipp = ipprob(ixlam)
        else if (ivarl .eq. iepap)   then
          ipp = ipprob(iepa)
        else if (ivarl .eq. iepsp)   then
          ipp = ipprob(ieps)
        else if (ivarl .eq. ifnetp)  then
          ipp = ipprob(ifnet)
        else if (ivarl .eq. ifconp) then
          ipp = ipprob(ifconv)
        else if (ivarl .eq. ihconp) then
          ipp = ipprob(ihconv)
        endif

        do iloc = 1, nfbrps
          ifac = lstfbr(iloc)
          trafbr(iloc) = propfb(ifac,ipp)
        enddo

        idimt  = 1
        ientla = 0
        ivarpr = 0

        call psteva(nummai, namevr, idimt, ientla, ivarpr,        &
                    ntcabs, ttcabs, rbid, rbid, trafbr)

      endif
    enddo

    name80 = 'radiative_boundary_zones'
    namevr = name80(1:32)

    do iloc = 1, nfbrps
      ifac = lstfbr(iloc)
      trafbr(iloc) = ia(iizfrd+ifac-1)
    enddo

    idimt  = 1
    ientla = 0
    ivarpr = 0
!
    call psteva(nummai, namevr, idimt, ientla, ivarpr,            &
         ntcabs, ttcabs, rbid, rbid, trafbr)

  endif
endif

!===============================================================================
!     2.3. VARIABLES ELECTRIQUES
!===============================================================================

if (     ippmod(ieljou).ge.1                                      &
    .or. ippmod(ielarc).ge.1                                      &
    .or. ippmod(ielion).ge.1) then

  ifinia = idebia

  iw1    = idebra
  iw2    = iw1    + ncelet*3
  ifinra = iw2    + ncelet*3

  call rasize ('dvvpst', ifinra)
  !==========

  call uselen                                                     &
  !==========
 ( ifinia , ifinra , nummai ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncelps , nfacps , nfbrps ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   lstcel , lstfac , lstfbr ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , rtp    , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ra(iw1) , ra(iw2) ,                                            &
   tracel , trafac , trafbr , rdevel , rtuser , ra     )

endif

!--------
! FORMATS
!--------

!----
! FIN
!----

return
end subroutine
