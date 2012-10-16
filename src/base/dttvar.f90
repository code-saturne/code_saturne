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

subroutine dttvar &
!================

 ( idbia0 , idbra0 ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iwarnp ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  , dam    , cofbdt , w1     , w2     , w3     , &
   coefbr , grarox , graroy , graroz , wcf    ,                   &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! CALCUL DU PAS DE TEMPS LOCAL
! AFFICHAGE DES NOMBRES DE COURANT + FOURIER MINIMUM, MAXIMUM
! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)

! Sous programme utilise dans le cas une seule phase (ou
! si seule la phase 1 pilote le pas de temps)
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
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iwarnp           ! i  ! <-- ! verbosity                                      !
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
! (ncesmp,nvar)    !    !     !  source de masse                               !
!                  !    !     ! pour ivar=ipr, smacel=flux de masse            !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! cofbdt(nfabor    ! tr ! --- ! condition limite pas de temps                  !
! w1,2,3(ncelet    ! tr ! --- ! tableaux de travail                            !
! graro.(ncelet    ! tr ! --- ! tableaux de travail (iptlro=1)                 !
! coefbr(nfabor    ! tr ! --- ! tableau de travail (iptlro=1)                  !
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

include "dimfbr.h"
include "paramx.h"
include "numvar.h"
include "cstnum.h"
include "cstphy.h"
include "optcal.h"
include "entsor.h"
include "parall.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse , iwarnp

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet ), cofbdt(nfabor)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
!   Attention, COEFBR n'est defini  que pour IPTLRO = 1
double precision coefbr(nfabor)
!   Attention, GRAROX, GRAROY, GRAROZ ne sont
!   definis que pour IPTLRO = 1 ou en compressible
double precision grarox(ncelet),graroy(ncelet),graroz(ncelet)
!   Attention, WCF n'est defini qu'en compressible
double precision wcf(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

character*8      cnom
integer          idebia, idebra
integer          ifac, iel, icfmax, icfmin, idiff0, iconv0, isym
integer          modntl
integer          iphas, iuiph, ipcvis, ipcvst
integer          iflmas, iflmab
integer          icou, ifou , icoucf
integer          inc, iccocg
integer          nswrgp, imligp , iphydp
integer          ipcrom, ipbrom, iivar
integer          nbrval
integer          ipccou, ipcfou
double precision epsrgp, climgp, extrap
double precision cfmax,cfmin, coufou, w1min, w2min, w3min
double precision unpvdt, rom
double precision xyzmax(3), xyzmin(3)
double precision dtsdtm,dtsdt0

!===============================================================================

!===============================================================================
! 0.  INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

iphas   = 1
iuiph   = iu(iphas)
iflmas  = ipprof(ifluma(iuiph))
iflmab  = ipprob(ifluma(iuiph))
ipcvis  = ipproc(iviscl(iphas))
ipcvst  = ipproc(ivisct(iphas))
ipcrom  = ipproc(irom  (iphas))
ipbrom  = ipprob(irom  (iphas))
ipccou  = ipproc(icour (iphas))
ipcfou  = ipproc(ifour (iphas))

if(ntlist.gt.0) then
  modntl = mod(ntcabs,ntlist)
elseif(ntlist.eq.-1.and.ntcabs.eq.ntmabs) then
  modntl = 0
else
  modntl = 1
endif

if (                                                              &
   .not. ( iconv(iuiph).ge.1.and.                                 &
           (iwarnp.ge.2.or.modntl.eq.0) ) .and.                   &
   .not. ( idiff(iuiph).ge.1.and.                                 &
           (iwarnp.ge.2.or.modntl.eq.0) ) .and.                   &
   .not. ( ippmod(icompf).ge.0.and.                               &
           (iwarnp.ge.2.or.modntl.eq.0) ) .and.                   &
   .not. ( idtvar.eq.-1.or.idtvar.eq.1.or.idtvar.eq.2.or.         &
           ( (iwarnp.ge.2.or.modntl.eq.0).and.                    &
             (idiff(iuiph).ge.1.or.iconv(iuiph).ge.1              &
                               .or.ippmod(icompf).ge.0)  ) )      &
   ) then

  return

endif

!===============================================================================
! 1.  CONDITION LIMITE POUR MATRDT
!===============================================================================


do ifac = 1, nfabor
  if(propfb(ifac,iflmab).lt.0.d0) then
    cofbdt(ifac) = 0.d0
  else
    cofbdt(ifac) = 1.d0
  endif
enddo

!===============================================================================
! 2.  CALCUL DE LA LIMITATION EN COMPRESSIBLE
!===============================================================================

!     On commence par cela afin de disposer de VISCF VISCB comme
!       tableaux de travail.

  if(ippmod(icompf).ge.0) then

    call cfdttv                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iwarnp ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   wcf    ,                                                       &
!        ---
   viscf  , viscb  , cofbdt , w1     , w2     , dam    ,          &
   grarox , graroy , graroz ,                                     &
   rdevel , rtuser , ra     )

  endif


!===============================================================================
! 3.  CALCUL DE LA VISCOSITE FACETTES
!===============================================================================


!     On s'en sert dans les divers matrdt suivants

!     "VITESSE" DE DIFFUSION FACETTE

if( idiff(iuiph).ge. 1 ) then
  do iel = 1, ncel
    w1    (iel) = propce(iel,ipcvis)                              &
                                +idifft(iuiph)*propce(iel,ipcvst)
  enddo
  call viscfa                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse , imvisf ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     ,                                                       &
   viscf  , viscb  ,                                              &
   rdevel , rtuser , ra     )

else
  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo
endif

!===============================================================================
! 4.  ALGORITHME INSTATIONNAIRE
!===============================================================================

if (idtvar.ge.0) then

!===============================================================================
! 4.1  PAS DE TEMPS VARIABLE A PARTIR DE COURANT ET FOURIER IMPOSES
!===============================================================================

!     On calcule le pas de temps thermique max (meme en IDTVAR=0, pour affichage)
!     DTTMAX = 1/SQRT(MAX(0+,gradRO.g/RO) -> W3

  if (iptlro.eq.1) then

    do ifac = 1, nfabor
      coefbr(ifac) = 0.d0
    enddo

    nswrgp = nswrgr(ipr(iphas))
    imligp = imligr(ipr(iphas))
    iwarnp = iwarni(ipr(iphas))
    epsrgp = epsrgr(ipr(iphas))
    climgp = climgr(ipr(iphas))
    extrap = 0.d0
    iphydp = 0

    iivar = 0
    inc   = 1
    iccocg = 1

    call grdcel                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iivar  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom), coefbr ,                   &
   grarox , graroy , graroz ,                                     &
!        ------   ------   ------
   w1     , w2     , dam    ,                                     &
   rdevel , rtuser , ra     )

    do iel = 1, ncel
      w3(iel) = (grarox(iel)*gx + graroy(iel)*gy + graroz(iel)*gz)&
           /propce(iel,ipcrom)
      w3(iel) = 1.d0/sqrt(max(epzero,w3(iel)))

    enddo

!     On met le nombre de clippings a 0 (il le restera pour IDTVAR=0)
    nclptr = 0

  endif


  if (idtvar.eq.1.or.idtvar.eq.2) then

    icou = 0
    ifou = 0

! 4.1.1 LIMITATION PAR LE COURANT
! =============================

    if ( coumax.gt.0.d0.and.iconv(iuiph).ge.1 ) then

!     ICOU = 1 marque l'existence d'une limitation par le COURANT
      icou = 1

! ---> CONSTRUCTION DE U/DX           (COURANT         ) =W1

      idiff0 = 0

!     Matrice a priori non symetrique
      isym = 2

      call matrdt                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iconv(iuiph)    , idiff0          , isym   ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   cofbdt , propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  , &
   dam    ,                                                       &
   rdevel , rtuser , ra     )

      do iel = 1, ncel
        rom = propce(iel,ipcrom)
        w1    (iel) = dam(iel)/(rom*volume(iel))
      enddo

! ---> CALCUL DE W1     = PAS DE TEMPS VARIABLE VERIFIANT
!       LE NOMBRE DE COURANT         MAXIMUM PRESCRIT PAR L'UTILISATEUR

      do iel = 1, ncel
        w1    (iel) = coumax/max( w1    (iel), epzero)
      enddo

! ---> PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE

      if (idtvar.eq.1) then
        w1min = grand
        do iel = 1, ncel
          w1min = min(w1min,w1(iel))
        enddo
        if (irangp.ge.0) then
          call parmin (w1min)
          !==========
        endif
        do iel = 1, ncel
          w1(iel) = w1min
        enddo
      endif

    endif

! 4.1.2 LIMITATION PAR LE FOURIER
! =============================

    if ( foumax.gt.0.d0.and.idiff(iuiph).ge.1 ) then

!     IFOU = 1 marque l'existence d'une limitation par le FOURIER
      ifou = 1

      iconv0 = 0
!                                   2
! ---> CONSTRUCTION DE      +2.NU/DX  (         FOURIER) =W2

!     Matrice a priori symetrique
      isym = 1

      call matrdt                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iconv0          , idiff(iuiph)    , isym   ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   cofbdt , propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  , &
   dam    ,                                                       &
   rdevel , rtuser , ra     )

      do iel = 1, ncel
        rom = propce(iel,ipcrom)
        w2    (iel) = dam(iel)/(rom*volume(iel))
      enddo

! ---> CALCUL DE W2     = PAS DE TEMPS VARIABLE VERIFIANT
!       LE NOMBRE DE         FOURIER MAXIMUM PRESCRIT PAR L'UTILISATEUR

      do iel = 1, ncel
        w2    (iel) = foumax/max( w2    (iel), epzero)
      enddo

! ---> PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE

      if (idtvar.eq.1) then
        w2min = grand
        do iel = 1, ncel
          w2min = min(w2min,w2(iel))
        enddo
        if (irangp.ge.0) then
          call parmin (w2min)
          !==========
        endif
        do iel = 1, ncel
          w2(iel) = w2min
        enddo
      endif

    endif

! 4.1.3 LIMITATION POUR L'ALGORITHME COMPRESSIBLE
! =============================================
!     Il est important de conserver WCF intact : on le reutilise
!     plus bas pour l'affichage

    icoucf = 0
    if ( coumax.gt.0.d0.and.ippmod(icompf).ge.0 ) then

      icoucf = 1

! ---> CALCUL DE DAM     = PAS DE TEMPS VARIABLE VERIFIANT
!       LA CONTRAINTE CFL MAXIMUM PRESCRITE PAR L'UTILISATEUR

      do iel = 1, ncel
        dam(iel) = coumax/max( wcf(iel), epzero)
      enddo

! ---> PAS DE TEMPS UNIFORME : ON PREND LE MINIMUM DE LA CONTRAINTE

      if (idtvar.eq.1) then
        w3min = grand
        do iel = 1, ncel
          w3min = min(w3min,dam(iel))
        enddo
        if (irangp.ge.0) then
          call parmin (w3min)
          !==========
        endif
        do iel = 1, ncel
          dam(iel) = w3min
        enddo
      endif

    endif

! 4.1.4 ON PREND LA PLUS CONTRAIGNANTE DES LIMITATIONS
! ==================================================
!    (le minimum des deux si elles existent et
!     celle qui existe s'il n'en existe qu'une)

    if(icou.eq.1.and.ifou.eq.1) then
      do iel = 1, ncel
        w1(iel) = min(w1(iel),w2(iel))
      enddo
    elseif(icou.eq.0.and.ifou.eq.1) then
      do iel = 1, ncel
        w1(iel) = w2(iel)
      enddo
    endif


!     En compressible, on prend obligatoirement
!     en compte la limitation associee à la masse volumique.

    if(icoucf.eq.1) then
      do iel = 1, ncel
        w1(iel) = min(w1(iel),dam(iel))
      enddo
    endif

! 4.1.5 ON CALCULE EFFECTIVEMENT LE PAS DE TEMPS
! ============================================

! --->  MONTEE           PROGRESSIVE DU PAS DE TEMPS
!              DESCENTE  IMMEDIATE   DU PAS DE TEMPS

    do iel = 1, ncel
      if( w1    (iel).ge.dt(iel) ) then
        unpvdt = 1.d0+varrdt
        dt(iel) = min( unpvdt*dt(iel), w1    (iel) )
      else
        dt(iel) =                      w1    (iel)
      endif
    enddo


! 4.1.6 ON LIMITE PAR LE PAS DE TEMPS "THERMIQUE" MAX
! =================================================
!     DTTMAX = W3 = 1/SQRT(MAX(0+,gradRO.g/RO)
!     on limite le pas de temps a DTTMAX

    if (iptlro.eq.1) then


!  On clippe le pas de temps a DTTMAX
!     (affiche dans ecrlis)

      nclptr = 0

      do iel = 1, ncel
        if ( dt(iel).gt.w3(iel) ) then
          nclptr = nclptr +1
          dt(iel) = w3(iel)
        endif
      enddo

      if (irangp.ge.0) then
        call parcpt (nclptr)
        !==========
      endif

! ---> PAS DE TEMPS UNIFORME : on reuniformise le pas de temps

      if (idtvar.eq.1) then
        w3min = grand
        do iel = 1, ncel
          w3min = min(w3min,dt(iel))
        enddo
        if (irangp.ge.0) then
          call parmin (w3min)
          !==========
        endif
        do iel = 1, ncel
          dt(iel) = w3min
        enddo
      endif

    endif

! 4.1.7 ON CLIPPE LE PAS DE TEMPS PAR RAPPORT A DTMIN ET DTMAX
! ==========================================================

    icfmin = 0
    icfmax = 0

    do iel = 1, ncel

      if( dt(iel).gt.dtmax      ) then
        icfmax = icfmax +1
        dt(iel) = dtmax
      endif
      if( dt(iel).lt.dtmin      ) then
        icfmin = icfmin +1
        dt(iel) = dtmin
      endif

    enddo

    if (irangp.ge.0) then
      call parcpt (icfmin)
      !==========
      call parcpt (icfmax)
      !==========
    endif

    iclpmx(ippdt) = icfmax
    iclpmn(ippdt) = icfmin

    if( iwarnp.ge.2) then
      write (nfecra,1003) icfmin,dtmin,icfmax,dtmax
    endif

  endif

!     Rapport DT sur DTmax lie aux effets de densite
!       (affichage dans ecrlis)
  if (iptlro.eq.1) then

    dtsdtm = 0.d0
    do iel = 1, ncel
      dtsdt0 = dt(iel)/w3(iel)
      if ( dtsdt0 .gt. dtsdtm ) then
        dtsdtm = dtsdt0
        icfmax = iel
      endif
    enddo
    xyzmax(1) = xyzcen(1,icfmax)
    xyzmax(2) = xyzcen(2,icfmax)
    xyzmax(3) = xyzcen(3,icfmax)

    if (irangp.ge.0) then
      nbrval = 3
      call parmxl (nbrval, dtsdtm, xyzmax)
      !==========
    endif
    rpdtro(1) = dtsdtm
    rpdtro(2) = xyzmax(1)
    rpdtro(3) = xyzmax(2)
    rpdtro(4) = xyzmax(3)

  endif

!===============================================================================
! 4.2  CALCUL DU NOMBRE DE COURANT POUR AFFICHAGE
!===============================================================================

  if ( iconv(iuiph).ge.1.and.                                     &
       (iwarnp.ge.2.or.modntl.eq.0) ) then

    idiff0 = 0
    CNOM   =' COURANT'

!     CONSTRUCTION DE U/DX           (COURANT         ) =W1

! MATRICE A PRIORI NON SYMETRIQUE

    isym = 2

    call matrdt                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iconv(iuiph)    , idiff0          , isym   ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   cofbdt , propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  , &
   dam    ,                                                       &
   rdevel , rtuser , ra     )

    do iel = 1, ncel
      rom = propce(iel,ipcrom)
      w1    (iel) = dam(iel)/(rom*volume(iel))
    enddo

!     CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax = -grand
    cfmin =  grand
    icfmax= 0
    icfmin= 0

    do iel = 1, ncel

      coufou = w1(iel)*dt(iel)
      propce(iel,ipccou) = coufou

      if( coufou.le.cfmin ) then
        cfmin  = coufou
        icfmin = iel
      endif

      if( coufou.ge.cfmax ) then
        cfmax  = coufou
        icfmax = iel
      endif

    enddo

    xyzmin(1) = xyzcen(1,icfmin)
    xyzmin(2) = xyzcen(2,icfmin)
    xyzmin(3) = xyzcen(3,icfmin)
    xyzmax(1) = xyzcen(1,icfmax)
    xyzmax(2) = xyzcen(2,icfmax)
    xyzmax(3) = xyzcen(3,icfmax)

    if (irangp.ge.0) then
      nbrval = 3
      call parmnl (nbrval, cfmin, xyzmin)
      !==========
      call parmxl (nbrval, cfmax, xyzmax)
      !==========
    endif

    if(iwarnp.ge.2) then
      write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)
    endif

!       -> pour listing
    ptploc(1,1) = cfmin
    ptploc(1,2) = xyzmin(1)
    ptploc(1,3) = xyzmin(2)
    ptploc(1,4) = xyzmin(3)
    ptploc(2,1) = cfmax
    ptploc(2,2) = xyzmax(1)
    ptploc(2,3) = xyzmax(2)
    ptploc(2,4) = xyzmax(3)

  endif

!===============================================================================
! 4.3  CALCUL DU NOMBRE DE FOURIER POUR AFFICHAGE
!===============================================================================

  if ( idiff(iuiph).ge.1.and.                                     &
       (iwarnp.ge.2.or.modntl.eq.0) ) then

    iconv0 = 0
    CNOM   =' FOURIER'
!                                   2
!     CONSTRUCTION DE      +2.NU/DX  (         FOURIER) =W1

! MATRICE A PRIORI SYMETRIQUE

    isym = 1

    call matrdt                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iconv0          , idiff(iuiph)    , isym   ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   cofbdt , propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  , &
   dam    ,                                                       &
   rdevel , rtuser , ra     )

    do iel = 1, ncel
      rom = propce(iel,ipcrom)
      w1    (iel) = dam(iel)/(rom*volume(iel))
    enddo

!     CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax  = -grand
    cfmin  =  grand
    icfmax = 0
    icfmin = 0

    do iel = 1, ncel

      coufou = w1(iel)*dt(iel)
      propce(iel,ipcfou) = coufou

      if( coufou.le.cfmin ) then
        cfmin  = coufou
        icfmin = iel
      endif

      if( coufou.ge.cfmax ) then
        cfmax  = coufou
        icfmax = iel
      endif

    enddo

    xyzmin(1) = xyzcen(1,icfmin)
    xyzmin(2) = xyzcen(2,icfmin)
    xyzmin(3) = xyzcen(3,icfmin)
    xyzmax(1) = xyzcen(1,icfmax)
    xyzmax(2) = xyzcen(2,icfmax)
    xyzmax(3) = xyzcen(3,icfmax)

    if (irangp.ge.0) then
      nbrval = 3
      call parmnl (nbrval, cfmin, xyzmin)
      !==========
      call parmxl (nbrval, cfmax, xyzmax)
      !==========
    endif

    if(iwarnp.ge.2) then
      write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)
    endif

!       -> pour listing
    ptploc(3,1) = cfmin
    ptploc(3,2) = xyzmin(1)
    ptploc(3,3) = xyzmin(2)
    ptploc(3,4) = xyzmin(3)
    ptploc(4,1) = cfmax
    ptploc(4,2) = xyzmax(1)
    ptploc(4,3) = xyzmax(2)
    ptploc(4,4) = xyzmax(3)

  endif

!===============================================================================
! 4.4  CALCUL DU NOMBRE DE COURANT/FOURIER POUR AFFICHAGE
!===============================================================================

!     En incompressible uniquement (en compressible, on preferera
!       afficher la contrainte liee a la masse volumique)

  if ( (iwarnp.ge.2.or.modntl.eq.0).and.                          &
      (idiff(iuiph).ge.1.or.iconv(iuiph).ge.1)                    &
      .and.(ippmod(icompf).lt.0)               ) then

    CNOM   =' COU/FOU'
!                                   2
!     CONSTRUCTION DE U/DX +2.NU/DX  (COURANT +FOURIER) =W1

! MATRICE A PRIORI NON SYMETRIQUE

    isym = 1
    if (iconv(iuiph).gt.0) isym = 2

    call matrdt                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iconv(iuiph)    , idiff(iuiph)    , isym   ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   cofbdt , propfa(1,iflmas), propfb(1,iflmab), viscf  , viscb  , &
   dam    ,                                                       &
   rdevel , rtuser , ra     )

    do iel = 1, ncel
      rom = propce(iel,ipcrom)
      w1    (iel) = dam(iel)/(rom*volume(iel))
    enddo

!     CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax  = -grand
    cfmin  =  grand
    icfmax = 0
    icfmin = 0

    do iel = 1, ncel

      coufou = w1(iel)*dt(iel)

      if( coufou.le.cfmin ) then
        cfmin  = coufou
        icfmin = iel
      endif

      if( coufou.ge.cfmax ) then
        cfmax  = coufou
        icfmax = iel
      endif

    enddo

    xyzmin(1) = xyzcen(1,icfmin)
    xyzmin(2) = xyzcen(2,icfmin)
    xyzmin(3) = xyzcen(3,icfmin)
    xyzmax(1) = xyzcen(1,icfmax)
    xyzmax(2) = xyzcen(2,icfmax)
    xyzmax(3) = xyzcen(3,icfmax)

    if (irangp.ge.0) then
      nbrval = 3
      call parmnl (nbrval, cfmin, xyzmin)
      !==========
      call parmxl (nbrval, cfmax, xyzmax)
      !==========
    endif

    if(iwarnp.ge.2) then
      write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)
    endif

!       -> pour listing
    ptploc(5,1) = cfmin
    ptploc(5,2) = xyzmin(1)
    ptploc(5,3) = xyzmin(2)
    ptploc(5,4) = xyzmin(3)
    ptploc(6,1) = cfmax
    ptploc(6,2) = xyzmax(1)
    ptploc(6,3) = xyzmax(2)
    ptploc(6,4) = xyzmax(3)

  endif

!===============================================================================
! 4.5  CALCUL DE LA CONTRAINTE CFL DE LA MASSE VOL. POUR AFFICHAGE
!===============================================================================

! En Compressible uniquement

  if ( (iwarnp.ge.2.or.modntl.eq.0).and.                          &
       (ippmod(icompf).ge.0)                        ) then

    CNOM   =' CFL/MAS'


!     CALCUL DU NOMBRE DE COURANT/FOURIER MAXIMUM ET MINIMUM

    cfmax  = -grand
    cfmin  =  grand
    icfmax = 0
    icfmin = 0

    do iel = 1, ncel

      coufou = wcf(iel)*dt(iel)

      if( coufou.le.cfmin ) then
        cfmin  = coufou
        icfmin = iel
      endif

      if( coufou.ge.cfmax ) then
        cfmax  = coufou
        icfmax = iel
      endif

    enddo

    xyzmin(1) = xyzcen(1,icfmin)
    xyzmin(2) = xyzcen(2,icfmin)
    xyzmin(3) = xyzcen(3,icfmin)
    xyzmax(1) = xyzcen(1,icfmax)
    xyzmax(2) = xyzcen(2,icfmax)
    xyzmax(3) = xyzcen(3,icfmax)

    if (irangp.ge.0) then
      nbrval = 3
      call parmnl (nbrval, cfmin, xyzmin)
      !==========
      call parmxl (nbrval, cfmax, xyzmax)
      !==========
    endif

    if(iwarnp.ge.2) then
      write(nfecra,1001) cnom,cfmax,xyzmax(1),xyzmax(2),xyzmax(3)
      write(nfecra,1002) cnom,cfmin,xyzmin(1),xyzmin(2),xyzmin(3)
    endif

!       -> pour listing
    ptploc(5,1) = cfmin
    ptploc(5,2) = xyzmin(1)
    ptploc(5,3) = xyzmin(2)
    ptploc(5,4) = xyzmin(3)
    ptploc(6,1) = cfmax
    ptploc(6,2) = xyzmax(1)
    ptploc(6,3) = xyzmax(2)
    ptploc(6,4) = xyzmax(3)

  endif

!===============================================================================
! 5.   ALGORITHME STATIONNAIRE
!===============================================================================
else

  isym = 1
  if (iconv(iuiph).gt.0) isym = 2

  call matrdt                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iconv(iuiph)    , idiff(iuiph)    , isym,                      &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   coefb(1,iuiph)  , propfa(1,iflmas), propfb(1,iflmab),          &
                                                viscf  , viscb  , &
   dt     , rdevel , rtuser , ra )

  do iel = 1, ncel
    dt(iel) = relaxv(iuiph)*propce(iel,ipcrom)                    &
         *volume(iel)/max(dt(iel),epzero)
  enddo

endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1001 FORMAT ( /,A8,' MAX= ',E11.4,                                     &
 ' EN ',E11.4,' ',E11.4,' ',E11.4)
 1002 FORMAT (   A8,' MIN= ',E11.4,                                     &
 ' EN ',E11.4,' ',E11.4,' ',E11.4)
 1003 FORMAT ( /,'CLIPPINGS DE DT : ',                                  &
                             I10,' A ',E11.4,', ',I10,' A ',E11.4)

#else

 1001 FORMAT ( /,A8,' MAX= ',E11.4,                                     &
 ' IN ',E11.4,' ',E11.4,' ',E11.4)
 1002 FORMAT (   A8,' MIN= ',E11.4,                                     &
 ' IN ',E11.4,' ',E11.4,' ',E11.4)
 1003 FORMAT ( /,'DT CLIPPING : ',                                      &
                             I10,' A ',E11.4,', ',I10,' A ',E11.4)

#endif

!----
! FIN
!----

return

end subroutine
