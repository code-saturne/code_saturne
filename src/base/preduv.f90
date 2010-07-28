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

subroutine preduv &
!================

 ( idbia0 , idbra0 , iappel ,                                     &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , iterns , ncepdp , ncesmp ,          &
   nideve , nrdeve , nituse , nrtuse , iphas  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   flumas , flumab ,                                              &
   tslagr , coefa  , coefb  ,                                     &
   ckupdc , smacel , frcxt  ,                                     &
   trava  , ximpa  , uvwk   , dfrcxt , tpucou , trav   ,          &
   viscf  , viscb  , viscfi , viscbi ,                            &
   dam    , xam    ,                                              &
   drtp   , smbr   , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , xnormp , coefu  ,                   &
   rdevel , rtuser , ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS N-S 1 PHASE INCOMPRESSIBLE OU RO VARIABLE
! SUR UN PAS DE TEMPS (CONVECTION/DIFFUSION - PRESSION /CONTINUITE)

! AU PREMIER APPEL,  ON EFFECTUE LA PREDICITION DES VITESSES
!               ET  ON CALCULE UN ESTIMATEUR SUR LA VITESSE PREDITE

! AU DEUXIEME APPEL, ON CALCULE UN ESTIMATEUR GLOBAL
!               POUR NAVIER-STOKES :
!   ON UTILISE TRAV, SMBR ET LES TABLEAUX DE TRAVAIL
!   ON APPELLE BILSC2 AU LIEU DE CODITS
!   ON REMPLIT LE PROPCE ESTIMATEUR IESTOT
!   CE DEUXIEME APPEL INTERVIENT DANS NAVSTO APRES RESOLP
!   LORS DE CE DEUXIEME APPEL
!     RTPA ET RTP SONT UN UNIQUE TABLEAU (= RTP)
!     LE FLUX DE MASSE EST LE FLUX DE MASSE DEDUIT DE LA VITESSE
!      AU CENTRE CONTENUE DANS RTP

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! iappel           ! e  ! <-- ! numero d'appel du sous programme               !
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
! iterns           ! e  ! <-- ! numero d'iteration sur navsto                  !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! iphas            ! i  ! <-- ! phase number                                   !
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
! flumas           ! tr ! <-- ! flux de masse aux faces internes               !
!  (nfac  )        !    !     !   (distinction iappel=1 ou 2)                  !
! flumab           ! tr ! <-- ! flux de masse aux faces de bord                !
!  (nfabor  )      !    !     !    (distinction iappel=1 ou 2)                 !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
! frcxt(ncelet,    ! tr ! <-- ! force exterieure generant la pression          !
!   3,nphas)       !    !     !  hydrostatique                                 !
! trava,ximpa      ! tr ! <-- ! tableau de travail pour couplage               !
! uvwk             ! tr ! <-- ! tableau de travail pour couplage u/p           !
!                  !    !     ! sert a stocker la vitesse de                   !
!                  !    !     ! l'iteration precedente                         !
!dfrcxt(ncelet,    ! tr ! --> ! variation de force exterieure                  !
!   3,nphas)       !    !     !  generant la pression hydrostatique            !
! tpucou           ! tr ! --> ! couplage vitesse pression                      !
! (ncelel,ndim)    !    !     !                                                !
! trav(ncelet,3    ! tr ! --> ! smb qui servira pour normalisation             !
!                  !    !     !  dans resolp                                   !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! viscfi(nfac)     ! tr ! --- ! idem viscf pour increments                     !
! viscbi(nfabor    ! tr ! --- ! idem viscb pour increments                     !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr  (ncelet    ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! xnormp(ncelet    ! tr ! <-- ! tableau reel pour norme resolp                 !
! coefu(nfab,3)    ! tr ! --- ! tableau de travail                             !
! tslagr(ncelet    ! tr ! <-- ! terme de couplage retour du                    !
!   ntersl)        !    !     !   lagrangien                                   !
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
include "pointe.h"
include "entsor.h"
include "cstphy.h"
include "cstnum.h"
include "optcal.h"
include "period.h"
include "parall.h"
include "lagpar.h"
include "lagran.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "matiss.h"
include "cplsat.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0 , iappel
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas  , iterns
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse , iphas

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
double precision flumas(nfac), flumab(nfabor)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision frcxt(ncelet,3,nphas), dfrcxt(ncelet,3,nphas)
double precision trava(ncelet,ndim,nphas)
double precision ximpa(ncelet,ndim,nphas),uvwk(ncelet,ndim,nphas)
double precision tpucou(ncelet,ndim), trav(ncelet,3)
double precision viscf(nfac), viscb(nfabor)
double precision viscfi(nfac), viscbi(nfabor)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet)
double precision smbr(ncelet), rovsdt(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision xnormp(ncelet)
double precision coefu(nfabor,3)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra, ifinia
integer          iel   , ielpdc, ifac  , ivar  , isou
integer          iccocg, inc   , init  , iphydp, ii    , isqrt
integer          ireslp, nswrgp, imligp, iwarnp, ippt  , ipp
integer          ipriph, ikiph , iuiph , iviph , iwiph
integer                  iclipr, icliup, iclivp, icliwp
integer                          iclik , iclvar, iclvaf
integer          ipcrom, ipcroa, ipcroo , ipcvis, ipcvst
integer          iconvp, idiffp, ndircp, nitmap, nswrsp
integer          ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          iesprp, iestop
integer          iptsna
integer          iflmb0, nswrp , imaspe, iismph, ipbrom
integer          idiaex, idtva0
integer          maxelt, ils
double precision rnorm , vitnor
double precision romvom, ro0iph, drom
double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision vit1  , vit2  , vit3, xkb, pip, pfac, pfac1
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23
double precision d2s3  , thetap, thetp1, thets , dtsrom
double precision diipbx, diipby, diipbz
double precision cx    , cy    , cz

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

ipriph = ipr(iphas)
iuiph  = iu(iphas)
iviph  = iv(iphas)
iwiph  = iw(iphas)
if(itytur(iphas).eq.2 .or. iturb(iphas).eq.50                     &
     .or. iturb(iphas).eq.60) then
  ikiph  = ik(iphas)
endif

iclipr = iclrtp(ipriph,icoef)
icliup = iclrtp(iuiph ,icoef)
iclivp = iclrtp(iviph ,icoef)
icliwp = iclrtp(iwiph ,icoef)

if(itytur(iphas).eq.2 .or. iturb(iphas).eq.50                     &
     .or. iturb(iphas).eq.60) then
  iclik  = iclrtp(ikiph ,icoef)
else
  iclik = 0
endif

ro0iph = ro0(iphas)

!     Reperage de rho au bord
ipbrom = ipprob(irom  (iphas))
!     Reperage de rho courant (ie en cas d'extrapolation rho^n+1/2)
ipcrom = ipproc(irom  (iphas))
!     Reperage de rho^n en cas d'extrapolation
if(iroext(iphas).gt.0) then
  ipcroa = ipproc(iroma(iphas))
else
  ipcroa = 0
endif

ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))
!     Theta relatif aux termes sources explicites
thets  = thetsn(iphas)
if(isno2t(iphas).gt.0) then
  iptsna = ipproc(itsnsa(iphas))
else
  iptsna = 0
endif


!===============================================================================
! 2.  GRADIENT DE PRESSION ET GRAVITE
!===============================================================================

! ---> PRISE EN COMPTE DE LA PRESSION HYDROSTATIQUE :
!       UNIQUEMENT AU PREMIER APPEL (I.E. POUR LE CALCUL NORMAL,
!       LE DEUXIEME APPEL EST DESTINE A UN CALCUL D'ESTIMATEUR)

if (iappel.eq.1.and.iphydr.eq.1) then

! force ext au pas de temps precedent :
!     FRCXT a ete initialise a zero
!     (est deja utilise dans typecl, et est mis a jour a la fin
!     de navsto)

  do iel = 1, ncel

! variation de force (utilise dans resolp)
    drom = (propce(iel,ipcrom)-ro0iph)
    dfrcxt(iel,1,iphas) = drom*gx - frcxt(iel,1,iphas)
    dfrcxt(iel,2,iphas) = drom*gy - frcxt(iel,2,iphas)
    dfrcxt(iel,3,iphas) = drom*gz - frcxt(iel,3,iphas)
  enddo
!     Ajout eventuel des pertes de charges
  if (ncepdp.gt.0) then
    do ielpdc = 1, ncepdp
      iel=icepdc(ielpdc)
      vit1   = rtpa(iel,iuiph)
      vit2   = rtpa(iel,iviph)
      vit3   = rtpa(iel,iwiph)
      cpdc11 = ckupdc(ielpdc,1)
      cpdc22 = ckupdc(ielpdc,2)
      cpdc33 = ckupdc(ielpdc,3)
      cpdc12 = ckupdc(ielpdc,4)
      cpdc13 = ckupdc(ielpdc,5)
      cpdc23 = ckupdc(ielpdc,6)
      dfrcxt(iel,1,iphas) = dfrcxt(iel,1,iphas)                   &
           -propce(iel,ipcrom)*(                                  &
           cpdc11*vit1+cpdc12*vit2+cpdc13*vit3)
      dfrcxt(iel,2,iphas) = dfrcxt(iel,2,iphas)                   &
           -propce(iel,ipcrom)*(                                  &
           cpdc12*vit1+cpdc22*vit2+cpdc23*vit3)
      dfrcxt(iel,3,iphas) = dfrcxt(iel,3,iphas)                   &
           -propce(iel,ipcrom)*(                                  &
           cpdc13*vit1+cpdc23*vit2+cpdc33*vit3)
    enddo
  endif
!     Ajout eventuel de la force de Coriolis
  if (icorio.eq.1) then
    do iel = 1, ncel
      cx = omegay*rtpa(iel,iwiph) - omegaz*rtpa(iel,iviph)
      cy = omegaz*rtpa(iel,iuiph) - omegax*rtpa(iel,iwiph)
      cz = omegax*rtpa(iel,iviph) - omegay*rtpa(iel,iuiph)
      dfrcxt(iel,1,iphas) = dfrcxt(iel,1,iphas) - 2.d0*propce(iel,ipcrom)*cx
      dfrcxt(iel,2,iphas) = dfrcxt(iel,2,iphas) - 2.d0*propce(iel,ipcrom)*cy
      dfrcxt(iel,3,iphas) = dfrcxt(iel,3,iphas) - 2.d0*propce(iel,ipcrom)*cz
    enddo
  endif

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(dfrcxt(1,1,iphas), dfrcxt(1,2,iphas), dfrcxt(1,3,iphas))
    !==========
  endif

endif



! ---> PRISE EN COMPTE DU GRADIENT DE PRESSION

iccocg = 1
inc    = 1
nswrgp = nswrgr(ipriph)
imligp = imligr(ipriph)
iwarnp = iwarni(ipriph)
epsrgp = epsrgr(ipriph)
climgp = climgr(ipriph)
extrap = extrag(ipriph)


call grdcel                                                       &
!==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ipriph , imrgra , inc    , iccocg , nswrgp , imligp , iphydr , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   frcxt(1,1,iphas), frcxt(1,2,iphas), frcxt(1,3,iphas),          &
   rtpa(1,ipriph)  , coefa(1,iclipr) , coefb(1,iclipr) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )


!    Calcul des efforts aux parois (partie 2/5), si demande
!    La pression a la face est calculee comme dans gradrc/gradmc
!    et on la transforme en pression totale
!    On se limite a la premiere iteration (pour faire simple par
!      rapport a la partie issue de condli, hors boucle)
if (ineedf.eq.1 .and. iterns.eq.1) then
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    ii = idiipb-1+3*(ifac-1)
    diipbx = ra(ii+1)
    diipby = ra(ii+2)
    diipbz = ra(ii+3)
    pip =  rtpa(iel,ipriph)                                       &
         +diipbx*w1(iel) +diipby*w2(iel)                          &
         +diipbz*w3(iel)
    pfac = coefa(ifac,iclipr) +coefb(ifac,iclipr)*pip
    pfac1= rtpa(iel,ipriph)                                       &
         +(cdgfbo(1,ifac)-xyzcen(1,iel))*w1(iel)                  &
         +(cdgfbo(2,ifac)-xyzcen(2,iel))*w2(iel)                  &
         +(cdgfbo(3,ifac)-xyzcen(3,iel))*w3(iel)
    pfac = coefb(ifac,iclipr)*(extrag(ipriph)*pfac1               &
         +(1.d0-extrag(ipriph))*pfac)                             &
         +(1.d0-coefb(ifac,iclipr))*pfac                          &
         + ro0(iphas)*(gx*(cdgfbo(1,ifac)-xyzp0(1,iphas))         &
         + gy*(cdgfbo(2,ifac)-xyzp0(2,iphas))                     &
         + gz*(cdgfbo(3,ifac)-xyzp0(3,iphas)) )                   &
         - pred0(iphas)
! on ne rajoute pas P0, pour garder un maximum de precision
!     &         + P0(IPHAS)
    do isou = 1, 3
      ra(iforbr+(ifac-1)*ndim + isou-1) =                         &
           ra(iforbr+(ifac-1)*ndim + isou-1)                      &
           + pfac*surfbo(isou,ifac)
    enddo
  enddo
endif

! ---> RESIDU DE NORMALISATION POUR RESOLP
!     Test d'un residu de normalisation de l'etape de pression
!       plus comprehensible = div(rho u* + dt gradP^(n))-Gamma
!       i.e. second membre du systeme en pression hormis la partie
!       pression (sinon a convergence, on tend vers 0)
!       Represente les termes que la pression doit equilibrer
!     On calcule ici div(rho dt/rho gradP^(n)) et on complete a la fin
!       avec  div(rho u*)
!     Pour grad P^(n) on suppose que des CL de Neumann homogenes
!       s'appliquent partout : on peut donc utiliser les CL de la
!       vitesse pour u*+dt/rho gradP^(n). Comme on calcule en deux fois,
!       on utilise les CL de vitesse homogenes pour dt/rho gradP^(n)
!       ci-dessous et les CL de vitesse completes pour u* a la fin.

if(iappel.eq.1.and.irnpnw.eq.1) then

!     Calcul de dt/rho*grad P
  do iel = 1, ncel
    dtsrom = dt(iel)/propce(iel,ipcrom)
    trav(iel,1) = w1(iel)*dtsrom
    trav(iel,2) = w2(iel)*dtsrom
    trav(iel,3) = w3(iel)*dtsrom
  enddo

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(trav(1,1), trav(1,2), trav(1,3))
    !==========
  endif

!     Calcul de rho dt/rho*grad P.n aux faces
!       Pour gagner du temps, on ne reconstruit pas.
  init   = 1
  inc    = 0
  iccocg = 1
  iflmb0 = 1
  nswrp  = 1
  imligp = imligr(iuiph )
  iwarnp = iwarni(ipriph)
  epsrgp = epsrgr(iuiph )
  climgp = climgr(iuiph )
  extrap = extrag(iuiph )

  imaspe = 1

  iismph = iisymp+nfabor*(iphas-1)

  call inimas                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   iuiph  , iviph  , iwiph  , imaspe , iphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrp  , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(iismph) ,               &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   trav(1,1) , trav(1,2) , trav(1,3) ,                            &
   coefa(1,icliup), coefa(1,iclivp), coefa(1,icliwp),             &
   coefb(1,icliup), coefb(1,iclivp), coefb(1,icliwp),             &
   viscf  , viscb  ,                                              &
   w4     , w5     , w6     , w7     , w8     , w9     ,          &
   smbr   , drtp   , rovsdt , coefu  ,                            &
   rdevel , rtuser , ra     )

!     Calcul de div(rho dt/rho*grad P)
  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
                                ifacel,ifabor,viscf,viscb,xnormp)

!     Ajout de -Gamma
  if (ncesmp.gt.0) then
    do ii = 1, ncesmp
      iel = icetsm(ii)
      xnormp(iel) = xnormp(iel)-volume(iel)*smacel(ii,ipriph)
    enddo
  endif

!     On conserve XNORMP, on complete avec u* a la fin et
!       on le transfere a resolp

endif


!     Au premier appel, TRAV est construit directement ici.
!     Au second  appel (estimateurs), TRAV contient deja
!       l'increment temporel.
!     On pourrait fusionner en initialisant TRAV a zero
!       avant le premier appel, mais ca fait des operations en plus.

!     Remarques :
!       rho g sera a l'ordre 2 s'il est extrapole.
!       si on itere sur navsto, ca ne sert a rien de recalculer rho g a
!         chaque fois (ie on pourrait le passer dans trava) mais ce n'est
!         pas cher.
if(iappel.eq.1) then

  if (iphydr.eq.1) then
    do iel = 1, ncel
      trav(iel,1) = (frcxt(iel,1,iphas) - w1(iel)) * volume(iel)
      trav(iel,2) = (frcxt(iel,2,iphas) - w2(iel)) * volume(iel)
      trav(iel,3) = (frcxt(iel,3,iphas) - w3(iel)) * volume(iel)
    enddo
  else
    do iel = 1, ncel
      drom = (propce(iel,ipcrom)-ro0iph)
      trav(iel,1) = ( drom*gx - w1(iel) ) * volume(iel)
      trav(iel,2) = ( drom*gy - w2(iel) ) * volume(iel)
      trav(iel,3) = ( drom*gz - w3(iel) ) * volume(iel)
    enddo
  endif

elseif(iappel.eq.2) then

  if (iphydr.eq.1) then
    do iel = 1, ncel
      trav(iel,1) =                                               &
           trav(iel,1) + ( frcxt(iel,1,iphas) -                   &
                           w1(iel) )*volume(iel)
      trav(iel,2) =                                               &
           trav(iel,2) + ( frcxt(iel,2,iphas) -                   &
                           w2(iel) )*volume(iel)
      trav(iel,3) =                                               &
           trav(iel,3) + ( frcxt(iel,3,iphas) -                   &
                           w3(iel) )*volume(iel)
    enddo
  else
    do iel = 1, ncel
      drom = (propce(iel,ipcrom)-ro0iph)
      trav(iel,1) =                                               &
           trav(iel,1) + ( drom*gx - w1(iel) )*volume(iel)
      trav(iel,2) =                                               &
           trav(iel,2) + ( drom*gy - w2(iel) )*volume(iel)
      trav(iel,3) =                                               &
           trav(iel,3) + ( drom*gz - w3(iel) )*volume(iel)
    enddo
  endif

endif

!   Pour IAPPEL = 1 (ie appel standard sans les estimateurs)
!     TRAV rassemble les termes sources  qui seront recalcules
!       a toutes les iterations sur navsto
!     Si on n'itere pas sur navsto et qu'on n'extrapole pas les
!       termes sources, TRAV contient tous les termes sources
!       jusqu'au basculement dans SMBR
!     A ce niveau, TRAV contient -grad P et rho g
!       P est suppose pris a n+1/2
!       rho est eventuellement interpole a n+1/2


! ---> INITIALISATION DU TABLEAU TRAVA et PROPCE AU PREMIER PASSAGE
!     (A LA PREMIERE ITER SUR NAVSTO)

!     TRAVA rassemble les termes sources qu'il suffit de calculer
!       a la premiere iteration sur navsto quand il y a plusieurs iter.
!     Quand il n'y a qu'une iter, on cumule directement dans TRAV
!       ce qui serait autrement alle dans TRAVA
!     PROPCE rassemble les termes sources explicites qui serviront
!       pour le pas de temps suivant en cas d'extrapolation (plusieurs
!       iter sur navsto ou pas)

!     A la premiere iter sur navsto
if(iterns.eq.1) then

!       Si on   extrapole     les T.S. : -theta*valeur precedente
    if(isno2t(iphas).gt.0) then
!         S'il n'y a qu'une    iter : TRAV  incremente
      if(nterup.eq.1) then
        do ii = 1, ndim
          do iel = 1, ncel
            trav (iel,ii)        = trav (iel,ii)                  &
                 - thets*propce(iel,iptsna+ii-1)
          enddo
        enddo
!         S'il   y a plusieurs iter : TRAVA initialise
      else
        do ii = 1, ndim
          do iel = 1, ncel
            trava(iel,ii,iphas)  =                                &
                 - thets*propce(iel,iptsna+ii-1)
          enddo
        enddo
      endif
!         Et on initialise PROPCE pour le remplir ensuite
      do ii = 1, ndim
        do iel = 1, ncel
          propce(iel,iptsna+ii-1) = 0.d0
        enddo
      enddo

!       Si on n'extrapole pas les T.S. : pas de PROPCE
    else
!         S'il   y a plusieurs iter : TRAVA initialise
!           sinon TRAVA n'existe pas
      if(nterup.gt.1) then
        do ii = 1, ndim
          do iel = 1, ncel
            trava(iel,ii,iphas)  = 0.d0
          enddo
        enddo
      endif
    endif

  endif

! ---> 2/3 RHO * GRADIENT DE K SI k-epsilon ou k-omega
!      NB : ON NE PREND PAS LE GRADIENT DE (RHO K), MAIS
!           CA COMPLIQUERAIT LA GESTION DES CL ...
!     On peut se demander si l'extrapolation en temps sert a
!       quelquechose

!     Ce terme explicite est calcule une seule fois,
!       a la premiere iter sur navsto : il va dans PROPCE si on
!       doit l'extrapoler en temps ; il va dans TRAVA si on n'extrapole
!       pas mais qu'on itere sur navsto. Il va dans TRAV si on
!       n'extrapole pas et qu'on n'itere pas sur navsto.
if( (itytur(iphas).eq.2 .or. iturb(iphas).eq.50                   &
     .or. iturb(iphas).eq.60)                                     &
       .and.igrhok(iphas).eq.1.and.iterns.eq.1)then
  iccocg = 1
  inc    = 1
  nswrgp = nswrgr(ikiph)
  imligp = imligr(ikiph)
  epsrgp = epsrgr(ikiph)
  climgp = climgr(ikiph)
  extrap = extrag(ikiph)

  iwarnp = iwarni(iuiph)
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ikiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w6     , w6     , w6     ,                                     &
   rtpa(1,ikiph)   , coefa(1,iclik)  , coefb(1,iclik)  ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

  d2s3 = 2.d0/3.d0

!     Si on extrapole les termes source en temps : PROPCE
  if(isno2t(iphas).gt.0) then
!       Calcul de rho^n grad k^n      si rho non extrapole
!                 rho^n grad k^n si rho     extrapole
    ipcroo = ipcrom
    if(iroext(iphas).gt.0) ipcroo = ipcroa
    do iel = 1, ncel
      romvom = -propce(iel,ipcroo)*volume(iel)*d2s3
      propce(iel,iptsna  )=propce(iel,iptsna  )+w1(iel)*romvom
      propce(iel,iptsna+1)=propce(iel,iptsna+1)+w2(iel)*romvom
      propce(iel,iptsna+2)=propce(iel,iptsna+2)+w3(iel)*romvom
    enddo
!     Si on n'extrapole pas les termes sources en temps : TRAV ou TRAVA
  else
    if(nterup.eq.1) then
      do iel = 1, ncel
        romvom = -propce(iel,ipcrom)*volume(iel)*d2s3
        trav (iel,1)       =                                      &
        trav (iel,1)       + w1(iel) * romvom
        trav (iel,2)       =                                      &
        trav (iel,2)       + w2(iel) * romvom
        trav (iel,3)       =                                      &
        trav (iel,3)       + w3(iel) * romvom
      enddo
    else
      do iel = 1, ncel
        romvom = -propce(iel,ipcrom)*volume(iel)*d2s3
        trava(iel,1,iphas) =                                      &
        trava(iel,1,iphas) + w1(iel) * romvom
        trava(iel,2,iphas) =                                      &
        trava(iel,2,iphas) + w2(iel) * romvom
        trava(iel,3,iphas) =                                      &
        trava(iel,3,iphas) + w3(iel) * romvom
      enddo
    endif
  endif

!    Calcul des efforts aux parois (partie 3/5), si demande
  if (ineedf.eq.1) then
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      ii = idiipb-1+3*(ifac-1)
      diipbx = ra(ii+1)
      diipby = ra(ii+2)
      diipbz = ra(ii+3)
      xkb = rtpa(iel,ikiph) + diipbx*w1(iel)                      &
           + diipby*w2(iel) + diipbz*w3(iel)
      xkb = coefa(ifac,iclik)+coefb(ifac,iclik)*xkb
      xkb = d2s3*propce(iel,ipcrom)*xkb
      do isou = 1, 3
        ra(iforbr+(ifac-1)*ndim + isou-1) =                       &
             ra(iforbr+(ifac-1)*ndim + isou-1)                    &
             + xkb*surfbo(isou,ifac)
      enddo
    enddo
  endif
endif


! ---> TERMES DE GRADIENT TRANSPOSE

!     Ce terme explicite est calcule une seule fois,
!       a la premiere iter sur navsto :
!       si on extrapole il va dans PROPCE,
!       sinon si on itere sur navsto dans TRAVA
!             sinon                  dans TRAV
if (ivisse(iphas).eq.1.and.iterns.eq.1) then

!     On utilise temporairement TRAV comme tableau de travail.
!     Son contenu est stocke dans W7, W8 et W9 jusqu'apres vissec
  do iel = 1,ncel
    w7(iel) = trav(iel,1)
    w8(iel) = trav(iel,2)
    w9(iel) = trav(iel,3)
    trav(iel,1) = 0.d0
    trav(iel,2) = 0.d0
    trav(iel,3) = 0.d0
  enddo

  call vissec                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iphas  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , icepdc , icetsm , itypsm , &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   trav   ,                                                       &
!        ------
   viscf  , viscb  , rovsdt ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes de gradient transpose et
!       TRAV retrouve sa valeur
  if(isno2t(iphas).gt.0) then
    do iel = 1, ncel
      propce(iel,iptsna  ) = propce(iel,iptsna  ) + trav(iel,1)
      propce(iel,iptsna+1) = propce(iel,iptsna+1) + trav(iel,2)
      propce(iel,iptsna+2) = propce(iel,iptsna+2) + trav(iel,3)
      trav(iel,1)  = w7(iel)
      trav(iel,2)  = w8(iel)
      trav(iel,3)  = w9(iel)
    enddo

!     Si on n'extrapole pas les termes source en temps :
!       Si on itere sur navsto
!         TRAVA recoit les termes de gradient transpose et
!         TRAV retrouve sa valeur
  else
    if(nterup.gt.1) then
      do iel = 1, ncel
        trava(iel,1,iphas) = trava(iel,1,iphas) + trav(iel,1)
        trava(iel,2,iphas) = trava(iel,2,iphas) + trav(iel,2)
        trava(iel,3,iphas) = trava(iel,3,iphas) + trav(iel,3)
        trav(iel,1)  = w7(iel)
        trav(iel,2)  = w8(iel)
        trav(iel,3)  = w9(iel)
      enddo
!       Si on n'itere pas sur navsto
!         TRAV retrouve sa valeur augmentee des termes de gradient transpose
    else
      do iel = 1, ncel
        trav(iel,1)  = w7(iel) + trav(iel,1)
        trav(iel,2)  = w8(iel) + trav(iel,2)
        trav(iel,3)  = w9(iel) + trav(iel,3)
      enddo
    endif
  endif

endif


! ---> TERMES DE PERTES DE CHARGE
!     SI IPHYDR=1 LE TERME A DEJA ETE PRIS EN COMPTE AVANT

if((ncepdp.gt.0).and.(iphydr.eq.0)) then

!     Les termes diagonaux sont places dans TRAV ou TRAVA,
!       La prise en compte de UVWK a partir de la seconde iteration
!       est faite directement dans codits.
  if(iterns.eq.1) then

!       On utilise temporairement TRAV comme tableau de travail.
!       Son contenu est stocke dans W7, W8 et W9 jusqu'apres tsepdc
    do iel = 1,ncel
      w7(iel) = trav(iel,1)
      w8(iel) = trav(iel,2)
      w9(iel) = trav(iel,3)
      trav(iel,1) = 0.d0
      trav(iel,2) = 0.d0
      trav(iel,3) = 0.d0
    enddo

    idiaex = 1
    call tsepdc                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdp ,                                                       &
   nideve , nrdeve , nituse , nrtuse , iphas  , idiaex ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , trav   ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

!     Si on itere sur navsto, on utilise TRAVA ; sinon TRAV
    if(nterup.gt.1) then
      do iel = 1, ncel
        trava(iel,1,iphas) = trava(iel,1,iphas) + trav(iel,1)
        trava(iel,2,iphas) = trava(iel,2,iphas) + trav(iel,2)
        trava(iel,3,iphas) = trava(iel,3,iphas) + trav(iel,3)
        trav(iel,1)  = w7(iel)
        trav(iel,2)  = w8(iel)
        trav(iel,3)  = w9(iel)
      enddo
    else
      do iel = 1, ncel
        trav(iel,1)  = w7(iel) + trav(iel,1)
        trav(iel,2)  = w8(iel) + trav(iel,2)
        trav(iel,3)  = w9(iel) + trav(iel,3)
      enddo
    endif
  endif

!     Les termes extradiagonaux ne sont calcules qu'au premier passage
!       si on extrapole il va dans PROPCE,
!       sinon si on itere sur navsto dans TRAVA
!             sinon                  dans TRAV

  if(iterns.eq.1) then

!       On utilise temporairement TRAV comme tableau de travail.
!       Son contenu est stocke dans W7, W8 et W9 jusqu'apres tsepdc
    do iel = 1,ncel
      w7(iel) = trav(iel,1)
      w8(iel) = trav(iel,2)
      w9(iel) = trav(iel,3)
      trav(iel,1) = 0.d0
      trav(iel,2) = 0.d0
      trav(iel,3) = 0.d0
    enddo

    idiaex = 2
    call tsepdc                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   ncepdp ,                                                       &
   nideve , nrdeve , nituse , nrtuse , iphas  , idiaex ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc ,                                                       &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , trav   ,                            &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )


!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes extradiagonaux et
!       TRAV retrouve sa valeur
    if(isno2t(iphas).gt.0) then
      do iel = 1, ncel
        propce(iel,iptsna  ) = propce(iel,iptsna  ) + trav(iel,1)
        propce(iel,iptsna+1) = propce(iel,iptsna+1) + trav(iel,2)
        propce(iel,iptsna+2) = propce(iel,iptsna+2) + trav(iel,3)
        trav(iel,1)  = w7(iel)
        trav(iel,2)  = w8(iel)
        trav(iel,3)  = w9(iel)
      enddo

!     Si on n'extrapole pas les termes source en temps :
!       Si on itere sur navsto
!         TRAVA recoit les termes extradiagonaux et
!         TRAV retrouve sa valeur
    else
      if(nterup.gt.1) then
        do iel = 1, ncel
          trava(iel,1,iphas) = trava(iel,1,iphas) + trav(iel,1)
          trava(iel,2,iphas) = trava(iel,2,iphas) + trav(iel,2)
          trava(iel,3,iphas) = trava(iel,3,iphas) + trav(iel,3)
          trav(iel,1)  = w7(iel)
          trav(iel,2)  = w8(iel)
          trav(iel,3)  = w9(iel)
        enddo
!       Si on n'itere pas sur navsto
!         TRAV retrouve sa valeur augmentee des termes extradiagonaux
      else
        do iel = 1, ncel
          trav(iel,1)  = w7(iel) + trav(iel,1)
          trav(iel,2)  = w8(iel) + trav(iel,2)
          trav(iel,3)  = w9(iel) + trav(iel,3)
        enddo
      endif
    endif

  endif

endif


! ---> TERMES DE CORIOLIS
!     SI IPHYDR=1 LE TERME A DEJA ETE PRIS EN COMPTE AVANT

if (icorio.eq.1.and.iphydr.eq.0) then

  ! A la premiere iter sur navsto, on ajoute la partie issue des
  ! termes explicites
  if (iterns.eq.1) then

    ! Si on extrapole les termes source en temps :
    if(isno2t(iphas).gt.0) then

      do iel = 1, ncel
        cx = omegay*rtpa(iel,iwiph) - omegaz*rtpa(iel,iviph)
        cy = omegaz*rtpa(iel,iuiph) - omegax*rtpa(iel,iwiph)
        cz = omegax*rtpa(iel,iviph) - omegay*rtpa(iel,iuiph)
        romvom = -2.d0*propce(iel,ipcrom)*volume(iel)
        propce(iel,iptsna  ) = propce(iel,iptsna  ) + romvom*cx
        propce(iel,iptsna+1) = propce(iel,iptsna+1) + romvom*cy
        propce(iel,iptsna+2) = propce(iel,iptsna+2) + romvom*cz
      enddo

    ! Si on n'extrapole pas les termes source en temps :
    else

      ! Si on n'itere pas sur navsto : TRAV
      if (nterup.eq.1) then

        do iel = 1, ncel
          cx = omegay*rtpa(iel,iwiph) - omegaz*rtpa(iel,iviph)
          cy = omegaz*rtpa(iel,iuiph) - omegax*rtpa(iel,iwiph)
          cz = omegax*rtpa(iel,iviph) - omegay*rtpa(iel,iuiph)
          romvom = -2.d0*propce(iel,ipcrom)*volume(iel)
          trav(iel,1) = trav(iel,1) + romvom*cx
          trav(iel,2) = trav(iel,2) + romvom*cy
          trav(iel,3) = trav(iel,3) + romvom*cz
        enddo

      ! Si on itere sur navsto : TRAVA
      else

        do iel = 1, ncel
          cx = omegay*rtpa(iel,iwiph) - omegaz*rtpa(iel,iviph)
          cy = omegaz*rtpa(iel,iuiph) - omegax*rtpa(iel,iwiph)
          cz = omegax*rtpa(iel,iviph) - omegay*rtpa(iel,iuiph)
          romvom = -2.d0*propce(iel,ipcrom)*volume(iel)
          trava(iel,1,iphas) = trava(iel,1,iphas) + romvom*cx
          trava(iel,2,iphas) = trava(iel,2,iphas) + romvom*cy
          trava(iel,3,iphas) = trava(iel,3,iphas) + romvom*cz
        enddo

      endif

    endif
  endif
endif


! ---> - DIVERGENCE DE RIJ

if(itytur(iphas).eq.3.and.iterns.eq.1) then

  do isou = 1, 3

    if(isou.eq.1) ivar = iuiph
    if(isou.eq.2) ivar = iviph
    if(isou.eq.3) ivar = iwiph

    call divrij                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse , isou   , ivar   , iphas  , &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  ,                                              &
   viscf  , viscb  ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   w7     , w8     , w9     , coefu  ,                            &
   rdevel , rtuser , ra     )

    init = 1
    call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,              &
                                   ifacel,ifabor,viscf,viscb,w1)

!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes de divergence
    if(isno2t(iphas).gt.0) then
      do iel = 1, ncel
        propce(iel,iptsna+isou-1 ) =                              &
        propce(iel,iptsna+isou-1 ) - w1(iel)
      enddo
!     Si on n'extrapole pas les termes source en temps :
    else
!       si on n'itere pas sur navsto : TRAV
      if(nterup.eq.1) then
        do iel = 1, ncel
          trav(iel,isou) = trav(iel,isou) - w1(iel)
        enddo
!       si on itere sur navsto       : TRAVA
      else
        do iel = 1, ncel
          trava(iel,isou,iphas) = trava(iel,isou,iphas) - w1(iel)
        enddo
      endif
    endif

  enddo

endif


! ---> "VITESSE" DE DIFFUSION FACETTE
!      SI ON FAIT AUTRE CHOSE QUE DU K EPS, IL FAUDRA LA METTRE
!        DANS LA BOUCLE

if( idiff(iuiph).ge. 1 ) then

! --- Si la vitesse doit etre diffusee, on calcule la viscosite
!       pour le second membre (selon Rij ou non)

  if (itytur(iphas).eq.3) then
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis)                                &
                            + idifft(iuiph)*propce(iel,ipcvst)
    enddo
  endif

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

!     Quand on n'est pas en Rij ou que irijnu = 0, les tableaux
!       VISCFI, VISCBI se trouvent remplis par la meme occasion
!       (ils sont confondus avec VISCF, VISCB)
!     En Rij avec irijnu = 1, on calcule la viscosite increment
!       de la matrice dans VISCFI, VISCBI

  if(itytur(iphas).eq.3.and.irijnu(iphas).eq.1) then
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis)                                &
                            + idifft(iuiph)*propce(iel,ipcvst)
    enddo
    call viscfa                                                   &
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
   viscfi , viscbi ,                                              &
   rdevel , rtuser , ra     )
  endif

else

! --- Si la vitesse n'a pas de diffusion, on annule la viscosite
!      (matrice et second membre sont dans le meme tableau,
!       sauf en Rij avec IRIJNU = 1)

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

  if(itytur(iphas).eq.3.and.irijnu(iphas).eq.1) then
    do ifac = 1, nfac
      viscfi(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscbi(ifac) = 0.d0
    enddo
  endif

endif


! 2.2  RESOLUTION IMPLICITE NON COUPLEE DES 3 COMPO. DE VITESSES
! ==============================================================


! ---> AU PREMIER APPEL,
!      MISE A ZERO DE L'ESTIMATEUR POUR LA VITESSE PREDITE
!      S'IL DOIT ETRE CALCULE

if (iappel.eq.1) then
  if(iescal(iespre,iphas).gt.0) then
    iesprp = ipproc(iestim(iespre,iphas))
    do iel = 1, ncel
      propce(iel,iesprp) =  0.d0
    enddo
  endif
endif

! ---> AU DEUXIEME APPEL,
!      MISE A ZERO DE L'ESTIMATEUR TOTAL POUR NAVIER-STOKES
!      (SI ON FAIT UN DEUXIEME APPEL, ALORS IL DOIT ETRE CALCULE)

if(iappel.eq.2) then
  iestop = ipproc(iestim(iestot,iphas))
  do iel = 1, ncel
    propce(iel,iestop) = 0.d0
  enddo
endif


! ---> BOUCLE SUR LES DIRECTIONS DE L'ESPACE (U, V, W)


! Remarque : On suppose que le couplage vitesse pression
!  n'est valable que pour une seule phase.

do isou = 1, 3

  if(isou.eq.1) then
    ivar = iuiph
    ippt = ipptx
  endif
  if(isou.eq.2) then
    ivar = iviph
    ippt = ippty
  endif
  if(isou.eq.3) then
    ivar = iwiph
    ippt = ipptz
  endif
  ipp  = ipprtp(ivar)

  iclvar = iclrtp(ivar,icoef)
  iclvaf = iclrtp(ivar,icoeff)


! ---> TERMES SOURCES UTILISATEURS

  do iel = 1, ncel
    w7    (iel) = 0.d0
    drtp  (iel) = 0.d0
  enddo

!     Le calcul des parties implicite et explicite des termes sources
!       utilisateurs est faite uniquement a la premiere iter sur navsto.
  if(iterns.eq.1) then

    if (imatis.eq.1) then

      call mttsns                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , iphas  ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w7     , drtp   ,                                              &
!        ------   ------
   dam    , xam    ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

    else

      maxelt = max(ncelet, nfac, nfabor)
      ils    = idebia
      ifinia = ils + maxelt
      CALL IASIZE('PREDUV',IFINIA)

      call ustsns                                                 &
      !==========
 ( ifinia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , iphas  ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel ,                            &
   w7     , drtp   ,                                              &
!        ------   ------
   dam    , xam    ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

    endif

    if (nbrcpl.gt.0) then
      call csccel                                                 &
      !==========
 ( ifinia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ivar   , iphas  ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   w7     , drtp   ,                                              &
!        ------   ------
   dam    , xam    ,                                              &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )
    endif

  endif

!     On conserve la partie implicite pour les autres iter sur navsto
  if(iterns.eq.1.and.nterup.gt.1) then
    do iel = 1, ncel
      ximpa(iel,isou,iphas) = drtp(iel)
    enddo
  endif

!     On ajoute a TRAV ou TRAVA la partie issue des termes implicites
!       en utilisant DRTP
!       La prise en compte de UVWK a partir de la seconde iteration
!       est faite directement dans codits.
!     En schema std en temps, on continue a mettre MAX(-DRTP,0) dans la matrice
!     Avec termes sources a l'ordre 2, on implicite DRTP quel que soit son signe
!       (si on le met dans la matrice ou non selon son signe, on risque de ne pas
!        avoir le meme traitement d'un pas de temps au suivant)
  if(iterns.eq.1) then
    if(nterup.gt.1) then
      do iel = 1, ncel
        trava(iel,isou,iphas) = trava(iel,isou,iphas)             &
               + drtp(iel)*rtpa(iel,ivar)
      enddo
    else
      do iel = 1, ncel
        trav(iel,isou) = trav(iel,isou)                           &
               + drtp(iel)*rtpa(iel,ivar)
      enddo
    endif
  endif

!     A la premiere iter sur navsto, on ajoute la partie issue des
!       termes explicites
  if(iterns.eq.1) then
!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes explicites
    if(isno2t(iphas).gt.0) then
      do iel = 1, ncel
        propce(iel,iptsna+isou-1 ) =                              &
        propce(iel,iptsna+isou-1 ) + w7(iel)
      enddo
!     Si on n'extrapole pas les termes source en temps :
    else
!       si on n'itere pas sur navsto : TRAV
      if(nterup.eq.1) then
        do iel = 1, ncel
          trav(iel,isou) = trav(iel,isou) + w7(iel)
        enddo
!       si on itere sur navsto : TRAVA
      else
        do iel = 1, ncel
          trava(iel,isou,iphas) =                                 &
          trava(iel,isou,iphas) + w7(iel)
        enddo
      endif
    endif
  endif


! ---> TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume

  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
                           ifacel,ifabor,flumas,flumab,w1)


!     On ajoute a TRAV ou TRAVA la partie issue des termes implicites
  if(iterns.eq.1) then
    if(nterup.gt.1) then
      do iel = 1, ncel
        trava(iel,isou,iphas) = trava(iel,isou,iphas)             &
               +iconv(ivar)*w1(iel)*rtpa(iel,ivar)
      enddo
    else
      do iel = 1, ncel
        trav(iel,isou) = trav(iel,isou)                           &
               +iconv(ivar)*w1(iel)*rtpa(iel,ivar)
      enddo
    endif
  endif

  if(iappel.eq.1) then
!     Extrapolation ou non, meme forme par coherence avec bilsc2
    do iel = 1, ncel
      rovsdt(iel) =                                               &
         istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel)     &
         -iconv(ivar)*w1(iel)*thetav(ivar)
    enddo

!     Le remplissage de ROVSDT est toujours indispensable,
!       meme si on peut se contenter de n'importe quoi pour IAPPEL=2.
  else
    do iel = 1, ncel
      rovsdt(iel) = 0.d0
    enddo
  endif

! ---> TERMES SOURCES UTILISATEUR

  if(iappel.eq.1) then
    if(isno2t(iphas).gt.0) then
      thetap = thetav(ivar)
      if(iterns.gt.1) then
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel) -ximpa(iel,isou,iphas)*thetap
        enddo
      else
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel) -drtp(iel)*thetap
        enddo
      endif
    else
      if(iterns.gt.1) then
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel)                               &
               + max(-ximpa(iel,isou,iphas),zero)
        enddo
      else
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel)                               &
               + max(-drtp(iel),zero)
        enddo
      endif
    endif
  endif


! ---> PERTES DE CHARGE

!  Au second appel, on n'a pas besoin de rovsdt
  if(iappel.eq.1) then
    if (ncepdp.gt.0) then
      if(isno2t(iphas).gt.0) then
        thetap = thetav(ivar)
      else
        thetap = 1.d0
      endif
      do ielpdc = 1, ncepdp
        iel = icepdc(ielpdc)
        rovsdt(iel) = rovsdt(iel) +                               &
        propce(iel,ipcrom)*volume(iel)*ckupdc(ielpdc,isou)*thetap
      enddo
    endif
  endif


! --->  TERMES DE SOURCE DE MASSE

  if (ncesmp.gt.0) then

!     On calcule les termes Gamma (uinj - u)
!       -Gamma u a la premiere iteration est mis dans
!          TRAV ou TRAVA selon qu'on itere ou non sur navsto
!       Gamma uinj a la premiere iteration est placee dans W1
!       ROVSDT a chaque iteration recoit Gamma
    if(nterup.eq.1) then
      call catsma                                                 &
      !==========
  ( ncelet , ncel , ncesmp , iterns , isno2t(iphas), thetav(ivar),&
    icetsm , itypsm(1,ivar) ,                                     &
    volume , rtpa(1,ivar) , smacel(1,ivar) ,smacel(1,ipr(iphas)) ,&
    trav(1,isou)        , rovsdt , w1 )
    else
      call catsma                                                 &
      !==========
  ( ncelet , ncel , ncesmp , iterns , isno2t(iphas), thetav(ivar),&
    icetsm , itypsm(1,ivar) ,                                     &
    volume , rtpa(1,ivar) , smacel(1,ivar) ,smacel(1,ipr(iphas)) ,&
    trava(1,isou,iphas) , rovsdt , w1 )
    endif

!     A la premiere iter sur navsto, on ajoute la partie Gamma uinj
    if(iterns.eq.1) then
!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes explicites
      if(isno2t(iphas).gt.0) then
        do iel = 1,ncel
          propce(iel,iptsna+isou-1 ) =                            &
          propce(iel,iptsna+isou-1 ) + w1(iel)
        enddo
!     Si on n'extrapole pas les termes source en temps :
      else
!       si on n'itere pas sur navsto : TRAV
        if(nterup.eq.1) then
          do iel = 1,ncel
            trav(iel,isou)  = trav(iel,isou) + w1(iel)
          enddo
!       si on itere sur navsto : TRAVA
        else
          do iel = 1,ncel
            trava(iel,isou,iphas) =                               &
            trava(iel,isou,iphas) + w1(iel)
          enddo
        endif
      endif
    endif

  endif


! ---> INITIALISATION DU SECOND MEMBRE

!     Si on extrapole les TS
  if(isno2t(iphas).gt.0) then
    thetp1 = 1.d0 + thets
!       Si on n'itere pas sur navsto : TRAVA n'existe pas
    if(nterup.eq.1) then
      do iel = 1, ncel
        smbr(iel) =  trav(iel,isou)                               &
             + thetp1*propce(iel,iptsna+isou-1)
      enddo
!       Si on   itere     sur navsto : tout existe
    else
      do iel = 1, ncel
        smbr(iel) =  trav(iel,isou) + trava(iel,isou,iphas)       &
             + thetp1*propce(iel,iptsna+isou-1)
      enddo
    endif
!     Si on n'extrapole pas les TS : PROPCE n'existe pas
  else
!       Si on n'itere pas sur navsto : TRAVA n'existe pas
    if(nterup.eq.1) then
      do iel = 1, ncel
        smbr(iel) =  trav(iel,isou)
      enddo
!       Si on   itere     sur navsto : TRAVA existe
    else
      do iel = 1, ncel
        smbr(iel) =  trav(iel,isou) + trava(iel,isou,iphas)
      enddo
    endif
  endif


! ---> LAGRANGIEN : COUPLAGE RETOUR

!     L'ordre 2 sur les termes issus du lagrangien necessiterait de
!       decomposer TSLAGR(IEL,ISOU) en partie implicite et
!       explicite, comme c'est fait dans ustsns.
!     Pour le moment, on n'y touche pas.
  if (iilagr.eq.2 .and. ltsdyn.eq.1 .and. iphas.eq.ilphas)  then

    do iel = 1, ncel
      smbr(iel)   = smbr(iel) + tslagr(iel,itsvx+isou-1)
    enddo
!  Au second appel, on n'a pas besoin de rovsdt
    if(iappel.eq.1) then
      do iel = 1, ncel
        rovsdt(iel) = rovsdt(iel) + max(-tslagr(iel,itsli),zero)
      enddo
    endif

  endif

! ---> VERSIONS ELECTRIQUES : Arc Electrique (Force de Laplace)
!     Pour le moment, pas d'ordre 2 en temps.

  if ( ippmod(ielarc) .ge. 1 ) then
    do iel = 1,ncel
      smbr(iel) = smbr(iel)                                       &
                 + volume(iel)*propce(iel,ipproc(ilapla(isou)))
    enddo
  endif


! ---> PARAMETRES POUR LA RESOLUTION DU SYSTEME OU LE CALCUL DE l'ESTIMATEUR

  iconvp = iconv (ivar)
  idiffp = idiff (ivar)
  ireslp = iresol(ivar)
  ndircp = ndircl(ivar)
  nitmap = nitmax(ivar)
!MO        IMRGRA
  nswrsp = nswrsm(ivar)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  imgrp  = imgr  (ivar)
  ncymxp = ncymax(ivar)
  nitmfp = nitmgf(ivar)
!MO        IPP
  iwarnp = iwarni(ivar)
  blencp = blencv(ivar)
  epsilp = epsilo(ivar)
  epsrsp = epsrsm(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)
  relaxp = relaxv(ivar)
  thetap = thetav(ivar)

  if(iappel.eq.1) then

    iescap = iescal(iespre,iphas)

! ---> FIN DE LA CONSTRUCTION ET DE LA RESOLUTION DU SYSTEME

    if(iterns.eq.1) then

!  Attention, dans le cas des estimateurs, DAM fournit l'estimateur
!     des vitesses predites
      call codits                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     flumas , flumab ,                            &
   viscfi , viscbi , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

    elseif(iterns.gt.1) then

      call codits                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa(1,ivar)    , uvwk(1,isou,iphas) ,                         &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     flumas , flumab ,                            &
   viscfi , viscbi , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

    endif

!     DANS LE CAS DE PERTES DE CHARGE, ON UTILISE LES TABLEAUX
!     TPUCOU POUR LA PHASE D'IMPLICITATION
!     Attention, il faut regarder s'il y a des pdc sur un proc quelconque,
!       pas uniquement en local.
    if((ncpdct(iphas).gt.0).and.(ipucou.eq.0)) then
      do iel = 1,ncel
        tpucou(iel,isou) = dt(iel)
      enddo
      do ielpdc = 1, ncepdp
        iel=icepdc(ielpdc)
        tpucou(iel,isou) = 1.d0/(                                 &
             1.d0/dt(iel)+ckupdc(ielpdc,isou))
      enddo
    endif

!     COUPLAGE P/U RENFORCE : CALCUL DU VECTEUR T, STOCKE DANS TPUCOU
!       ON PASSE DANS CODITS, EN NE FAISANT QU'UN SEUL SWEEP, ET EN
!       INITIALISANT TPUCOU A 0 POUR QUE LA PARTIE CV/DIFF AJOUTEE
!       PAR BILSC2 SOIT NULLE
!     NSWRSP = -1 INDIQUERA A CODITS QU'IL NE FAUT FAIRE QU'UN SWEEP
!       ET QU'IL FAUT METTRE INC A 0 (POUR OTER LES DIRICHLETS DANS LES
!       CL DES MATRICES POIDS)

    if (ipucou.eq.1) then
      nswrsp = -1
      do iel = 1, ncel
        smbr(iel) = volume(iel)
      enddo
      do iel = 1, ncelet
        tpucou(iel,isou) = 0.d0
      enddo
      iescap = 0

      call codits                                                 &
      !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ippt   , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   tpucou(1,isou)  , tpucou(1,isou)  ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     flumas , flumab ,                            &
   viscfi , viscbi , viscf  , viscb  ,                            &
   rovsdt , smbr   , tpucou(1,isou)  ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

      do iel = 1, ncelet
        tpucou(iel,isou) = propce(iel,ipcrom)*tpucou(iel,isou)
      enddo

    endif


! --->  ESTIMATEUR SUR LA VITESSE PREDITE : ON SOMME SUR LES COMPOSANTES

    if(iescal(iespre,iphas).gt.0) then
      iesprp = ipproc(iestim(iespre,iphas))
      do iel = 1, ncel
        propce(iel,iesprp) =  propce(iel,iesprp) + dam(iel)
      enddo
    endif

  elseif(iappel.eq.2) then

! ---> FIN DE LA CONSTRUCTION DE L'ESTIMATEUR
!        RESIDU SECOND MEMBRE(Un+1,Pn+1) + RHO*VOLUME*( Un+1 - Un )/DT

    inc = 1
    iccocg = 1
!     Pas de relaxation en stationnaire
    idtva0 = 0

    call bilsc2                                                   &
    !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idtva0 , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp ,                                              &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtp(1,ivar)     , rtp(1,ivar)     ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   flumas , flumab , viscf  , viscb  ,                            &
   smbr   ,                                                       &
   w1     , w2     , w3     , w4     , w5     , w6     ,          &
   rdevel , rtuser , ra     )

    iestop = ipproc(iestim(iestot,iphas))
    do iel = 1, ncel
      propce(iel,iestop) =                                        &
           propce(iel,iestop)+ (smbr(iel)/volume(iel))**2
    enddo

  endif

enddo


! --->  APRES LA BOUCLE SUR U, V, W,
!        FIN DU CALCUL DE LA NORME POUR RESOLP

if(iappel.eq.1.and.irnpnw.eq.1) then

!     Calcul de div(rho u*)

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(rtp(1,iuiph), rtp(1,iviph), rtp(1,iwiph))
    !==========
  endif

!       Pour gagner du temps, on ne reconstruit pas.
  init   = 1
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  nswrp  = 1
  imligp = imligr(iuiph )
  iwarnp = iwarni(ipriph)
  epsrgp = epsrgr(iuiph )
  climgp = climgr(iuiph )
  extrap = extrag(iuiph )

  imaspe = 1

  iismph = iisymp+nfabor*(iphas-1)

  call inimas                                                     &
  !==========
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   iuiph  , iviph  , iwiph  , imaspe , iphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrp  , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(iismph) ,               &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   rtp(1,iuiph) , rtp(1,iviph) , rtp(1,iwiph) ,                   &
   coefa(1,icliup), coefa(1,iclivp), coefa(1,icliwp),             &
   coefb(1,icliup), coefb(1,iclivp), coefb(1,icliwp),             &
   viscf  , viscb  ,                                              &
   w4     , w5     , w6     , w7     , w8     , w9     ,          &
   smbr   , drtp   , rovsdt , coefu  ,                            &
   rdevel , rtuser , ra     )

  init = 0
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
                                ifacel,ifabor,viscf,viscb,xnormp)

!     Calcul de la norme
!       RNORMP qui servira dans resolp
  isqrt = 1
  call prodsc(ncelet,ncel,isqrt,xnormp,xnormp,rnormp(iphas))

endif

! --->  APRES LA BOUCLE SUR U, V, W,
!        FIN DU CALCUL DES ESTIMATEURS ET IMPRESSION

if(iappel.eq.1) then

! --->  ESTIMATEUR SUR LA VITESSE PREDITE : ON PREND LA RACINE (NORME)
!         SANS OU AVEC VOLUME (ET DANS CE CAS C'EST LA NORME L2)

  if(iescal(iespre,iphas).gt.0) then
    iesprp = ipproc(iestim(iespre,iphas))
    if(iescal(iespre,iphas).eq.1) then
      do iel = 1, ncel
        propce(iel,iesprp) =  sqrt(propce(iel,iesprp)            )
      enddo
    elseif(iescal(iespre,iphas).eq.2) then
      do iel = 1, ncel
        propce(iel,iesprp) =  sqrt(propce(iel,iesprp)*volume(iel))
      enddo
    endif
  endif

! ---> IMPRESSION DE NORME

  if (iwarni(iuiph).ge.2) then
    rnorm = -1.d0
    do iel = 1, ncel
      vitnor =                                                    &
       sqrt(rtp(iel,iuiph)**2+rtp(iel,iviph)**2+rtp(iel,iwiph)**2)
      rnorm = max(rnorm,vitnor)
    enddo
    if (irangp.ge.0) call parmax (rnorm)
                               !==========
    write(nfecra,1100) iphas, rnorm
  endif

elseif (iappel.eq.2) then

! --->  ESTIMATEUR SUR NAVIER-STOKES TOTAL : ON PREND LA RACINE (NORME)
!         SANS OU AVEC VOLUME (ET DANS CE CAS C'EST LA NORME L2)

  iestop = ipproc(iestim(iestot,iphas))
  if(iescal(iestot,iphas).eq.1) then
    do iel = 1, ncel
      propce(iel,iestop) = sqrt(propce(iel,iestop)            )
    enddo
  elseif(iescal(iestot,iphas).eq.2) then
    do iel = 1, ncel
      propce(iel,iestop) = sqrt(propce(iel,iestop)*volume(iel))
    enddo
  endif

endif

!--------
! FORMATS
!--------
#if defined(_CS_LANG_FR)

 1100 format(/,                                                   &
 1X,'Phase ',I4,' : Vitesse maximale apres prediction ',E12.4)

#else

 1100 format(/,                                                   &
 1X,'Phase ',I4,' : Maximum velocity after prediction ',E12.4)

#endif

!----
! FIN
!----

return

end subroutine
