!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine preduv &
!================

 ( iappel ,                                                       &
   nvar   , nscal  , iterns , ncepdp , ncesmp ,                   &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   flumas , flumab ,                                              &
   tslagr , coefa  , coefb  ,                                     &
   ckupdc , smacel , frcxt  , grdphd ,                            &
   trava  , ximpa  , uvwk   , dfrcxt , tpucou , trav   ,          &
   viscf  , viscb  , viscfi , viscbi ,                            &
   drtp   , smbr   , rovsdt ,                                     &
   w1     , w7     , w8     , w9     , xnormp )

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
! iappel           ! e  ! <-- ! numero d'appel du sous programme               !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! iterns           ! e  ! <-- ! numero d'iteration sur navsto                  !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
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
! frcxt(ncelet,3)  ! tr ! <-- ! force exterieure generant la pression          !
!                  !    !     !  hydrostatique                                 !
! grdphd(ncelet,3) ! tr ! <-- ! hydrostatic pressure gradient                  !
! trava,ximpa      ! tr ! <-- ! tableau de travail pour couplage               !
! uvwk             ! tr ! <-- ! tableau de travail pour couplage u/p           !
!                  !    !     ! sert a stocker la vitesse de                   !
!                  !    !     ! l'iteration precedente                         !
!dfrcxt(ncelet,3)  ! tr ! --> ! variation de force exterieure                  !
!                  !    !     !  generant la pression hydrostatique            !
! tpucou           ! tr ! --> ! couplage vitesse pression                      !
! (ncelel,ndim)    !    !     !                                                !
! trav(ncelet,3    ! tr ! --> ! smb qui servira pour normalisation             !
!                  !    !     !  dans resolp                                   !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! viscfi(nfac)     ! tr ! --- ! idem viscf pour increments                     !
! viscbi(nfabor    ! tr ! --- ! idem viscb pour increments                     !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr  (ncelet    ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! w1...9(ncelet    ! tr ! --- ! tableau de travail                             !
! xnormp(ncelet    ! tr ! <-- ! tableau reel pour norme resolp                 !
! tslagr(ncelet    ! tr ! <-- ! terme de couplage retour du                    !
!   ntersl)        !    !     !   lagrangien                                   !
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
use dimens, only: ndimfb
use numvar
use pointe, only: forbr, porosi
use entsor
use cstphy
use cstnum
use optcal
use parall
use period
use lagpar
use lagran
use ppppar
use ppthch
use ppincl
use cplsat
use mesh

!===============================================================================

implicit none

! Arguments

integer          iappel
integer          nvar   , nscal  , iterns
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision flumas(nfac), flumab(nfabor)
double precision tslagr(ncelet,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision frcxt(ncelet,3), dfrcxt(ncelet,3)
double precision grdphd(ncelet,3)
double precision trava(ncelet,ndim)
double precision ximpa(ncelet,ndim),uvwk(ncelet,ndim)
double precision tpucou(ncelet,ndim), trav(ncelet,3)
double precision viscf(nfac), viscb(nfabor)
double precision viscfi(nfac), viscbi(nfabor)
double precision drtp(ncelet)
double precision smbr(ncelet), rovsdt(ncelet)
double precision w1(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision xnormp(ncelet)

! Local variables

integer          iel   , ielpdc, ifac  , ivar  , isou
integer          iccocg, inc   , init  , ii    , isqrt
integer          ireslp, nswrgp, imligp, iwarnp, ippt  , ipp
integer          iclipr, iclu  , iclv  , iclw
integer          iclik , iclvar, iclvaf
integer          ipcrom, ipcroa, ipcroo, ipcrho, ipcvis, ipcvst
integer          iconvp, idiffp, ndircp, nitmap, nswrsp
integer          ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          iesprp, iestop
integer          iptsna
integer          iflmb0, nswrp , imaspe, ipbrom, itypfl
integer          idiaex, idtva0
integer          imucpp, idftnp, iswdyp

double precision rnorm , vitnor
double precision romvom, rom, drom
double precision epsrgp, climgp, extrap, relaxp, blencp, epsilp
double precision epsrsp
double precision vit1  , vit2  , vit3, xkb, pip, pfac, pfac1
double precision cpdc11, cpdc22, cpdc33, cpdc12, cpdc13, cpdc23
double precision d2s3  , thetap, thetp1, thets , dtsrom
double precision diipbx, diipby, diipbz
double precision cx    , cy    , cz

double precision rvoid(1)

double precision, allocatable, dimension(:) :: eswork
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: dpvar

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate a temporary array for the prediction-stage error estimator
if (iescal(iespre).gt.0) then
  allocate(eswork(ncelet))
endif
allocate(dpvar(ncelet))


iclipr = iclrtp(ipr,icoef)
iclu = iclrtp(iu ,icoef)
iclv = iclrtp(iv ,icoef)
iclw = iclrtp(iw ,icoef)

if(itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) then
  iclik  = iclrtp(ik ,icoef)
else
  iclik = 0
endif

!     Reperage de rho au bord
ipbrom = ipprob(irom  )
!     Reperage de rho courant (ie en cas d'extrapolation rho^n+1/2)
ipcrom = ipproc(irom  )
!     Reperage de rho^n en cas d'extrapolation
if(iroext.gt.0.or.idilat.gt.1) then
  ipcroa = ipproc(iroma)
else
  ipcroa = 0
endif

ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
!     Theta relatif aux termes sources explicites
thets  = thetsn
if(isno2t.gt.0) then
  iptsna = ipproc(itsnsa)
else
  iptsna = 0
endif

! Compute the porosity if needed
if (iterns.eq.1.and.iporos.eq.1) then
  call usporo
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
    drom = (propce(iel,ipcrom)-ro0)
    dfrcxt(iel,1) = drom*gx - frcxt(iel,1)
    dfrcxt(iel,2) = drom*gy - frcxt(iel,2)
    dfrcxt(iel,3) = drom*gz - frcxt(iel,3)
  enddo
!     Ajout eventuel des pertes de charges
  if (ncepdp.gt.0) then
    do ielpdc = 1, ncepdp
      iel=icepdc(ielpdc)
      vit1   = rtpa(iel,iu)
      vit2   = rtpa(iel,iv)
      vit3   = rtpa(iel,iw)
      cpdc11 = ckupdc(ielpdc,1)
      cpdc22 = ckupdc(ielpdc,2)
      cpdc33 = ckupdc(ielpdc,3)
      cpdc12 = ckupdc(ielpdc,4)
      cpdc23 = ckupdc(ielpdc,5)
      cpdc13 = ckupdc(ielpdc,6)
      dfrcxt(iel,1) = dfrcxt(iel,1)                   &
           -propce(iel,ipcrom)*(                                  &
           cpdc11*vit1+cpdc12*vit2+cpdc13*vit3)
      dfrcxt(iel,2) = dfrcxt(iel,2)                   &
           -propce(iel,ipcrom)*(                                  &
           cpdc12*vit1+cpdc22*vit2+cpdc23*vit3)
      dfrcxt(iel,3) = dfrcxt(iel,3)                   &
           -propce(iel,ipcrom)*(                                  &
           cpdc13*vit1+cpdc23*vit2+cpdc33*vit3)
    enddo
  endif
!     Ajout eventuel de la force de Coriolis
  if (icorio.eq.1) then
    do iel = 1, ncel
      cx = omegay*rtpa(iel,iw) - omegaz*rtpa(iel,iv)
      cy = omegaz*rtpa(iel,iu) - omegax*rtpa(iel,iw)
      cz = omegax*rtpa(iel,iv) - omegay*rtpa(iel,iu)
      dfrcxt(iel,1) = dfrcxt(iel,1) - 2.d0*propce(iel,ipcrom)*cx
      dfrcxt(iel,2) = dfrcxt(iel,2) - 2.d0*propce(iel,ipcrom)*cy
      dfrcxt(iel,3) = dfrcxt(iel,3) - 2.d0*propce(iel,ipcrom)*cz
    enddo
  endif

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(dfrcxt(1,1), dfrcxt(1,2), dfrcxt(1,3))
    !==========
  endif

endif



! ---> PRISE EN COMPTE DU GRADIENT DE PRESSION

! Allocate a work array for the gradient calculation
allocate(grad(ncelet,3))

iccocg = 1
inc    = 1
nswrgp = nswrgr(ipr)
imligp = imligr(ipr)
iwarnp = iwarni(ipr)
epsrgp = epsrgr(ipr)
climgp = climgr(ipr)
extrap = extrag(ipr)


call grdpot &
!==========
 ( ipr , imrgra , inc    , iccocg , nswrgp , imligp , iphydr ,    &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rvoid  ,                                                       &
   frcxt(1,1), frcxt(1,2), frcxt(1,3),                            &
   rtpa(1,ipr)  , coefa(1,iclipr) , coefb(1,iclipr) ,             &
   grad   )

! With porosity
if (iporos.eq.1) then
  do iel = 1, ncel
    grad(iel,1) = grad(iel,1)*porosi(iel)
    grad(iel,2) = grad(iel,2)*porosi(iel)
    grad(iel,3) = grad(iel,3)*porosi(iel)
  enddo
endif

!    Calcul des efforts aux parois (partie 2/5), si demande
!    La pression a la face est calculee comme dans gradrc/gradmc
!    et on la transforme en pression totale
!    On se limite a la premiere iteration (pour faire simple par
!      rapport a la partie issue de condli, hors boucle)
if (ineedf.eq.1 .and. iterns.eq.1) then
  do ifac = 1, nfabor
    iel = ifabor(ifac)
    diipbx = diipb(1,ifac)
    diipby = diipb(2,ifac)
    diipbz = diipb(3,ifac)
    pip = rtpa(iel,ipr) &
        + diipbx*grad(iel,1) + diipby*grad(iel,2) + diipbz*grad(iel,3)
    pfac = coefa(ifac,iclipr) +coefb(ifac,iclipr)*pip
    pfac1= rtpa(iel,ipr)                                          &
         +(cdgfbo(1,ifac)-xyzcen(1,iel))*grad(iel,1)              &
         +(cdgfbo(2,ifac)-xyzcen(2,iel))*grad(iel,2)              &
         +(cdgfbo(3,ifac)-xyzcen(3,iel))*grad(iel,3)
    pfac = coefb(ifac,iclipr)*(extrag(ipr)*pfac1                  &
         +(1.d0-extrag(ipr))*pfac)                                &
         +(1.d0-coefb(ifac,iclipr))*pfac                          &
         + ro0*(gx*(cdgfbo(1,ifac)-xyzp0(1))                      &
         + gy*(cdgfbo(2,ifac)-xyzp0(2))                           &
         + gz*(cdgfbo(3,ifac)-xyzp0(3)) )                         &
         - pred0
! on ne rajoute pas P0, pour garder un maximum de precision
!     &         + P0
    do isou = 1, 3
      forbr(isou,ifac) = forbr(isou,ifac) + pfac*surfbo(isou,ifac)
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
    trav(iel,1) = grad(iel,1)*dtsrom
    trav(iel,2) = grad(iel,2)*dtsrom
    trav(iel,3) = grad(iel,3)*dtsrom
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
  imligp = imligr(iu )
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu )
  climgp = climgr(iu )
  extrap = extrag(iu )

  imaspe = 1
  itypfl = 1

  call inimas &
  !==========
 ( nvar   , nscal  ,                                              &
   iu  , iv  , iw  , imaspe , itypfl ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrp  , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   trav(1,1) , trav(1,2) , trav(1,3) ,                            &
   coefa(1,iclu), coefa(1,iclv), coefa(1,iclw),                   &
   coefb(1,iclu), coefb(1,iclv), coefb(1,iclw),                   &
   viscf  , viscb  )

!     Calcul de div(rho dt/rho*grad P)
  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
                                ifacel,ifabor,viscf,viscb,xnormp)

!     Ajout de -Gamma
  if (ncesmp.gt.0) then
    do ii = 1, ncesmp
      iel = icetsm(ii)
      xnormp(iel) = xnormp(iel)-volume(iel)*smacel(ii,ipr)
    enddo
  endif
  ! Semi-analytic weakly compressible algorithm add + 1/rho Drho/Dt
  if (idilat.eq.4)then
    do iel = 1, ncel
      xnormp(iel) = xnormp(iel) + propce(iel,ipproc(iustdy(itsrho)))
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
      trav(iel,1) = (frcxt(iel,1) - grad(iel,1)) * volume(iel)
      trav(iel,2) = (frcxt(iel,2) - grad(iel,2)) * volume(iel)
      trav(iel,3) = (frcxt(iel,3) - grad(iel,3)) * volume(iel)
    enddo

  elseif(iphydr.eq.2) then
    do iel = 1, ncel
      rom = propce(iel,ipcrom)
      trav(iel,1) = (rom*gx - grdphd(iel,1) - grad(iel,1)) * volume(iel)
      trav(iel,2) = (rom*gy - grdphd(iel,2) - grad(iel,2)) * volume(iel)
      trav(iel,3) = (rom*gz - grdphd(iel,3) - grad(iel,3)) * volume(iel)
    enddo

  else
    do iel = 1, ncel
      drom = (propce(iel,ipcrom)-ro0)
      trav(iel,1) = ( drom*gx - grad(iel,1) ) * volume(iel)
      trav(iel,2) = ( drom*gy - grad(iel,2) ) * volume(iel)
      trav(iel,3) = ( drom*gz - grad(iel,3) ) * volume(iel)
    enddo
  endif

elseif(iappel.eq.2) then

  if (iphydr.eq.1) then
    do iel = 1, ncel
      trav(iel,1) = trav(iel,1) + (frcxt(iel,1) - grad(iel,1))*volume(iel)
      trav(iel,2) = trav(iel,2) + (frcxt(iel,2) - grad(iel,2))*volume(iel)
      trav(iel,3) = trav(iel,3) + (frcxt(iel,3) - grad(iel,3))*volume(iel)
    enddo

  elseif (iphydr.eq.2) then
    do iel = 1, ncel
      rom = propce(iel,ipcrom)
      trav(iel,1) = trav(iel,1) + (rom*gx - grdphd(iel,1) - grad(iel,1))*volume(iel)
      trav(iel,2) = trav(iel,2) + (rom*gy - grdphd(iel,2) - grad(iel,2))*volume(iel)
      trav(iel,3) = trav(iel,3) + (rom*gz - grdphd(iel,3) - grad(iel,3))*volume(iel)
    enddo

  else
    do iel = 1, ncel
      drom = (propce(iel,ipcrom)-ro0)
      trav(iel,1) = trav(iel,1) + (drom*gx - grad(iel,1))*volume(iel)
      trav(iel,2) = trav(iel,2) + (drom*gy - grad(iel,2))*volume(iel)
      trav(iel,3) = trav(iel,3) + (drom*gz - grad(iel,3))*volume(iel)
    enddo
  endif

endif

! Free memory
deallocate(grad)


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
    if(isno2t.gt.0) then
!         S'il n'y a qu'une    iter : TRAV  incremente
      if(nterup.eq.1) then
        do ii = 1, ndim
          do iel = 1, ncel
            trav (iel,ii) = trav (iel,ii) - thets*propce(iel,iptsna+ii-1)
          enddo
        enddo
!         S'il   y a plusieurs iter : TRAVA initialise
      else
        do ii = 1, ndim
          do iel = 1, ncel
            trava(iel,ii) = - thets*propce(iel,iptsna+ii-1)
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
            trava(iel,ii)  = 0.d0
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
if(     (itytur.eq.2 .or. itytur.eq.5 .or. iturb.eq.60) &
   .and. igrhok.eq.1 .and. iterns.eq.1) then

  ! Allocate a work array for the gradient calculation
  allocate(grad(ncelet,3))

  iccocg = 1
  inc    = 1
  nswrgp = nswrgr(ik)
  imligp = imligr(ik)
  epsrgp = epsrgr(ik)
  climgp = climgr(ik)
  extrap = extrag(ik)

  iwarnp = iwarni(iu)

  call grdcel &
  !==========
 ( ik  , imrgra , inc    , iccocg , nswrgp , imligp ,             &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,ik)   , coefa(1,iclik)  , coefb(1,iclik)  ,             &
   grad   )

  ! With porosity
  if (iporos.eq.1) then
    do iel = 1, ncel
      grad(iel,1) = grad(iel,1)*porosi(iel)
      grad(iel,2) = grad(iel,2)*porosi(iel)
      grad(iel,3) = grad(iel,3)*porosi(iel)
    enddo
  endif

  d2s3 = 2.d0/3.d0

!     Si on extrapole les termes source en temps : PROPCE
  if(isno2t.gt.0) then
!       Calcul de rho^n grad k^n      si rho non extrapole
!                 rho^n grad k^n si rho     extrapole
    ipcroo = ipcrom
    if(iroext.gt.0) ipcroo = ipcroa
    do iel = 1, ncel
      romvom = -propce(iel,ipcroo)*volume(iel)*d2s3
      propce(iel,iptsna  )=propce(iel,iptsna  )+grad(iel,1)*romvom
      propce(iel,iptsna+1)=propce(iel,iptsna+1)+grad(iel,2)*romvom
      propce(iel,iptsna+2)=propce(iel,iptsna+2)+grad(iel,3)*romvom
    enddo
!     Si on n'extrapole pas les termes sources en temps : TRAV ou TRAVA
  else
    if(nterup.eq.1) then
      do iel = 1, ncel
        romvom = -propce(iel,ipcrom)*volume(iel)*d2s3
        trav (iel,1) = trav (iel,1) + grad(iel,1) * romvom
        trav (iel,2) = trav (iel,2) + grad(iel,2) * romvom
        trav (iel,3) = trav (iel,3) + grad(iel,3) * romvom
      enddo
    else
      do iel = 1, ncel
        romvom = -propce(iel,ipcrom)*volume(iel)*d2s3
        trava(iel,1) = trava(iel,1) + grad(iel,1) * romvom
        trava(iel,2) = trava(iel,2) + grad(iel,2) * romvom
        trava(iel,3) = trava(iel,3) + grad(iel,3) * romvom
      enddo
    endif
  endif

!    Calcul des efforts aux parois (partie 3/5), si demande
  if (ineedf.eq.1) then
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      diipbx = diipb(1,ifac)
      diipby = diipb(2,ifac)
      diipbz = diipb(3,ifac)
      xkb = rtpa(iel,ik) + diipbx*grad(iel,1)                      &
           + diipby*grad(iel,2) + diipbz*grad(iel,3)
      xkb = coefa(ifac,iclik)+coefb(ifac,iclik)*xkb
      xkb = d2s3*propce(iel,ipcrom)*xkb
      do isou = 1, 3
        forbr(isou,ifac) = forbr(isou,ifac) + xkb*surfbo(isou,ifac)
      enddo
    enddo
  endif

  ! Free memory
  deallocate(grad)

endif


! ---> TERMES DE GRADIENT TRANSPOSE

!     Ce terme explicite est calcule une seule fois,
!       a la premiere iter sur navsto :
!       si on extrapole il va dans PROPCE,
!       sinon si on itere sur navsto dans TRAVA
!             sinon                  dans TRAV
if (ivisse.eq.1.and.iterns.eq.1) then

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
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   trav   ,                                                       &
!  ------
   viscf  , viscb  )

!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes de gradient transpose et
!       TRAV retrouve sa valeur
  if(isno2t.gt.0) then
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
        trava(iel,1) = trava(iel,1) + trav(iel,1)
        trava(iel,2) = trava(iel,2) + trav(iel,2)
        trava(iel,3) = trava(iel,3) + trav(iel,3)
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
 ( nvar   , nscal  ,                                              &
   ncepdp ,                                                       &
   idiaex ,                                                       &
   icepdc ,                                                       &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , trav   )

    ! With porosity
    if (iporos.eq.1) then
      do iel = 1, ncel
        trav(iel,1) = trav(iel,1)*porosi(iel)
        trav(iel,2) = trav(iel,2)*porosi(iel)
        trav(iel,3) = trav(iel,3)*porosi(iel)
      enddo
    endif

!     Si on itere sur navsto, on utilise TRAVA ; sinon TRAV
    if(nterup.gt.1) then
      do iel = 1, ncel
        trava(iel,1) = trava(iel,1) + trav(iel,1)
        trava(iel,2) = trava(iel,2) + trav(iel,2)
        trava(iel,3) = trava(iel,3) + trav(iel,3)
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
 ( nvar   , nscal  ,                                              &
   ncepdp ,                                                       &
   idiaex ,                                                       &
   icepdc ,                                                       &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , trav   )

    ! With porosity
    if (iporos.eq.1) then
      do iel = 1, ncel
        trav(iel,1) = trav(iel,1)*porosi(iel)
        trav(iel,2) = trav(iel,2)*porosi(iel)
        trav(iel,3) = trav(iel,3)*porosi(iel)
      enddo
    endif

!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes extradiagonaux et
!       TRAV retrouve sa valeur
    if(isno2t.gt.0) then
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
          trava(iel,1) = trava(iel,1) + trav(iel,1)
          trava(iel,2) = trava(iel,2) + trav(iel,2)
          trava(iel,3) = trava(iel,3) + trav(iel,3)
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


! ---> TERMES DE CORIOLIS ! FIXME with porosity
!     SI IPHYDR=1 LE TERME A DEJA ETE PRIS EN COMPTE AVANT

if (icorio.eq.1.and.iphydr.eq.0) then

  ! A la premiere iter sur navsto, on ajoute la partie issue des
  ! termes explicites
  if (iterns.eq.1) then

    ! Si on extrapole les termes source en temps :
    if(isno2t.gt.0) then

      do iel = 1, ncel
        cx = omegay*rtpa(iel,iw) - omegaz*rtpa(iel,iv)
        cy = omegaz*rtpa(iel,iu) - omegax*rtpa(iel,iw)
        cz = omegax*rtpa(iel,iv) - omegay*rtpa(iel,iu)
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
          cx = omegay*rtpa(iel,iw) - omegaz*rtpa(iel,iv)
          cy = omegaz*rtpa(iel,iu) - omegax*rtpa(iel,iw)
          cz = omegax*rtpa(iel,iv) - omegay*rtpa(iel,iu)
          romvom = -2.d0*propce(iel,ipcrom)*volume(iel)
          trav(iel,1) = trav(iel,1) + romvom*cx
          trav(iel,2) = trav(iel,2) + romvom*cy
          trav(iel,3) = trav(iel,3) + romvom*cz
        enddo

      ! Si on itere sur navsto : TRAVA
      else

        do iel = 1, ncel
          cx = omegay*rtpa(iel,iw) - omegaz*rtpa(iel,iv)
          cy = omegaz*rtpa(iel,iu) - omegax*rtpa(iel,iw)
          cz = omegax*rtpa(iel,iv) - omegay*rtpa(iel,iu)
          romvom = -2.d0*propce(iel,ipcrom)*volume(iel)
          trava(iel,1) = trava(iel,1) + romvom*cx
          trava(iel,2) = trava(iel,2) + romvom*cy
          trava(iel,3) = trava(iel,3) + romvom*cz
        enddo

      endif

    endif
  endif
endif


! ---> - DIVERGENCE DE RIJ

if(itytur.eq.3.and.iterns.eq.1) then

  do isou = 1, 3

    if(isou.eq.1) ivar = iu
    if(isou.eq.2) ivar = iv
    if(isou.eq.3) ivar = iw

    call divrij                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   isou   , ivar   ,                                              &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  ,                                              &
   viscf  , viscb  )

    init = 1
    call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,              &
                                   ifacel,ifabor,viscf,viscb,w1)

!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes de divergence
    if(isno2t.gt.0) then
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
          trava(iel,isou) = trava(iel,isou) - w1(iel)
        enddo
      endif
    endif

  enddo

endif


! ---> "VITESSE" DE DIFFUSION FACETTE
!      SI ON FAIT AUTRE CHOSE QUE DU K EPS, IL FAUDRA LA METTRE
!        DANS LA BOUCLE

if( idiff(iu).ge. 1 ) then

! --- Si la vitesse doit etre diffusee, on calcule la viscosite
!       pour le second membre (selon Rij ou non)

  if (itytur.eq.3) then
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis) + idifft(iu)*propce(iel,ipcvst)
    enddo
  endif

  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

!     Quand on n'est pas en Rij ou que irijnu = 0, les tableaux
!       VISCFI, VISCBI se trouvent remplis par la meme occasion
!       (ils sont confondus avec VISCF, VISCB)
!     En Rij avec irijnu = 1, on calcule la viscosite increment
!       de la matrice dans VISCFI, VISCBI

  if(itytur.eq.3.and.irijnu.eq.1) then
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis)                                &
                            + idifft(iu)*propce(iel,ipcvst)
    enddo
    call viscfa                                                   &
    !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscfi , viscbi )
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

  if(itytur.eq.3.and.irijnu.eq.1) then
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
  if(iescal(iespre).gt.0) then
    iesprp = ipproc(iestim(iespre))
    do iel = 1, ncel
      propce(iel,iesprp) =  0.d0
    enddo
  endif
endif

! ---> AU DEUXIEME APPEL,
!      MISE A ZERO DE L'ESTIMATEUR TOTAL POUR NAVIER-STOKES
!      (SI ON FAIT UN DEUXIEME APPEL, ALORS IL DOIT ETRE CALCULE)

if(iappel.eq.2) then
  iestop = ipproc(iestim(iestot))
  do iel = 1, ncel
    propce(iel,iestop) = 0.d0
  enddo
endif


! ---> BOUCLE SUR LES DIRECTIONS DE L'ESPACE (U, V, W)


! Remarque : On suppose que le couplage vitesse pression
!  n'est valable que pour une seule phase.

do isou = 1, 3

  if(isou.eq.1) then
    ivar = iu
    ippt = ipptx
  endif
  if(isou.eq.2) then
    ivar = iv
    ippt = ippty
  endif
  if(isou.eq.3) then
    ivar = iw
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
! FIXME with porosity
  if(iterns.eq.1) then

    call ustsns                                                    &
    !==========
  ( nvar   , nscal  , ncepdp , ncesmp ,                            &
    ivar   ,                                                       &
    icepdc , icetsm , itypsm ,                                     &
    dt     , rtpa   , propce , propfa , propfb ,                   &
    ckupdc , smacel ,                                              &
    w7     , drtp   )

    if (nbrcpl.gt.0) then
      call csccel                                                 &
      !==========
 ( nvar   , nscal  ,                                              &
   ivar   ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  ,                                              &
   w7     , drtp   )
    endif

  endif

!     On conserve la partie implicite pour les autres iter sur navsto
  if(iterns.eq.1.and.nterup.gt.1) then
    do iel = 1, ncel
      ximpa(iel,isou) = drtp(iel)
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
        trava(iel,isou) = trava(iel,isou) + drtp(iel)*rtpa(iel,ivar)
      enddo
    else
      do iel = 1, ncel
        trav(iel,isou) = trav(iel,isou) + drtp(iel)*rtpa(iel,ivar)
      enddo
    endif
  endif

!     A la premiere iter sur navsto, on ajoute la partie issue des
!       termes explicites
  if(iterns.eq.1) then
!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes explicites
    if(isno2t.gt.0) then
      do iel = 1, ncel
        propce(iel,iptsna+isou-1 ) = propce(iel,iptsna+isou-1 ) + w7(iel)
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
          trava(iel,isou) = trava(iel,isou) + w7(iel)
        enddo
      endif
    endif
  endif


! ---> TERME D'ACCUMULATION DE MASSE -(dRO/dt)*Volume

  if(iappel.eq.1) then

    ! Low Mach compressible Algos
    if (idilat.gt.1) then
      ipcrho = ipcroa

    ! Standard algo
    else
      ipcrho = ipcrom
    endif

    do iel = 1, ncel
      rovsdt(iel) = istat(ivar)*propce(iel,ipcrho)/dt(iel)*volume(iel)
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
    if(isno2t.gt.0) then
      thetap = thetav(ivar)
      if(iterns.gt.1) then
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel) -ximpa(iel,isou)*thetap
        enddo
      else
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel) -drtp(iel)*thetap
        enddo
      endif
    else
      if(iterns.gt.1) then
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel) + max(-ximpa(iel,isou),zero)
        enddo
      else
        do iel = 1, ncel
          rovsdt(iel) = rovsdt(iel) + max(-drtp(iel),zero)
        enddo
      endif
    endif
  endif


! ---> PERTES DE CHARGE

!  Au second appel, on n'a pas besoin de rovsdt
  if(iappel.eq.1) then
    if (ncepdp.gt.0) then
      if(isno2t.gt.0) then
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
  ( ncelet , ncel , ncesmp , iterns , isno2t, thetav(ivar),&
    icetsm , itypsm(1,ivar) ,                                     &
    volume , rtpa(1,ivar) , smacel(1,ivar) ,smacel(1,ipr) ,&
    trav(1,isou)        , rovsdt , w1 )
    else
      call catsma                                                 &
      !==========
  ( ncelet , ncel , ncesmp , iterns , isno2t, thetav(ivar),&
    icetsm , itypsm(1,ivar) ,                                     &
    volume , rtpa(1,ivar) , smacel(1,ivar) ,smacel(1,ipr) ,&
    trava(1,isou) , rovsdt , w1 )
    endif

!     A la premiere iter sur navsto, on ajoute la partie Gamma uinj
    if(iterns.eq.1) then
!     Si on extrapole les termes source en temps :
!       PROPCE recoit les termes explicites
      if(isno2t.gt.0) then
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
            trava(iel,isou) =                               &
            trava(iel,isou) + w1(iel)
          enddo
        endif
      endif
    endif

  endif


! ---> INITIALISATION DU SECOND MEMBRE

!     Si on extrapole les TS
  if(isno2t.gt.0) then
    thetp1 = 1.d0 + thets
!       Si on n'itere pas sur navsto : TRAVA n'existe pas
    if(nterup.eq.1) then
      do iel = 1, ncel
        smbr(iel) =  trav(iel,isou) + thetp1*propce(iel,iptsna+isou-1)
      enddo
!       Si on   itere     sur navsto : tout existe
    else
      do iel = 1, ncel
        smbr(iel) =  trav(iel,isou) + trava(iel,isou)       &
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
        smbr(iel) =  trav(iel,isou) + trava(iel,isou)
      enddo
    endif
  endif


! ---> LAGRANGIEN : COUPLAGE RETOUR

!     L'ordre 2 sur les termes issus du lagrangien necessiterait de
!       decomposer TSLAGR(IEL,ISOU) en partie implicite et
!       explicite, comme c'est fait dans ustsns.
!     Pour le moment, on n'y touche pas.
  if (iilagr.eq.2 .and. ltsdyn.eq.1)  then

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
      smbr(iel) = smbr(iel) + volume(iel)*propce(iel,ipproc(ilapla(isou)))
    enddo
  endif

  ! With porosity
  if (iporos.eq.1) then
    do iel = 1, ncel
      rovsdt(iel) = rovsdt(iel)*porosi(iel)
    enddo
  endif

! ---> PARAMETRES POUR LA RESOLUTION DU SYSTEME OU LE CALCUL DE l'ESTIMATEUR

  iconvp = iconv (ivar)
  idiffp = idiff (ivar)
  ireslp = iresol(ivar)
  ndircp = ndircl(ivar)
  nitmap = nitmax(ivar)
  nswrsp = nswrsm(ivar)
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  ircflp = ircflu(ivar)
  ischcp = ischcv(ivar)
  isstpp = isstpc(ivar)
  imucpp = 0
  idftnp = 1 !no tensorial diffusion avalaible for the velocity
  iswdyp = iswdyn(ivar)
  imgrp  = imgr  (ivar)
  ncymxp = ncymax(ivar)
  nitmfp = nitmgf(ivar)
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

    iescap = iescal(iespre)

! ---> FIN DE LA CONSTRUCTION ET DE LA RESOLUTION DU SYSTEME

    if(iterns.eq.1) then

!  Attention, dans le cas des estimateurs, eswork fournit l'estimateur
!     des vitesses predites
      call codits &
      !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   flumas , flumab ,                                              &
   viscfi , viscbi , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   rovsdt , smbr   , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , eswork )

    elseif(iterns.gt.1) then

      call codits &
      !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   rtpa(1,ivar)    , uvwk(1,isou) ,                               &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   flumas , flumab ,                                              &
   viscfi , viscbi , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   rovsdt , smbr   , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , eswork )

    endif

!     DANS LE CAS DE PERTES DE CHARGE, ON UTILISE LES TABLEAUX
!     TPUCOU POUR LA PHASE D'IMPLICITATION
!     Attention, il faut regarder s'il y a des pdc sur un proc quelconque,
!       pas uniquement en local.
    if((ncpdct.gt.0).and.(ipucou.eq.0)) then
      do iel = 1,ncel
        tpucou(iel,isou) = dt(iel)
      enddo
      do ielpdc = 1, ncepdp
        iel=icepdc(ielpdc)
        tpucou(iel,isou) = 1.d0 / (1.d0/dt(iel) + ckupdc(ielpdc,isou))
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

      call codits &
      !==========
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ippt   , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetap ,                                              &
   tpucou(1,isou)  , tpucou(1,isou)  ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   flumas , flumab ,                                              &
   viscfi , viscbi , rvoid  , viscf  , viscb  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   rovsdt , smbr   , tpucou(1,isou)  , dpvar  ,                   &
   rvoid  , rvoid  )

      do iel = 1, ncelet
        tpucou(iel,isou) = propce(iel,ipcrom)*tpucou(iel,isou)
      enddo

    endif


! --->  ESTIMATEUR SUR LA VITESSE PREDITE : ON SOMME SUR LES COMPOSANTES

    if(iescal(iespre).gt.0) then
      iesprp = ipproc(iestim(iespre))
      do iel = 1, ncel
        propce(iel,iesprp) =  propce(iel,iesprp) + eswork(iel)
      enddo
    endif

  elseif(iappel.eq.2) then

! ---> FIN DE LA CONSTRUCTION DE L'ESTIMATEUR
!        RESIDU SECOND MEMBRE(Un+1,Pn+1) + RHO*VOLUME*( Un+1 - Un )/DT

    inc = 1
    iccocg = 1
!     Pas de relaxation en stationnaire
    idtva0 = 0

    call bilsca &
    !==========
 ( nvar   , nscal  ,                                              &
   idtva0 , ivar   , iconvp , idiffp , nswrgp , imligp , ircflp , &
   ischcp , isstpp , inc    , imrgra , iccocg ,                   &
   ipp    , iwarnp , imucpp , idftnp ,                            &
   blencp , epsrgp , climgp , extrap , relaxp , thetap ,          &
   rtp(1,ivar)     , rtp(1,ivar)     ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   flumas , flumab , viscf  , viscb  , rvoid  , rvoid  ,          &
   rvoid  , rvoid  ,                                              &
   smbr   )

    iestop = ipproc(iestim(iestot))
    do iel = 1, ncel
      propce(iel,iestop) = propce(iel,iestop)+ (smbr(iel)/volume(iel))**2
    enddo

  endif

enddo


! --->  APRES LA BOUCLE SUR U, V, W,
!        FIN DU CALCUL DE LA NORME POUR RESOLP

if(iappel.eq.1.and.irnpnw.eq.1) then

!     Calcul de div(rho u*)

  if (irangp.ge.0.or.iperio.eq.1) then
    call synvec(rtp(1,iu), rtp(1,iv), rtp(1,iw))
    !==========
  endif

!       Pour gagner du temps, on ne reconstruit pas.
  init   = 1
  inc    = 1
  iccocg = 1
  iflmb0 = 1
  nswrp  = 1
  imligp = imligr(iu )
  iwarnp = iwarni(ipr)
  epsrgp = epsrgr(iu )
  climgp = climgr(iu )
  extrap = extrag(iu )

  imaspe = 1
  itypfl = 1

  call inimas &
  !==========
 ( nvar   , nscal  ,                                              &
   iu  , iv  , iw  , imaspe , itypfl ,                            &
   iflmb0 , init   , inc    , imrgra , iccocg , nswrp  , imligp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   propce(1,ipcrom), propfb(1,ipbrom),                            &
   rtp(1,iu) , rtp(1,iv) , rtp(1,iw) ,                            &
   coefa(1,iclu), coefa(1,iclv), coefa(1,iclw),                   &
   coefb(1,iclu), coefb(1,iclv), coefb(1,iclw),                   &
   viscf  , viscb  )

  init = 0
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
                                ifacel,ifabor,viscf,viscb,xnormp)

!     Calcul de la norme
!       RNORMP qui servira dans resolp
  isqrt = 1
  call prodsc(ncel,isqrt,xnormp,xnormp,rnormp)

endif

! --->  APRES LA BOUCLE SUR U, V, W,
!        FIN DU CALCUL DES ESTIMATEURS ET IMPRESSION

if(iappel.eq.1) then

! --->  ESTIMATEUR SUR LA VITESSE PREDITE : ON PREND LA RACINE (NORME)
!         SANS OU AVEC VOLUME (ET DANS CE CAS C'EST LA NORME L2)

  if(iescal(iespre).gt.0) then
    iesprp = ipproc(iestim(iespre))
    if(iescal(iespre).eq.1) then
      do iel = 1, ncel
        propce(iel,iesprp) =  sqrt(propce(iel,iesprp)            )
      enddo
    elseif(iescal(iespre).eq.2) then
      do iel = 1, ncel
        propce(iel,iesprp) =  sqrt(propce(iel,iesprp)*volume(iel))
      enddo
    endif
  endif

! ---> IMPRESSION DE NORME

  if (iwarni(iu).ge.2) then
    rnorm = -1.d0
    do iel = 1, ncel
      vitnor = sqrt(rtp(iel,iu)**2+rtp(iel,iv)**2+rtp(iel,iw)**2)
      rnorm = max(rnorm,vitnor)
    enddo
    if (irangp.ge.0) call parmax (rnorm)
                               !==========
    write(nfecra,1100) rnorm
  endif

elseif (iappel.eq.2) then

! --->  ESTIMATEUR SUR NAVIER-STOKES TOTAL : ON PREND LA RACINE (NORME)
!         SANS OU AVEC VOLUME (ET DANS CE CAS C'EST LA NORME L2)

  iestop = ipproc(iestim(iestot))
  if(iescal(iestot).eq.1) then
    do iel = 1, ncel
      propce(iel,iestop) = sqrt(propce(iel,iestop)            )
    enddo
  elseif(iescal(iestot).eq.2) then
    do iel = 1, ncel
      propce(iel,iestop) = sqrt(propce(iel,iestop)*volume(iel))
    enddo
  endif

endif

! Free memory
deallocate(dpvar)

!--------
! FORMATS
!--------
#if defined(_CS_LANG_FR)

 1100 format(/,                                                   &
 1X,'Vitesse maximale apres prediction ',E12.4)

#else

 1100 format(/,                                                   &
 1X,'Maximum velocity after prediction ',E12.4)

#endif

!----
! FIN
!----

return

end subroutine
