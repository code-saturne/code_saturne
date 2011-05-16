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

subroutine resv2f &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  , prdv2f ,                                     &
   dam    , xam    , drtp   , smbr   , rovsdt ,                   &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     , w10    ,          &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR PHI ET DE DIFFUSION POUR F_BARRE DANS LE CADRE DU
!   MODELE V2F PHI-MODEL

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! ivar             ! i  ! <-- ! variable number                                !
! isou             ! e  ! <-- ! numero de passage                              !
! ipp              ! e  ! <-- ! numero de variable pour sorties post           !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! ia(*)            ! ia ! --- ! main integer work array                        !
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
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! prdv2f(ncelet    ! tr ! <-- ! tableau de stockage du terme de                !
!                  !    !     ! prod de turbulence pour le v2f                 !
! dam(ncelet       ! tr ! --- ! tableau de travail pour matrice                !
! xam(nfac,*)      ! tr ! --- ! tableau de travail pour matrice                !
! drtp(ncelet      ! tr ! --- ! tableau de travail pour increment              !
! smbr(ncelet      ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
! w1..10(ncelet    ! tr ! --- ! tableau de travail                             !
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
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstnum
use cstphy
use pointe
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision prdv2f(ncelet)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbr(ncelet), rovsdt(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra, ifinia
integer          init  , ifac  , iel   , inc   , iccocg
integer          ivar, ipp
integer          iiun
integer          ipriph,iuiph
integer          ikiph, ieiph, iphiph, ifbiph
integer          iclvar, iclvaf
integer          iclikp, iclphi, iclfbp
integer          ipcrom, ipcroo, ipcvis, ipcvlo, ipcvst, ipcvso
integer          iflmas, iflmab
integer          nswrgp, imligp, iwarnp, iphydp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          iptsta
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision tuexpe, thets , thetv , thetap, thetp1
double precision d2s3, d1s4, d3s2
double precision xk, xe, xnu, xrom, ttke, ttmin, llke, llmin

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

ipriph = ipr
iuiph  = iu
ikiph  = ik
ieiph  = iep
iphiph = iphi
ifbiph = ifb

ipcrom = ipproc(irom  )
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
iflmas = ipprof(ifluma(iuiph))
iflmab = ipprob(ifluma(iuiph))

iclikp = iclrtp(ikiph,icoef)
iclphi = iclrtp(iphiph,icoef)
iclfbp = iclrtp(ifbiph,icoef)

if(isto2t.gt.0) then
  iptsta = ipproc(itstua)
else
  iptsta = 0
endif

d2s3 = 2.0d0/3.0d0
d1s4 = 1.0d0/4.0d0
d3s2 = 3.0d0/2.0d0

if(iwarni(iphiph).ge.1) then
  write(nfecra,1000)
endif

!===============================================================================
! 2. CALCUL DU TERME EN GRAD PHI.GRAD K
!===============================================================================

  iccocg = 1
  inc = 1
  ivar = iphiph

  nswrgp = nswrgr(ivar )
  imligp = imligr(ivar )
  iwarnp = iwarni(ivar )
  epsrgp = epsrgr(ivar )
  climgp = climgr(ivar )
  extrap = extrag(ivar )
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( iphiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iphiph ) , coefa(1,iclphi) , coefb(1,iclphi) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   ra     )

  iccocg = 1
  inc = 1
  ivar = ikiph

  nswrgp = nswrgr(ivar )
  imligp = imligr(ivar )
  iwarnp = iwarni(ivar )
  epsrgp = epsrgr(ivar )
  climgp = climgr(ivar )
  extrap = extrag(ivar )
  iphydp = 0

  call grdcel                                                     &
  !==========
 ( ikiph  , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ia     ,                                                       &
   w4     , w4     , w4     ,                                     &
   rtpa(1,ikiph )  , coefa(1,iclikp) , coefb(1,iclikp) ,          &
   w4     , w5     , w6     ,                                     &
!        ------   ------   ------
   w7     , w8     , w9     ,                                     &
   ra     )

  do iel = 1, ncel
    w1(iel) = w1(iel)*w4(iel) + w2(iel)*w5(iel) + w3(iel)*w6(iel)
  enddo

!===============================================================================
! 3. RESOLUTION DE L'EQUATION DE F_BARRE
!===============================================================================

ivar = ifbiph
iclvar = iclfbp
iclvaf = iclfbp
ipp    = ipprtp(ivar)

if(iwarni(ivar).ge.1) then
  write(nfecra,1100) nomvar(ipp)
endif

!     S pour Source, V pour Variable
thets  = thetst
thetv  = thetav(ivar )

ipcroo = ipcrom
ipcvlo = ipcvis
if(isto2t.gt.0) then
  if (iroext.gt.0) then
    ipcroo = ipproc(iroma)
  endif
  if(iviext.gt.0) then
    ipcvlo = ipproc(ivisla)
  endif
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 3.1 TERMES SOURCES  UTILISATEURS
!===============================================================================

call ustsv2                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , prdv2f , w1     ,          &
   smbr   , rovsdt ,                                              &
!        ------   ------
   viscf  , viscb  , xam    ,                                     &
   w2     , w3     , w4     , w5     , w6     , w7     ,          &
   w8     , w9     , w10    , dam    , drtp   ,                   &
   ra     )

!     Si on extrapole les T.S.
if(isto2t.gt.0) then
  do iel = 1, ncel
!       Sauvegarde pour echange
    tuexpe = propce(iel,iptsta+2)
!       Pour la suite et le pas de temps suivant
!       On met un signe "-" car on résout en fait "-div(grad fb) = ..."
    propce(iel,iptsta+2) = - smbr(iel)
!       Second membre du pas de temps precedent
!       on implicite le terme source utilisateur (le reste)
    smbr(iel) = - rovsdt(iel)*rtpa(iel,ivar) - thets*tuexpe
!       Diagonale
    rovsdt(iel) = thetv*rovsdt(iel)
  enddo
else
  do iel = 1, ncel
!       On met un signe "-" car on résout en fait "-div(grad fb) = ..."
!       On resout par gradient conjugue, donc on n'impose pas le signe
!          de ROVSDT
    smbr(iel)   = -rovsdt(iel)*rtpa(iel,ivar) - smbr(iel)
!          ROVSDT(IEL) =  ROVSDT(IEL)
  enddo
endif


!===============================================================================
! 3.2 TERME SOURCE DE F_BARRE
!     SMBR=1/L^2*(f_b + 1/T(C1-1)(phi-2/3) - C2*Pk/k/rho
!     -2*nu/k*grad_phi*grad_k -nu*div(grad(phi)) )
!  En fait on met un signe "-" car l'eq resolue est
!    -div(grad f_b) = SMBR
!===============================================================================

!     On calcule le terme en -VOLUME*div(grad(phi)) par itrgrp,
!     et on le stocke dans W2
!     Attention, les VISCF et VISCB calcules ici servent a ITRGRP mais
!     aussi a CODITS qui suit

do iel = 1, ncel
  w3(iel) = 1.d0
enddo
call viscfa                                                       &
!==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w3     ,                                                       &
   viscf  , viscb  ,                                              &
   ra     )


iccocg = 1
inc = 1
init = 1

nswrgp = nswrgr(iphiph)
imligp = imligr(iphiph)
iwarnp = iwarni(iphiph)
epsrgp = epsrgr(iphiph)
climgp = climgr(iphiph)
extrap = extrag(iphiph)
iphydp = 0

call itrgrp                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ia     ,                                                       &
   w2     , w2     , w2     ,                                     &
   rtpa(1,iphiph)   , coefa(1,iclphi) , coefb(1,iclphi) ,         &
   viscf  , viscb  ,                                              &
   w3     , w3     , w3     ,                                     &
   w2     ,                                                       &
!        --
   w4     , w5     , w6     , w7     , w8     , w9     ,          &
   ra     )

!      On stocke T dans W3 et L^2 dans W4
!      Dans le cas de l'ordre 2 en temps, T est calcule en n
!      (il sera extrapole) et L^2 en n+theta (meme si k et eps restent en n)
do iel=1,ncel
  xk = rtpa(iel,ikiph)
  xe = rtpa(iel,ieiph)
  xnu  = propce(iel,ipcvlo)/propce(iel,ipcroo)
  ttke = xk / xe
  ttmin = cv2fct*sqrt(xnu/xe)
  w3(iel) = max(ttke,ttmin)

  xnu  = propce(iel,ipcvis)/propce(iel,ipcrom)
  llke = xk**d3s2/xe
  llmin = cv2fet*(xnu**3/xe)**d1s4
  w4(iel) = ( cv2fcl*max(llke,llmin) )**2
enddo

!     Terme explicite, stocke temporairement dans W5
!     W2 est deja multiplie par le volume et contient deja
!     un signe "-" (issu de ITRGRP)
do iel = 1, ncel
    xrom = propce(iel,ipcroo)
    xnu  = propce(iel,ipcvlo)/xrom
    xk = rtpa(iel,ikiph)
    xe = rtpa(iel,ieiph)
    w5(iel) = - volume(iel)*                                      &
         ( (cv2fc1-1.d0)*(rtpa(iel,iphiph)-d2s3)/w3(iel)          &
           -cv2fc2*prdv2f(iel)/xrom/xk                            &
           -2.0d0*xnu/xe/w3(iel)*w1(iel) ) - xnu*w2(iel)
enddo
!     Si on extrapole les T.S : PROPCE
if(isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    propce(iel,iptsta+2) =                                   &
    propce(iel,iptsta+2) + w5(iel)
    smbr(iel) = smbr(iel) + thetp1*propce(iel,iptsta+2)
  enddo
!     Sinon : SMBR
else
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w5(iel)
  enddo
endif

!     Terme implicite
do iel = 1, ncel
  smbr(iel) = ( - volume(iel)*rtpa(iel,ifbiph) + smbr(iel) )      &
       /w4(iel)
enddo

! ---> Matrice

if(isto2t.gt.0) then
  thetap = thetv
else
  thetap = 1.d0
endif
do iel = 1, ncel
  rovsdt(iel) = (rovsdt(iel) + volume(iel)*thetap)/w4(iel)
enddo



!===============================================================================
! 3.3 RESOLUTION EFFECTIVE DE L'EQUATION DE F_BARRE
!===============================================================================


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
iescap = 0
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

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w2     , w3     , w4     , w5     , w6     ,                   &
   w7     , w8     , w9     , w10    ,                            &
   ra     )


!===============================================================================
! 4. RESOLUTION DE L'EQUATION DE PHI
!===============================================================================

ivar = iphiph
iclvar = iclphi
iclvaf = iclphi
ipp    = ipprtp(ivar)

if(iwarni(ivar).ge.1) then
  write(nfecra,1100) nomvar(ipp)
endif

!     S pour Source, V pour Variable
thets  = thetst
thetv  = thetav(ivar )

ipcroo = ipcrom
ipcvso = ipcvst
if(isto2t.gt.0) then
  if (iroext.gt.0) then
    ipcroo = ipproc(iroma)
  endif
  if(iviext.gt.0) then
    ipcvso = ipproc(ivista)
  endif
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 4.1 TERMES SOURCES  UTILISATEURS
!===============================================================================

call ustsv2                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   ia     ,                                                       &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , prdv2f , w1     ,          &
   smbr   , rovsdt ,                                              &
!        ------   ------
   viscf  , viscb  , xam    ,                                     &
   w2     , w3     , w4     , w5     , w6     , w7     ,          &
   w8     , w9     , w10    , dam    , drtp   ,                   &
   ra     )

!     Si on extrapole les T.S.
if(isto2t.gt.0) then
  do iel = 1, ncel
!       Sauvegarde pour echange
    tuexpe = propce(iel,iptsta+3)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta+3) = smbr(iel)
!       Second membre du pas de temps precedent
!       On suppose -ROVSDT > 0 : on implicite
!          le terme source utilisateur (le reste)
    smbr(iel) = rovsdt(iel)*rtpa(iel,ivar) - thets*tuexpe
!       Diagonale
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo
else
  do iel = 1, ncel
    smbr(iel)   = rovsdt(iel)*rtpa(iel,ivar) + smbr(iel)
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif

!===============================================================================
! 4.2 TERME SOURCE DE MASSE
!===============================================================================


if (ncesmp.gt.0) then

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBR par -Gamma RTPA et ROVSDT par Gamma (*theta)
  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isto2t , thetv ,    &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipriph) ,    &
   smbr   ,  rovsdt , w2 )

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+3) =                                 &
      propce(iel,iptsta+3) + w2(iel)
    enddo
!       Sinon on le met directement dans SMBR
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w2(iel)
    enddo
  endif

endif


!===============================================================================
! 4.3 TERME D'ACCUMULATION DE MASSE -(dRO/dt)*VOLUME
!    ET TERME INSTATIONNAIRE
!===============================================================================

! ---> Calcul de mij

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w2)

! ---> Ajout au second membre

do iel = 1, ncel
  smbr(iel) = smbr(iel)                                           &
              + iconv(ivar)*w2(iel)*rtpa(iel,ivar)
enddo

! ---> Ajout dans la diagonale de la matrice
!     Extrapolation ou non, meme forme par coherence avec bilsc2

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                       &
           + istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel) &
           - iconv(ivar)*w2(iel)*thetv
enddo

!===============================================================================
! 4.4 TERME SOURCE DE PHI
!     SMBR=rho*f_barre - phi/k*Pk +2/k*mu_t/sigmak*grad_phi*grad_k
!===============================================================================

!     Terme explicite, stocke temporairement dans W2

do iel = 1, ncel
!   Le terme en f_barre est pris en RTP et pas en RTPA ... a priori meilleur
!  Rq : si on reste en RTP, il faut modifier le cas de l'ordre 2 (qui
!       necessite RTPA pour l'extrapolation).
  w2(iel)   =  volume(iel)*                                       &
       ( propce(iel,ipcroo)*rtp(iel,ifbiph)                       &
         +2.d0/rtpa(iel,ikiph)*propce(iel,ipcvso)/sigmak*w1(iel) )
enddo

!     Si on extrapole les T.S : PROPCE
if(isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    propce(iel,iptsta+3) =                                   &
    propce(iel,iptsta+3) + w2(iel)
    smbr(iel) = smbr(iel) + thetp1*propce(iel,iptsta+3)
  enddo
!     Sinon : SMBR
else
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w2(iel)
  enddo
endif

!     Terme implicite
do iel = 1, ncel
  smbr(iel) = smbr(iel)                                           &
       - volume(iel)*prdv2f(iel)*rtpa(iel,iphiph)/rtpa(iel,ikiph)
enddo

! ---> Matrice

if(isto2t.gt.0) then
  thetap = thetv
else
  thetap = 1.d0
endif
do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                       &
       + volume(iel)*prdv2f(iel)/rtpa(iel,ikiph)*thetap
enddo



!===============================================================================
! 4.5 TERMES DE DIFFUSION
!===============================================================================
! ---> Viscosite
! Normalement, dans les equations du phi-model, seul la viscosite
!  turbulente intervient dans la diffusion de phi (le terme en mu
!  a disparu passant de f a f_barre). Mais tel
!  quel, cela rend le calcul instable (car mu_t tend vers 0 a la paroi
!  ce qui decouple phi de sa condition a la limite et le terme de diffusion
!  moleculaire etant integre dans f_barre, c'est comme s'il etait traite
!  en explicite).
!  -> on rajoute artificiellement de la diffusion (sachant que comme k=0 a
!  la paroi, on se moque de la valeur de phi).

  if( idiff(ivar).ge. 1 ) then
    do iel = 1, ncel
      w2(iel) = propce(iel,ipcvis) + propce(iel,ipcvst)/sigmak
    enddo

    call viscfa                                                   &
   !==========
 ( idebia , idebra ,                                              &
   imvisf ,                                                       &
   ia     ,                                                       &
   w2     ,                                                       &
   viscf  , viscb  ,                                              &
   ra     )

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

  endif



!===============================================================================
! 4.6 RESOLUTION EFFECTIVE DE L'EQUATION DE PHI
!===============================================================================

if(isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + thetp1*propce(iel,iptsta+3)
  enddo
endif

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
iescap = 0
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

call codits                                                       &
!==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ia     ,                                                       &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   ra     )

!===============================================================================
! 10. CLIPPING
!===============================================================================

   call clpv2f                                                    &
   !==========
 ( ncelet , ncel   , nvar   ,                                     &
   iwarni(iphi) ,                                                 &
   propce , rtp    )

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000    format(/,                                                &
'   ** RESOLUTION DU V2F (PHI ET F_BARRE)        ',/,&
'      ----------------------------------        ',/)
 1100    format(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)

#else

 1000    format(/,                                                &
'   ** SOLVING V2F (PHI AND F_BAR)'               ,/,&
'      ---------------------------'               ,/)
 1100    format(/,'           SOLVING VARIABLE ',A8                  ,/)

#endif

!12345678 : MAX: 12345678901234 MIN: 12345678901234 NORM: 12345678901234
!----
! FIN
!----

return

end subroutine
