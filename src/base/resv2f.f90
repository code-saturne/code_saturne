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
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse , iphas  ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , ckupdc , smacel ,                            &
   viscf  , viscb  , prdv2f ,                                     &
   dam    , xam    , drtp   , smbr   , rovsdt ,                   &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     , w9     , w10    ,          &
   rdevel , rtuser , ra     )

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
! iphas            ! i  ! <-- ! phase number                                   !
! ivar             ! i  ! <-- ! variable number                                !
! isou             ! e  ! <-- ! numero de passage                              !
! ipp              ! e  ! <-- ! numero de variable pour sorties post           !
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
include "cstnum.h"
include "cstphy.h"
include "pointe.h"
include "period.h"
include "parall.h"

!===============================================================================

! Arguments

integer          idbia0 , idbra0
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nvar   , nscal  , nphas
integer          ncepdp , ncesmp
integer          nideve , nrdeve , nituse , nrtuse
integer          iphas

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
double precision prdv2f(ncelet)
double precision dam(ncelet), xam(nfac,2)
double precision drtp(ncelet), smbr(ncelet), rovsdt(ncelet)
double precision w1(ncelet), w2(ncelet), w3(ncelet)
double precision w4(ncelet), w5(ncelet), w6(ncelet)
double precision w7(ncelet), w8(ncelet), w9(ncelet)
double precision w10(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

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
integer          maxelt, ils
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

ipriph = ipr (iphas)
iuiph  = iu  (iphas)
ikiph  = ik(iphas)
ieiph  = iep(iphas)
iphiph = iphi(iphas)
ifbiph = ifb(iphas)

ipcrom = ipproc(irom  (iphas))
ipcvis = ipproc(iviscl(iphas))
ipcvst = ipproc(ivisct(iphas))
iflmas = ipprof(ifluma(iuiph))
iflmab = ipprob(ifluma(iuiph))

iclikp = iclrtp(ikiph,icoef)
iclphi = iclrtp(iphiph,icoef)
iclfbp = iclrtp(ifbiph,icoef)

if(isto2t(iphas).gt.0) then
  iptsta = ipproc(itstua(iphas))
else
  iptsta = 0
endif

d2s3 = 2.0d0/3.0d0
d1s4 = 1.0d0/4.0d0
d3s2 = 3.0d0/2.0d0

if(iwarni(iphiph).ge.1) then
  write(nfecra,1000)iphas
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
 ( idebia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr , nphas  ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphiph , imrgra , inc    , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w1     , w1     , w1     ,                                     &
   rtpa(1,iphiph ) , coefa(1,iclphi) , coefb(1,iclphi) ,          &
   w1     , w2     , w3     ,                                     &
!        ------   ------   ------
   w4     , w5     , w6     ,                                     &
   rdevel , rtuser , ra     )

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
   w4     , w4     , w4     ,                                     &
   rtpa(1,ikiph )  , coefa(1,iclikp) , coefb(1,iclikp) ,          &
   w4     , w5     , w6     ,                                     &
!        ------   ------   ------
   w7     , w8     , w9     ,                                     &
   rdevel , rtuser , ra     )

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
thets  = thetst(iphas)
thetv  = thetav(ivar )

ipcroo = ipcrom
ipcvlo = ipcvis
if(isto2t(iphas).gt.0) then
  if (iroext(iphas).gt.0) then
    ipcroo = ipproc(iroma(iphas))
  endif
  if(iviext(iphas).gt.0) then
    ipcvlo = ipproc(ivisla(iphas))
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

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
CALL IASIZE('RESV2F',IFINIA)

call ustsv2                                                       &
!==========
 ( ifinia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , ivar   ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , prdv2f , w1     ,          &
   smbr   , rovsdt ,                                              &
!        ------   ------
   viscf  , viscb  , xam    ,                                     &
   w2     , w3     , w4     , w5     , w6     , w7     ,          &
   w8     , w9     , w10    , dam    , drtp   ,                   &
   rdevel , rtuser , ra     )

!     Si on extrapole les T.S.
if(isto2t(iphas).gt.0) then
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
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse , imvisf ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w3     ,                                                       &
   viscf  , viscb  ,                                              &
   rdevel , rtuser , ra     )


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
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp , &
   iwarnp , nfecra ,                                              &
   epsrgp , climgp , extrap ,                                     &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w2     , w2     , w2     ,                                     &
   rtpa(1,iphiph)   , coefa(1,iclphi) , coefb(1,iclphi) ,         &
   viscf  , viscb  ,                                              &
   w3     , w3     , w3     ,                                     &
   w2     ,                                                       &
!        --
   w4     , w5     , w6     , w7     , w8     , w9     ,          &
   rdevel , rtuser , ra     )

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
if(isto2t(iphas).gt.0) then
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

if(isto2t(iphas).gt.0) then
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
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w2     , w3     , w4     , w5     , w6     ,                   &
   w7     , w8     , w9     , w10    ,                            &
   rdevel , rtuser , ra     )


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
thets  = thetst(iphas)
thetv  = thetav(ivar )

ipcroo = ipcrom
ipcvso = ipcvst
if(isto2t(iphas).gt.0) then
  if (iroext(iphas).gt.0) then
    ipcroo = ipproc(iroma(iphas))
  endif
  if(iviext(iphas).gt.0) then
    ipcvso = ipproc(ivista(iphas))
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

maxelt = max(ncelet, nfac, nfabor)
ils    = idebia
ifinia = ils + maxelt
CALL IASIZE('RESV2F',IFINIA)

call ustsv2                                                       &
!==========
 ( ifinia , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  , ncepdp , ncesmp ,                   &
   nideve , nrdeve , nituse , nrtuse ,                            &
   iphas  , ivar   ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   icepdc , icetsm , itypsm ,                                     &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smacel , prdv2f , w1     ,          &
   smbr   , rovsdt ,                                              &
!        ------   ------
   viscf  , viscb  , xam    ,                                     &
   w2     , w3     , w4     , w5     , w6     , w7     ,          &
   w8     , w9     , w10    , dam    , drtp   ,                   &
   rdevel , rtuser , ra     )

!     Si on extrapole les T.S.
if(isto2t(iphas).gt.0) then
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
 ( ncelet , ncel   , ncesmp , iiun   , isto2t(iphas) , thetv ,    &
   icetsm , itypsm(1,ivar) ,                                      &
   volume , rtpa(1,ivar) , smacel(1,ivar) , smacel(1,ipriph) ,    &
   smbr   ,  rovsdt , w2 )

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isto2t(iphas).gt.0) then
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
if(isto2t(iphas).gt.0) then
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

if(isto2t(iphas).gt.0) then
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
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse , imvisf ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   w2     ,                                                       &
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
! 4.6 RESOLUTION EFFECTIVE DE L'EQUATION DE PHI
!===============================================================================

if(isto2t(iphas).gt.0) then
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
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nvar   , nscal  , nphas  ,                                     &
   nideve , nrdeve , nituse , nrtuse ,                            &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   dam    , xam    , drtp   ,                                     &
   w1     , w2     , w3     , w4     , w5     ,                   &
   w6     , w7     , w8     , w9     ,                            &
   rdevel , rtuser , ra     )

!===============================================================================
! 10. CLIPPING
!===============================================================================

   call clpv2f                                                    &
   !==========
 ( ncelet , ncel   , nvar   , nphas  ,                            &
   iphas  , iwarni(iphi(iphas)) ,                                 &
   propce , rtp    )

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000    format(/,                                                &
'   ** PHASE ',I4,' RESOLUTION DU V2F (PHI ET F_BARRE)        ',/,&
'      -----------------------------------------------        ',/)
 1100    format(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)

#else

 1000    format(/,                                                &
'   ** PHASE ',I4,' SOLVING V2F (PHI AND F_BAR)'               ,/,&
'      ----------------------------------------'               ,/)
 1100    format(/,'           SOLVING VARIABLE ',A8                  ,/)

#endif

!12345678 : MAX: 12345678901234 MIN: 12345678901234 NORM: 12345678901234
!----
! FIN
!----

return

end subroutine
