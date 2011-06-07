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

subroutine reseps &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   , isou   , ipp    ,                                     &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , grdvit , produc , gradro ,                   &
   ckupdc , smcelp , gamma  ,                                     &
   viscf  , viscb  ,                                              &
   tslagr ,                                                       &
   smbr   , rovsdt )

!===============================================================================
! FONCTION :
! ----------

! RESOLUTION DES EQUATIONS CONVECTION DIFFUSION TERME SOURCE
!   POUR EPSILON DANS LE CAS DU MODELE Rij-epsilon

!   On a ici ISOU   = 7 (7 ieme variable du Rij-epsilon)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! ivar             ! i  ! <-- ! variable number                                !
! isou             ! e  ! <-- ! numero de passage                              !
! ipp              ! e  ! <-- ! numero de variable pour sorties post           !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itpsmp           ! te ! <-- ! type de source de masse pour la                !
! (ncesmp)         !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! grdvit           ! tr ! --- ! tableau de travail pour terme grad             !
!  (ncelet,3,3)    !    !     !    de vitesse     uniqt pour iturb=31          !
! produc           ! tr ! <-- ! tableau de travail pour production             !
!  (6,ncelet)      !    !     ! (sans rho volume) uniqt pour iturb=30          !
! gradro(ncelet,3) ! tr ! <-- ! tableau de travail pour grad rom               !
! ckupdc           ! tr ! <-- ! tableau de travail pour pdc                    !
!  (ncepdp,6)      !    !     !                                                !
! smcelp(ncesmp    ! tr ! <-- ! valeur de la variable associee a la            !
!                  !    !     !  source de masse                               !
! gamma(ncesmp)    ! tr ! <-- ! valeur du flux de masse                        !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
! tslagr           ! tr ! <-- ! terme de couplage retour du                    !
!  (ncelet,*)      !    !     !   lagrangien                                   !
! smbr(ncelet      ! tr ! --- ! tableau de travail pour sec mem                !
! rovsdt(ncelet    ! tr ! --- ! tableau de travail pour terme instat           !
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
use parall
use period
use lagran
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar   , isou   , ipp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itpsmp(ncesmp)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision produc(6,ncelet), grdvit(ncelet,3,3)
double precision gradro(ncelet,3)
double precision ckupdc(ncepdp,6)
double precision smcelp(ncesmp), gamma(ncesmp)
double precision viscf(nfac), viscb(nfabor)
double precision tslagr(ncelet,*)
double precision smbr(ncelet), rovsdt(ncelet)

! Local variables

integer          init  , ifac  , iel   , inc   , iccocg
integer          ii    ,jj     , iiun
integer          iclvar, iclvaf
integer          ipcrom, ipcroo, ipcvis, ipcvst, iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          iptsta

double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision trprod , trrij ,csteps, rctse
double precision grdpx,grdpy,grdpz,grdsn
double precision surfn2
double precision tseps , kseps , ceps2
double precision tuexpe, thets , thetv , thetap, thetp1

double precision rvoid(1)

double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1, w2, w3
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7, w8, w9

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet), w3(ncelet))
allocate(w4(ncelet), w5(ncelet), w6(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))


if(iwarni(ivar).ge.1) then
  write(nfecra,1000) nomvar(ipp)
endif

ipcrom = ipproc(irom  )
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
iflmas = ipprof(ifluma(iu))
iflmab = ipprob(ifluma(iu))

iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

!     Constante Ce2, qui vaut CE2 pour ITURB=30 et CSSGE2 pour ITRUB=31
if (iturb.eq.30) then
  ceps2 = ce2
else
  ceps2 = cssge2
endif

!     S pour Source, V pour Variable
thets  = thetst
thetv  = thetav(ivar )

ipcroo = ipcrom
if(isto2t.gt.0.and.iroext.gt.0) then
  ipcroo = ipproc(iroma)
endif
if(isto2t.gt.0) then
  iptsta = ipproc(itstua)
else
  iptsta = 0
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 2. TERMES SOURCES  UTILISATEURS
!===============================================================================

call ustsri                                                       &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , ckupdc , smcelp , gamma  , grdvit , produc , &
   smbr   , rovsdt )

!     Si on extrapole les T.S.
if(isto2t.gt.0) then
  do iel = 1, ncel
!       Sauvegarde pour echange
    tuexpe = propce(iel,iptsta+isou-1)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta+isou-1) = smbr(iel)
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
! 3.  TERMES SOURCES LAGRANGIEN : COUPLAGE RETOUR
!===============================================================================

!     Ordre 2 non pris en compte
 if (iilagr.eq.2 .and. ltsdyn.eq.1) then

   do iel = 1,ncel
! Ts sur eps
     tseps = -0.5d0 * ( tslagr(iel,itsr11)                        &
                      + tslagr(iel,itsr22)                        &
                      + tslagr(iel,itsr33) )
! rapport k/eps
     kseps = 0.5d0 * ( rtpa(iel,ir11)                           &
                     + rtpa(iel,ir22)                           &
                     + rtpa(iel,ir33) )                         &
                     / rtpa(iel,iep)

     smbr(iel)   = smbr(iel) + ce4 *tseps *rtpa(iel,iep) /kseps
     rovsdt(iel) = rovsdt(iel) + max( (-ce4*tseps/kseps) , zero)
   enddo

 endif

!===============================================================================
! 4. TERME SOURCE DE MASSE
!===============================================================================


if (ncesmp.gt.0) then

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBR par -Gamma RTPA et ROVSDT par Gamma (*theta)
  call catsma                                                     &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isto2t , thetv ,    &
   icetsm , itpsmp ,                                              &
   volume , rtpa(1,ivar) , smcelp , gamma  ,                      &
   smbr   ,  rovsdt , w1 )

!       Si on extrapole les TS on met Gamma Pinj dans PROPCE
  if(isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) =                                 &
      propce(iel,iptsta+isou-1) + w1(iel)
    enddo
!       Sinon on le met directement dans SMBR
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w1(iel)
    enddo
  endif

endif

!===============================================================================
! 5. TERME D'ACCUMULATION DE MASSE -(dRO/dt)*VOLUME
!    ET TERME INSTATIONNAIRE
!===============================================================================

! ---> Calcul de mij

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                  &
               ifacel,ifabor,propfa(1,iflmas),propfb(1,iflmab),w1)

! ---> Ajout au second membre

do iel = 1, ncel
  smbr(iel) = smbr(iel)                                           &
              + iconv(ivar)*w1(iel)*rtpa(iel,ivar)
enddo

! ---> Ajout dans la diagonale de la matrice
!     Extrapolation ou non, meme forme par coherence avec bilsc2

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                       &
           + istat(ivar)*(propce(iel,ipcrom)/dt(iel))*volume(iel) &
           - iconv(ivar)*w1(iel)*thetv
enddo


!===============================================================================
! 6. PRODUCTION RHO * Ce1 * epsilon / k * P
!    DISSIPATION RHO*Ce2.epsilon/k*epsilon
!===============================================================================


! ---> Calcul de k pour la suite du sous-programme
!       on utilise un tableau de travail puisqu'il y en a...
do iel = 1, ncel
  w8(iel) = 0.5d0 * (rtpa(iel,ir11) + rtpa(iel,ir22) + rtpa(iel,ir33))
enddo
! ---> Calcul de la trace de la production, suivant qu'on est en
!     Rij standard ou en SSG (utilisation de PRODUC ou GRDVIT)
if (iturb.eq.30) then
  do iel = 1, ncel
    w9(iel) = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
  enddo
else
  do iel = 1, ncel
    w9(iel) = -( rtpa(iel,ir11)*grdvit(iel,1,1) +               &
                 rtpa(iel,ir12)*grdvit(iel,1,2) +               &
                 rtpa(iel,ir13)*grdvit(iel,1,3) +               &
                 rtpa(iel,ir12)*grdvit(iel,2,1) +               &
                 rtpa(iel,ir22)*grdvit(iel,2,2) +               &
                 rtpa(iel,ir23)*grdvit(iel,2,3) +               &
                 rtpa(iel,ir13)*grdvit(iel,3,1) +               &
                 rtpa(iel,ir23)*grdvit(iel,3,2) +               &
                 rtpa(iel,ir33)*grdvit(iel,3,3) )
  enddo
endif


!     Terme explicite (Production)

do iel = 1, ncel
!     Demi-traces
  trprod = w9(iel)
  trrij  = w8(iel)
  w1(iel)   =             propce(iel,ipcroo)*volume(iel)*         &
       ce1*rtpa(iel,iep)/trrij*trprod
enddo

!     Si on extrapole les T.S : PROPCE
if(isto2t.gt.0) then
  do iel = 1, ncel
    propce(iel,iptsta+isou-1) =                                   &
    propce(iel,iptsta+isou-1) + w1(iel)
  enddo
!     Sinon : SMBR
else
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w1(iel)
  enddo
endif

!     Terme implicite (Dissipation)
do iel = 1, ncel
  trrij  = w8(iel)
  smbr(iel) = smbr(iel) - propce(iel,ipcrom)*volume(iel)*         &
                         ceps2*rtpa(iel,iep)**2/trrij
enddo

! ---> Matrice

if(isto2t.gt.0) then
  thetap = thetv
else
  thetap = 1.d0
endif
do iel = 1, ncel
  trrij  = w8(iel)
  rovsdt(iel) = rovsdt(iel) + ceps2*propce(iel,ipcrom)*volume(iel)&
                     *rtpa(iel,iep)/trrij*thetap
enddo

!===============================================================================
! 7. TERMES DE GRAVITE
!===============================================================================

if(igrari.eq.1) then

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  call rijthe                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   ivar   , isou   , ipp    ,                                     &
   rtp    , rtpa   , propce , propfa , propfb ,                   &
   coefa  , coefb  , gradro , w7     )

!     Si on extrapole les T.S. : PROPCE
if(isto2t.gt.0) then
  do iel = 1, ncel
     propce(iel,iptsta+isou-1) =                                  &
     propce(iel,iptsta+isou-1) + w7(iel)
   enddo
!     Sinon SMBR
 else
   do iel = 1, ncel
     smbr(iel) = smbr(iel) + w7(iel)
   enddo
 endif

endif


!===============================================================================
! 8.a TERMES DE DIFFUSION  A.grad(Eps) : PARTIE EXTRADIAGONALE EXPLICITE
!     RIJ STANDARD
!===============================================================================

if (iturb.eq.30) then

! ---> Calcul du grad(Eps)

  ! Allocate a temporary array for the gradient calculation
  allocate(grad(ncelet,3))

  iccocg = 1
  inc = 1

  nswrgp = nswrgr(ivar )
  imligp = imligr(ivar )
  iwarnp = iwarni(ivar )
  epsrgp = epsrgr(ivar )
  climgp = climgr(ivar )
  extrap = extrag(ivar )

  call grdcel                                                     &
  !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,ivar )   , coefa(1,iclvar) , coefb(1,iclvar) ,          &
   grad   )

! ---> Calcul des termes extradiagonaux de A.grad(Eps)

  do iel = 1 , ncel
    trrij  = w8(iel)
    csteps  = propce(iel,ipcroo) * crijep *trrij / rtpa(iel,iep)
    w4(iel) = csteps * (rtpa(iel,ir12)*grad(iel,2) + rtpa(iel,ir13)*grad(iel,3))
    w5(iel) = csteps * (rtpa(iel,ir12)*grad(iel,1) + rtpa(iel,ir23)*grad(iel,3))
    w6(iel) = csteps * (rtpa(iel,ir13)*grad(iel,1) + rtpa(iel,ir23)*grad(iel,2))
  enddo

! ---> Assemblage de { A.grad(Eps) } .S aux faces

  call vectds                                                     &
  !==========
( w4     , w5     , w6     ,                                      &
  viscf  , viscb  )

  init = 1
  call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,                &
                                   ifacel,ifabor,viscf,viscb,w4)

!     Si on extrapole les termes sources
  if(isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) =                                 &
           propce(iel,iptsta+isou-1) + w4(iel)
    enddo
!     Sinon
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w4(iel)
    enddo
  endif


!===============================================================================
! 8.b TERMES DE DIFFUSION  A.grad(Eps) : PARTIE DIAGONALE
!     RIJ STANDARD
!===============================================================================
!     Implicitation de (grad(eps).n)n en gradient facette
!     Si IDIFRE=1, terme correctif explicite
!        grad(eps)-(grad(eps).n)n calcule en gradient cellule
!     Les termes de bord sont uniquement pris en compte dans la partie
!        en (grad(eps).n)n
!     (W1,W2,W3) contient toujours le gradient de la variable traitee

!     La synchronisation des halos du gradient de epsilon a ete faite dans
!       grdcel. Pas utile de recommencer.

  if (idifre.eq.1) then

    do iel = 1, ncel
      trrij  = w8(iel)
      csteps = propce(iel,ipcroo) * crijep *trrij/rtpa(iel,iep)
      w4(iel)=csteps*rtpa(iel,ir11)
      w5(iel)=csteps*rtpa(iel,ir22)
      w6(iel)=csteps*rtpa(iel,ir33)
    enddo

! --->  TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

    if (irangp.ge.0.or.iperio.eq.1) then
      call syndia(w4, w5, w6)
      !==========
    endif


    do ifac = 1, nfac

      ii=ifacel(1,ifac)
      jj=ifacel(2,ifac)

      surfn2 = surfan(ifac)**2

      grdpx=0.5d0*(grad(ii,1)+grad(jj,1))
      grdpy=0.5d0*(grad(ii,2)+grad(jj,2))
      grdpz=0.5d0*(grad(ii,3)+grad(jj,3))
      grdsn=grdpx*surfac(1,ifac)+grdpy*surfac(2,ifac)             &
           +grdpz*surfac(3,ifac)
      grdpx=grdpx-grdsn*surfac(1,ifac)/surfn2
      grdpy=grdpy-grdsn*surfac(2,ifac)/surfn2
      grdpz=grdpz-grdsn*surfac(3,ifac)/surfn2

      viscf(ifac)= 0.5d0*(                                        &
           (w4(ii)+w4(jj))*grdpx*surfac(1,ifac)                   &
           +(w5(ii)+w5(jj))*grdpy*surfac(2,ifac)                  &
           +(w6(ii)+w6(jj))*grdpz*surfac(3,ifac))

    enddo

    ! Free memory
    deallocate(grad)

    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

    init = 1
    call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,              &
       ifacel,ifabor,viscf,viscb,w1)

!     Si on extrapole les termes sources
    if(isto2t.gt.0) then
      do iel = 1, ncel
        propce(iel,iptsta+isou-1) =                               &
             propce(iel,iptsta+isou-1) + w1(iel)
      enddo
!     Sinon
    else
      do iel = 1, ncel
        smbr(iel) = smbr(iel) + w1(iel)
      enddo
    endif

  endif

! ---> Viscosite orthotrope pour partie implicite

  if( idiff(ivar).ge. 1 ) then
    do iel = 1, ncel
      trrij  = w8(iel)
      rctse = propce(iel,ipcrom) * crijep * trrij/rtpa(iel,iep)
      w1(iel) = propce(iel,ipcvis) + idifft(ivar)*rctse*rtpa(iel,ir11)
      w2(iel) = propce(iel,ipcvis) + idifft(ivar)*rctse*rtpa(iel,ir22)
      w3(iel) = propce(iel,ipcvis) + idifft(ivar)*rctse*rtpa(iel,ir33)
    enddo

    call visort                                                   &
    !==========
 ( imvisf ,                                                       &
   w1     , w2     , w3     ,                                     &
   viscf  , viscb  )

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

  endif

else
!===============================================================================
! 8.c TERMES DE DIFFUSION
!     RIJ SSG
!===============================================================================
! ---> Viscosite

  if( idiff(ivar).ge. 1 ) then
    do iel = 1, ncel
      w1(iel) = propce(iel,ipcvis)                                &
           + idifft(ivar)*propce(iel,ipcvst)/sigmae
    enddo

    call viscfa                                                   &
   !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )

  else

    do ifac = 1, nfac
      viscf(ifac) = 0.d0
    enddo
    do ifac = 1, nfabor
      viscb(ifac) = 0.d0
    enddo

  endif

endif


!===============================================================================
! 9. RESOLUTION
!===============================================================================

if(isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + thetp1*propce(iel,iptsta+isou-1)
  enddo
endif

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ndircp = ndircl(ivar)
ireslp = iresol(ivar)
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
!MO      IPP    =
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
 ( nvar   , nscal  ,                                              &
   idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap ,                                     &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
                     coefa(1,iclvar) , coefb(1,iclvar) ,          &
                     coefa(1,iclvaf) , coefb(1,iclvaf) ,          &
                     propfa(1,iflmas), propfb(1,iflmab),          &
   viscf  , viscb  , viscf  , viscb  ,                            &
   rovsdt , smbr   , rtp(1,ivar)     ,                            &
   rvoid  )

!===============================================================================
! 10. IMPRESSIONS
!===============================================================================

! Free memory
deallocate(w1, w2, w3)
deallocate(w4, w5, w6)
deallocate(w7, w8, w9)

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)

#else

 1000 format(/,'           SOLVING VARIABLE ',A8           ,/)

#endif

!12345678 : MAX: 12345678901234 MIN: 12345678901234 NORM: 12345678901234
!----
! FIN
!----

return

end subroutine
