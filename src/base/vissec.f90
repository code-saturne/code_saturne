!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine vissec &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   icepdc , icetsm , itypsm ,                                     &
   rtpa   , propce , propfa , propfb ,                            &
   coefa  , coefb  , ckupdc , smacel ,                            &
   trav   ,                                                       &
   viscf  , viscb  )

!===============================================================================
! FONCTION :
! ----------

! AJOUT AU SECOND MEMBRE DES TERMES

!         GRAD( (K -2/3 MU) DIV(U) ) + DIV( MU (GRAD_TRANSPOSE(U)) )

! AVEC MU = MU_LAMINAIRE + MU_TURBULENT
!  ET  K = VISCOSITE EN VOLUME (NULLE EN GENERAL)

! DIV(U)            EST CALCULE A PARTIR DE FLUMAS/RHO
! GRAD_TRANSPOSE(U) EST UN GRADIENT CELLULE

! REMARQUES :
!  - Theoriquement le terme en div(u) devrait plutot etre calcule
!      par un gradient cellule, pour correspondre exactement au
!      terme en dUj/dxi. Mais comme la partie en dUi/dxj est
!      calculee completement autrement (gradient facette et implicitation)
!      de toute facon on n'aura jamais Trace(tau_ij)=0 exactement.
!  - Pour la meme raison, comme le terme en dUi/dxj est calcule sur les
!      elements de bord et pas celui en dUj/dxi, il est difficile de
!      traiter le terme en div(u) de maniere rigoureuse. Il est donc
!      conserve sur les elements de bord.
!  - En LES, le tenseur <(u-<u>)(u-<u>)> est modelise par mut <S>
!      et non pas par mut <S> - 2/3 mut Tr(<S>) Id + 2/3 k Id
!      de sorte que il n'apparait pas ici de mut div<u>
!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! rtpa             ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant prec)                     !
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
! trav(ncelet,3    ! tr ! <-- ! tableau de travail pour sec mem                !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist aux faces internes           !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist aux faces de bord            !
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
use cstphy
use entsor
use numvar
use optcal
use pointe, only: forbr, porosi
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(ndimfb,*)
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision trav(ncelet,3)
double precision viscf(nfac), viscb(nfabor)

! Local variables

integer          iccocg, inc, iel, ifac, ivar, isou, ii, jj, init
integer          idim, ig, it
integer          iclvaf
integer          nswrgp, imligp, iwarnp
integer          ipcrom, ipbrom, ipcvis, ipcvst, iflmas, iflmab
integer          ipcvsv

double precision epsrgp, climgp, extrap
double precision romf, d2s3m, vecfac
double precision hint

double precision, allocatable, dimension(:) :: vistot
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w4, w6
double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, dimension(:), pointer :: imasfl, bmasfl

!===============================================================================

!===============================================================================
! 1.  INITIALISATION
!===============================================================================

! Allocate temporary arrays
allocate(vistot(ncelet))


ipcrom = ipproc(irom  )
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)


if (ippmod(icompf).ge.0) then
  if (iviscv.gt.0) then
    ipcvsv = ipproc(iviscv)
  else
    ipcvsv = 0
  endif
else
  ipcvsv = -1
endif

ipbrom = ipprob(irom  )

call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

!     Si on extrapole les termes sources, on prend les prop a l'instant n
if (isno2t.gt.0) then
  if (iroext.gt.0) then
    ipcrom = ipproc(iroma )
    ipbrom = ipprob(iroma )
  endif
  if (iviext.gt.0) then
    ipcvis = ipproc(ivisla)
    ipcvst = ipproc(ivista)
  endif
!     Il faudrait aussi faire quelque chose pour le flux de masse, non ?
endif

! --- Calcul de la viscosite totale

if (itytur.eq.3) then
  !$omp parallel do firstprivate(ipcvis)
  do iel = 1, ncel
    vistot(iel) = propce(iel,ipcvis)
  enddo
else
  !$omp parallel do firstprivate(ipcvis, ipcvst)
  do iel = 1, ncel
    vistot(iel) = propce(iel,ipcvis) + propce(iel,ipcvst)
  enddo
endif

!===============================================================================
! 2.  CALCUL DES TERMES EN GRAD_TRANSPOSE
!===============================================================================

! Allocate a temporary array for the gradient calculation
allocate(grad(ncelet,3))
allocate(w4(ncelet), w6(ncelet))
allocate(coefap(nfabor), coefbp(nfabor))
do isou = 1, 3

  if (isou.eq.1) ivar = iu
  if (isou.eq.2) ivar = iv
  if (isou.eq.3) ivar = iw

  ! FIXME: the previous version of CS took the flux BCS to compute the gradient
  iclvaf = iclrtp(ivar,icoeff)

  do ifac = 1, nfabor
    iel = ifabor(ifac)

    if (itytur.eq.3) then
      hint = propce(iel,ipcvis)/distb(ifac)
    else
      hint = (propce(iel,ipcvis) + propce(iel,ipcvst))/distb(ifac)
    endif

    coefap(ifac) = -coefa(ifac,iclvaf)/hint
    coefbp(ifac) = 1-coefb(ifac,iclvaf)/hint
  enddo


! --- Calcul du gradient de la vitesse

  iccocg = 1
  inc    = 1
  nswrgp = nswrgr(ivar)
  imligp = imligr(ivar)
  iwarnp = iwarni(ivar)
  epsrgp = epsrgr(ivar)
  climgp = climgr(ivar)
  extrap = extrag(ivar)

  call grdcel &
  !==========
 ( ivar   , imrgra , inc    , iccocg , nswrgp , imligp ,          &
   iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
   rtpa(1,ivar)    , coefap , coefbp ,                            &
   grad   )

  !$omp parallel do
  do iel = 1, ncel
    w6(iel) = 1.d0
  enddo

  do ig = 1, ngrpb
    !$omp parallel do private(ifac) if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)
        w6(ifabor(ifac)) = 0.d0
      enddo
    enddo
  enddo

  ! Ghost cells must be updated before the loop on faces.

  if (irangp.ge.0.or.iperio.eq.1) then
    call synsca(w6)
    !==========
  endif

! --- Assemblage sur les faces internes

  do idim = 1, 3

    if (idim.eq.1) then
      !$omp parallel do
      do iel = 1, ncel
        w4(iel) = vistot(iel)*grad(iel,1)
      enddo
    elseif (idim.eq.2) then
      !$omp parallel do
      do iel = 1, ncel
        w4(iel) = vistot(iel)*grad(iel,2)
      enddo
    elseif (idim.eq.3) then
      !$omp parallel do
      do iel = 1, ncel
        w4(iel) = vistot(iel)*grad(iel,3)
      enddo
    endif

    ! With porosity
    if (iporos.eq.1) then
      !$omp parallel do
      do iel = 1, ncel
        w4(iel) = w4(iel)*porosi(iel)
      enddo
    endif

! On initialise TRAV(NCEL+1, NCELET)
!     (valeur bidon, mais pas NaN : les calculs sur le halo sont
!      par principe denue de sens, sauf exception)
    if (ncelet.gt.ncel) then
      !$omp parallel do if (ncelet - ncel > thr_n_min)
      do iel = ncel+1, ncelet
        trav(iel,idim) = 0.d0
      enddo
    endif

    ! If ghost cells are present, they must be updated before the
    ! loop on faces.

    if (irangp.ge.0.or.iperio.eq.1) then
      call synsca(w4)
      !==========
    endif

    do ig = 1, ngrpi
      !$omp parallel do firstprivate(isou) private(ifac, ii, jj, vecfac)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)
          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)
          vecfac = surfac(isou,ifac)*(w4(ii)+w4(jj))*0.5d0
          trav(ii,idim) = trav(ii,idim) + vecfac*w6(ii)
          trav(jj,idim) = trav(jj,idim) - vecfac*w6(jj)
        enddo
      enddo
    enddo

  enddo

enddo

! Free memory
deallocate(grad)
deallocate(w6)
deallocate(coefap, coefbp)

!===============================================================================
! 3.  CALCUL DES TERMES EN DIV
!===============================================================================
!  Pour periodicite et parallelisme, ROM est echange dans phyvar.
!     ou apres avoir ete calcule dans cfmsvl en compressible

!  Ici pour l'ordre 2 en temps, il faudrait tout prendre en n...

! Allocate a temporary array
allocate(w1(ncelet))

!$omp parallel do firstprivate(ipcrom) private(ii, jj, romf)
do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
  romf = (propce(ii,ipcrom)+propce(jj,ipcrom))*0.5d0
  viscf(ifac) = imasfl(ifac)/romf
enddo

!$omp parallel do firstprivate(ipbrom) if(nfabor > thr_n_min)
do ifac = 1, nfabor
  viscb(ifac) = bmasfl(ifac)/propfb(ifac,ipbrom)
enddo

init = 1
call divmas(ncelet,ncel,nfac,nfabor,init,nfecra,            &
!==========
            ifacel,ifabor,viscf ,viscb ,w1)

d2s3m = -2.d0/3.d0


if (ipcvsv.gt.0) then
  !$omp parallel do
  do iel = 1, ncel
    w4(iel) = ( propce(iel,ipcvsv) + d2s3m*vistot(iel) )   &
              * w1(iel)/volume(iel)
  enddo
else if (ipcvsv.eq.0) then
  !$omp parallel do
  do iel = 1, ncel
    w4(iel) = ( viscv0      + d2s3m*vistot(iel) )          &
              * w1(iel)/volume(iel)
  enddo
else
  if (itytur.eq.4) then
    !$omp parallel do firstprivate(d2s3m, ipcvis)
    do iel = 1, ncel
      w4(iel) = d2s3m*propce(iel,ipcvis)*w1(iel)/volume(iel)
    enddo
  else
    !$omp parallel do firstprivate(d2s3m)
    do iel = 1, ncel
      w4(iel) = d2s3m*vistot(iel)*w1(iel)/volume(iel)
    enddo
  endif
endif

! Free memory
deallocate(vistot)
deallocate(w1)

! ---> TRAITEMENT DU PARALLELISME ET DE LA PERIODICITE

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w4)
  !==========
endif


!$omp parallel do private(ii, jj)
do ifac = 1, nfac
  ii = ifacel(1,ifac)
  jj = ifacel(2,ifac)
  viscf (ifac) = (w4(ii)+w4(jj))*0.5d0
enddo

do isou = 1, 3

  idim = isou

  ! --- Assemblage sur les faces internes

  do ig = 1, ngrpi
    !$omp parallel do firstprivate(isou) private(ifac, ii, jj, vecfac) &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)
        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        vecfac = surfac(isou,ifac)*viscf(ifac)
        trav(ii,idim) = trav(ii,idim) + vecfac
        trav(jj,idim) = trav(jj,idim) - vecfac
      enddo
    enddo
  enddo

  ! --- Assemblage sur les faces de bord

  do ig = 1, ngrpb
    !$omp parallel do firstprivate(idim, isou) private(ifac, ii) &
    !$omp          if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)
        ii = ifabor(ifac)
        trav(ii,idim) = trav(ii,idim) + surfbo(isou,ifac)*w4(ii)
      enddo
    enddo
  enddo


  ! --- Calcul des efforts aux bords (partie 4/5)

  if (ineedf.eq.1) then

    do ig = 1, ngrpb
      !$omp parallel do firstprivate(isou) private(ifac, ii) &
      !$omp          if(nfabor > thr_n_min)
      do it = 1, nthrdb
        do ifac = iomplb(1,ig,it), iomplb(2,ig,it)
          ii = ifabor(ifac)
          forbr(isou,ifac) = forbr(isou,ifac) + surfbo(isou,ifac)*w4(ii)
        enddo
      enddo
    enddo

  endif

enddo

! Free memory
deallocate(w4)

return

end subroutine
