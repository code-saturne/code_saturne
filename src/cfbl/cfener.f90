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

subroutine cfener &
!================

 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce ,                            &
   ckupdc , smacel ,                                              &
   viscf  , viscb  ,                                              &
   smbrs  , rovsdt )

!===============================================================================
! FONCTION :
! ----------

! SOLVING OF A CONVECTION-DIFFUSION EQUATION WITH SOURCE TERMS
!   FOR TOTAL ENERGY ON ONE TIME-STEP
!   (COMPRESSIBLE ALGORITHM IN P,U,E)

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ncepdp           ! i  ! <-- ! number of cells with head loss                 !
! ncesmp           ! i  ! <-- ! number of cells with mass source term          !
! iscal            ! i  ! <-- ! scalar number                                  !
! icepdc(ncelet    ! te ! <-- ! numero des ncepdp cellules avec pdc            !
! icetsm(ncesmp    ! te ! <-- ! numero des cellules a source de masse          !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. ustsma)                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! ckupdc           ! tr ! <-- ! work array for the head loss                   !
!  (ncepdp,6)      !    !     !                                                !
! smacel           ! tr ! <-- ! variable value associated to the mass source   !
! (ncesmp,*   )    !    !     ! term (for ivar=ipr, smacel is the mass flux    !
!                  !    !     ! \f$ \Gamma^n \f$)                              !
! viscf(nfac)      ! tr ! --- ! visc*surface/dist at internal faces            !
! viscb(nfabor     ! tr ! --- ! visc*surface/dist at boundary faces            !
! smbrs(ncelet     ! tr ! --- ! work array for second member                   !
! rovsdt(ncelet    ! tr ! --- ! work array for unsteady term                   !
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
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use cfpoin
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          iscal

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision smbrs(ncelet)
double precision rovsdt(ncelet)

! Local variables

character*80     chaine
integer          ivar
integer          ifac  , iel
integer          init  , isqrt , iii
integer          ipcvst, ipcvsl, iflmas, iflmab
integer          ippvar, ipp   , icvflb
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp, nitmap
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
double precision epsrgp, climgp, extrap, blencp, epsilp
double precision sclnor, thetap, epsrsp, relaxp

integer          inc    , iccocg , imucpp , idftnp , iswdyp
integer          ivar0  , ii , jj
integer          iccfth , imodif
integer          iel1  , iel2
integer          iterns

double precision flux
double precision dijpfx, dijpfy, dijpfz, pnd  , pip   , pjp
double precision diipfx, diipfy, diipfz, djjpfx, djjpfy, djjpfz

double precision rvoid(1)

double precision, allocatable, dimension(:) :: wb
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w4, w5, w6
double precision, allocatable, dimension(:) :: w7, w8, w9
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: brom, crom, cromo
double precision, dimension(:,:), pointer :: coefau
double precision, dimension(:,:,:), pointer :: coefbu
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: coefa_p, coefb_p

!===============================================================================

!===============================================================================
! 1. INITIALIZATION
!===============================================================================

! Allocate a temporary array
allocate(wb(nfabor))

! Allocate work arrays
allocate(grad(ncelet,3))
allocate(w1(ncelet))
allocate(w4(ncelet), w5(ncelet), w6(ncelet))
allocate(w7(ncelet), w8(ncelet), w9(ncelet))
allocate(dpvar(ncelet))

! Computation number and post-treatment number of the scalar total energy
ivar   = isca(iscal)
ippvar = ipprtp(ivar)

! Physical property numbers
call field_get_val_s(icrom, crom)
call field_get_val_s(ibrom, brom)
call field_get_val_prev_s(icrom, cromo)

ipcvst = ipproc(ivisct)

call field_get_key_int(ivarfl(ivar), kimasf, iflmas)
call field_get_key_int(ivarfl(ivar), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

if(ivisls(iscal).gt.0) then
  ipcvsl = ipproc(ivisls(iscal))
else
  ipcvsl = 0
endif

! Prints
call field_get_label(ivarfl(ivar), chaine)

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) chaine(1:8)
endif

!===============================================================================
! 2. SOURCE TERMS
!===============================================================================

! Theta-scheme:
! for now, theta=1 is assumed and the theta-scheme is not implemented

! --> Initialization

do iel = 1, ncel
  smbrs(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

! HEAT VOLUMIC SOURCE TERM: RHO * PHI *VOLUME
! =================================

call ustssc                                                                     &
!==========
( nvar   , nscal  , ncepdp , ncesmp ,                                           &
  iscal  ,                                                                      &
  icepdc , icetsm , itypsm ,                                                    &
  dt     , rtpa   , rtp    , propce ,                                           &
  ckupdc , smacel , smbrs  , rovsdt )

do iel = 1, ncel
  smbrs(iel) = smbrs(iel) + rovsdt(iel)*rtp(iel,ivar)
  rovsdt(iel) = max(-rovsdt(iel),zero)
enddo


! MASS SOURCE TERMS
! =================

! GAMMA(IEL) = SMACEL(IEL,IPR)

! Implicit term : GAMMA*VOLUME
!                                                        n
! Explicit term : GAMMA*VOLUME*e   - GAMMA*VOLUME*e
!                                     inj
if (ncesmp.gt.0) then
  iterns = 1
  call catsma ( ncelet , ncel , ncesmp , iterns ,                               &
                isno2t, thetav(ivar),                                           &
                icetsm , itypsm(1,ivar) ,                                       &
                volume , rtpa(1,ivar) , smacel(1,ivar) ,                        &
                smacel(1,ipr) , smbrs , rovsdt , w1    )
endif


!                          RHO*VOLUME
! UNSTEADY IMPLICIT TERM : ----------
! ======================       DT

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                                     &
                + istat(ivar)*(cromo(iel)/dt(iel))*volume(iel)
enddo

!                                       __        v
!     TERME DE DISSIPATION VISQUEUS   : >  ((SIGMA *U).n)  *S
!     ==============================    --               ij  ij

if (idiff(iu).ge. 1) then
  call cfdivs                                                     &
  !==========
 ( rtpa   , propce ,                                              &
   smbrs  , rtp(1,iu), rtp(1,iv), rtp(1,iw) )
endif

!                              __   P        n+1
! PRESSURE TRANSPORT TERM  : - >  (---)  *(Q    .n)  *S
! =======================      --  RHO ij   pr     ij  ij

do iel = 1, ncel
  w9(iel) = crom(iel)
enddo

!     Avec Reconstruction : ca pose probleme pour l'instant

!   Calcul du gradient de P/RHO

!      do iel = 1, ncel
!        w7(iel) = rtp(iel,ipr)/w9(iel)
!      enddo

! Rq : A defaut de connaitre les parametres pour P/RHO on prend ceux de P

!      iii = ipr
!      inc = 1
!      iccocg = 1
!      nswrgp = nswrgr(iii)
!      imligp = imligr(iii)
!      iwarnp = iwarni(iii)
!      epsrgp = epsrgr(iii)
!      climgp = climgr(iii)
!      extrap = extrag(iii)

!       On alloue localement 2 tableaux de NFABOR pour le calcul
!       de COEFA et COEFB de P/RHO

!      allocate(coefap(nfabor))
!      allocate(coefbp(nfabor))

!      do ifac = 1, nfabor
!        coefap(ifac) = zero
!        coefbp(ifac) = 1.d0
!      enddo

!  IVAR0 = 0 (indique pour la periodicite de rotation que la variable
!     n'est pas la vitesse ni Rij)
!      ivar0 = 0
!      call grdcel
!      !==========
!     & ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,
!     &   iwarnp , nfecra , epsrgp , climgp , extrap ,
!     &   w7     , coefap , coefbp ,
!     &   grad   )

!     Faces internes
!      do ifac = 1, nfac

!        ii = ifacel(1,ifac)
!        jj = ifacel(2,ifac)

!        dijpfx = dijpf(1,ifac)
!        dijpfy = dijpf(2,ifac)
!        dijpfz = dijpf(3,ifac)

!        pnd   = pond(ifac)

!        Calcul II' et JJ'

!        diipfx = cdgfac(1,ifac) - (xyzcen(1,ii)+ (1.d0-pnd) * dijpfx)
!        diipfy = cdgfac(2,ifac) - (xyzcen(2,ii)+ (1.d0-pnd) * dijpfy)
!        diipfz = cdgfac(3,ifac) - (xyzcen(3,ii)+ (1.d0-pnd) * dijpfz)
!        djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj)+       pnd  * dijpfx
!        djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj)+       pnd  * dijpfy
!        djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj)+       pnd  * dijpfz

!        pip = w7(ii) +grad(ii,1)*diipfx+grad(ii,2)*diipfy+grad(ii,3)*diipfz
!        pjp = w7(jj) +grad(jj,1)*djjpfx+grad(jj,2)*djjpfy+grad(jj,3)*djjpfz

!        flui = (imasfl(ifac)+abs(imasfl(ifac)))
!        fluj = (imasfl(ifac)-abs(imasfl(ifac)))

!        viscf(ifac) = -(pnd*pip*flui+pnd*pjp*fluj)

!      enddo

!     Sans Reconstruction

!     En periodique et parallele, echange avant utilisation
!       des valeurs aux faces
if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(w9)
  !==========
endif

! Internal faces
do ifac = 1, nfac
  iel1 = ifacel(1,ifac)
  iel2 = ifacel(2,ifac)
  viscf(ifac) =                                                                 &
     - rtp(iel1,ipr)/w9(iel1) * 0.5d0*(imasfl(ifac) +abs(imasfl(ifac)))         &
     - rtp(iel2,ipr)/w9(iel2) * 0.5d0*(imasfl(ifac) -abs(imasfl(ifac)))
enddo

! Boundary faces: for the faces where a flux (Rusanov or analytical) has been
! computed, the standard contribution is replaced by this flux in bilsc2.

call field_get_coefa_s(ivarfl(ipr), coefa_p)
call field_get_coefb_s(ivarfl(ipr), coefb_p)

do ifac = 1, nfabor
  if (icvfli(ifac).eq.0) then
    iel = ifabor(ifac)
    viscb(ifac) = - bmasfl(ifac)                                                &
                    * (coefa_p(ifac) + coefb_p(ifac)*rtp(iel,ipr))              &
                    / brom(ifac)
  else
    viscb(ifac) = 0.d0
  endif
enddo

!     Divergence
init = 0
call divmas(init, viscf, viscb, smbrs)


! GRAVITATION FORCE TERM : RHO*g.U *VOLUME
! ======================

do iel = 1, ncel
  smbrs(iel) = smbrs(iel) + w9(iel)*volume(iel)                                 &
                           *( gx*rtp(iel,iu)                                    &
                            + gy*rtp(iel,iv)                                    &
                            + gz*rtp(iel,iw) )
enddo

!                                  Kij*Sij           LAMBDA   Cp   MUT
!     FACE DIFFUSION "VELOCITY" : --------- avec K = ------ + -- .------
!     =========================    IJ.nij              Cv     Cv  SIGMAS

if( idiff(ivar).ge. 1 ) then

!     MUT/SIGMAS
  do iel = 1, ncel
    w1(iel) = propce(iel,ipcvst)/sigmas(iscal)
  enddo
!     CP*MUT/SIGMAS
  if(icp.gt.0) then
    do iel = 1, ncel
      w1(iel) = w1(iel)*propce(iel,ipproc(icp))
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel)*cp0
    enddo
  endif
!     (CP/CV)*MUT/SIGMAS
  if(icv.gt.0) then
    do iel = 1, ncel
      w1(iel) = w1(iel)/propce(iel,ipproc(icv))
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel)/cv0
    enddo
  endif
!     (CP/CV)*MUT/SIGMAS+LAMBDA/CV
  if(ipcvsl.eq.0)then
    do iel = 1, ncel
      w1(iel) = w1(iel) + visls0(iscal)
    enddo
  else
    do iel = 1, ncel
      w1(iel) = w1(iel) + propce(iel,ipcvsl)
    enddo
  endif

  call viscfa                                                     &
  !==========
 ( imvisf ,                                                       &
   w1     ,                                                       &
   viscf  , viscb  )


!     COMPLEMENTARY DIFFUSIVE TERM : - div( K grad ( epsilon - Cv.T ) )
!     ============================                   1  2
!                                    - div( K grad ( -.u  ) )
!                                                    2

! Complementary term at cell centers
  iccfth = 7
  imodif = 0
  call cfther                                                                   &
  !==========
 ( nvar   ,                                                                     &
   iccfth , imodif ,                                                            &
   rtp    ,                                                                     &
   w9     , wb     , w8     , rvoid  , rvoid )

! Divergence computation with reconstruction


! Computation of the gradient of (0.5*u*u+EPSILONsup)

  do iel = 1, ncel
    w7(iel) =0.5d0*( rtp(iel,iu)**2                                             &
                    +rtp(iel,iv)**2                                             &
                    +rtp(iel,iw)**2 ) + w9(iel)
  enddo

! Note : by default, since the parameters are unknowns, the velocity parameters
! are taken

  iii = iu
  inc = 1
  iccocg = 1
  nswrgp = nswrgr(iii)
  imligp = imligr(iii)
  iwarnp = iwarni(iii)
  epsrgp = epsrgr(iii)
  climgp = climgr(iii)
  extrap = extrag(iii)

! Allocate temporary arrays
  allocate(coefap(nfabor))
  allocate(coefbp(nfabor))

  do ifac = 1, nfabor
    coefap(ifac) = zero
    coefbp(ifac) = 1.d0
  enddo

!  IVAR0 = 0 (indicates, for the rotation periodicity,
!  that the variable is not Rij)
  ivar0 = 0
  call grdcel                                                       &
  !==========
   ( ivar0  , imrgra , inc    , iccocg , nswrgp , imligp ,          &
     iwarnp , nfecra , epsrgp , climgp , extrap ,                   &
     w7     , coefap , coefbp ,                                     &
     grad   )

! Free memory
  deallocate(coefap, coefbp)

! Internal faces

  do ifac = 1, nfac

    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)

    dijpfx = dijpf(1,ifac)
    dijpfy = dijpf(2,ifac)
    dijpfz = dijpf(3,ifac)

    pnd   = pond(ifac)

! Computation of II' and JJ'

    diipfx = cdgfac(1,ifac) - (xyzcen(1,ii) + (1.d0-pnd) * dijpfx)
    diipfy = cdgfac(2,ifac) - (xyzcen(2,ii) + (1.d0-pnd) * dijpfy)
    diipfz = cdgfac(3,ifac) - (xyzcen(3,ii) + (1.d0-pnd) * dijpfz)
    djjpfx = cdgfac(1,ifac) -  xyzcen(1,jj) +  pnd  * dijpfx
    djjpfy = cdgfac(2,ifac) -  xyzcen(2,jj) +  pnd  * dijpfy
    djjpfz = cdgfac(3,ifac) -  xyzcen(3,jj) +  pnd  * dijpfz

    pip = w7(ii) + grad(ii,1)*diipfx+grad(ii,2)*diipfy+grad(ii,3)*diipfz
    pjp = w7(jj) + grad(jj,1)*djjpfx+grad(jj,2)*djjpfy+grad(jj,3)*djjpfz

    flux = viscf(ifac)*(pip-pjp)

    smbrs(ii) = smbrs(ii) + flux
    smbrs(jj) = smbrs(jj) - flux

  enddo

  ! Assembling based on boundary faces
  ! for the faces where a flux or a temperature is imposed,
  ! all is taken into account by the energy diffusion term.
  ! Hence the contribution of the terms in u2 and e-CvT shouldn't be taken into
  ! account when ifbet(ifac).ne.0

  call field_get_coefa_v(ivarfl(iu), coefau)
  call field_get_coefb_v(ivarfl(iu), coefbu)

  do ifac = 1, nfabor

    if (ifbet(ifac).eq.0) then

      iel = ifabor(ifac)

      flux = viscb(ifac)*(w1(iel)/distb(ifac))*                  &
            ( w9(iel) - wb(ifac)                                 &
            + 0.5d0*( rtp(iel,iu)**2 -                           &
              ( coefau(1, ifac) + coefbu(1, 1, ifac)*rtp(iel,iu) &
                                + coefbu(1, 2, ifac)*rtp(iel,iv) &
                                + coefbu(1, 3, ifac)*rtp(iel,iw) &
              )**2                                               &
                   + rtp(iel,iv)**2 -                            &
              ( coefau(2, ifac) + coefbu(2, 1, ifac)*rtp(iel,iu) &
                                + coefbu(2, 2, ifac)*rtp(iel,iv) &
                                + coefbu(2, 3, ifac)*rtp(iel,iw) &
              )**2                                               &
                   + rtp(iel,iw)**2 -                            &
              ( coefau(3, ifac) + coefbu(3, 1, ifac)*rtp(iel,iu) &
                                + coefbu(3, 2, ifac)*rtp(iel,iv) &
                                + coefbu(3, 3, ifac)*rtp(iel,iw) &
              )**2                                               &
            ))


      smbrs(iel) = smbrs(iel) + flux
    endif

  enddo

else

  do ifac = 1, nfac
    viscf(ifac) = 0.d0
  enddo
  do ifac = 1, nfabor
    viscb(ifac) = 0.d0
  enddo

endif

!===============================================================================
! 4. Solving
!===============================================================================

iconvp = iconv (ivar)
idiffp = idiff (ivar)
ireslp = iresol(ivar)
nitmap = nitmax(ivar)
ndircp = ndircl(ivar)
nswrsp = nswrsm(ivar)
nswrgp = nswrgr(ivar)
imligp = imligr(ivar)
ircflp = ircflu(ivar)
ischcp = ischcv(ivar)
isstpp = isstpc(ivar)
imgrp  = imgr  (ivar)
ncymxp = ncymax(ivar)
nitmfp = nitmgf(ivar)
ipp    = ippvar
iwarnp = iwarni(ivar)
blencp = blencv(ivar)
epsilp = epsilo(ivar)
epsrsp = epsrsm(ivar)
epsrgp = epsrgr(ivar)
climgp = climgr(ivar)
extrap = extrag(ivar)
relaxp = relaxv(ivar)
thetap = thetav(ivar)
iescap = 0

!  idtvar = 1  => unsteady
imucpp = 0  ! not a thermal scalar
idftnp = 1  ! scalar viscosity
iswdyp = 0  ! no dynamic relaxation

! impose boundary convective at some faces (face indicator icvfli)
icvflb = 1

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

call codits                                                      &
!==========
( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
  imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
  ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
  imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
  blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
  relaxp , thetap ,                                              &
  rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
  coefap , coefbp , cofafp , cofbfp ,                            &
  imasfl, bmasfl,                                                &
  viscf  , viscb  , rvoid  , viscf  , viscb  , rvoid  ,          &
  rvoid  , rvoid  ,                                              &
  icvflb , icvfli ,                                              &
  rovsdt , smbrs  , rtp(1,ivar)     , dpvar  ,                   &
  rvoid  , rvoid  )


!===============================================================================
! 5. PRINTINGS AND CLIPPINGS
!===============================================================================

! dummy value
iii = 1

call clpsca(ncelet, ncel, iscal, rtp(1,iii), rtp)
!==========

! --- Traitement utilisateur pour gestion plus fine des bornes
!       et actions correctives éventuelles.
  iccfth = -4
  imodif = 0
  call cfther                                                     &
  !==========
 ( nvar   ,                                                       &
   iccfth , imodif ,                                              &
   rtp    ,                                                       &
   w6     , w7     , w8     , rvoid  , rvoid )


! Explicit balance (see codits : the increment is removed)

if (iwarni(ivar).ge.2) then
  do iel = 1, ncel
    smbrs(iel) = smbrs(iel)                                                     &
            - istat(ivar)*(crom(iel)/dt(iel))*volume(iel)              &
                *(rtp(iel,ivar)-rtpa(iel,ivar))                                 &
                * max(0,min(nswrsm(ivar)-2,1))
  enddo
  isqrt = 1
  call prodsc(ncel,isqrt,smbrs,smbrs,sclnor)
  write(nfecra,1200)chaine(1:8) ,sclnor
endif

!===============================================================================
! 6. FINAL UPDATING OF THE PRESSURE (AND TEMPERATURE)
!===============================================================================
!                             n+1      n+1  n+1
! The state equation is used P   =P(RHO   ,H   )

! Computation of P and T at cell centers
iccfth = 24
imodif = 0
call cfther                                                                     &
!==========
( nvar   ,                                                                      &
  iccfth , imodif ,                                                             &
  rtp    ,                                                                      &
  rtp(1,ipr) , rtp(1,isca(itempk)) , w8     ,                                   &
  rvoid  , rvoid )

!===============================================================================
! 7. COMMUNICATION OF PRESSURE, ENERGY AND TEMPERATURE
!===============================================================================

if (irangp.ge.0.or.iperio.eq.1) then
  call synsca(rtp(1,ipr))
  !==========
  call synsca(rtp(1,ivar))
  !==========
  call synsca(rtp(1,isca(itempk)))
  !==========
endif

! Free memory
deallocate(wb)
deallocate(grad)
deallocate(w1)
deallocate(w4, w5, w6)
deallocate(w7, w8, w9)

!--------
! FORMATS
!--------

 1000 format(/,                                                   &
'   ** RESOLUTION FOR THE VARIABLE ',A8                        ,/,&
'      ---------------------------                            ',/)
 1200 format(1X,A8,' : EXPLICIT BALANCE = ',E14.5)

!----
! FIN
!----

return

end subroutine
