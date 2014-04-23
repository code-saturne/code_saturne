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

!===============================================================================
! Function:
! ---------

!> \file resssg.f90
!>
!> \brief This subroutine performs the solving of the Reynolds stress components
!> in \f$ R_{ij} - \varepsilon \f$ RANS (SSG) turbulence model.
!>
!> Remark:
!> - isou=1 for \f$ R_{11} \f$
!> - isou=2 for \f$ R_{22} \f$
!> - isou=3 for \f$ R_{33} \f$
!> - isou=4 for \f$ R_{12} \f$
!> - isou=5 for \f$ R_{13} \f$
!> - isou=6 for \f$ R_{23} \f$
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     ivar          variable number
!> \param[in]     isou          local variable number (7 here)
!> \param[in]     ipp           index for writing
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        type of mass source term for each variable
!>                               (see \ref ustsma)
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     gradv         tableau de travail pour terme grad
!>                                 de vitesse     uniqt pour iturb=31
!> \param[in]     produc        tableau de travail pour production
!> \param[in]     gradro        tableau de travail pour grad rom
!>                              (sans rho volume) uniqt pour iturb=30
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        value associated to each variable in the mass
!>                               source terms or mass rate (see \ref ustsma)
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     tslagr        coupling term for lagrangian
!> \param[in]     tslage        explicit source terms for the Lagrangian module
!> \param[in]     tslagi        implicit source terms for the Lagrangian module
!> \param[in]     smbr          working array
!> \param[in]     rovsdt        working array
!_______________________________________________________________________________

subroutine resssg &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   , isou   , ipp    ,                                     &
   icepdc , icetsm , itypsm ,                                     &
   dt     , rtp    , rtpa   , propce ,                            &
   gradv  , gradro ,                                              &
   ckupdc , smacel ,                                              &
   viscf  , viscb  ,                                              &
   tslage , tslagi ,                                              &
   smbr   , rovsdt )

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
use lagran
use pointe, only: visten
use mesh
use field
use field_operator
use cs_f_interfaces
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar   , isou   , ipp

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet), rtp(ncelet,nflown:nvar), rtpa(ncelet,nflown:nvar)
double precision propce(ncelet,*)
double precision gradv(3, 3, ncelet)
double precision gradro(ncelet,3)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision tslage(ncelet),tslagi(ncelet)
double precision smbr(ncelet), rovsdt(ncelet)

! Local variables

integer          iel
integer          ii    , jj    , kk    , iiun  , iii   , jjj
integer          ipcvis, iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          iptsta
integer          iprev , inc, iccocg, ll
integer          imucpp, idftnp, iswdyp
integer          indrey(3,3)
integer          icvflb
integer          ivoid(1)

double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision trprod, trrij , deltij
double precision tuexpr, thets , thetv , thetp1
double precision aiksjk, aikrjk, aii ,aklskl, aikakj
double precision xaniso(3,3), xstrai(3,3), xrotac(3,3), xprod(3,3), matrot(3,3)
double precision xrij(3,3), xnal(3), xnoral, xnnd
double precision d1s2, d1s3, d2s3
double precision alpha3
double precision pij, phiij1, phiij2, epsij
double precision phiijw, epsijw
double precision ccorio

double precision rvoid(1)

character(len=80) :: label
double precision, allocatable, dimension(:,:) :: grad
double precision, allocatable, dimension(:) :: w1, w2
double precision, allocatable, dimension(:) :: w7
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, cromo
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet), w2(ncelet))
allocate(dpvar(ncelet))
allocate(viscce(6,ncelet))
allocate(weighf(2,nfac))
allocate(weighb(nfabor))

! Initialize variables to avoid compiler warnings
iii = 0
jjj = 0

if (iwarni(ivar).ge.1) then
  call field_get_label(ivarfl(ivar), label)
  write(nfecra,1000) label
endif

call field_get_val_s(icrom, crom)
ipcvis = ipproc(iviscl)
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

d1s2   = 1.d0/2.d0
d1s3   = 1.d0/3.d0
d2s3   = 2.d0/3.d0

deltij = 1.0d0
if(isou.gt.3) then
  deltij = 0.0d0
endif

!     S pour Source, V pour Variable
thets  = thetst
thetv  = thetav(ivar )

if (isto2t.gt.0.and.iroext.gt.0) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
endif
iptsta = 0
if (isto2t.gt.0) then
  iptsta = ipproc(itstua)
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

! Coefficient of the "Coriolis-type" term
if (icorio.eq.1) then
  ! Relative velocity formulation
  ccorio = 2.d0
elseif (iturbo.eq.1) then
  ! Mixed relative/absolute velocity formulation
  ccorio = 1.d0
else
  ccorio = 0.d0
endif

!===============================================================================
! 2. User source terms
!===============================================================================

call cs_user_turbulence_source_terms &
!===================================
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivarfl(ivar)    ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   smbr   , rovsdt )

!     Si on extrapole les T.S.
if(isto2t.gt.0) then
  do iel = 1, ncel
!       Sauvegarde pour echange
    tuexpr = propce(iel,iptsta+isou-1)
!       Pour la suite et le pas de temps suivant
    propce(iel,iptsta+isou-1) = smbr(iel)
!       Second membre du pas de temps precedent
!       On suppose -ROVSDT > 0 : on implicite
!          le terme source utilisateur (le reste)
    smbr(iel) = rovsdt(iel)*rtpa(iel,ivar)  - thets*tuexpr
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
! 3. Lagrangian source terms
!===============================================================================

!     Ordre 2 non pris en compte
 if (iilagr.eq.2 .and. ltsdyn.eq.1) then
   do iel = 1,ncel
     smbr(iel)   = smbr(iel)   + tslage(iel)
     rovsdt(iel) = rovsdt(iel) + max(-tslagi(iel),zero)
   enddo
 endif

!===============================================================================
! 4. Mass source term
!===============================================================================

if (ncesmp.gt.0) then

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

! On incremente SMBR par -Gamma RTPA et ROVSDT par Gamma (*theta)
  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isto2t , thetv  ,          &
   icetsm , itypsm(:,ivar)  ,                                     &
   volume , rtpa(:,ivar)    , smacel(:,ivar)   , smacel(:,ipr) ,  &
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
! 5. Non-stationary term
!===============================================================================

! ---> Ajout dans la diagonale de la matrice

do iel=1,ncel
  rovsdt(iel) = rovsdt(iel)                                       &
            + istat(ivar)*(crom(iel)/dt(iel))*volume(iel)
enddo


!===============================================================================
! 6. Production, Pressure-Strain correlation, dissipation, Coriolis
!===============================================================================

! ---> Terme source
!     -rho*epsilon*( Cs1*aij + Cs2*(aikajk -1/3*aijaij*deltaij))
!     -Cr1*P*aij + Cr2*rho*k*sij - Cr3*rho*k*sij*sqrt(aijaij)
!     +Cr4*rho*k(aik*sjk+ajk*sik-2/3*akl*skl*deltaij)
!     +Cr5*rho*k*(aik*rjk + ajk*rik)
!     -2/3*epsilon*deltaij

if(isou.eq.1)then
  iii = 1
  jjj = 1
elseif(isou.eq.2)then
  iii = 2
  jjj = 2
elseif(isou.eq.3)then
  iii = 3
  jjj = 3
elseif(isou.eq.4)then
  iii = 1
  jjj = 2
elseif(isou.eq.5)then
  iii = 1
  jjj = 3
elseif(isou.eq.6)then
  iii = 2
  jjj = 3
endif

! EBRSM
if (iturb.eq.32) then
  allocate(grad(3,ncelet))

  ! Compute the gradient of Alpha
  iprev  = 1
  inc    = 1
  iccocg = 1

  call field_gradient_scalar(ivarfl(ial), iprev, imrgra, inc,     &
                             iccocg,                              &
                             grad)

endif

if (icorio.eq.1 .or. iturbo.eq.1) then

  ! Compute the rotation matrix (dual matrix of the rotation vector)
  matrot(1,2) = -rotax(3)
  matrot(1,3) =  rotax(2)
  matrot(2,3) = -rotax(1)

  do ii = 1, 3
    matrot(ii,ii) = 0.d0
    do jj = ii+1, 3
      matrot(jj,ii) = -matrot(ii,jj)
    enddo
  enddo

else
  do ii = 1, 3
    do jj = 1, 3
      matrot(ii,jj) = 0.d0
    enddo
  enddo
endif

! Index of the Reynolds stress variables in rtpa array
indrey(1,1) = ir11
indrey(2,2) = ir22
indrey(3,3) = ir33
indrey(1,2) = ir12
indrey(1,3) = ir13
indrey(2,3) = ir23
indrey(2,1) = indrey(1,2)
indrey(3,1) = indrey(1,3)
indrey(3,2) = indrey(2,3)

do iel=1,ncel

  ! EBRSM
  if (iturb.eq.32) then
    ! Compute the magnitude of the Alpha gradient
    xnoral = ( grad(1,iel)*grad(1,iel)          &
           +   grad(2,iel)*grad(2,iel)          &
           +   grad(3,iel)*grad(3,iel) )
    xnoral = sqrt(xnoral)
   ! Compute the unitary vector of Alpha
    if (xnoral.le.epzero) then
      xnal(1) = 0.d0
      xnal(2) = 0.d0
      xnal(3) = 0.d0
    else
      xnal(1) = grad(1,iel)/xnoral
      xnal(2) = grad(2,iel)/xnoral
      xnal(3) = grad(3,iel)/xnoral
    endif
  endif

  ! Pij
  xprod(1,1) = -2.0d0*(rtpa(iel,ir11)*gradv(1, 1, iel) +         &
                       rtpa(iel,ir12)*gradv(2, 1, iel) +         &
                       rtpa(iel,ir13)*gradv(3, 1, iel) )
  xprod(1,2) = -(      rtpa(iel,ir11)*gradv(1, 2, iel) +         &
                       rtpa(iel,ir12)*gradv(2, 2, iel) +         &
                       rtpa(iel,ir13)*gradv(3, 2, iel) )         &
               -(      rtpa(iel,ir12)*gradv(1, 1, iel) +         &
                       rtpa(iel,ir22)*gradv(2, 1, iel) +         &
                       rtpa(iel,ir23)*gradv(3, 1, iel) )
  xprod(1,3) = -(      rtpa(iel,ir11)*gradv(1, 3, iel) +         &
                       rtpa(iel,ir12)*gradv(2, 3, iel) +         &
                       rtpa(iel,ir13)*gradv(3, 3, iel) )         &
               -(      rtpa(iel,ir13)*gradv(1, 1, iel) +         &
                       rtpa(iel,ir23)*gradv(2, 1, iel) +         &
                       rtpa(iel,ir33)*gradv(3, 1, iel) )
  xprod(2,2) = -2.0d0*(rtpa(iel,ir12)*gradv(1, 2, iel) +         &
                       rtpa(iel,ir22)*gradv(2, 2, iel) +         &
                       rtpa(iel,ir23)*gradv(3, 2, iel) )
  xprod(2,3) = -(      rtpa(iel,ir12)*gradv(1, 3, iel) +         &
                       rtpa(iel,ir22)*gradv(2, 3, iel) +         &
                       rtpa(iel,ir23)*gradv(3, 3, iel) )         &
               -(      rtpa(iel,ir13)*gradv(1, 2, iel) +         &
                       rtpa(iel,ir23)*gradv(2, 2, iel) +         &
                       rtpa(iel,ir33)*gradv(3, 2, iel) )
  xprod(3,3) = -2.0d0*(rtpa(iel,ir13)*gradv(1, 3, iel) +         &
                       rtpa(iel,ir23)*gradv(2, 3, iel) +         &
                       rtpa(iel,ir33)*gradv(3, 3, iel) )

  ! Rotating frame of reference => "Coriolis production" term
  if (icorio.eq.1 .or. iturbo.eq.1) then
    if (irotce(iel).gt.0) then
      do ii = 1, 3
        do jj = ii, 3
          do kk = 1, 3
            xprod(ii,jj) = xprod(ii,jj)                             &
                 - ccorio*( matrot(ii,kk)*rtpa(iel,indrey(jj,kk))   &
                 + matrot(jj,kk)*rtpa(iel,indrey(ii,kk)) )
          enddo
        enddo
      enddo
    endif
  endif

  xprod(2,1) = xprod(1,2)
  xprod(3,1) = xprod(1,3)
  xprod(3,2) = xprod(2,3)

  trprod = d1s2 * (xprod(1,1) + xprod(2,2) + xprod(3,3) )
  trrij  = d1s2 * (rtpa(iel,ir11) + rtpa(iel,ir22) + rtpa(iel,ir33))
!-----> aII = aijaij
  aii    = 0.d0
  aklskl = 0.d0
  aiksjk = 0.d0
  aikrjk = 0.d0
  aikakj = 0.d0
  ! aij
  xaniso(1,1) = rtpa(iel,ir11)/trrij - d2s3
  xaniso(2,2) = rtpa(iel,ir22)/trrij - d2s3
  xaniso(3,3) = rtpa(iel,ir33)/trrij - d2s3
  xaniso(1,2) = rtpa(iel,ir12)/trrij
  xaniso(1,3) = rtpa(iel,ir13)/trrij
  xaniso(2,3) = rtpa(iel,ir23)/trrij
  xaniso(2,1) = xaniso(1,2)
  xaniso(3,1) = xaniso(1,3)
  xaniso(3,2) = xaniso(2,3)
  ! Sij
  xstrai(1,1) = gradv(1, 1, iel)
  xstrai(1,2) = d1s2*(gradv(2, 1, iel)+gradv(1, 2, iel))
  xstrai(1,3) = d1s2*(gradv(3, 1, iel)+gradv(1, 3, iel))
  xstrai(2,1) = xstrai(1,2)
  xstrai(2,2) = gradv(2, 2, iel)
  xstrai(2,3) = d1s2*(gradv(3, 2, iel)+gradv(2, 3, iel))
  xstrai(3,1) = xstrai(1,3)
  xstrai(3,2) = xstrai(2,3)
  xstrai(3,3) = gradv(3, 3, iel)
  ! omegaij
  xrotac(1,1) = 0.d0
  xrotac(1,2) = d1s2*(gradv(2, 1, iel)-gradv(1, 2, iel))
  xrotac(1,3) = d1s2*(gradv(3, 1, iel)-gradv(1, 3, iel))
  xrotac(2,1) = -xrotac(1,2)
  xrotac(2,2) = 0.d0
  xrotac(2,3) = d1s2*(gradv(3, 2, iel)-gradv(2, 3, iel))
  xrotac(3,1) = -xrotac(1,3)
  xrotac(3,2) = -xrotac(2,3)
  xrotac(3,3) = 0.d0

  ! Rotating frame of reference => "absolute" vorticity
  if (icorio.eq.1) then
    do ii = 1, 3
      do jj = 1, 3
        xrotac(ii,jj) = xrotac(ii,jj) + matrot(ii,jj)
      enddo
    enddo
  endif

  do ii=1,3
    do jj = 1,3
      ! aii = aij.aij
      aii    = aii+xaniso(ii,jj)*xaniso(ii,jj)
      ! aklskl = aij.Sij
      aklskl = aklskl + xaniso(ii,jj)*xstrai(ii,jj)
    enddo
  enddo

  do kk = 1,3
    ! aiksjk = aik.Sjk+ajk.Sik
    aiksjk = aiksjk + xaniso(iii,kk)*xstrai(jjj,kk)              &
              +xaniso(jjj,kk)*xstrai(iii,kk)
    ! aikrjk = aik.Omega_jk + ajk.omega_ik
    aikrjk = aikrjk + xaniso(iii,kk)*xrotac(jjj,kk)              &
              +xaniso(jjj,kk)*xrotac(iii,kk)
    ! aikakj = aik*akj
    aikakj = aikakj + xaniso(iii,kk)*xaniso(kk,jjj)
  enddo

!     Si on extrapole les TS (rarissime), on met tout dans PROPCE.
!     On n'implicite pas le terme en Cs1*aij ni le terme en Cr1*P*aij.
!     Sinon, on met tout dans SMBR et on peut impliciter Cs1*aij
!     et Cr1*P*aij. Ici on stocke le second membre et le terme implicite
!     dans W1 et W2, pour eviter d'avoir un test IF(ISTO2T.GT.0)
!     dans la boucle NCEL
!     Dans le terme en W1, qui a vocation a etre extrapole, on utilise
!     naturellement CROMO.
!     L'implicitation des deux termes pourrait se faire aussi en cas
!     d'extrapolation, en isolant ces deux termes et les mettant dans
!     SMBR et pas PROPCE et en utilisant IPCROM ... a modifier si le
!     besoin s'en fait vraiment sentir           !

  if (iturb.eq.31) then

    pij = xprod(iii,jjj)
    phiij1 = -rtpa(iel,iep)* &
       (cssgs1*xaniso(iii,jjj)+cssgs2*(aikakj-d1s3*deltij*aii))
    phiij2 = - cssgr1*trprod*xaniso(iii,jjj)                             &
           +   trrij*xstrai(iii,jjj)*(cssgr2-cssgr3*sqrt(aii))           &
           +   cssgr4*trrij*(aiksjk-d2s3*deltij*aklskl)                  &
           +   cssgr5*trrij* aikrjk
    epsij = -d2s3*rtpa(iel,iep)*deltij

    w1(iel) = cromo(iel)*volume(iel)*(pij+phiij1+phiij2+epsij)

    w2(iel) = volume(iel)/trrij*crom(iel)*(                &
           cssgs1*rtpa(iel,iep) + cssgr1*max(trprod,0.d0) )

  ! EBRSM
  else

    xrij(1,1) = rtpa(iel,ir11)
    xrij(2,2) = rtpa(iel,ir22)
    xrij(3,3) = rtpa(iel,ir33)
    xrij(1,2) = rtpa(iel,ir12)
    xrij(1,3) = rtpa(iel,ir13)
    xrij(2,3) = rtpa(iel,ir23)
    xrij(2,1) = xrij(1,2)
    xrij(3,1) = xrij(1,3)
    xrij(3,2) = xrij(2,3)

    ! Compute the explicit term

    ! Calcul des termes de proches parois et quasi-homgene de phi et
    ! epsilon

    ! Calcul du terme de proche paroi \Phi_{ij}^w --> W3
    phiijw = 0.d0
    xnnd = d1s2*( xnal(iii)*xnal(jjj) + deltij )
    do kk = 1, 3
      phiijw = phiijw + xrij(iii,kk)*xnal(jjj)*xnal(kk)
      phiijw = phiijw + xrij(jjj,kk)*xnal(iii)*xnal(kk)
      do ll = 1, 3
        phiijw = phiijw - xrij(kk,ll)*xnal(kk)*xnal(ll)*xnnd
      enddo
    enddo
    phiijw = -5.d0*rtpa(iel,iep)/trrij * phiijw

    ! Calcul du terme quasi-homogene \Phi_{ij}^h --> W4
    phiij1 = -rtpa(iel,iep)*cebms1*xaniso(iii,jjj)
    phiij2 = -cebmr1*trprod*xaniso(iii,jjj)                       &
               +trrij*xstrai(iii,jjj)*(cebmr2-cebmr3*sqrt(aii))   &
               +cebmr4*trrij   *(aiksjk-d2s3*deltij*aklskl)       &
               +cebmr5*trrij   * aikrjk

    ! Calcul de \e_{ij}^w --> W5 (Rotta model)
    ! Rij/k*epsilon
    epsijw =  xrij(iii,jjj)/trrij   *rtpa(iel,iep)

    ! Calcul de \e_{ij}^h --> W6
    epsij =  d2s3*rtpa(iel,iep)*deltij

    ! Calcul du terme source explicite de l'equation des Rij
    !   [ P_{ij} + (1-\alpha^3)\Phi_{ij}^w + \alpha^3\Phi_{ij}^h
    !            - (1-\alpha^3)\e_{ij}^w   - \alpha^3\e_{ij}^h  ] --> W1
    alpha3 = rtp(iel,ial)**3

    w1(iel) = volume(iel)*crom(iel)*(                    &
               xprod(iii,jjj)                                     &
            + (1.d0-alpha3)*phiijw + alpha3*(phiij1+phiij2)       &
            - (1.d0-alpha3)*epsijw - alpha3*epsij)

    !  Implicite term

    ! le terme ci-dessous correspond a la partie implicitee du SSG
    ! dans le cadre de la ponderation elliptique, il est multiplie par
    ! \alpha^3
    w2(iel) = volume(iel)*crom(iel)*(                    &
              cebms1*rtpa(iel,iep)/trrij*alpha3                   &
             +cebmr1*max(trprod/trrij,0.d0)*alpha3                &
    ! Implicitation de epsijw
    ! (le facteur 5 apparait lorsqu'on fait Phi_{ij}^w - epsijw)
            + 5.d0 * (1.d0-alpha3)*rtpa(iel,iep)/trrij            &
            +        (1.d0-alpha3)*rtpa(iel,iep)/trrij)
  endif

enddo

if (iturb.eq.32) then
  deallocate(grad)
endif

if (isto2t.gt.0) then

  do iel = 1, ncel
    propce(iel,iptsta+isou-1) = propce(iel,iptsta+isou-1) + w1(iel)
  enddo

else

  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w1(iel)
    rovsdt(iel) = rovsdt(iel) + w2(iel)
  enddo

endif

!===============================================================================
! 7. Buoyancy source term
!===============================================================================

if (igrari.eq.1) then

  ! Allocate a work array
  allocate(w7(ncelet))

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  call rijthe(nscal, ivar, rtpa, gradro, w7)
  !==========

  ! Si on extrapole les T.S. : PROPCE
  if(isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) = propce(iel,iptsta+isou-1) + w7(iel)
    enddo
  ! Sinon SMBR
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w7(iel)
    enddo
  endif

  ! Free memory
  deallocate(w7)

endif

!===============================================================================
! 8. Diffusion term (Daly Harlow: generalized gradient hypothesis method)
!===============================================================================

! Symmetric tensor diffusivity (GGDH)
if (idften(ivar).eq.6) then

  do iel = 1, ncel
    viscce(1,iel) = visten(1,iel) + propce(iel,ipcvis)
    viscce(2,iel) = visten(2,iel) + propce(iel,ipcvis)
    viscce(3,iel) = visten(3,iel) + propce(iel,ipcvis)
    viscce(4,iel) = visten(4,iel)
    viscce(5,iel) = visten(5,iel)
    viscce(6,iel) = visten(6,iel)
  enddo

  iwarnp = iwarni(ivar)

  call vitens &
  !==========
 ( viscce , iwarnp ,             &
   weighf , weighb ,             &
   viscf  , viscb  )

else
  call csexit(1)
endif

!===============================================================================
! 9. Solving
!===============================================================================

if (isto2t.gt.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + thetp1*propce(iel,iptsta+isou-1)
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
imucpp = 0
idftnp = idften(ivar)
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
! all boundary convective flux with upwind
icvflb = 0

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

call codits &
!==========
 ( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscce , viscf  , viscb  , viscce ,          &
   weighf , weighb ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbr   , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , rvoid  )

! Free memory
deallocate(w1, w2)
deallocate(dpvar)
deallocate(viscce)
deallocate(weighf, weighb)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(/,'           RESOLUTION POUR LA VARIABLE ',A8,/)

#else

 1000 format(/,'           SOLVING VARIABLE ',A8           ,/)

#endif

!----
! End
!----

return

end subroutine
