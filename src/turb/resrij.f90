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

!===============================================================================
! Function:
! ---------

!> \file resrij.f90
!>
!> \brief This subroutine performs the solving of the Reynolds stress components
!> in \f$ R_{ij} - \varepsilon \f$ RANS (LRR) turbulence model.
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
!> \param[in]     itpsmp        type of mass source term for the variables
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rtp, rtpa     calculated variables at cell centers
!>                               (at current and previous time steps)
!> \param[in]     propce        physical properties at cell centers
!> \param[in]     coefa, coefb  boundary conditions
!> \param[in]     grdvit        tableau de travail pour terme grad
!>                                 de vitesse     uniqt pour iturb=31
!> \param[in]     produc        tableau de travail pour production
!> \param[in]     gradro        tableau de travail pour grad rom
!>                              (sans rho volume) uniqt pour iturb=30
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smcelp        variable value associated to the mass source
!>                               term
!> \param[in]     gamma         valeur du flux de masse
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     tslagr        coupling term for lagrangian
!> \param[in]     tslage        explicit source terms for the Lagrangian module
!> \param[in]     tslagi        implicit source terms for the Lagrangian module
!> \param[in]     smbr          working array
!> \param[in]     rovsdt        working array
!_______________________________________________________________________________

subroutine resrij &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   , isou   , ipp    ,                                     &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtp    , rtpa   , propce ,                            &
   coefa  , coefb  , produc , gradro ,                            &
   ckupdc , smcelp , gamma  ,                                     &
   viscf  , viscb  ,                                              &
   tslage , tslagi ,                                              &
   smbr   , rovsdt )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use numvar
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use lagran
use pointe, only:visten
use mesh
use field
use cs_f_interfaces

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
double precision coefa(ndimfb,*), coefb(ndimfb,*)
double precision produc(6,ncelet)
double precision gradro(ncelet,3)
double precision ckupdc(ncepdp,6)
double precision smcelp(ncesmp), gamma(ncesmp)
double precision viscf(nfac), viscb(nfabor)
double precision tslage(ncelet),tslagi(ncelet)
double precision smbr(ncelet), rovsdt(ncelet)

! Local variables

integer          iel
integer          ii    , jj    , kk    , iiun
integer          ipcvis, iflmas, iflmab
integer          iclvar, iclvaf
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          iptsta
integer          isoluc
integer          imucpp, idftnp, iswdyp
integer          indrey(3,3)
integer          icvflb
integer          ivoid(1)

double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision trprod, trrij , deltij
double precision tuexpr, thets , thetv , thetp1
double precision d1s3  , d2s3
double precision matrot(3,3)

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w7, w8
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer ::  crom, cromo

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet))
allocate(w7(ncelet), w8(ncelet))
allocate(dpvar(ncelet))
allocate(viscce(6,ncelet))
allocate(weighf(2,nfac))
allocate(weighb(nfabor))

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) nomvar(ipp)
endif

call field_get_val_s(icrom, crom)
ipcvis = ipproc(iviscl)
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

iclvar = iclrtp(ivar,icoef)
iclvaf = iclrtp(ivar,icoeff)

deltij = 1.0d0
if(isou.gt.3) then
  deltij = 0.0d0
endif
d1s3 = 1.d0/3.d0
d2s3 = 2.d0/3.d0

!     S pour Source, V pour Variable
thets  = thetst
thetv  = thetav(ivar )

if (isto2t.gt.0.and.iroext.gt.0) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
endif
iptsta = 0
if(isto2t.gt.0) then
  iptsta = ipproc(itstua)
endif

do iel = 1, ncel
  smbr(iel) = 0.d0
enddo
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

!===============================================================================
! 2. User source terms
!===============================================================================
!(le premier argument PRODUC est lu en GRDVIT dans ustsri, mais ce
! tableau n'est dimensionne et utilise qu'en modele Rij SSG)

call ustsri &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smcelp , gamma  , produc , produc ,                   &
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

!       On incremente SMBR par -Gamma RTPA et ROVSDT par Gamma (*theta)
  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isto2t , thetv  ,          &
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
! 5. Non-stationary term
!===============================================================================

do iel=1,ncel
  rovsdt(iel) = rovsdt(iel)                                          &
              + istat(ivar)*(crom(iel)/dt(iel))*volume(iel)
enddo


!===============================================================================
! 6. Production, Pressure-Strain correlation, dissipation
!===============================================================================

! ---> Calcul de k pour la suite du sous-programme
!       on utilise un tableau de travail puisqu'il y en a...
do iel = 1, ncel
  w8(iel) = 0.5d0 * (rtpa(iel,ir11) + rtpa(iel,ir22) + rtpa(iel,ir33))
enddo

! ---> Terme source

!      (1-CRIJ2) Pij (pour toutes les composantes de Rij)

!      DELTAIJ*(2/3.CRIJ2.P+2/3.CRIJ1.EPSILON)
!                    (termes diagonaux pour R11, R22 et R33)

!      -DELTAIJ*2/3*EPSILON

!     Si on extrapole les TS
!       On modifie la partie implicite :
!         Dans PHI1, on ne prendra que RHO CRIJ1 EPSILON/K et non pas
!                                  RHO CRIJ1 EPSILON/K (1-2/3 DELTAIJ)
!         Cela permet de conserver k^n au lieu de (R11^(n+1)+R22^n+R33^n)
!         Ce choix est discutable. C'est la solution ISOLUC = 1
!       Si on veut tout prendre en implicite (comme c'est fait
!         en ordre 1 std), c'est la solution ISOLUC = 2
!       -> a tester plus precisement si necessaire


!     Si on extrapole les TS
if (isto2t.gt.0) then

  isoluc = 1

  do iel = 1, ncel

!     Demi-traces de Prod et R
    trprod = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
    trrij  = w8(iel)

!     Calcul de Prod+Phi1+Phi2-Eps
!       = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
!       Dans PROPCE :
!       = rhoPij-C1rho eps/k(   -2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
!       = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij           }
    propce(iel,iptsta+isou-1) = propce(iel,iptsta+isou-1)         &
                          + cromo(iel) * volume(iel)              &
      *(   deltij*d2s3*                                           &
           (  crij2*trprod                                        &
            +(crij1-1.d0)* rtpa(iel,iep)  )                     &
         +(1.0d0-crij2)*produc(isou,iel)               )
!       Dans SMBR
!       =       -C1rho eps/k(Rij         )
!       = rho{                                     -C1eps/kRij}
    smbr(iel) = smbr(iel) + crom(iel) * volume(iel)      &
      *( -crij1*rtpa(iel,iep)/trrij * rtpa(iel,ivar)  )

!     Calcul de la partie implicite issue de Phi1
!       = C1rho eps/k(1        )
    rovsdt(iel) = rovsdt(iel) + crom(iel) * volume(iel)  &
                            *crij1*rtpa(iel,iep)/trrij*thetv

  enddo

!     Si on veut impliciter un bout de -C1rho eps/k(   -2/3k dij)
  if(isoluc.eq.2) then

    do iel = 1, ncel

      trrij  = w8(iel)

!    On enleve a CROMO
!       =       -C1rho eps/k(   -1/3Rij dij)
      propce(iel,iptsta+isou-1) = propce(iel,iptsta+isou-1)       &
                          - cromo(iel) * volume(iel)              &
      *(deltij*d1s3*crij1*rtpa(iel,iep)/trrij * rtpa(iel,ivar))
!    On ajoute a SMBR (avec CROM)
!       =       -C1rho eps/k(   -1/3Rij dij)
      smbr(iel)                 = smbr(iel)                       &
                          + crom(iel) * volume(iel)      &
      *(deltij*d1s3*crij1*rtpa(iel,iep)/trrij * rtpa(iel,ivar))
!    On ajoute a ROVSDT (avec CROM)
!       =        C1rho eps/k(   -1/3    dij)
      rovsdt(iel) = rovsdt(iel) + crom(iel) * volume(iel)&
      *(deltij*d1s3*crij1*rtpa(iel,iep)/trrij                 )
    enddo

  endif

! Si on n'extrapole pas les termes sources
else

  do iel = 1, ncel

!     Demi-traces de Prod et R
    trprod = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
    trrij  = w8(iel)

!     Calcul de Prod+Phi1+Phi2-Eps
!       = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
!       = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij-C1eps/kRij}
    smbr(iel) = smbr(iel) + crom(iel) * volume(iel)      &
      *(   deltij*d2s3*                                           &
           (  crij2*trprod                                        &
            +(crij1-1.d0)* rtpa(iel,iep)  )                     &
         +(1.0d0-crij2)*produc(isou,iel)                          &
         -crij1*rtpa(iel,iep)/trrij * rtpa(iel,ivar)  )

!     Calcul de la partie implicite issue de Phi1
!       = C1rho eps/k(1-1/3 dij)
    rovsdt(iel) = rovsdt(iel) + crom(iel) * volume(iel)  &
         *(1.d0-d1s3*deltij)*crij1*rtpa(iel,iep)/trrij
  enddo

endif

!===============================================================================
! 6-bis. Coriolis terms in the Phi1 and production
!===============================================================================

if (icorio.eq.1) then

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  ! Rotation matrix: dual antisymmetric matrix of the rotation vector omega
  matrot(1,2) = -omegaz
  matrot(1,3) =  omegay
  matrot(2,3) = -omegax

  do ii = 1, 3
    matrot(ii,ii) = 0.d0
    do jj = ii+1, 3
      matrot(jj,ii) = -matrot(ii,jj)
    enddo
  enddo

  ! Index Connectivity
  indrey(1,1) = ir11
  indrey(2,2) = ir22
  indrey(3,3) = ir33
  indrey(1,2) = ir12
  indrey(1,3) = ir13
  indrey(2,3) = ir23
  indrey(2,1) = indrey(1,2)
  indrey(3,1) = indrey(1,3)
  indrey(3,2) = indrey(2,3)

  if (isou.eq.1) then
    ii = 1
    jj = 1
  elseif (isou.eq.2) then
    ii = 2
    jj = 2
  elseif (isou.eq.3) then
    ii = 3
    jj = 3
  elseif (isou.eq.4) then
    ii = 1
    jj = 2
  elseif (isou.eq.5) then
    ii = 1
    jj = 3
  elseif (isou.eq.6) then
    ii = 2
    jj = 3
  endif

  do iel = 1, ncel
    ! Compute Gij: (i,j) component of the Coriolis production
    do kk = 1, 3
      w7(iel) = w7(iel) - 2.d0*( matrot(ii,kk)*rtpa(iel,indrey(jj,kk)) &
                               + matrot(jj,kk)*rtpa(iel,indrey(ii,kk)) )
    enddo
    ! Coriolis contribution in the Phi1 term:
    ! (1-C2/2)Gij
    w7(iel) = crom(iel) * volume(iel)                   &
            * (1.d0 - 0.5d0*crij2)*w7(iel)
  enddo

  ! If source terms are extrapolated
  if (isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) =                                 &
      propce(iel,iptsta+isou-1) + w7(iel)
    enddo
  ! Otherwise, directly in smbr
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w7(iel)
    enddo
  endif

endif

!===============================================================================
! 7. Wall echo terms
!===============================================================================

if (irijec.eq.1) then

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  call rijech(isou, rtpa, produc, w7)
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

endif


!===============================================================================
! 8. Buoyancy source term
!===============================================================================

if (igrari.eq.1) then

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  call rijthe(nscal, ivar, rtpa, gradro, w7)
  !==========

  ! If source terms are extrapolated
  if (isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) =                                  &
      propce(iel,iptsta+isou-1) + w7(iel)
    enddo
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w7(iel)
    enddo
  endif

endif

!===============================================================================
! 9. Diffusion term (Daly Harlow: generalized gradient hypothesis method)
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
! 10. Solving
!===============================================================================

if (isto2t.gt.0) then
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

call codits &
!==========
 ( idtvar , ivar   , iconvp , idiffp , ireslp , ndircp , nitmap , &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   imgrp  , ncymxp , nitmfp , ipp    , iwarnp ,                   &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   rtpa(1,ivar)    , rtpa(1,ivar)    ,                            &
   coefa(1,iclvar) , coefb(1,iclvar) ,                            &
   coefa(1,iclvaf) , coefb(1,iclvaf) ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscce , viscf  , viscb  , viscce ,          &
   weighf , weighb ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbr   , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , rvoid  )

! Free memory
deallocate(w1)
deallocate(w7, w8)
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
