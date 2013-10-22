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

!> \file reseps.f90
!>
!> \brief This subroutine performs the solving of epsilon in
!> \f$ R_{ij} - \varepsilon \f$ RANS turbulence model.
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
!> \param[in]     gradv         tableau de travail pour terme grad
!>                                 de vitesse     uniqt pour iturb=31
!> \param[in]     produc        tableau de travail pour production
!>                              (sans rho volume) uniqt pour iturb=30
!> \param[in]     gradro        tableau de travail pour grad rom
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smcelp        variable value associated to the mass source
!>                               term
!> \param[in]     gamma         valeur du flux de masse
!> \param[in]     viscf         visc*surface/dist aux faces internes
!> \param[in]     viscb         visc*surface/dist aux faces de bord
!> \param[in]     tslagr        coupling term for lagrangian
!> \param[in]     smbr          working array
!> \param[in]     rovsdt        working array
!_______________________________________________________________________________

subroutine reseps &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   , isou   , ipp    ,                                     &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtp    , rtpa   , propce ,                            &
   gradv  , produc , gradro ,                                     &
   ckupdc , smcelp , gamma  ,                                     &
   viscf  , viscb  ,                                              &
   tslagr ,                                                       &
   smbr   , rovsdt )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstnum
use cstphy
use parall
use period
use lagran
use pointe, only:visten
use mesh
use field
use field_operator
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
double precision produc(6,ncelet), gradv(3, 3, ncelet)
double precision gradro(ncelet,3)
double precision ckupdc(ncepdp,6)
double precision smcelp(ncesmp), gamma(ncesmp)
double precision viscf(nfac), viscb(nfabor)
double precision tslagr(ncelet,*)
double precision smbr(ncelet), rovsdt(ncelet)

! Local variables

integer          iel
integer          iiun
integer          ipcvis, ipcvst, iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp, ireslp
integer          nitmap, nswrsp, ircflp, ischcp, isstpp, iescap
integer          imgrp , ncymxp, nitmfp
integer          iptsta
integer          imucpp, idftnp, iswdyp
integer          icvflb
integer          ivoid(1)
double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp, alpha3
double precision trprod , trrij
double precision tseps , kseps , ceps2
double precision tuexpe, thets , thetv , thetap, thetp1
double precision prdeps, xttdrb, xttke , xttkmg

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w7, w9
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer ::  crom, cromo
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp

!===============================================================================

!===============================================================================
! 1. Initialisation
!===============================================================================

! Allocate work arrays
allocate(w1(ncelet))
allocate(w9(ncelet))
allocate(dpvar(ncelet))
allocate(viscce(6,ncelet))
allocate(weighf(2,nfac))
allocate(weighb(nfabor))

if(iwarni(ivar).ge.1) then
  write(nfecra,1000) nomvar(ipp)
endif

call field_get_val_s(icrom, crom)
ipcvis = ipproc(iviscl)
ipcvst = ipproc(ivisct)
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

! Constante Ce2, qui vaut CE2 pour ITURB=30 et CSSGE2 pour ITRUB=31
if (iturb.eq.30) then
  ceps2 = ce2
elseif (iturb.eq.31) then
  ceps2 = cssge2
else
  ceps2 = cebme2
endif

! S pour Source, V pour Variable
thets  = thetst
thetv  = thetav(ivar )

if (isto2t.gt.0.and.iroext.gt.0) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
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
! 2. User source terms
!===============================================================================

call ustsri                                                       &
!==========
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   ,                                                       &
   icepdc , icetsm , itpsmp ,                                     &
   dt     , rtpa   , propce ,                                     &
   ckupdc , smcelp , gamma  , gradv  , produc ,                   &
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
! 3. Lagrangian source terms
!===============================================================================

!     Ordre 2 non pris en compte
if (iilagr.eq.2 .and. ltsdyn.eq.1) then

  do iel = 1, ncel
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
! 4. Mass source term
!===============================================================================

if (ncesmp.gt.0) then

!       Entier egal a 1 (pour navsto : nb de sur-iter)
  iiun = 1

!       On incremente SMBR par -Gamma RTPA et ROVSDT par Gamma (*theta)
  call catsma &
  !==========
 ( ncelet , ncel   , ncesmp , iiun   , isto2t , thetv ,           &
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

do iel = 1, ncel
  rovsdt(iel) = rovsdt(iel)                                          &
              + istat(ivar)*(crom(iel)/dt(iel))*volume(iel)
enddo

!===============================================================================
! 6. Production (rho * Ce1 * epsilon / k * P)
!    Dissipation (rho*Ce2.epsilon/k*epsilon)
!===============================================================================

if (isto2t.gt.0) then
  thetap = thetv
else
  thetap = 1.d0
endif

! ---> Calcul de la trace de la production, suivant qu'on est en
!     Rij standard ou en SSG (utilisation de PRODUC ou GRDVIT)
if (iturb.eq.30) then
  do iel = 1, ncel
    w9(iel) = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
  enddo
else
  do iel = 1, ncel
    w9(iel) = -( rtpa(iel,ir11)*gradv(1, 1, iel) +               &
                 rtpa(iel,ir12)*gradv(2, 1, iel) +               &
                 rtpa(iel,ir13)*gradv(3, 1, iel) +               &
                 rtpa(iel,ir12)*gradv(1, 2, iel) +               &
                 rtpa(iel,ir22)*gradv(2, 2, iel) +               &
                 rtpa(iel,ir23)*gradv(3, 2, iel) +               &
                 rtpa(iel,ir13)*gradv(1, 3, iel) +               &
                 rtpa(iel,ir23)*gradv(2, 3, iel) +               &
                 rtpa(iel,ir33)*gradv(3, 3, iel) )
  enddo
endif


! EBRSM
if (iturb.eq.32) then

  do iel = 1, ncel
    ! Demi-traces
    trprod = w9(iel)
    trrij  = 0.5d0 * (rtpa(iel,ir11) + rtpa(iel,ir22) + rtpa(iel,ir33))
    ! Calcul de l echelle de temps de Durbin
    xttke  = trrij/rtpa(iel,iep)
    xttkmg = xct*sqrt(propce(iel,ipcvis)/crom(iel)           &
                                        /rtpa(iel,iep))
    xttdrb = max(xttke,xttkmg)

    prdeps = trprod/rtpa(iel,iep)
    alpha3 = rtp(iel,ial)**3

    ! Production (explicit)
    ! Compute of C_eps_1'
    w1(iel) = cromo(iel)*volume(iel)*                         &
              ce1*(1.d0+xa1*(1.d0-alpha3)*prdeps)*trprod/xttdrb


    ! Dissipation (implicit)
    smbr(iel) = smbr(iel) - crom(iel)*volume(iel)*           &
                             ceps2*rtpa(iel,iep)/xttdrb

    rovsdt(iel) = rovsdt(iel)                                         &
                + ceps2*crom(iel)*volume(iel)*thetap/xttdrb
  enddo

! SSG and LRR
else

  do iel = 1, ncel
    ! Demi-traces
    trprod = w9(iel)
    trrij  = 0.5d0 * (rtpa(iel,ir11) + rtpa(iel,ir22) + rtpa(iel,ir33))
    xttke  = trrij/rtpa(iel,iep)
    ! Production (explicit)
    w1(iel) = cromo(iel)*volume(iel)*ce1/xttke*trprod

    ! Dissipation (implicit)
    smbr(iel) = smbr(iel)                              &
              - crom(iel)*volume(iel)*ceps2*rtpa(iel,iep)**2/trrij
    rovsdt(iel) = rovsdt(iel)                                        &
                + ceps2*crom(iel)*volume(iel)/xttke*thetap
  enddo

endif

! Extrapolation of source terms (2nd order in time)
if (isto2t.gt.0) then
  do iel = 1, ncel
    propce(iel,iptsta+isou-1) = propce(iel,iptsta+isou-1) + w1(iel)
  enddo
else
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + w1(iel)
  enddo
endif

!===============================================================================
! 7. Buoyancy term
!===============================================================================

if (igrari.eq.1) then

  ! Allocate a work array
  allocate(w7(ncelet))

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  call rijthe(nscal, ivar, rtpa, gradro, w7)
  !==========

  ! Extrapolation of source terms (2nd order in time)
  if (isto2t.gt.0) then
    do iel = 1, ncel
      propce(iel,iptsta+isou-1) = propce(iel,iptsta+isou-1) + w7(iel)
    enddo
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
    viscce(1,iel) = visten(1,iel)/sigmae + propce(iel,ipcvis)
    viscce(2,iel) = visten(2,iel)/sigmae + propce(iel,ipcvis)
    viscce(3,iel) = visten(3,iel)/sigmae + propce(iel,ipcvis)
    viscce(4,iel) = visten(4,iel)/sigmae
    viscce(5,iel) = visten(5,iel)/sigmae
    viscce(6,iel) = visten(6,iel)/sigmae
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
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscce  , viscf  , viscb  , viscce  ,        &
   weighf , weighb ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbr   , rtp(1,ivar)     , dpvar  ,                   &
   rvoid  , rvoid  )

! Free memory
deallocate(w1)
deallocate(w9)
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
