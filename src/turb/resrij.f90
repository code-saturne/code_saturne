!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2017 EDF S.A.
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
!> \remark
!> - isou=1 for \f$ R_{11} \f$
!> - isou=2 for \f$ R_{22} \f$
!> - isou=3 for \f$ R_{33} \f$
!> - isou=4 for \f$ R_{12} \f$
!> - isou=5 for \f$ R_{13} \f$
!> - isou=6 for \f$ R_{23} \f$
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     ncepdp        number of cells with head loss
!> \param[in]     ncesmp        number of cells with mass source term
!> \param[in]     ivar          variable number
!> \param[in]     isou          local variable number (7 here)
!> \param[in]     icepdc        index of cells with head loss
!> \param[in]     icetsm        index of cells with mass source term
!> \param[in]     itypsm        type of mass source term for each variable
!>                               (see \ref cs_user_mass_source_terms)
!> \param[in]     dt            time step (per cell)
!> \param[in]     produc        work array for production
!> \param[in]     gradro        work array for grad rom
!>                              (without rho volume) only for iturb=30
!> \param[in]     ckupdc        work array for the head loss
!> \param[in]     smacel        value associated to each variable in the mass
!>                              source terms or mass rate (see
!>                              \ref cs_user_mass_source_terms)
!> \param[in]     viscf         visc*surface/dist at internal faces
!> \param[in]     viscb         visc*surface/dist at edge faces
!> \param[in]     tslagi        implicit source terms for the Lagrangian module
!> \param[in]     smbr          working array
!> \param[in]     rovsdt        working array
!______________________________________________________________________________!

subroutine resrij &
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivar   , isou   ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   produc , gradro ,                                              &
   ckupdc , smacel ,                                              &
   viscf  , viscb  ,                                              &
   tslagi ,                                                       &
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
use mesh
use field
use cs_f_interfaces
use rotation
use turbomachinery
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal
integer          ncepdp , ncesmp
integer          ivar   , isou

integer          icepdc(ncepdp)
integer          icetsm(ncesmp), itypsm(ncesmp,nvar)

double precision dt(ncelet)
double precision produc(6,ncelet)
double precision gradro(3,ncelet)
double precision ckupdc(ncepdp,6), smacel(ncesmp,nvar)
double precision viscf(nfac), viscb(nfabor)
double precision tslagi(ncelet)
double precision smbr(ncelet), rovsdt(ncelet)

! Local variables

integer          iel
integer          ii    , jj    , kk    , iiun, comp_id
integer          iflmas, iflmab
integer          nswrgp, imligp, iwarnp
integer          iconvp, idiffp, ndircp
integer          nswrsp, ircflp, ischcp, isstpp, iescap
integer          st_prv_id
integer          isoluc
integer          imucpp, idftnp, iswdyp
integer          ivar_r(3,3)
integer          icvflb
integer          ivoid(1)

double precision blencp, epsilp, epsrgp, climgp, extrap, relaxp
double precision epsrsp
double precision trprod, trrij , deltij
double precision tuexpr, thets , thetv , thetp1
double precision d1s3  , d2s3
double precision ccorio, matrot(3,3)
double precision rctse

double precision rvoid(1)

character(len=80) :: label
double precision, allocatable, dimension(:) :: w1
double precision, allocatable, dimension(:) :: w7, w8
double precision, allocatable, dimension(:) :: dpvar
double precision, allocatable, dimension(:,:) :: viscce
double precision, allocatable, dimension(:,:) :: weighf
double precision, allocatable, dimension(:) :: weighb
double precision, dimension(:), pointer :: imasfl, bmasfl
double precision, dimension(:), pointer :: crom, cromo
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:,:), pointer :: visten, lagr_st_rij
double precision, dimension(:), pointer :: cvara_ep
double precision, dimension(:), pointer :: cvara_r11, cvara_r22, cvara_r33
double precision, dimension(:), pointer :: cvar_var, cvara_var
double precision, dimension(:), pointer :: cvara_rik, cvara_rjk
double precision, dimension(:), pointer :: viscl, c_st_prv

type(var_cal_opt) :: vcopt

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

call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)

if (vcopt%iwarni.ge.1) then
  call field_get_label(ivarfl(ivar), label)
  write(nfecra,1000) label
endif

call field_get_val_s(icrom, crom)
call field_get_val_s(iviscl, viscl)
call field_get_key_int(ivarfl(iu), kimasf, iflmas)
call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
call field_get_val_s(iflmas, imasfl)
call field_get_val_s(iflmab, bmasfl)

call field_get_val_prev_s(ivarfl(iep), cvara_ep)

call field_get_val_prev_s(ivarfl(ir11), cvara_r11)
call field_get_val_prev_s(ivarfl(ir22), cvara_r22)
call field_get_val_prev_s(ivarfl(ir33), cvara_r33)

call field_get_val_s(ivarfl(ivar), cvar_var)
call field_get_val_prev_s(ivarfl(ivar), cvara_var)

call field_get_coefa_s(ivarfl(ivar), coefap)
call field_get_coefb_s(ivarfl(ivar), coefbp)
call field_get_coefaf_s(ivarfl(ivar), cofafp)
call field_get_coefbf_s(ivarfl(ivar), cofbfp)

deltij = 1.0d0
if (isou.gt.3) then
  deltij = 0.0d0
endif
d1s3 = 1.d0/3.d0
d2s3 = 2.d0/3.d0

!     S as Source, V as Variable
thets  = thetst
thetv  = vcopt%thetav

call field_get_key_int(ivarfl(ivar), kstprv, st_prv_id)
if (st_prv_id.ge.0) then
  call field_get_val_s(st_prv_id, c_st_prv)
else
  c_st_prv=> null()
endif

if (st_prv_id.ge.0.and.iroext.gt.0) then
  call field_get_val_prev_s(icrom, cromo)
else
  call field_get_val_s(icrom, cromo)
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
 ( nvar   , nscal  , ncepdp , ncesmp ,                            &
   ivarfl(ivar)    ,                                              &
   icepdc , icetsm , itypsm ,                                     &
   ckupdc , smacel ,                                              &
   smbr   , rovsdt )

!     If we extrapolate the source terms
if (st_prv_id.ge.0) then
  do iel = 1, ncel
!       Save for exchange
    tuexpr = c_st_prv(iel)
!       For continuation and next time step
    c_st_prv(iel) = smbr(iel)
!       Second member of previous time step
!       We suppose -rovsdt > 0: we implicite
!          the user source term (the rest)
    smbr(iel) = rovsdt(iel)*cvara_var(iel)  - thets*tuexpr
!       Diagonal
    rovsdt(iel) = - thetv*rovsdt(iel)
  enddo
else
  do iel = 1, ncel
    smbr(iel)   = rovsdt(iel)*cvara_var(iel) + smbr(iel)
    rovsdt(iel) = max(-rovsdt(iel),zero)
  enddo
endif

!===============================================================================
! 3. Lagrangian source terms
!===============================================================================

!     Second order is not taken into account
if (iilagr.eq.2 .and. ltsdyn.eq.1) then
  call field_get_val_v_by_name('rij_st_lagr', lagr_st_rij)
  comp_id = isou
  if (isou .eq. 5) then
    comp_id = 6
  else if (isou .eq.6) then
    comp_id = 5
  endif
  do iel = 1,ncel
    smbr(iel)   = smbr(iel)   + lagr_st_rij(comp_id,iel)
    rovsdt(iel) = rovsdt(iel) + max(-tslagi(iel),zero)
  enddo
endif

!===============================================================================
! 4. Mass source term
!===============================================================================

if (ncesmp.gt.0) then

!       Integer equal to 1 (for navsto: nb of sur-iter)
  iiun = 1

!       We increment smbr with -Gamma.var_prev. and rovsdt with Gamma
  call catsma &
 ( ncelet , ncel   , ncesmp , iiun   , isto2t ,                   &
   icetsm , itypsm(:,ivar)  ,                                     &
   cell_f_vol , cvara_var   , smacel(:,ivar)   , smacel(:,ipr) ,  &
   smbr   ,  rovsdt , w1 )

!       If we extrapolate the source terms we put Gamma Pinj in c_st_prv
  if (st_prv_id.ge.0) then
    do iel = 1, ncel
      c_st_prv(iel) = c_st_prv(iel) + w1(iel)
    enddo
!       Otherwise we put it directly in smbr
  else
    do iel = 1, ncel
      smbr(iel) = smbr(iel) + w1(iel)
    enddo
  endif

endif

!===============================================================================
! 5. Unsteady term
!===============================================================================

do iel=1,ncel
  rovsdt(iel) = rovsdt(iel)                                          &
              + vcopt%istat*(crom(iel)/dt(iel))*cell_f_vol(iel)
enddo


!===============================================================================
! 6. Production, Pressure-Strain correlation, dissipation
!===============================================================================

! ---> Calculation of k for the sub-routine continuation
!       we use a work array
do iel = 1, ncel
  w8(iel) = 0.5d0 * (cvara_r11(iel) + cvara_r22(iel) + cvara_r33(iel))
enddo

! ---> Source term

!      (1-CRIJ2) Pij (for all components of Rij)

!      DELTAIJ*(2/3.CRIJ2.P+2/3.CRIJ1.EPSILON)
!                    (diagonal terms for R11, R22 et R33)

!      -DELTAIJ*2/3*EPSILON

!     If we extrapolate the source terms
!     We modify the implicit part:
!     In PHI1, we will only take rho CRIJ1 epsilon/k and not
!                                rho CRIJ1 epsilon/k (1-2/3 DELTAIJ)
!     It allow to keep  k^n instead of (R11^(n+1)+R22^n+R33^n)
!     This choice is questionable. It is the solution isoluc = 1
!     If we want to take all as implicit (like it is done in
!     standard first order), it is the solution isoluc = 2
!     -> to  be tested more precisely if necessary


!     If we extrapolate the source terms
if (st_prv_id.ge.0) then

  isoluc = 1

  do iel = 1, ncel

    !     Half-traces of Prod and R
    trprod = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
    trrij  = w8(iel)

    !     Calculation of Prod+Phi1+Phi2-Eps
    !       = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
    !       In c_st_prv:
    !       = rhoPij-C1rho eps/k(   -2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
    !       = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij           }
    c_st_prv(iel) = c_st_prv(iel) + cromo(iel) * cell_f_vol(iel)  &
      *(   deltij*d2s3*                                           &
           (  crij2*trprod                                        &
            +(crij1-1.d0)* cvara_ep(iel)  )                       &
         +(1.0d0-crij2)*produc(isou,iel)               )
    !       In smbr
    !       =       -C1rho eps/k(Rij         )
    !       = rho{                                     -C1eps/kRij}
    smbr(iel) = smbr(iel) + crom(iel) * cell_f_vol(iel)               &
      *( -crij1*cvara_ep(iel)/trrij * cvara_var(iel) )

    !     Calculation of the implicit part coming from Phil
    !       = C1rho eps/k(1        )
    rovsdt(iel) = rovsdt(iel) + crom(iel) * cell_f_vol(iel)           &
                            *crij1*cvara_ep(iel)/trrij*thetv

  enddo

  !     If we want to implicit a part of -C1rho eps/k(   -2/3k dij)
  if (isoluc.eq.2) then

    do iel = 1, ncel

      trrij  = w8(iel)

     !    We remove of cromo
     !       =       -C1rho eps/k(   -1/3Rij dij)
      c_st_prv(iel) = c_st_prv(iel) - cromo(iel) * cell_f_vol(iel)    &
      *(deltij*d1s3*crij1*cvara_ep(iel)/trrij * cvara_var(iel))
      !    We add to smbr (with crom)
      !       =       -C1rho eps/k(   -1/3Rij dij)
      smbr(iel)                 = smbr(iel)                       &
                          + crom(iel) * cell_f_vol(iel)           &
      *(deltij*d1s3*crij1*cvara_ep(iel)/trrij * cvara_var(iel))
      !    We add to rovsdt (woth crom)
      !       =        C1rho eps/k(   -1/3    dij)
      rovsdt(iel) = rovsdt(iel) + crom(iel) * cell_f_vol(iel)         &
      *(deltij*d1s3*crij1*cvara_ep(iel)/trrij                 )
    enddo

  endif

! If we do not extrapolate the source terms
else

  do iel = 1, ncel

    !     Half-traces of Prod and R
    trprod = 0.5d0*(produc(1,iel)+produc(2,iel)+produc(3,iel))
    trrij  = w8(iel)

    !     Calculation of Prod+Phi1+Phi2-Eps
    !       = rhoPij-C1rho eps/k(Rij-2/3k dij)-C2rho(Pij-1/3Pkk dij)-2/3rho eps dij
    !       = rho{2/3dij[C2 Pkk/2+(C1-1)eps)]+(1-C2)Pij-C1eps/kRij}
    smbr(iel) = smbr(iel) + crom(iel) * cell_f_vol(iel)           &
      *(   deltij*d2s3*                                           &
           (  crij2*trprod                                        &
            +(crij1-1.d0)* cvara_ep(iel)  )                           &
         +(1.0d0-crij2)*produc(isou,iel)                          &
         -crij1*cvara_ep(iel)/trrij * cvara_var(iel)  )

    !     Calculation of the implicit part coming from Phi1
    !       = C1rho eps/k(1-1/3 dij)
    rovsdt(iel) = rovsdt(iel) + crom(iel) * cell_f_vol(iel)           &
         *(1.d0-d1s3*deltij)*crij1*cvara_ep(iel)/trrij
  enddo

endif

!===============================================================================
! 6-bis. Coriolis terms in the Phi1 and production
!===============================================================================

if (icorio.eq.1 .or. iturbo.eq.1) then

  do iel = 1, ncel
    w7(iel) = 0.d0
  enddo

  ! Index connectivity (i,j) <-> ivar
  ivar_r(1,1) = ir11
  ivar_r(2,2) = ir22
  ivar_r(3,3) = ir33
  ivar_r(1,2) = ir12
  ivar_r(1,3) = ir13
  ivar_r(2,3) = ir23
  ivar_r(2,1) = ivar_r(1,2)
  ivar_r(3,1) = ivar_r(1,3)
  ivar_r(3,2) = ivar_r(2,3)

  if (ivar.eq.ir11) then
    ii = 1
    jj = 1
  else if (ivar.eq.ir22) then
    ii = 2
    jj = 2
  else if (ivar.eq.ir33) then
    ii = 3
    jj = 3
  else if (ivar.eq.ir12) then
    ii = 1
    jj = 2
  else if (ivar.eq.ir13) then
    ii = 1
    jj = 3
  else if (ivar.eq.ir23) then
    ii = 2
    jj = 3
  end if

  ! Compute Gij: (i,j) component of the Coriolis production
  do kk = 1, 3
    call field_get_val_prev_s(ivarfl(ivar_r(ii,kk)), cvara_rik)
    call field_get_val_prev_s(ivarfl(ivar_r(jj,kk)), cvara_rjk)

    do iel = 1, ncel
      call coriolis_t(irotce(iel), 1.d0, matrot)

      w7(iel) = w7(iel) - ccorio*(  matrot(ii,kk)*cvara_rjk(iel) &
                                  + matrot(jj,kk)*cvara_rik(iel) )
    enddo
  enddo

  ! Coriolis contribution in the Phi1 term: (1-C2/2)Gij
  if (icorio.eq.1) then
    do iel = 1, ncel
      w7(iel) = crom(iel)*cell_f_vol(iel)*(1.d0 - 0.5d0*crij2)*w7(iel)
    enddo
  endif

  ! If source terms are extrapolated
  if (st_prv_id.ge.0) then
    do iel = 1, ncel
      c_st_prv(iel) = c_st_prv(iel) + w7(iel)
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

  call rijech(isou, produc, w7)

  ! If we extrapolate the source terms: c_st_prv
  if (st_prv_id.ge.0) then
    do iel = 1, ncel
       c_st_prv(iel) = c_st_prv(iel) + w7(iel)
     enddo
  ! Otherwise smbr
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

  call rijthe(nscal, ivar, gradro, w7)

  ! If source terms are extrapolated
  if (st_prv_id.ge.0) then
    do iel = 1, ncel
      c_st_prv(iel) = c_st_prv(iel) + w7(iel)
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
if (vcopt%idften.eq.6) then

  call field_get_val_v(ivsten, visten)

  do iel = 1, ncel
    viscce(1,iel) = visten(1,iel) + viscl(iel)
    viscce(2,iel) = visten(2,iel) + viscl(iel)
    viscce(3,iel) = visten(3,iel) + viscl(iel)
    viscce(4,iel) = visten(4,iel)
    viscce(5,iel) = visten(5,iel)
    viscce(6,iel) = visten(6,iel)
  enddo

  iwarnp = vcopt%iwarni

  call vitens &
 ( viscce , iwarnp ,             &
   weighf , weighb ,             &
   viscf  , viscb  )

! Scalar diffusivity
else

  do iel = 1, ncel
    trrij = 0.5d0 * (cvara_r11(iel) + cvara_r22(iel) + cvara_r33(iel))
    rctse = crom(iel) * csrij * trrij**2 / cvara_ep(iel)
    w1(iel) = viscl(iel) + vcopt%idifft*rctse
  enddo

  call viscfa                    &
 ( imvisf ,                      &
   w1     ,                      &
   viscf  , viscb  )

endif

!===============================================================================
! 10. Solving
!===============================================================================

if (st_prv_id.ge.0) then
  thetp1 = 1.d0 + thets
  do iel = 1, ncel
    smbr(iel) = smbr(iel) + thetp1*c_st_prv(iel)
  enddo
endif

iconvp = vcopt%iconv
idiffp = vcopt%idiff
ndircp = ndircl(ivar)
nswrsp = vcopt%nswrsm
nswrgp = vcopt%nswrgr
imligp = vcopt%imligr
ircflp = vcopt%ircflu
ischcp = vcopt%ischcv
isstpp = vcopt%isstpc
iescap = 0
imucpp = 0
idftnp = vcopt%idften
iswdyp = vcopt%iswdyn
iwarnp = vcopt%iwarni
blencp = vcopt%blencv
epsilp = vcopt%epsilo
epsrsp = vcopt%epsrsm
epsrgp = vcopt%epsrgr
climgp = vcopt%climgr
extrap = vcopt%extrag
relaxp = vcopt%relaxv
! all boundary convective flux with upwind
icvflb = 0

call codits &
 ( idtvar , ivarfl(ivar)    , iconvp , idiffp , ndircp ,          &
   imrgra , nswrsp , nswrgp , imligp , ircflp ,                   &
   ischcp , isstpp , iescap , imucpp , idftnp , iswdyp ,          &
   iwarnp ,                                                       &
   blencp , epsilp , epsrsp , epsrgp , climgp , extrap ,          &
   relaxp , thetv  ,                                              &
   cvara_var       , cvara_var       ,                            &
   coefap , coefbp , cofafp , cofbfp ,                            &
   imasfl , bmasfl ,                                              &
   viscf  , viscb  , viscf  , viscb  , viscce ,                   &
   weighf , weighb ,                                              &
   icvflb , ivoid  ,                                              &
   rovsdt , smbr   , cvar_var        , dpvar  ,                   &
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

 1000 format(/,'           Resolution pour la variable ',A8,/)

#else

 1000 format(/,'           Solving variable ',A8           ,/)

#endif

!----
! End
!----

return

end subroutine
