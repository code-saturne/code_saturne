!-------------------------------------------------------------------------------

!VERS

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
! Purpose:
! -------

!> \file cs_user_physical_properties-scalar-drift.f90
!>
!> brief Definition of physical variable laws for scalars with a drift.
!>
!> See \subpage physical_properties for examples.
!>

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     mbrom         indicator of filling of romb array
!> \param[in]     dt            time step (per cell)
!_______________________________________________________________________________

subroutine usphyv &
 ( nvar   , nscal  ,                                              &
   mbrom  ,                                                       &
   dt     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use field
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          mbrom

double precision dt(ncelet)

! Local variables

!< [loc_var_dec]
integer          ivart, iel, ifac
integer          iscal, iflid, iscdri, ifcvsl
integer          f_id, keydri, nfld, keysca
double precision rho, viscl
double precision diamp, rhop, cuning
double precision xvart, xk, xeps, beta1
double precision turb_schmidt

character*80     fname

double precision, dimension(:), pointer :: cpro_taup
double precision, dimension(:), pointer :: cpro_taufpt
double precision, dimension(:), pointer :: cpro_rom
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_omg
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cpro_viscl, cpro_vscal
double precision, dimension(:), pointer :: cvar_scalt

!< [loc_var_dec]

!===============================================================================

!===============================================================================
! 0. Initializations to keep
!===============================================================================

!< [init]
call field_get_val_s(iviscl, cpro_viscl)
call field_get_val_s(icrom, cpro_rom)

! Key id for drift scalar
call field_get_key_id("drift_scalar_model", keydri)

! Key id for scalar id
call field_get_key_id("scalar_id", keysca)

! Number of fields
call field_get_n_fields(nfld)
!< [init]

!===============================================================================

!   The following examples should be adapted by the user
!   ====================================================

!===============================================================================
!  Example: If thermophorese is required, one MUST set the diffusivity
!  =======
!  (Brownian motion) to be variable in space and set the proper relation
!  between the molecular diffusivity and T:
!  ex: Kb x T x cuning /(3*pi*diamp(iscal)*cpro_viscl(iel))
!===============================================================================

!    Excluding:
!      - temperature, enthalpy (handled above)
!      - fluctuation variances (property equal to that of the associated scalar)
!
!    Below, we define the same diffusivity law for all scalars (except the
!      ones excluded above).
!    Values of this property must be defined at cell centers
!  ===================================================================

!< [example_1]

! Loop over fields which are scalar with a drift
do iflid = 0, nfld-1

  call field_get_key_int(iflid, keydri, iscdri)

  ! We only handle here scalar with a drift
  if (btest(iscdri, DRIFT_SCALAR_ADD_DRIFT_FLUX)) then

    ! Position of variables, coefficients
    ! -----------------------------------

    ! Index of the scalar
    call field_get_key_int(iflid, keysca, iscal)

    ! --- Number of the thermal variable
    if (iscalt.gt.0) then
      ivart = isca(iscalt)
      call field_get_val_s(ivarfl(ivart), cvar_scalt)
    else
      write(nfecra,9010) iscalt
      call csexit (1)
    endif

    ! --- Scalar's diffusivity (Brownian motion)

    call field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
    if (ifcvsl.ge.0) then
      call field_get_val_s(ifcvsl, cpro_vscal)
    else
      cpro_vscal => NULL()
    endif

    ! --- Coefficients of drift scalar CHOSEN BY THE USER
    !       Values given here are fictitious

    ! diamp: is the diameter of the particle class
    ! cuning: is the Cuningham correction factor
    ! rhop: particle density
    diamp = 1.d-4
    cuning = 1.d0
    rhop = 1.d4

    ! Name of the scalar with a drift
    call field_get_name(iflid, fname)

    ! Index of the corresponding relaxation time (cpro_taup)
    call field_get_id('drift_tau_'//trim(fname), f_id)
    call field_get_val_s(f_id, cpro_taup)

    ! Index of the corresponding interaction time particle--eddies (cpro_taufpt)
    if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS)) then
      call field_get_id('drift_turb_tau_'//trim(fname), f_id)
      call field_get_val_s(f_id, cpro_taufpt)
    endif

    ! Computation of the relaxation time of the particles
    !----------------------------------------------------

    if (diamp.le.1.d-6) then
      ! Cuningham's correction for submicronic particules
      do iel = 1, ncel
        cpro_taup(iel) = cuning*diamp**2*rhop/(18.d0*cpro_viscl(iel))
      enddo
    else
      do iel = 1, ncel
        cpro_taup(iel) = diamp**2*rhop/(18.d0*cpro_viscl(iel))
      enddo
    endif

    ! Compute the interaction time particle--eddies (tau_fpt)
    !--------------------------------------------------------

    if (btest(iscdri, DRIFT_SCALAR_TURBOPHORESIS)) then

      ! k-epsilon or v2-f models
      if (itytur.eq.2 .or. itytur.eq.5) then
        call field_get_val_s(ivarfl(ik), cvar_k)
        call field_get_val_s(ivarfl(iep), cvar_ep)
        call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)
        do iel = 1, ncel
          xk = cvar_k(iel)
          xeps = cvar_ep(iel)
          cpro_taufpt(iel) = (3.d0/2.d0)*(cmu/turb_schmidt)*xk/xeps
        enddo

      ! Rij-epsilon models
      else if (itytur.eq.3) then
        call field_get_val_s(ivarfl(ir11), cvar_r11)
        call field_get_val_s(ivarfl(ir22), cvar_r22)
        call field_get_val_s(ivarfl(ir33), cvar_r33)
        call field_get_val_s(ivarfl(iep), cvar_ep)
        beta1  = 0.5d0+3.d0/(4.d0*xkappa)
        do iel = 1, ncel
          xk = 0.5d0*( cvar_r11(iel)                    &
                      +cvar_r22(iel)                    &
                      +cvar_r33(iel) )
          xeps = cvar_ep(iel)
          cpro_taufpt(iel) = xk/xeps/beta1
        enddo

      ! k-omega models
      else if (iturb.eq.60) then
        call field_get_val_s(ivarfl(ik), cvar_k)
        call field_get_val_s(ivarfl(iomg), cvar_omg)
        call field_get_key_double(ivarfl(isca(iscal)), ksigmas, turb_schmidt)
        do iel = 1, ncel
          xk = cvar_k(iel)
          xeps = cmu*xk*cvar_omg(iel)
          cpro_taufpt(iel) = (3.d0/2.d0)*(cmu/turb_schmidt)*xk/xeps
        enddo
      endif

    endif

    ! Brownian diffusion at cell centers
    !-----------------------------------

    ! --- Stop if the diffusivity is not variable
    call field_get_key_int(ivarfl(isca(iscal)), kivisl, ifcvsl)
    if (ifcvsl.lt.0) then
      write(nfecra,1010) iscal
      call csexit (1)
    endif

    ! Homogeneous to a dynamic viscosity
    do iel = 1, ncel
      xvart = cvar_scalt(iel)
      rho = cpro_rom(iel)
      viscl = cpro_viscl(iel)
      cpro_vscal(iel) = rho*kboltz*xvart*cuning/(3.d0*pi*diamp*viscl)
    enddo

  endif ! --- Tests on drift scalar
enddo
!< [example_1]

!===============================================================================

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    DONNEES DE CALCUL INCOHERENTES                          ',/,&
'@                                                            ',/,&
'@    Pour le scalaire ',I10                                   ,/,&
'@      la diffusivite est uniforme alors que                 ',/,&
'@      usphyv impose une diffusivite variable.               ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier usipsu ou usphyv.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9010 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DU CALCUL DES GRANDEURS PHYSIQUES',/,&
'@    =========                                               ',/,&
'@    APPEL A csexit DANS LE SOUS PROGRAMME usphyv            ',/,&
'@                                                            ',/,&
'@    La variable dont dependent les proprietes physiques ne  ',/,&
'@      semble pas etre une variable de calcul.               ',/,&
'@    En effet, on cherche a utiliser la temperature alors que',/,&
'@      ISCALT = ',I10                                         ,/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier le codage de usphyv (et le test lors de la     ',/,&
'@      definition de IVART).                                 ',/,&
'@    Verifier la definition des variables de calcul dans     ',/,&
'@      usipsu. Si un scalaire doit jouer le role de la       ',/,&
'@      temperature, verifier que ISCALT a ete renseigne.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@    Inconsistent calculation data',/,                           &
'@',/,                                                            &
'@    For scalar', i10,/,                                         &
'@      the diffusivity is uniform while',/,                      &
'@      usphyv prescribes a variable diffusivity.',/,             &
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Modify usipsu or usphyv.',/,                                &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 9010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:  stop when computing physical quantities',/,       &
'@    =======',/,                                                 &
'@',/,                                                            &
'@    The variable on which physical properties depend does',/,   &
'@      seem to be a calculation variable.',/,                    &
'@    Indeed, we are trying to use the temperature while',/,      &
'@      iscalt = ',i10                                  ,/,&
'@',/,                                                            &
'@    The calculation will not be run.',/,                        &
'@',/,                                                            &
'@    Check the programming in usphyv (and the test when',/,      &
'@      defining ivart).',/,                                      &
'@    Check the definition of calculation variables in',/,        &
'@      usipsu. If a scalar should represent the,',/,             &
'@      temperature, check that iscalt has been defined',/,       &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)

#endif

!----
! End
!----

return
end subroutine usphyv
