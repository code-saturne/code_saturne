!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \file inivar.f90
!> \brief Initialization of calculation variables, time step
!> and table that stores distance to the wall
!> by the user (after reading a restart file).
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!______________________________________________________________________________

subroutine inivar &
 ( nvar   , nscal )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use pointe
use entsor
use parall
use period
use ppppar
use ppthch
use ppincl
use mesh
use field
use cfpoin, only:ithvar
use cs_c_bindings
use cs_cf_bindings
use cs_f_interfaces
use vof

use, intrinsic :: iso_c_binding

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

! Local variables

character(len=80) :: chaine
integer          ivar  , iscal
integer          iel
integer          iok   , ii, iclvfl
integer          kscmin, kscmax, keyvar, n_fields, kclvfl
integer          f_id, f_dim
integer          iflid, iflidp
integer          idimf
integer          ivoid, uprtot

double precision valmax, valmin, vfmin , vfmax
double precision gravn
double precision xxp0, xyp0, xzp0
double precision scmaxp, scminp

double precision rvoid(1)
double precision vvoid(3)

double precision, dimension(:), pointer :: dt
double precision, dimension(:), pointer :: field_s_v
double precision, dimension(:,:), pointer :: field_v_v
double precision, dimension(:), pointer :: cvar_pr
double precision, dimension(:), pointer :: cvar_var
double precision, dimension(:), pointer :: cpro_prtot

type(var_cal_opt) :: vcopt

!===============================================================================

procedure() :: ppiniv0, ppiniv1,  cs_user_f_initialization

interface

  subroutine cs_gui_initial_conditions()  &
      bind(C, name='cs_gui_initial_conditions')
    use, intrinsic :: iso_c_binding
    implicit none
  end subroutine cs_gui_initial_conditions

  function cs_turbulence_init_clip_and_verify() result(n_errors) &
      bind(C, name='cs_turbulence_init_clip_and_verify')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int) :: n_errors
  end function cs_turbulence_init_clip_and_verify

  subroutine cs_navstv_total_pressure() &
   bind(C, name='cs_f_navier_stokes_total_pressure')
   use,intrinsic :: iso_c_binding
   implicit none
 end subroutine cs_navstv_total_pressure

 subroutine vof_compute_linear_rho_mu() &
      bind(C, name='cs_f_vof_compute_linear_rho_mu')
   use, intrinsic :: iso_c_binding
   implicit none
 end subroutine vof_compute_linear_rho_mu

end interface

!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_val_s_by_name('dt', dt)

call field_get_n_fields(n_fields)

call field_get_key_id("variable_id", keyvar)

! Key ids for clipping
call field_get_key_id("min_scalar_clipping", kscmin)
call field_get_key_id("max_scalar_clipping", kscmax)
call field_get_key_id("variance_clipping", kclvfl)

iok = 0

gravn = sqrt(gx**2+gy**2+gz**2)

!===============================================================================
! 2. ON REPASSE LA MAIN A L'UTILISATEUR POUR LA PROGRAMMATION DES
!    INITIALISATIONS QUI LUI SONT PROPRES
!===============================================================================

iflidp = -1

do ivar = 1, nvar
  call field_get_key_struct_var_cal_opt(ivarfl(ivar), vcopt)
  if (vcopt%iwgrec.eq.1) then

    if (vcopt%idiff.lt.1) cycle
    iflid = ivarfl(ivar)
    if (iflid.eq.iflidp) cycle
    iflidp = iflid

    call field_get_key_int(iflid, kwgrec, f_id)
    call field_get_dim(f_id, idimf)

    if (idimf.eq.6) then
      call field_get_val_v(f_id, field_v_v)
      do iel = 1, ncelet
        field_v_v(1,iel) = 1.d0
        field_v_v(2,iel) = 1.d0
        field_v_v(3,iel) = 1.d0
        field_v_v(4,iel) = 0.d0
        field_v_v(5,iel) = 0.d0
        field_v_v(6,iel) = 0.d0
      enddo
    else if (idimf.eq.1) then
      call field_get_val_s(f_id, field_s_v)
      do iel = 1, ncelet
        field_s_v(iel) = 1.d0
      enddo
    endif

  endif
enddo

! First pre-initialization for specific physic modules
! before the GUI
if (ippmod(iphpar).gt.0) then
  call ppiniv0
endif

! GUI definitions
! ===============

call cs_gui_initial_conditions

! User subroutine
! ===============

call cs_user_f_initialization(nvar, nscal, dt)

call user_initialization()

! Second stage of initialization for specific physic modules
! after the user
if (ippmod(iphpar).gt.0) then
  call ppiniv1
endif

! VoF algorithm
if (ivofmt.gt.0) then
  call vof_compute_linear_rho_mu
  ! density is stored at the two previous time steps
  call field_current_to_previous(icrom)
  call field_current_to_previous(ibrom)
  call field_current_to_previous(icrom)
  call field_current_to_previous(ibrom)
endif

if (ippmod(icompf).ge.0.and.(    isuite.eq.0                 &
                             .or.isuite.eq.1.and.ileaux.eq.0)) then

  if (     ithvar.ne. 60000.and.ithvar.ne.100000                    &
      .and.ithvar.ne.140000.and.ithvar.ne.150000.and.ithvar.ne.210000) then
      write(nfecra,1000) ithvar
      iok = iok + 1
  endif

  ivoid = -1
  call cs_cf_thermo(ithvar, ivoid,  rvoid, rvoid, rvoid, vvoid)

endif

! Pressure / Total pressure initialisation

! Standard:
! If the user has initialized the total pressure Ptot, P* is initialized
! accordingly, only if the user has speficied the reference point.
! (all values of the total pressure have to be initialized).
! Otherwise, the total pressure is initialized using P*,
! Ptot = P* + P0 + rho.g.r

! In case of restart without auxiliary, Ptot is recomputed with P*.
! (For EVM models, the shift by 2/3*rho*k is missing)
! In case of restart with auxiliary, nothing is done.

! Compressible:
! The total pressure field does not need to be defined. The solved pressure is
! the total pressure.

if (ippmod(icompf).lt.0) then

  call field_get_val_s(ivarfl(ipr), cvar_pr)
  call field_get_val_s(iprtot, cpro_prtot)

  uprtot = 0

  if (ixyzp0.gt.-1.and.(isuite.eq.0.or.ileaux.eq.0)) then
    uprtot = 1
    do iel = 1, ncel
      if (cpro_prtot(iel).le.-0.5d0*rinfin) then
        uprtot = 0
        exit
      endif
    enddo
  endif

  xxp0   = xyzp0(1)
  xyp0   = xyzp0(2)
  xzp0   = xyzp0(3)

  if (uprtot.gt.0) then
    do iel = 1, ncel
      cvar_pr(iel) =  cpro_prtot(iel)               &
                    - ro0*( gx*(xyzcen(1,iel)-xxp0) &
                    + gy*(xyzcen(2,iel)-xyp0)       &
                    + gz*(xyzcen(3,iel)-xzp0) )     &
                    + pred0 - p0
    enddo
  elseif (isuite.eq.0.or.ileaux.eq.0) then
    call cs_navstv_total_pressure
  endif

endif

!===============================================================================
! 3.  CLIPPING DES GRANDEURS TURBULENTES (UTILISATEUR OU SUITE)
!===============================================================================

if (cs_turbulence_init_clip_and_verify() .gt. 0) then
  iok = iok + 1
endif

!===============================================================================
! 4.  CLIPPING DES SCALAIRES (UTILISATEUR OU SUITE)
!     Si l'utilisateur est intervenu dans USINIV ou PPINIV et
!       a impose des valeurs "correctes" (au sens comprises dans des bornes
!         simplifiees a base de 0, scamin, scamax)
!         on considere qu'il s'agit d'une initialisation admissible,
!         on la clippe pour la rendre coherente avec le clipping du code
!         et on continue le calcul
!       si l'utilisateur a impose des valeurs visiblement erronees
!         (au sens comprises dans des bornes simplifiees a base de 0, scamin,
!          scamax), on s'arrete (il s'est sans doute trompe).
!     On adopte le meme traitement en suite de calcul
!       pour assurer un comportement identique en suite entre un calcul
!       ou l'utilisateur modifie une variable avec usiniv (mais pas un
!       scalaire) et un calcul ou l'utilisateur ne modifie pas usiniv.
!     Sinon, les grandeurs ont deja ete clippees apres les init par defaut

!     Pour resumer :
!      -en   suite  avec des valeurs grossierement admissibles : on clippe
!      -avec usiniv ou ppiniv
!                   avec des valeurs grossierement admissibles : on clippe
!      -non suite sans usiniv ni ppiniv :
!                                      grandeurs par defaut (deja clippees)
!      -suite ou usiniv ou ppiniv
!                   avec une valeur grossierement non admissible : stop
!         (on souhaite indiquer a l'utilisateur que son fichier suite est
!          bizarre ou que son initialisation est fausse et qu'il a donc
!          fait au moins une erreur qui peut en cacher d'autres)
!===============================================================================

! On traite tous les scalaires d'abord, car ils peuvent etre necessaires
!     pour clipper les variances

if (nscal.gt.0) then

!     Scalaires non variance

  do ii = 1, nscal
    f_id = ivarfl(isca(ii))
    call field_get_dim(f_id, f_dim)
    if((iscavr(ii).le.0.or.iscavr(ii).gt.nscal).and.f_dim.eq.1) then

      ! Get the min clipping
      call field_get_key_double(f_id, kscmin, scminp)
      call field_get_key_double(f_id, kscmax, scmaxp)

      if (scminp.le.scmaxp) then
        ivar = isca(ii)
        call field_get_val_s(ivarfl(ivar), cvar_var)
        valmax = cvar_var(1)
        valmin = cvar_var(1)
        do iel = 1, ncel
          valmax = max(valmax,cvar_var(iel))
          valmin = min(valmin,cvar_var(iel))
        enddo
        if (irangp.ge.0) then
          call parmax (valmax)
          call parmin (valmin)
        endif

        ! Check coherence for clippings of non-variance scalars.
        if (valmin.ge.scminp.and.valmax.le.scmaxp) then
          iscal = ii
          call clpsca(iscal)
        else
          call field_get_label(ivarfl(isca(ii)), chaine)
          write(nfecra,3040) ii,chaine(1:16),                     &
                             valmin,scminp,valmax,scmaxp
          iok = iok + 1
        endif
      endif

    endif
  enddo

  ! Variances

  do ii = 1, nscal
    if(iscavr(ii).gt.0.and.iscavr(ii).le.nscal) then

      ! Get the min clipping
      f_id = ivarfl(isca(ii))
      call field_get_key_double(f_id, kscmin, scminp)
      call field_get_key_double(f_id, kscmax, scmaxp)


      if (scminp.le.scmaxp) then
        ivar = isca(ii)
        call field_get_val_s(ivarfl(ivar), cvar_var)
        valmax = cvar_var(1)
        valmin = cvar_var(1)
        do iel = 1, ncel
          valmax = max(valmax,cvar_var(iel))
          valmin = min(valmin,cvar_var(iel))
        enddo
        if (irangp.ge.0) then
          call parmax (valmax)
          call parmin (valmin)
        endif

        call field_get_key_int(ivarfl(ivar), kclvfl, iclvfl)

        ! Verification de la coherence pour les clippings de variance.
        ! Pour iclvfl = 1 on ne verifie que > 0 sinon ca va devenir difficile
        ! de faire une initialisation correcte.

        if (iclvfl.eq.0) then
          ! On pourrait clipper dans le cas ou VALMIN.GE.0, mais ca
          ! n'apporterait rien, par definition
          if(valmin.lt.0.d0) then
            call field_get_name(ivarfl(isca(ii)), chaine)
            write(nfecra,3050)ii,chaine(1:16),                     &
                              valmin,scminp,valmax,scmaxp
            iok = iok + 1
          endif
        elseif (iclvfl.eq.1) then
          ! Here we clip to be coherent with the scalar's value.
          if(valmin.ge.0.d0) then
            iscal = ii
            call clpsca(iscal)
          else
            call field_get_name(ivarfl(isca(ii)), chaine)
            write(nfecra,3050)ii,chaine(1:16),valmin,scminp,valmax,scmaxp
            iok = iok + 1
          endif
        elseif (iclvfl.eq.2) then
          vfmin = 0.d0
          vfmin = max(scminp, vfmin)
          vfmax = scmaxp
          ! We could clip when valmin >= vfmin and valmax <= vfmax
          ! but by definition, this would add nothing.
          if(valmin.lt.vfmin.or.valmax.gt.vfmax) then
            call field_get_name(ivarfl(isca(ii)), chaine)
            write(nfecra,3051)ii,chaine(1:16),                     &
                              valmin,scminp,valmax,scmaxp,         &
                              ii,iclvfl
            iok = iok + 1
          endif
        endif
      endif

    endif
  enddo

endif

call user_extra_operations_initialize()

!===============================================================================
! 5.  IMPRESSIONS DE CONTROLE POUR LES INCONNUES, LE PAS DE TEMPS
!        LE CUMUL DES DUREE POUR LES MOYENNES
!===============================================================================

write(nfecra,2000)

!===============================================================================
! 6.  ARRET GENERAL SI PB
!===============================================================================

if (iok.gt.0) then
  write(nfecra,3090) iok
  call csexit (1)
endif

!----
! Formats
!----

 1000 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING :     stop in compressible thermodynamics at    ',/,&
'@    =========     initialisation.                           ',/,&
'@                                                            ',/,&
'@    The computation will stop.                              ',/,&
'@                                                            ',/,&
'@    Unexpected value of the indicator ithvar (',i10,').     ',/,&
'@                                                            ',/,&
'@    Two and only two independant variables among            ',/,&
'@    P, rho, T and E (except T and E) should be imposed at   ',/,&
'@    the initialisation in the GUI or in the user subroutine ',/,&
'@    of initialization (cs_user_initialization.f90) or       ',/,&
'@    in both.                                                ',/,&
'@                                                            ',/,&
'@    Check if iccfth has not been partially set through the  ',/,&
'@    GUI in a non consistant manner with                     ',/,&
'@    cs_user_initialization.f90.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

 2000 format(                                                     &
                                                                /,&
' ** VARIABLES INITIALIZATION',                                 /,&
'    ------------------------',                                 /,&
''                                                             )

 3040 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     SCALAR QUANTITIES OUT OF BOUNDS',                        /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Scalar number ',i10,': ',a16,                               /,&
'@  Minimum value                = ',e14.5,                     /,&
'@    Desired clipping at SCAMIN = ',e14.5,                     /,&
'@  Maximum value                = ',e14.5,                     /,&
'@    Desired clipping at SCAMAX = ',e14.5,                     /,&
'@  The bounds are not coherent with the limits SCAMIN and',    /,&
'@    SCAMAX.',                                                 /,&
'@',                                                            /,&
'@  Verify the initialization and the restart file.',           /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@  Verify the clipping values.',                               /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3050 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     NEGATIVE VARIANCE',                                      /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Scalar number ',i10,': ',a16,                               /,&
'@  Minimum value               = ',e14.5,                      /,&
'@  This scalar is a variance (ISCAVR is positive)',            /,&
'@    but the initialization has some negative values.',        /,&
'@',                                                            /,&
'@  Verify the initialization and the restart file.',           /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@  Verify the variance definition.',                           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3051 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@     VARIANCE OUT OF BOUNDS',                                 /,&
'@',                                                            /,&
'@  The calculation will not be run.',                          /,&
'@',                                                            /,&
'@  Scalar number ',i10,': ',a16,                               /,&
'@  Minimum value                = ',e14.5,                     /,&
'@    Desired clipping at SCAMIN = ',e14.5,                     /,&
'@  Maximum value                = ',e14.5,                     /,&
'@    Desired clipping at SCAMAX = ',e14.5,                     /,&
'@  This scalar is a variance (ISCAVR is positive)',            /,&
'@    but the initialization has some values out',              /,&
'@    of the bounds SCAMIN, SCAMAX or lower than 0 and the',    /,&
'@    desired clipping mode is ICLVFL(',i10,') = ',i10,         /,&
'@',                                                            /,&
'@  Verify the initialization and the restart file.',           /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@  Verify the variance definition and the clipping mode.',     /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 3090 format(                                                     &
'@',                                                            /,&
'@',                                                            /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ WARNING: ABORT IN THE VARIABLES INITIALIZATION',          /,&
'@    ========',                                                /,&
'@',                                                            /,&
'@    THE VARIABLES INITIALIZATION IS INCOMPLETE OR',           /,&
'@    INCOHERENT WITH THE PARAMETERS VALUE OF THE CALCULATION', /,&
'@',                                                            /,&
'@  The calculation will not be run (',i10,' errors).',         /,&
'@',                                                            /,&
'@  Refer to the previous warnings for further information.',   /,&
'@  Pay attention to the initialization of',                    /,&
'@                                the time-step',               /,&
'@                                the turbulence',              /,&
'@                                the scalars and variances',   /,&
'@                                the time averages',           /,&
'@',                                                            /,&
'@  Verify the initialization and the restart file.',           /,&
'@  In the case where the values read in the restart file',     /,&
'@    are incorrect, they may be modified with',                /,&
'@    cs_user_initialization.f90 or with the interface.',       /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

!----
! End
!----

return
end subroutine
