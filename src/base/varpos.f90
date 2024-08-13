!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!> \file varpos.f90
!> \brief Variables location initialization, according to calculation type
!> selected by the user.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine varpos () &
 bind(C, name='cs_f_varpos')

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

integer          iscal , id, pflag
integer          ii
integer          iok
integer          key_buoyant_id, coupled_with_vel_p_fld, key_restart_file
integer          kturt, turb_flux_model, kisso2t
integer          isso2t
integer          ibeta

double precision gravn2

character(len=80) :: f_name, s_label, s_name

type(var_cal_opt) :: vcopt

procedure() :: add_property_field_1d, hide_property
procedure() :: add_source_term_prev_field

!===============================================================================
! Initialization
!===============================================================================

! Key id for buoyant field (inside the Navier Stokes loop)
call field_get_key_id("coupled_with_vel_p", key_buoyant_id)
call field_get_key_id('turbulent_flux_model', kturt)
call field_get_key_id("scalar_time_scheme", kisso2t)
call field_get_key_id("restart_file", key_restart_file)

pflag = POST_ON_LOCATION + POST_MONITOR

!===============================================================================
! Time-scheme related properties
!===============================================================================

! Dans le cas ou on a un schema en temps d'ordre 2, il faut aussi
!   prevoir les proprietes au temps n-1. Ce sera fait au dernier appel

! Dans le cas ou on calcule des moments, il faut en prevoir le nombre
!   et prevoir le nombre de tableaux necessaires pour stocker le
!   temps cumule de moyenne. On suppose que l'on peut s'aider
!   d'infos situees en tete de fichier suite (si on fait une
!   suite avec des moments non reinitialises).

! 1.1 PROPRIETES ADDITIONNELLES POUR LES ET SCHEMA EN TEMPS
! ---------------------------------------------------------

! Initialisations par defaut eventuelles et verifications
!   des options utilisees ci-dessous pour decider si l'on
!   reserve des tableaux supplementaires pour des grandeurs
!   au pas de temps precedent

iok = 0

!  Pression hydrostatique
if (iphydr.eq.1.and.icalhy.eq.-1) then
  gravn2 = gx**2+gy**2+gz**2
  if (gravn2.lt.epzero**2) then
    icalhy = 0
  else
    icalhy = 1
  endif
else
  icalhy = 0
endif

!   Global time stepping
!       for LES: 2nd order; 1st order otherwise
!       (2nd order forbidden for "coupled" k-epsilon)
if (ischtp.eq.-1) then
  if ((itytur.eq.4).or.(hybrid_turb.eq.4)) then
    ischtp = 2
  else
    ischtp = 1
  endif
endif

!   Schemas en temps : variables deduites
!   Schema pour le Flux de masse
if (istmpf.eq.-999) then
  if (ischtp.eq.1) then
    istmpf = 1
  else if (ischtp.eq.2) then
    istmpf = 2
  endif
endif

!   Collocated time scheme for gaz combustion
if (itpcol.eq.-1) then
  if (ischtp.eq.2.and.ippmod(islfm).ge.0) then
    itpcol = 1
  else
    itpcol = 0
  endif
endif

!     Termes sources NS,
if (isno2t.eq.-999) then
  if (ischtp.eq.1) then
    isno2t = 0
  else if (ischtp.eq.2) then
    !       Pour le moment par defaut on prend l'ordre 2
    isno2t = 1
  endif
endif
!     Termes sources turbulence (k-eps, Rij, v2f ou k-omega)
!     On n'autorise de changer ISTO2T qu'en Rij (sinon avec
!       le couplage k-eps/omega il y a pb)
if (isto2t.eq.-999) then
  if (ischtp.eq.1) then
    isto2t = 0
  else if (ischtp.eq.2) then
    !       Pour le moment par defaut on ne prend pas l'ordre 2
    !              ISTO2T = 1
    isto2t = 0
  endif
else if (itytur.eq.2.or.iturb.eq.50.or.iturb.ne.60) then
  write(nfecra,8132) iturb,isto2t
  iok = iok + 1
endif

do iscal = 1, nscal
  ! Termes sources Scalaires,
  call field_get_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
  if (isso2t.eq.-1) then
    if (ischtp.eq.1) then
      isso2t = 0
      call field_set_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
    else if (ischtp.eq.2) then
      ! Pour coherence avec Navier Stokes on prend l'ordre 2
      ! mais de toute facon qui dit ordre 2 dit LES et donc
      ! generalement pas de TS scalaire a interpoler.
      isso2t = 1
      call field_set_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
      if (iscal.eq.iscalt .and. iirayo.gt.0) then
        isso2t = 0
        call field_set_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
      end if
    endif
  endif

  call field_get_key_int(ivarfl(isca(iscal)), kturt, turb_flux_model)

  if (iscal.eq.iscalt) then
    call field_get_id_try("thermal_expansion", ibeta)
    if (turb_flux_model.gt.0.and.irovar.eq.1.and.ibeta.eq.-1) then
      call add_property_field_1d('thermal_expansion', 'Beta', ibeta)
    endif
  endif

enddo

! Schemas en temps

!     Schema en temps global.
if (ischtp.ne. 1.and.ischtp.ne.2) then
  write(nfecra,8101) 'ISCHTP',ischtp
  iok = iok + 1
endif
if (ischtp.eq. 2.and.idtvar.ne.0) then
  write(nfecra,8111) ischtp,idtvar
  iok = iok + 1
endif
if (ischtp.eq. 2.and.itytur.eq.2) then
  write(nfecra,8112) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq.1.and.itytur.eq.4) then
  write(nfecra,8113) ischtp,iturb
endif
if (ischtp.eq. 2.and.iturb.eq.50) then
  write(nfecra,8114) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq. 2.and.iturb.eq.51.and.hybrid_turb.ne.4) then
  write(nfecra,8117) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq. 2.and.iturb.eq.60.and.hybrid_turb.ne.4) then
  write(nfecra,8115) ischtp,iturb
  iok = iok + 1
endif
if (ischtp.eq. 2.and.iturb.eq.70) then
  write(nfecra,8116) ischtp,iturb
  iok = iok + 1
endif

! Schema en temps pour le flux de masse
if (istmpf.ne. 2.and.istmpf.ne.0.and.istmpf.ne. 1) then
  write(nfecra,8121) 'ISTMPF',istmpf
  iok = iok + 1
endif

! Schema en temps pour les termes sources de NS
if (isno2t.ne.0.and.isno2t.ne. 1.and.isno2t.ne.2) then
  write(nfecra,8131) 'ISNO2T',isno2t
  iok = iok + 1
endif
! Schema en temps pour les termes sources des grandeurs turbulentes
if (isto2t.ne.0.and.isto2t.ne. 1.and.isto2t.ne.2) then
  write(nfecra,8131) 'ISTO2T',isto2t
  iok = iok + 1
endif

do iscal = 1, nscal
  ! Schema en temps pour les termes sources des scalaires
  call field_get_key_int(ivarfl(isca(iscal)), kisso2t, isso2t)
  if (isso2t.ne.0.and.isso2t.ne. 1.and.isso2t.ne.2) then
    write(nfecra,8141) iscal,'ISSO2T',isso2t
    iok = iok + 1
  endif
enddo

! Stop si probleme
if (iok.gt.0) then
  call csexit(1)
endif

! add thermal expansion field for Boussinesq approximation
! if not already added
if (idilat.eq.0) then
  call field_get_id_try('thermal_expansion', ibeta)
  if (ibeta.lt.0) then
    call add_property_field_1d('thermal_expansion', 'Beta', ibeta)
  endif
endif

! Source term for weakly compressible algorithm (semi analytic scheme)
if (idilat.ge.4) then
  do iscal = 1, nscal
    id = ivarfl(isca(iscal))
    call field_get_name(id, s_name)
    call field_get_label(id, s_label)
    f_name  = trim(s_name) // '_dila_st'
    call add_property_field_1d(f_name, '', iustdy(iscal))
    id = iustdy(iscal)
    call field_set_key_int(id, keyvis, 0)
    ! Set restart file option for source terms
    call field_set_key_int(id, key_restart_file, RESTART_AUXILIARY)
  enddo
  itsrho = nscal + 1
  call add_property_field_1d('dila_st', '', iustdy(itsrho))
  id = iustdy(iscal)
  call field_set_key_int(id, keyvis, 0)
endif

! On a besoin d'un tableau pour les termes sources de Navier Stokes
!  a extrapoler. Ce tableau est NDIM dans le cas general et NDIM+1
!  si on extrapole aussi les termes sources de l equation sur le taux
!  de vide pour l'algo. VOF.
if (isno2t.gt.0) then
  call add_source_term_prev_field(ivarfl(iu))
  if (ivofmt.gt.0) then
    call add_source_term_prev_field(ivarfl(ivolf2))
  endif
endif

if (isto2t.gt.0) then
  ! The dimension of this array depends on turbulence model:
  if (itytur.eq.2) then
    call add_source_term_prev_field(ivarfl(ik))
    call add_source_term_prev_field(ivarfl(iep))
  else if (itytur.eq.3) then
    call add_source_term_prev_field(ivarfl(irij))
    call add_source_term_prev_field(ivarfl(iep))
    if (iturb.eq.32) then
      call add_source_term_prev_field(ivarfl(ial))
    endif
  else if (itytur.eq.5) then
    call add_source_term_prev_field(ivarfl(ik))
    call add_source_term_prev_field(ivarfl(iep))
    call add_source_term_prev_field(ivarfl(iphi))
    if (iturb.eq.50) then
      call add_source_term_prev_field(ivarfl(ifb))
    else if (iturb.eq.51) then
      call add_source_term_prev_field(ivarfl(ial))
    endif
  else if (iturb.eq.60) then
    call add_source_term_prev_field(ivarfl(ik))
    call add_source_term_prev_field(ivarfl(iomg))
  else if (iturb.eq.70) then
    call add_source_term_prev_field(ivarfl(inusa))
  endif
endif

! Proprietes des scalaires : termes sources pour theta schema
if (nscal.ge.1) then
  do ii = 1, nscal
    call field_get_key_int(ivarfl(isca(ii)), kisso2t, isso2t)
    if (isso2t.gt.0) then
      ! For buoyant scalars, save the current user source term
      call field_get_key_int(ivarfl(isca(ii)), key_buoyant_id, coupled_with_vel_p_fld)
      if (coupled_with_vel_p_fld.eq.1) then
        call add_source_term_field(ivarfl(isca(ii)))
      endif
      call add_source_term_prev_field(ivarfl(isca(ii)))
    endif
    ! Only usefull for Min/Max limiter
    call field_get_key_struct_var_cal_opt(ivarfl(isca(ii)), vcopt)
    if (vcopt%isstpc.eq.2) then
      call add_source_term_field(ivarfl(isca(ii)))
    endif
  enddo
endif

return

!===============================================================================
! Formats
!===============================================================================

 8101 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 1 OR 2               ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8111 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    WITH A SECOND ORDER SCHEME IN TIME: ISCHTP = ', I10      ,/,&
'@    IT IS NECESSARY TO USE A CONSTANT AND UNIFORM TIME STEP ',/,&
'@    BUT IDTVAR = ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8112 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2ND ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    WITH K-EPSILON (ITURB = ',I10,' )'                       ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8113 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   :      AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 1st ORDER SCHEME HAS BEEN IMPOSSED   (ISCHTP = ',I10   ,/,&
'@    FOR LES (ITURB = ',I10,' )'                              ,/,&
'@                                                            ',/,&
'@  The calculation will   be executed                        ',/,&
'@                                                            ',/,&
'@  It is recommended to verify  parameters.                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8114 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR PHI_FBAR (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8117 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA                    ',/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR BL-V2/K  (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8115 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR K-OMEGA   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-omega.                 ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8116 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR SPALART   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of Spalart-Allmaras.        ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8121 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8131 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8132 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  With the chosen turbulence model   , ITURB = ',I10         ,/,&
'@    the value of ISTO2T (extrapolation of the source terms  ',/,&
'@    for the turbulent variables) cannot be modified         ',/,&
'@    yet ISTO2T has been forced to ',I10                      ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8141 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITAL DATA FOR SCALARS    ',I10 ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 5. End
!===============================================================================

return
end subroutine varpos

!===============================================================================

!> \brief add field defining previous source term values for a given field
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     f_id          base field id
!_______________________________________________________________________________

subroutine add_source_term_prev_field &
 ( f_id )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: f_id

! Local variables

character(len=64) :: f_name

integer :: type_flag, location_id, st_id, f_dim
logical :: has_previous

!===============================================================================

type_flag = FIELD_EXTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells
has_previous = .false.

! Define associated field

call field_get_dim(f_id, f_dim)
call field_get_name (f_id, f_name)

call field_create(trim(f_name)//'_prev_st', type_flag,               &
                  location_id, f_dim, has_previous, st_id)

call field_set_key_int(f_id, kstprv, st_id)

call hide_property(st_id)

return

end subroutine add_source_term_prev_field


!===============================================================================

!> \brief add field defining current source term values for a given field
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     f_id          base field id
!_______________________________________________________________________________

subroutine add_source_term_field &
 ( f_id )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: f_id

! Local variables

character(len=64) :: f_name

integer :: type_flag, location_id, st_id, f_dim

!===============================================================================

type_flag = FIELD_EXTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells

! Define asscociated field

call field_get_dim(f_id, f_dim)
call field_get_name (f_id, f_name)

call field_find_or_create(trim(f_name)//'_st', type_flag,       &
                            location_id, f_dim, st_id)

call field_set_key_int(f_id, kst, st_id)

return

end subroutine add_source_term_field
