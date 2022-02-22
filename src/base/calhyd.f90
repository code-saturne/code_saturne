!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!> \file calhyd.f90
!> \brief Poisson equation resolution for hydrostatic pressure:
!>  \f$ \divs ( \grad P ) = \divs ( f ) \f$
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[out]    indhyd        indicateur de mise a jour de phydr
!> \param[in]     iterns        Navier-Stokes iteration number
!> \param[in]     fext          external force generating hydrostatic pressure
!> \param[in]     dfext         external force increment
!>                              generating hydrostatic pressure
!> \param[out]    phydr         hydrostatic pressure increment
!> \param[in]     flumas        work array
!> \param[in]     flumab        work array
!> \param[in,out] viscf         work array
!> \param[in,out] viscb         work array
!> \param[in,out] dam           work array
!> \param[in,out] xam           work array
!> \param[in,out] dpvar         work array
!> \param[in,out] rhs           work array
!______________________________________________________________________________

subroutine calhyd &
 ( indhyd , iterns ,                                              &
   fext   , dfext  ,                                              &
   phydr  , flumas , flumab ,                                     &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   dpvar  , rhs   )                                               &
  bind(C, name='cs_hydrostatic_pressure_compute')

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding
use atincl, only: iatmst
use paramx
use numvar
use entsor
use cstnum
use cstphy
use optcal
use parall
use period
use mesh
use field_operator
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer(c_int) :: indhyd
integer(c_int), value :: iterns

double precision fext(3,ncelet)
double precision dfext(3,ncelet)
double precision phydr(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac)
double precision dpvar(ncelet)
double precision rhs(ncelet)

! Local variables

character(len=80) :: chaine
integer          lchain
integer          f_id, iccocg, inc   , init  , isym
integer          iel   , ical, ifac
integer          nswmpr
integer          isweep, niterf
integer          iphydp
integer          imrgrp, nswrgp, imligp, iwarnp, iwgrp
integer          idiffp, iconvp, ndircp, imvisp
integer          ibsize, iesize
integer          imucpp, f_id0

double precision ressol, residu, rnorm , rnrmf , rnrmdf
double precision epsrgp, climgp, extrap, epsilp
double precision precre, precab, thetap
double precision qimp, hint

double precision rvoid(1)

double precision, allocatable, dimension(:) :: rovsdt, div_fext, viscce
double precision, allocatable, dimension(:,:) :: next_fext
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp

type(solving_info) sinfo
type(var_cal_opt) :: vcopt_pr

!===============================================================================

!===============================================================================
! 1.  Initialisations
!===============================================================================

! Allocate temporary arrays
allocate(rovsdt(ncelet), div_fext(ncelet), viscce(ncelet))
allocate(next_fext(3, ncelet))

! Field id
call field_get_id_try("hydrostatic_pressure", f_id)

! Get variables calculation options
call field_get_key_struct_var_cal_opt(f_id, vcopt_pr)
call field_get_key_struct_solving_info(f_id, sinfo)

if (iterns.le.1) then
  sinfo%nbivar = 0
endif

call field_get_coefa_s (f_id, coefap)
call field_get_coefb_s (f_id, coefbp)
call field_get_coefaf_s(f_id, cofafp)
call field_get_coefbf_s(f_id, cofbfp)

chaine = 'hydrostatic_p'
lchain = 16

! Symmetric
isym  = 1

! Matrix block size
ibsize = 1
iesize = 1

!     TEST DE VARIATION DE LA PRESSION HYDROSTATIQUE EN SORTIE

!     on regarde si les terme source ont varie
!     on ne passe dans calhyd que si on a des faces de sortie std
!     la precision pour les tests est a peu pres arbitraire.
precre = sqrt(epzero)
precab = 1.d2*epzero

ical = 0
do iel = 1, ncel
  rnrmf  = fext(1,iel)**2+fext(2,iel)**2+fext(3,iel)**2
  rnrmdf = dfext(1,iel)**2+dfext(2,iel)**2+dfext(3,iel)**2
  if ((rnrmdf.ge.precre*rnrmf).and.(rnrmdf.ge.precab)) then
    ical = 1
  endif
enddo
if (irangp.ge.0) then
  call parcpt (ical)
endif

if (iatmst.eq.0) then
  if (ical.eq.0) then
    indhyd = 0
    return
  endif
endif

if (mod(ntcabs,ntlist).eq.0.or.vcopt_pr%iwarni.ge.1) write(nfecra,1000)

 1000 format(                                                           &
'  Hydrostatic pressure computation: ',/,                   &
'         updating the Dirichlets at the end (CALHYD)',/)

indhyd = 1

f_id0 = -1

do iel = 1, ncel
  next_fext(1 ,iel) = fext(1 ,iel) * cell_is_active(iel) + dfext(1 ,iel)
  next_fext(2 ,iel) = fext(2 ,iel) * cell_is_active(iel) + dfext(2 ,iel)
  next_fext(3 ,iel) = fext(3 ,iel) * cell_is_active(iel) + dfext(3 ,iel)
enddo

! Parallel or periodicity synchronization
call synvin(next_fext)

!===============================================================================
! 2. Prepare matrix and boundary conditions
!===============================================================================

! Boundary conditions

!  On resout avec des CL de flux nul partout
ndircp = 0

do ifac = 1, nfabor
  ! LOCAL Neumann Boundary Conditions on the hydrostatic pressure
  !--------------------------------------------------------------

  hint = 1.d0/distb(ifac)
  qimp = 0.d0

  call set_neumann_scalar &
    !==================
  ( coefap(ifac), cofafp(ifac),             &
    coefbp(ifac), cofbfp(ifac),             &
    qimp        , hint )

enddo

! Unsteady term
do iel = 1, ncel
  rovsdt(iel) = 0.d0
enddo

! Face diffusivity
do iel = 1, ncel
  viscce(iel) = 1.d0
enddo

imvisp = vcopt_pr%imvisf

call viscfa &
 ( imvisp ,                                                       &
   viscce ,                                                       &
   viscf  , viscb  )

iconvp = vcopt_pr%iconv
idiffp = vcopt_pr%idiff
thetap = vcopt_pr%thetav
imucpp = 0

call matrix &
!==========
 ( iconvp , idiffp , ndircp , isym ,                              &
   thetap , imucpp ,                                              &
   coefbp , cofbfp , rovsdt ,                                     &
   flumas , flumab , viscf  , viscb  ,                            &
   rvoid  , dam    , xam    )

!===============================================================================
! 3. Compute right hand side
!===============================================================================

init   = 1
inc    = 0
iccocg = 1
nswrgp = vcopt_pr%nswrgr
imligp = vcopt_pr%imligr
iwarnp = vcopt_pr%iwarni
epsrgp = vcopt_pr%epsrgr
climgp = vcopt_pr%climgr

! Compute div(f_ext^n+1)
call projts &
!==========
 ( init   , nswrgp ,                                              &
   next_fext,                                                     &
   cofbfp ,                                                       &
   flumas, flumab ,                                               &
   viscf  , viscb  ,                                              &
   viscce , viscce , viscce    )

init = 1
call divmas(init,flumas,flumab, div_fext)
rnorm = sqrt(cs_gdot(ncel,div_fext,div_fext))

!===============================================================================
! 4.  BOUCLES SUR LES NON ORTHOGONALITES (RESOLUTION)
!===============================================================================

! --- Nombre de sweeps
nswmpr = vcopt_pr%nswrsm

! --- Mise a zero des variables
!       dpvar      sera l'increment d'increment a chaque sweep
!       div_fext         sera la divergence du flux de masse predit
do iel = 1, ncel
  dpvar(iel) = 0.d0
enddo

! Reconstruction loop (beginning)
!--------------------------------

do isweep = 1, nswmpr

  ! --- Update the right hand side and update the residual
  !      rhs^{k+1} = - div(fext^n+1) - D(1, p_h^{k+1})
  !-------------------------------------------------------------

    iccocg = 1
    init = 1 ! re-init rhs to 0 if init = 1
    inc  = 1
    imrgrp = vcopt_pr%imrgra
    nswrgp = vcopt_pr%nswrgr
    imligp = vcopt_pr%imligr
    iwgrp  = vcopt_pr%iwgrec
    iwarnp = vcopt_pr%iwarni
    epsrgp = vcopt_pr%epsrgr
    climgp = vcopt_pr%climgr
    extrap = 0.d0
    iphydp = 1

    call itrgrp &
    !==========
 ( f_id0  , init   , inc    , imrgrp , iccocg , nswrgp , imligp , iphydp ,     &
   iwgrp  , iwarnp ,                                                           &
   epsrgp , climgp , extrap ,                                                  &
   next_fext,                                                                  &
   phydr  ,                                                                    &
   coefap , coefbp ,                                                           &
   cofafp , cofbfp ,                                                           &
   viscf  , viscb  ,                                                           &
   viscce ,                                                                    &
   rhs   )


  do iel = 1, ncel
    rhs(iel) = - div_fext(iel) - rhs(iel)
  enddo

  ! --- Convergence test
  residu = sqrt(cs_gdot(ncel,rhs,rhs))
  if (vcopt_pr%iwarni.ge.2) then
    write(nfecra,1400) chaine(1:16), isweep,residu, rnorm
  endif

  if (isweep.eq.1) then
    sinfo%rnsmbr = residu
  endif

  if (residu.le.vcopt_pr%epsrsm*rnorm) then
    !     Si convergence,  sortie

    goto 101

  endif

  ! Solving on the increment
  !-------------------------
  do iel = 1, ncel
    dpvar(iel) = 0.d0
  enddo

  iwarnp = vcopt_pr%iwarni
  epsilp = vcopt_pr%epsilo

  ! Solver residual
  ressol = residu

  call sles_solve_native(f_id, '',                              &
                         isym, ibsize, iesize, dam, xam,          &
                         epsilp, rnorm, niterf, ressol, rhs, dpvar)

  ! Writing
  sinfo%nbivar = sinfo%nbivar + niterf

  ! Update the increment of pressure
  !---------------------------------

  do iel = 1, ncel
    phydr(iel) = phydr(iel) + dpvar(iel)
  enddo

enddo

! --- Reconstruction loop (end)

if (vcopt_pr%iwarni.ge.2) then
  write(nfecra,1600) chaine(1:16), nswmpr
endif

 101  continue

 ! For logging
if (abs(rnorm).gt.0.d0) then
  sinfo%resvar = residu/rnorm
else
  sinfo%resvar = 0.d0
endif

! Save convergence info
call field_set_key_struct_solving_info(f_id, sinfo)

! Free memory
deallocate(rovsdt, div_fext, viscce)
deallocate(next_fext)

!===============================================================================
! 5. Free solver setup
!===============================================================================

call sles_free_native(f_id, '')

!--------
! Formats
!--------

 1400 format(1X,A16,' : sweep = ',I5,' RHS residual = ',E14.6, ' NORM =',E14.6)
 1600 format(                                                           &
'@                                                            ',/,&
'@ @@ WARNING: ', A16,' HYDROSTATIC PRESSURE STEP             ',/,&
'@    ========                                                ',/,&
'@  Maximum number of iterations ',I10   ,' reached           ',/,&
'@                                                            '  )

!----
! End
!----

return

end subroutine
