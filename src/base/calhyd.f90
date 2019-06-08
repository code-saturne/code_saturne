!-------------------------------------------------------------------------------

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
!> \param[in,out] smbr          work array
!______________________________________________________________________________

subroutine calhyd &
 ( indhyd ,                                                       &
   fext   , dfext  ,                                              &
   phydr  , flumas , flumab ,                                     &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   dpvar  , smbr   )

!===============================================================================
! Module files
!===============================================================================

use atincl, only: iatmst
use dimens, only: ndimfb
use paramx
use numvar
use entsor
use cstnum
use optcal
use parall
use period
use mesh
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          indhyd

double precision fext(3,ncelet)
double precision dfext(3,ncelet)
double precision phydr(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac)
double precision dpvar(ncelet)
double precision smbr(ncelet)

! Local variables

character(len=80) :: chaine
integer          lchain
integer          f_id, iccocg, inc   , init  , isym
integer          iel   , ical, ifac
integer          nswmpr
integer          isweep, niterf
integer          iphydp
integer          nswrgp, imligp, iwarnp
integer          idiffp, iconvp, ndircp
integer          ibsize, iesize
integer          imucpp, f_id0

double precision residu, rnorm , rnrmf , rnrmdf
double precision epsrgp, climgp, extrap, epsilp
double precision precre, precab, thetap
double precision qimp, hint

double precision rvoid(1)

double precision, allocatable, dimension(:) :: rovsdt, div_dfext, viscce
double precision, allocatable, dimension(:) :: coefap, coefbp
double precision, allocatable, dimension(:) :: cofafp, cofbfp

type(var_cal_opt) :: vcopt_u, vcopt_pr

!===============================================================================

!===============================================================================
! 1.  Initialisations
!===============================================================================

! Allocate temporary arrays
allocate(rovsdt(ncelet), div_dfext(ncelet), viscce(ncelet))

! Get variables calculation options

call field_get_key_struct_var_cal_opt(ivarfl(iu), vcopt_u)
call field_get_key_struct_var_cal_opt(ivarfl(ipr), vcopt_pr)

! Variable name

f_id = -1
chaine = 'hydrostatic_p'
lchain = 16

! Symetric
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
    do iel = 1,ncel
      phydr(iel) = 0.d0
    enddo
    indhyd = 0
    return
  endif
endif

if (mod(ntcabs,ntlist).eq.0.or.vcopt_u%iwarni.ge.0) write(nfecra,1000)

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'  Calcul de la pression hydrostatique : ',/,               &
'         mise a jour des Dirichlets en sortie (CALHYD)',/)

#else

 1000 format(                                                           &
'  Hydrostatic pressure computation: ',/,                   &
'         updating the Dirichlets at the end (CALHYD)',/)

#endif

indhyd = 1

f_id0 = -1

!===============================================================================
! 2. Prepare matrix and boundary conditions
!===============================================================================

! Boundary conditions

! Work arrays for BCs
allocate(coefap(ndimfb), cofafp(ndimfb))
allocate(coefbp(ndimfb), cofbfp(ndimfb))

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

call viscfa &
 ( imvisf ,                                                       &
   viscce ,                                                       &
   viscf  , viscb  )


iconvp = 0
idiffp = 1
!  On resout avec des CL de flux nul partout
ndircp = 0

thetap = 1.d0
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

call projts &
!==========
 ( init   , nswrgp ,                                              &
   dfext  ,                                                       &
   cofbfp ,                                                       &
   flumas, flumab ,                                               &
   viscf  , viscb  ,                                              &
   viscce , viscce , viscce    )

init = 1
call divmas(init,flumas,flumab,div_dfext)
rnorm = sqrt(cs_gdot(ncel,div_dfext,div_dfext))

!===============================================================================
! 4.  BOUCLES SUR LES NON ORTHOGONALITES (RESOLUTION)
!===============================================================================

! --- Nombre de sweeps
nswmpr = vcopt_pr%nswrsm

! --- Mise a zero des variables
!       phydr      sera l'increment de pression cumule
!       dpvar      sera l'increment d'increment a chaque sweep
!       div_dfext         sera la divergence du flux de masse predit
do iel = 1,ncel
  phydr(iel) = 0.d0
  dpvar(iel) = 0.d0
  smbr(iel) = 0.d0
enddo


! --- Boucle de reconstruction : debut
do isweep = 1, nswmpr

! --- Mise a jour du second membre
!     (signe "-" a cause de celui qui est implicitement dans la matrice)
  do iel = 1, ncel
    smbr(iel) = - div_dfext(iel) - smbr(iel)
  enddo

! --- Test de convergence du calcul

  residu = sqrt(cs_gdot(ncel,smbr,smbr))
  if (vcopt_pr%iwarni.ge.2) then
    write(nfecra,1400)chaine(1:16),isweep,residu
  endif

  ! FIXME coherent with resopv!
  if (residu.le.10.d0*vcopt_pr%epsrsm*rnorm) then
!     Si convergence,  sortie

    goto 101

  endif

! --- Resolution implicite sur l'increment d'increment DPVAR
  do iel = 1, ncel
    dpvar(iel) = 0.d0
  enddo

  iwarnp = vcopt_pr%iwarni
  epsilp = vcopt_pr%epsilo
  ibsize = 1
  iesize = 1

  call sles_solve_native(f_id, chaine,                            &
                         isym, ibsize, iesize, dam, xam,          &
                         epsilp, rnorm, niterf, residu, smbr, dpvar)

  if (isweep.eq.nswmpr ) then
!     Mise a jour de l'increment de pression
    do iel = 1, ncel
      phydr(iel) = phydr(iel) + dpvar(iel)
    enddo


  else

! --- Si ce n'est pas le dernier sweep
!       Mise a jour de l'increment de pression et calcul direct de la
!       partie en gradient d'increment de pression du second membre
!       (avec reconstruction)

    do iel = 1, ncel
      phydr(iel) = phydr(iel) + dpvar(iel)
    enddo

    iccocg = 1
    init = 1
    inc  = 1
    nswrgp = vcopt_pr%nswrgr
    imligp = vcopt_pr%imligr
    iwarnp = vcopt_pr%iwarni
    epsrgp = vcopt_pr%epsrgr
    climgp = vcopt_pr%climgr
    extrap = 0.d0
    iphydp = 1

    call itrgrp &
    !==========
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp ,     &
   iwarnp ,                                                                    &
   epsrgp , climgp , extrap ,                                                  &
   dfext  ,                                                                    &
   phydr  ,                                                                    &
   coefap , coefbp ,                                                           &
   cofafp , cofbfp ,                                                           &
   viscf  , viscb  ,                                                           &
   viscce ,                                                                    &
   smbr   )

  endif

enddo

! --- Boucle de reconstruction : fin

if(vcopt_pr%iwarni.ge.2) then
   write( nfecra,1600)chaine(1:16),nswmpr
endif

 101  continue

! Free memory
deallocate(rovsdt, div_dfext, viscce)
deallocate(coefap, cofafp)
deallocate(coefbp, cofbfp)

!===============================================================================
! 5. Free solver setup
!===============================================================================

call sles_free_native(-1, chaine)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1400 format(1X,A16,' : SWEEP = ',I5,' NORME SECOND MEMBRE = ',E14.6)
 1600 format(                                                           &
'@                                                            ',/,&
'@ @@ ATTENTION : ', A16,' ETAPE DE PRESSION HYDROSTATIQUE    ',/,&
'@    =========                                               ',/,&
'@  Nombre d''iterations maximal ',I10   ,' atteint           ',/,&
'@                                                            '  )

#else

 1400 format(1X,A16,' : SWEEP = ',I5,' RIGHT HAND SIDE NORM = ',E14.6)
 1600 format(                                                           &
'@                                                            ',/,&
'@ @@ WARNING: ', A16,' HYDROSTATIC PRESSURE STEP             ',/,&
'@    ========                                                ',/,&
'@  Maximum number of iterations ',I10   ,' reached           ',/,&
'@                                                            '  )

#endif

!----
! End
!----

return

end subroutine
