!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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
!> \param[in]     coefap        boundary conditions coefficient
!> \param[in]     coefbp        boundary conditions coefficient
!> \param[in]     cofafp        boundary conditions coefficient
!> \param[in]     cofbfp        boundary conditions coefficient
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
   coefap , coefbp ,                                              &
   cofafp , cofbfp ,                                              &
   viscf  , viscb  ,                                              &
   dam    , xam    ,                                              &
   dpvar  , smbr   )

!===============================================================================
! Module files
!===============================================================================

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
double precision coefap(nfabor), coefbp(nfabor)
double precision cofafp(nfabor), cofbfp(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision dam(ncelet), xam(nfac)
double precision dpvar(ncelet)
double precision smbr(ncelet)

! Local variables

character(len=80) :: chaine
integer          lchain
integer          f_id, iccocg, inc   , init  , isym
integer          iel   , ical
integer          nswmpr
integer          isweep, niterf
integer          iphydp
integer          nswrgp, imligp, iwarnp
integer          iinvpe
integer          idiffp, iconvp, ndircp
integer          ibsize, iesize
integer          imucpp, f_id0

double precision residu, rnorm , rnrmf , rnrmdf
double precision epsrgp, climgp, extrap, epsilp
double precision precre, precab, thetap

double precision rvoid(1)

double precision, allocatable, dimension(:) :: w1, w7, w10

type(var_cal_opt) :: vcopt_u, vcopt_pr

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

! Allocate temporary arrays
allocate(w1(ncelet), w7(ncelet), w10(ncelet))

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
if (ical.eq.0) then
  do iel = 1,ncel
    phydr(iel) = 0.d0
  enddo
  indhyd = 0
  return
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
! 2.  PREPARATION DE LA MATRICE DU SYSTEME A RESOUDRE
!===============================================================================

! ---> TERME INSTATIONNAIRE

do iel = 1, ncel
  w1(iel) = 0.d0
enddo

! ---> "VITESSE" DE DIFFUSION FACETTE

do iel = 1, ncel
  w10(iel) = 1.d0
enddo

call viscfa                                                       &
!==========
 ( imvisf ,                                                       &
   w10    ,                                                       &
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
   coefbp , cofbfp , w1     ,                                     &
   flumas , flumab , viscf  , viscb  ,                            &
   rvoid  , dam    , xam    )

!===============================================================================
! 3.  INITIALISATION DU FLUX DE MASSE
!===============================================================================


!     PROJECTION AUX FACES DES TERMES SOURCES

init   = 1
inc    = 0
iccocg = 1
nswrgp = vcopt_pr%nswrgr
imligp = vcopt_pr%imligr
iwarnp = vcopt_pr%iwarni
epsrgp = vcopt_pr%epsrgr
climgp = vcopt_pr%climgr

call projts                                                       &
!==========
 ( init   , nswrgp ,                                              &
   dfext  ,                                                       &
   cofbfp ,                                                       &
   flumas, flumab ,                                               &
   viscf  , viscb  ,                                              &
   w10    , w10    , w10    )

init = 1
call divmas(init,flumas,flumab,w7)
rnorm = sqrt(cs_gdot(ncel,w7,w7))

!===============================================================================
! 4.  BOUCLES SUR LES NON ORTHOGONALITES (RESOLUTION)
!===============================================================================

! --- Nombre de sweeps
nswmpr = vcopt_pr%nswrsm

! --- Mise a zero des variables
!       RTP(.,IPR) sera l'increment de pression cumule
!       DPVAR      sera l'increment d'increment a chaque sweep
!       W7         sera la divergence du flux de masse predit
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
    smbr(iel) = - w7(iel) - smbr(iel)
  enddo

! --- Test de convergence du calcul

  residu = sqrt(cs_gdot(ncel,smbr,smbr))
  if (vcopt_pr%iwarni.ge.2) then
    write(nfecra,1400)chaine(1:16),isweep,residu
  endif

!MO IL FAUDRA VERIFIER LA PERTINENCE DU TEST

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
  iinvpe = 1
  ibsize = 1
  iesize = 1

  call sles_solve_native(f_id, chaine,                            &
                         isym, ibsize, iesize, dam, xam, iinvpe,  &
                         epsilp, rnorm, niterf, residu, smbr, dpvar)

  if( isweep.eq.nswmpr ) then
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
 ( f_id0  , init   , inc    , imrgra , iccocg , nswrgp , imligp , iphydp , &
   iwarnp ,                                                                    &
   epsrgp , climgp , extrap ,                                                  &
   dfext  ,                                                                    &
   phydr  ,                                                                    &
   coefap , coefbp ,                                                           &
   cofafp , cofbfp ,                                                           &
   viscf  , viscb  ,                                                           &
   w10         ,                                                               &
   smbr   )

  endif

enddo
! --- Boucle de reconstruction : fin

if(vcopt_pr%iwarni.ge.2) then
   write( nfecra,1600)chaine(1:16),nswmpr
endif

 101  continue

! Free memory
deallocate(w1, w7, w10)

!===============================================================================
! 5. Free solver setup
!===============================================================================

call sles_free_native(-1, chaine)

!--------
! FORMATS
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
! FIN
!----

return

end subroutine
