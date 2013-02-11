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

!> \file matrix.f90
!>
!> \brief This function builds the matrix of advection/diffusion for a scalar
!> field.
!>
!> The advection is upwind, the diffusion is not reconstructed.
!> The matrix is splitted into a diagonal block (number of cells)
!> and an extra diagonal part (of dimension 2 time the number of internal
!> faces).
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     nfac          number of interior faces
!> \param[in]     nfabor        number of boundary faces
!> \param[in]     iconvp        indicator
!>                               - 1 advection
!>                               - 0 otherwise
!> \param[in]     idiffp        indicator
!>                               - 1 diffusion
!>                               - 0 otherwise
!> \param[in]     ndircp        indicator
!>                               - 0 if the diagonal stepped aside
!> \param[in]     isym          indicator
!>                               - 1 symmetric matrix
!>                               - 2 non symmmetric matrix
!> \param[in]     thetap        weightening coefficient for the theta-schema,
!>                               - thetap = 0: explicit scheme
!>                               - thetap = 0.5: time-centred
!>                               scheme (mix between Crank-Nicolson and
!>                               Adams-Bashforth)
!>                               - thetap = 1: implicit scheme
!> \param[in]     imucpp        indicator
!>                               - 0 do not multiply the convectiv term by Cp
!>                               - 1 do multiply the convectiv term by Cp
!> \param[in]     ifacel        cell indexes of interior faces
!> \param[in]     ifabor        no de l'elt voisin d'une face de bord
!> \param[in]     coefbp        boundary condition array for the variable
!>                               (Impplicit part)
!> \param[in]     cofbfp        boundary condition array for the variable flux
!>                               (Impplicit part)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at border faces
!> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the matrix
!> \param[in]     viscb         \f$ S_\fib \f$
!>                               at border faces for the matrix
!> \param[in]     xcpp          array of specific heat (Cp)
!> \param[out]    da            diagonal part of the matrix
!> \param[out]    xa            extra interleaved diagonal part of the matrix
!_______________________________________________________________________________


subroutine matrix &
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp , isym   , nfecra ,                   &
   thetap , imucpp ,                                              &
   ifacel , ifabor ,                                              &
   coefbp , cofbfp ,rovsdt , flumas , flumab , viscf  , viscb  ,  &
   xcpp   , da     , xa     )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use parall

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , nfac   , nfabor
integer          iconvp , idiffp , ndircp , isym
integer          nfecra
integer          imucpp
double precision thetap

integer          ifacel(2,nfac), ifabor(nfabor)
double precision coefbp(nfabor), cofbfp(nfabor), rovsdt(ncelet)
double precision flumas(nfac), flumab(nfabor)
double precision viscf(nfac), viscb(nfabor)
double precision da(ncelet ), xa(nfac ,isym)
double precision xcpp(ncelet)

! Local variables

integer          ifac, ii, jj, iel, ig, it
double precision flui, fluj, epsi

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

if (isym.ne.1.and.isym.ne.2) then
  write(nfecra,1000) isym
  call csexit (1)
endif

epsi = 1.d-7

!$omp parallel do
do iel = 1, ncel
  da(iel) = rovsdt(iel)
enddo
if (ncelet.gt.ncel) then
  !$omp parallel do if (ncelet - ncel > thr_n_min)
  do iel = ncel+1, ncelet
    da(iel) = 0.d0
  enddo
endif

if (isym.eq.2) then
  !$omp parallel do
  do ifac = 1, nfac
    xa(ifac,1) = 0.d0
    xa(ifac,2) = 0.d0
  enddo
else
  !$omp parallel do
  do ifac = 1, nfac
    xa(ifac,1) = 0.d0
  enddo
endif

! When solving the temperature, the convective part is multiplied by Cp
if (imucpp.eq.0) then

!===============================================================================
! 2.    CALCUL DES TERMES EXTRADIAGONAUX
!===============================================================================

  if (isym.eq.2) then

    !$omp parallel do firstprivate(thetap, iconvp, idiffp) private(flui, fluj)
    do ifac = 1, nfac
      flui = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))
      fluj =-0.5d0*(flumas(ifac) +abs(flumas(ifac)))
      xa(ifac,1) = thetap*(iconvp*flui -idiffp*viscf(ifac))
      xa(ifac,2) = thetap*(iconvp*fluj -idiffp*viscf(ifac))
    enddo

  else

    !$omp parallel do firstprivate(thetap, iconvp, idiffp) private(flui)
    do ifac = 1, nfac
      flui = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )
      xa(ifac,1) = thetap*(iconvp*flui -idiffp*viscf(ifac))
    enddo

  endif

!===============================================================================
! 3.     CONTRIBUTION DES TERMES X-TRADIAGONAUX A LA DIAGONALE
!===============================================================================

  if (isym.eq.2) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)
          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)
          da(ii) = da(ii) - xa(ifac,1)
          da(jj) = da(jj) - xa(ifac,2)
        enddo
      enddo
    enddo

  else

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)
          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)
          da(ii) = da(ii) - xa(ifac,1)
          da(jj) = da(jj) - xa(ifac,1)
        enddo
      enddo
    enddo

  endif

!===============================================================================
! 4.     CONTRIBUTION DES FACETTES DE BORDS A LA DIAGONALE
!===============================================================================

  do ig = 1, ngrpb
    !$omp parallel do firstprivate(thetap, iconvp, idiffp) &
    !$omp          private(ifac, ii, flui) if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)
        ii = ifabor(ifac)
        flui = 0.5d0*(flumab(ifac) - abs(flumab(ifac)))
        da(ii) = da(ii) + thetap*( iconvp*flui*(coefbp(ifac)-1.d0)         &
                                 + idiffp*viscb(ifac)*cofbfp(ifac))
      enddo
    enddo
  enddo

! When solving the temperature, the convective part is multiplied by Cp
else

!===============================================================================
! 2.    CALCUL DES TERMES EXTRADIAGONAUX
!===============================================================================

  if (isym.eq.2) then

    !$omp parallel do firstprivate(thetap, iconvp, idiffp) private(flui, fluj)
    do ifac = 1, nfac
      flui = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))
      fluj =-0.5d0*(flumas(ifac) +abs(flumas(ifac)))
      xa(ifac,1) = thetap*(iconvp*xcpp(ifacel(1,ifac))*flui -idiffp*viscf(ifac))
      xa(ifac,2) = thetap*(iconvp*xcpp(ifacel(2,ifac))*fluj -idiffp*viscf(ifac))
    enddo

  else

    !$omp parallel do firstprivate(thetap, iconvp, idiffp) private(flui)
    do ifac = 1, nfac
      flui = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )
      xa(ifac,1) = thetap*(iconvp*xcpp(ifacel(1,ifac))*flui -idiffp*viscf(ifac))
    enddo

  endif

!===============================================================================
! 3.     CONTRIBUTION DES TERMES X-TRADIAGONAUX A LA DIAGONALE
!===============================================================================

  if (isym.eq.2) then

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)
          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)
          da(ii) = da(ii) - xa(ifac,1)
          da(jj) = da(jj) - xa(ifac,2)
        enddo
      enddo
    enddo

  else

    do ig = 1, ngrpi
      !$omp parallel do private(ifac, ii, jj)
      do it = 1, nthrdi
        do ifac = iompli(1,ig,it), iompli(2,ig,it)
          ii = ifacel(1,ifac)
          jj = ifacel(2,ifac)
          da(ii) = da(ii) - xa(ifac,1)
          da(jj) = da(jj) - xa(ifac,1)
        enddo
      enddo
    enddo

  endif

!===============================================================================
! 4.     CONTRIBUTION DES FACETTES DE BORDS A LA DIAGONALE
!===============================================================================

  do ig = 1, ngrpb
    !$omp parallel do firstprivate(thetap, iconvp, idiffp) &
    !$omp          private(ifac, ii, flui) if(nfabor > thr_n_min)
    do it = 1, nthrdb
      do ifac = iomplb(1,ig,it), iomplb(2,ig,it)
        ii = ifabor(ifac)
        flui = 0.5d0*(flumab(ifac) - abs(flumab(ifac)))
        da(ii) = da(ii) + thetap*( iconvp*flui*xcpp(ii)*(coefbp(ifac)-1.d0)   &
                                 + idiffp*viscb(ifac)*cofbfp(ifac))
      enddo
    enddo
  enddo

endif

!===============================================================================
! 5.  NON PRESENCE DE PTS DIRICHLET --> LEGER RENFORCEMENT DE LA
!     DIAGONALE POUR DECALER LE SPECTRE DES VALEURS PROPRES
!===============================================================================
!     (si IDIRCL=0, on a force NDIRCP a valoir au moins 1 pour ne pas
!      decaler la diagonale)

if (ndircp.le.0) then
  !$omp parallel do firstprivate(epsi)
  do iel=1,ncel
    da(iel) = (1.d0+epsi)*da(iel)
  enddo
endif

!--------
! FORMATS
!--------

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS matrix                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE matrix              AVEC ISYM   = ',I10        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ WARNING: ABORT IN matrix'                                ,/,&
'@    ========'                                                ,/,&
'@     matrix CALLED                WITH ISYM   = ',I10        ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run.'                          ,/,&
'@'                                                            ,/,&
'@  Contact support.'                                          ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

!----
! FIN
!----

return

end subroutine
