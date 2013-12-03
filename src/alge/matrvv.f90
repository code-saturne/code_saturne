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

!> \file matrvv.f90
!>
!> \brief This function builds the matrix of advection/diffusion for a vector
!> field with a tensorial diffusivity.
!>
!> The advection is upwind, the diffusion is not reconstructed.
!> The matrix is splitted into a diagonal block (3x3 times number of cells)
!> and an extra diagonal part (of dimension 2 times 3x3 the number of internal
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
!> \param[in]     ifacel        cell indexes of interior faces
!> \param[in]     ifabor        no de l'elt voisin d'une face de bord
!> \param[in]     coefbu        boundary condition array for the variable
!>                               (Impplicit part - 3x3 tensor array)
!> \param[in]     cofbfu        boundary condition array for the variable flux
!>                               (Impplicit part - 3x3 tensor array)
!> \param[in]     flumas        mass flux at interior faces
!> \param[in]     flumab        mass flux at border faces
!> \param[in]     viscf         \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$
!>                               at interior faces for the matrix
!> \param[in]     viscb         \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$
!>                               at border faces for the matrix
!> \param[out]    da            diagonal part of the matrix
!> \param[out]    xa            extra interleaved diagonal part of the matrix
!_______________________________________________________________________________

subroutine matrvv &
 ( ncelet , ncel   , nfac   , nfabor ,                            &
   iconvp , idiffp , ndircp , isym   , nfecra ,                   &
   thetap ,                                                       &
   ifacel , ifabor ,                                              &
   coefbu , cofbfu , fimp   , flumas , flumab , viscf  , viscb  , &
   da     , xa     )

!===============================================================================
! Module files
!===============================================================================

use parall
use mesh, only:surfbn

!===============================================================================

implicit none

! Arguments

integer          ncelet , ncel   , nfac   , nfabor
integer          iconvp , idiffp , ndircp , isym
integer          nfecra
double precision thetap

integer          ifacel(2,nfac), ifabor(nfabor)
double precision coefbu(3,3,nfabor), fimp(3,3,ncelet), cofbfu(3,3,nfabor)
double precision flumas(nfac), flumab(nfabor)
double precision viscf(3,3,nfac), viscb(nfabor)
double precision da(3,3,ncelet),xa(3,3,isym,nfac)

! Local variables

integer          ifac,ii,jj,iel, isou, jsou
double precision flui,fluj,epsi

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

if(isym.ne.1.and.isym.ne.2) then
  write(nfecra,1000) isym
  call csexit (1)
endif

epsi = 1.d-7

do iel = 1, ncel
  do isou = 1, 3
    do jsou = 1, 3
      da(isou,jsou,iel) = fimp(isou,jsou,iel)
    enddo
  enddo
enddo
if(ncelet.gt.ncel) then
  do iel = ncel+1, ncelet
    do isou = 1, 3
      do jsou = 1, 3
        da(isou,jsou,iel) = 0.d0
      enddo
    enddo
  enddo
endif

if(isym.eq.2) then
  do ifac = 1, nfac
    do isou = 1, 3
      do jsou = 1, 3
        xa(isou,jsou,1,ifac) = 0.d0
        xa(isou,jsou,2,ifac) = 0.d0
      enddo
    enddo
  enddo
else
  do ifac = 1, nfac
    do isou = 1, 3
      do jsou = 1, 3
        xa(isou,jsou,1,ifac) = 0.d0
      enddo
    enddo
  enddo
endif

!===============================================================================
! 2. Computation of extradiagonal terms
!===============================================================================

if(isym.eq.2) then

  do ifac = 1, nfac
    flui = 0.5d0*( flumas(ifac) -abs(flumas(ifac)) )
    fluj =-0.5d0*( flumas(ifac) +abs(flumas(ifac)) )
    do isou = 1, 3
      xa(isou,isou,1,ifac) = iconvp*flui
      xa(isou,isou,2,ifac) = iconvp*fluj
      do jsou = 1, 3
        xa(isou,jsou,1,ifac) = thetap*( xa(isou,jsou,1,ifac)          &
                                      - idiffp*viscf(isou,jsou,ifac))
        xa(isou,jsou,2,ifac) = thetap*( xa(isou,jsou,2,ifac)          &
                                      - idiffp*viscf(isou,jsou,ifac))
      enddo
    enddo
  enddo

else

  do ifac = 1, nfac
    flui = 0.5d0*(flumas(ifac) -abs(flumas(ifac)))
    do isou = 1, 3
      xa(isou,isou,1,ifac) = iconvp*flui
      do jsou = 1, 3
        xa(isou,jsou,1,ifac) = thetap*( xa(isou,jsou,1,ifac)          &
                                      - idiffp*viscf(isou,jsou,ifac))
      enddo
    enddo
  enddo

endif

!===============================================================================
! 3. Contribution of the extra-diagonal terms to the diagonal
!===============================================================================

if (isym.eq.2) then

  do ifac = 1, nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    do isou = 1, 3
      do jsou = 1, 3
        da(isou,jsou,ii) = da(isou,jsou,ii) - xa(isou,jsou,1,ifac)
        da(isou,jsou,jj) = da(isou,jsou,jj) - xa(isou,jsou,2,ifac)
      enddo
    enddo
  enddo

else

  do ifac = 1,nfac
    ii = ifacel(1,ifac)
    jj = ifacel(2,ifac)
    do isou = 1, 3
      do jsou = 1, 3
        da(isou,jsou,ii) = da(isou,jsou,ii) -xa(isou,jsou,1,ifac)
        da(isou,jsou,jj) = da(isou,jsou,jj) -xa(isou,jsou,1,ifac)
      enddo
    enddo
  enddo

endif

!===============================================================================
! 4. Contribution of border faces to the diagonal
!===============================================================================

do ifac = 1,nfabor
  ii = ifabor(ifac)
  flui = 0.5d0*( flumab(ifac) -abs(flumab(ifac)) )
  do isou = 1, 3
    do jsou = 1, 3
      if(isou.eq.jsou) then
        da(isou,jsou,ii) = da(isou,jsou,ii) + thetap*(                    &
                       iconvp*flui*(coefbu(isou,jsou,ifac)-1.d0)          &
                      +idiffp*viscb(ifac)*cofbfu(isou,jsou,ifac)          &
                             )
      else
        da(isou,jsou,ii) = da(isou,jsou,ii) + thetap*(                    &
                       iconvp*( flui*coefbu(isou,jsou,ifac) )             &
                      +idiffp*viscb(ifac)*cofbfu(isou,jsou,ifac)          &
                             )
      endif
    enddo
  enddo
enddo


!===============================================================================
! 5. If no Dirichlet condition, the diagonal is slightly increased so that
!    the eigenvalues are stepped aside
!===============================================================================
!     (si IDIRCL=0, on a force NDIRCP a valoir au moins 1 pour ne pas
!      decaler la diagonale)

if ( ndircp.le.0 ) then
  do iel=1,ncel
    do isou = 1, 3
      da(isou,isou,iel) = (1.d0+epsi)*da(isou,isou,iel)
    enddo
  enddo
endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS matrxv                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE matrxv              AVEC ISYM   = ',I10        ,/,&
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
'@ @@ WARNING: ABORT IN matrxv'                                ,/,&
'@    ========'                                                ,/,&
'@     matrxv CALLED                WITH ISYM   = ',I10        ,/,&
'@'                                                            ,/,&
'@  The calculation will not be run.'                          ,/,&
'@'                                                            ,/,&
'@  Contact support.'                                          ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

!----
! End
!----

return

end subroutine
