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

!> \file inimav.f90
!>
!> \brief This function adds \f$ \rho \vect{u} \cdot \vect{S}_\ij\f$ to the mass
!> flux \f$ \dot{m}_\ij \f$.
!>
!> For the reconstruction, \f$ \gradt \left(\rho \vect{u} \right) \f$ is
!> computed with the following approximated boundary conditions:
!>  - \f$ \vect{A}_{\rho u} = \rho_\fib \vect{A}_u \f$
!>  - \f$ \tens{B}_{\rho u} = \tens{B}_u \f$
!>
!> For the mass flux at the boundary we have:
!> \f[
!> \dot{m}_\ib = \left[ \rho_\fib \vect{A}_u  + \rho_\fib \tens{B}_u \vect{u}
!> + \tens{B}_u \left(\gradt \vect{u} \cdot \vect{\centi \centip}\right)\right]
!> \cdot \vect{S}_\ij
!> \f]
!> The last equation uses some approximations detailed in the theory guide.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     ivar          index of the current variable
!> \param[in]     itypfl        indicator (take rho into account or not)
!>                               - 1 compute \f$ \rho\vect{u}\cdot\vect{S} \f$
!>                               - 0 compute \f$ \vect{u}\cdot\vect{S} \f$
!> \param[in]     iflmb0        the mass flux is set to 0 on walls and
!>                               symmetries if = 1
!> \param[in]     init          the mass flux is initialize to 0 if > 0
!> \param[in]     inc           indicator
!>                               - 0 solve an increment
!>                               - 1 otherwise
!> \param[in]     imrgra        indicator
!>                               - 0 iterative gradient
!>                               - 1 least square gradient
!> \param[in]     nswrgu        number of sweeps for the reconstruction
!>                               of the gradients
!> \param[in]     imligu        clipping gradient method
!>                               - < 0 no clipping
!>                               - = 0 thank to neighbooring gradients
!>                               - = 1 thank to the mean gradient
!> \param[in]     iwarnu        verbosity
!> \param[in]     nfecra        unit of the standard output file
!> \param[in]     epsrgu        relative precision for the gradient
!>                               reconstruction
!> \param[in]     climgu        clipping coeffecient for the computation of
!>                               the gradient
!> \param[in]     rom           cell density
!> \param[in]     romb          border face density
!> \param[in]     vel           vector variable
!> \param[in]     coefav        boundary condition array for the variable
!>                               (Explicit part - vector array )
!> \param[in]     coefbv        boundary condition array for the variable
!>                               (Impplicit part - 3x3 tensor array)
!> \param[in,out] flumas        interior mass flux \f$ \dot{m}_\fij \f$
!> \param[in,out] flumab        border mass flux \f$ \dot{m}_\fib \f$
!_______________________________________________________________________________

subroutine inimav &
 ( ivar   , itypfl ,                                              &
   iflmb0 , init   , inc    , imrgra , nswrgu , imligu ,          &
   iwarnu , nfecra ,                                              &
   epsrgu , climgu ,                                              &
   rom    , romb   ,                                              &
   vel    ,                                                       &
   coefav , coefbv ,                                              &
   flumas , flumab )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only: ndimfb
use pointe
use optcal, only: iporos
use parall
use period
use mesh

!===============================================================================

implicit none

! Arguments

integer          ivar   , itypfl
integer          iflmb0 , init   , inc    , imrgra
integer          nswrgu , imligu
integer          iwarnu , nfecra
double precision epsrgu , climgu


double precision rom(ncelet), romb(nfabor)
double precision vel(3,ncelet)
double precision coefav(3,ndimfb)
double precision coefbv(3,3,nfabor)
double precision flumas(nfac), flumab(nfabor)

! Local variables

integer          ifac, ii, jj, iel, ig, it
integer          isou, jsou
double precision pfac, pip
double precision dofx,dofy,dofz,pnd
double precision diipbx, diipby, diipbz
logical          ilved

double precision, dimension(:,:), allocatable :: qdm, coefaq
double precision, dimension(:,:,:), allocatable :: grdqdm

allocate(qdm(3,ncelet))
allocate(coefaq(3,ndimfb))

!===============================================================================

!===============================================================================
! 1.  Initialization
!===============================================================================

! ---> Momentum computation

if (init.eq.1) then
  !$omp parallel do
  do ifac = 1, nfac
    flumas(ifac) = 0.d0
  enddo
  !$omp parallel do if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    flumab(ifac) = 0.d0
  enddo

elseif (init.ne.0) then
  write(nfecra,1000) init
  call csexit (1)
endif

! Standard mass flux
if (itypfl.eq.1) then

  ! Without porosity
  if (iporos.eq.0) then
    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
        qdm(isou,iel) = rom(iel)*vel(isou,iel)
      enddo
    enddo
  ! With porosity
  else if (iporos.eq.1) then
    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
        qdm(isou,iel) = rom(iel)*vel(isou,iel)*porosi(iel)
      enddo
    enddo
  ! With anisotropic porosity
  else if (iporos.eq.2) then
    !$omp parallel do private(isou)
    do iel = 1, ncel
      qdm(1, iel) = ( porosf(1, iel)*vel(1, iel)   &
                    + porosf(4, iel)*vel(2, iel)   &
                    + porosf(6, iel)*vel(3, iel) ) &
                  * rom(iel)
      qdm(2, iel) = ( porosf(4, iel)*vel(1, iel)   &
                    + porosf(2, iel)*vel(2, iel)   &
                    + porosf(5, iel)*vel(3, iel) ) &
                  * rom(iel)
      qdm(3, iel) = ( porosf(6, iel)*vel(1, iel)   &
                    + porosf(5, iel)*vel(2, iel)   &
                    + porosf(3, iel)*vel(3, iel) ) &
                  * rom(iel)
    enddo
  endif

! Velocity flux
else

  ! Without porosity
  if (iporos.eq.0) then
    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
        qdm(isou,iel) = vel(isou,iel)
      enddo
    enddo
  ! With porosity
  else if (iporos.eq.1) then
    !$omp parallel do private(isou)
    do iel = 1, ncel
      do isou = 1, 3
        qdm(isou,iel) = vel(isou,iel)*porosi(iel)
      enddo
    enddo
  ! With anisotropic porosity
  else if (iporos.eq.2) then
    !$omp parallel do private(isou)
    do iel = 1, ncel
      qdm(1, iel) = porosf(1, iel)*vel(1, iel)  &
                  + porosf(4, iel)*vel(2, iel)  &
                  + porosf(6, iel)*vel(3, iel)
      qdm(2, iel) = porosf(4, iel)*vel(1, iel)  &
                  + porosf(2, iel)*vel(2, iel)  &
                  + porosf(5, iel)*vel(3, iel)
      qdm(3, iel) = porosf(6, iel)*vel(1, iel)  &
                  + porosf(5, iel)*vel(2, iel)  &
                  + porosf(3, iel)*vel(3, iel)
    enddo
  endif
endif

! ---> Periodicity and parallelism treatment

if (irangp.ge.0.or.iperio.eq.1) then
  call synvin(qdm)
endif

! Standard mass flux
if (itypfl.eq.1) then

  ! Without porosity
  if (iporos.eq.0) then
    !$omp parallel do private(isou) if(nfabor > thr_n_min)
    do ifac =1, nfabor
      do isou = 1, 3
        coefaq(isou,ifac) = romb(ifac)*coefav(isou,ifac)
      enddo
    enddo
  ! With porosity
  else if (iporos.eq.1) then
    !$omp parallel do private(iel, isou) if(nfabor > thr_n_min)
    do ifac =1, nfabor
      iel = ifabor(ifac)
      do isou = 1, 3
        coefaq(isou,ifac) = romb(ifac)*coefav(isou,ifac)*porosi(iel)
      enddo
    enddo
  ! With anisotropic porosity
  else if (iporos.eq.2) then
    !$omp parallel do private(iel) if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      coefaq(1, ifac) = ( porosf(1, iel)*coefav(1, ifac)   &
                        + porosf(4, iel)*coefav(2, ifac)   &
                        + porosf(6, iel)*coefav(3, ifac) ) &
                      * romb(ifac)
      coefaq(2, ifac) = ( porosf(4, iel)*coefav(1, ifac)   &
                        + porosf(2, iel)*coefav(2, ifac)   &
                        + porosf(5, iel)*coefav(3, ifac) ) &
                      * romb(ifac)
      coefaq(3, ifac) = ( porosf(6, iel)*coefav(1, ifac)   &
                        + porosf(5, iel)*coefav(2, ifac)   &
                        + porosf(3, iel)*coefav(3, ifac) ) &
                      * romb(ifac)
    enddo
  endif

! Velocity flux
else

  ! Without porosity
  if (iporos.eq.0) then
    !$omp parallel do private(isou) if(nfabor > thr_n_min)
    do ifac =1, nfabor
      do isou = 1, 3
        coefaq(isou,ifac) = coefav(isou,ifac)
      enddo
    enddo
  ! With porosity
  else if (iporos.eq.1) then
    !$omp parallel do private(iel, isou) if(nfabor > thr_n_min)
    do ifac =1, nfabor
      iel = ifabor(ifac)
      do isou = 1, 3
        coefaq(isou,ifac) = coefav(isou,ifac)*porosi(iel)
      enddo
    enddo
  ! With anisotropic porosity
  else if (iporos.eq.2) then
    !$omp parallel do private(iel) if(nfabor > thr_n_min)
    do ifac = 1, nfabor
      iel = ifabor(ifac)
      coefaq(1, ifac) = porosf(1, iel)*coefav(1, ifac)  &
                      + porosf(4, iel)*coefav(2, ifac)  &
                      + porosf(6, iel)*coefav(3, ifac)
      coefaq(2, ifac) = porosf(4, iel)*coefav(1, ifac)  &
                      + porosf(2, iel)*coefav(2, ifac)  &
                      + porosf(5, iel)*coefav(3, ifac)
      coefaq(3, ifac) = porosf(6, iel)*coefav(1, ifac)  &
                      + porosf(5, iel)*coefav(2, ifac)  &
                      + porosf(3, iel)*coefav(3, ifac)
    enddo
  endif

endif

!===============================================================================
! 2. Compute mass flux without recontructions
!===============================================================================

if (nswrgu.le.1) then

  ! Interior faces

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, pnd, isou)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)
        pnd = pond(ifac)
        ! u, v, w Components
        do isou = 1, 3
          flumas(ifac) = flumas(ifac) +                                       &
               (pnd*qdm(isou,ii)+(1.d0-pnd)*qdm(isou,jj)) * surfac(isou,ifac)
        enddo

      enddo
    enddo
  enddo

  ! Boundary faces

  ! Standard mass flux
  if (itypfl.eq.1) then
    do ig = 1, ngrpb
      !$omp parallel do private(ifac, ii, isou, jsou, pfac) &
      !$omp          if(nfabor > thr_n_min)
      do it = 1, nthrdb
        do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

          ii = ifabor(ifac)
          ! u, v, w Components
          do isou = 1, 3
            pfac = inc*coefaq(isou,ifac)

            ! coefbv is a matrix
            do jsou = 1, 3
              pfac = pfac + romb(ifac)*coefbv(isou,jsou,ifac)*vel(jsou,ii)
            enddo

            flumab(ifac) = flumab(ifac) + pfac*surfbo(isou,ifac)
          enddo

        enddo
      enddo
    enddo

  ! Velocity flux
  else
    do ig = 1, ngrpb
      !$omp parallel do private(ifac, ii, isou, jsou, pfac) &
      !$omp          if(nfabor > thr_n_min)
      do it = 1, nthrdb
        do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

          ii = ifabor(ifac)
          ! u, v, w Components
          do isou = 1, 3
            pfac = inc*coefaq(isou,ifac)

            ! coefbv is a matrix
            do jsou = 1, 3
              pfac = pfac + coefbv(isou,jsou,ifac)*vel(jsou,ii)
            enddo

            flumab(ifac) = flumab(ifac) + pfac*surfbo(isou,ifac)
          enddo

        enddo
      enddo
    enddo

  endif

endif

!===============================================================================
! 4.  CALCUL DU FLUX DE MASSE AVEC TECHNIQUE DE RECONSTRUCTION
!        SI LE MAILLAGE EST NON ORTHOGONAL
!===============================================================================

if (nswrgu.gt.1) then

  allocate(grdqdm(3,3,ncelet))


!     CALCUL DU GRADIENT SUIVANT de QDM
!     =================================
  ! gradient vectoriel la periodicite est deja traitee
  ilved = .true.

  call grdvec &
  !==========
( ivar   , imrgra , inc    , nswrgu , imligu ,                   &
  iwarnu , epsrgu , climgu ,                                     &
  ilved  ,                                                       &
  qdm    , coefaq , coefbv ,                                     &
  grdqdm )


  ! Mass flow through interior faces

  do ig = 1, ngrpi
    !$omp parallel do private(ifac, ii, jj, isou, pnd, dofx, dofy, dofz)
    do it = 1, nthrdi
      do ifac = iompli(1,ig,it), iompli(2,ig,it)

        ii = ifacel(1,ifac)
        jj = ifacel(2,ifac)

        pnd = pond(ifac)

        dofx = dofij(1,ifac)
        dofy = dofij(2,ifac)
        dofz = dofij(3,ifac)

        ! Terms along U, V, W
        do isou = 1, 3

          flumas(ifac) = flumas(ifac)                                &
               ! Non-reconstructed term
               + (pnd*qdm(isou,ii) + (1.d0-pnd)*qdm(isou,jj)         &

               !  --->     ->    -->      ->
               ! (Grad(rho U ) . OFij ) . Sij
               + 0.5d0*(grdqdm(1,isou,ii) +grdqdm(1,isou,jj))*dofx   &
               + 0.5d0*(grdqdm(2,isou,ii) +grdqdm(2,isou,jj))*dofy   &
               + 0.5d0*(grdqdm(3,isou,ii) +grdqdm(3,isou,jj))*dofz   &
               )*surfac(isou,ifac)
        enddo

      enddo
    enddo

  enddo

  ! Mass flow though boundary faces

  ! Standard mass flux
  if (itypfl.eq.1) then
    do ig = 1, ngrpb
      !$omp parallel do private(ifac, ii, isou, jsou, diipbx, diipby, diipbz, &
      !$omp                     pfac, pip) if(nfabor > thr_n_min)
      do it = 1, nthrdb
        do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

          ii = ifabor(ifac)
          diipbx = diipb(1,ifac)
          diipby = diipb(2,ifac)
          diipbz = diipb(3,ifac)

          ! Terms along U, V, W
          do isou = 1, 3

            pfac = inc*coefaq(isou,ifac)

            ! coefu is a matrix
            do jsou = 1, 3

              pip =  romb(ifac)*vel(jsou,ii)                &
                   + grdqdm(1,jsou,ii)*diipbx               &
                   + grdqdm(2,jsou,ii)*diipby               &
                   + grdqdm(3,jsou,ii)*diipbz

              pfac = pfac + coefbv(isou,jsou,ifac)*pip

            enddo

            flumab(ifac) = flumab(ifac) + pfac*surfbo(isou,ifac)

          enddo

        enddo
      enddo
    enddo

  ! Velocity flux
  else
    do ig = 1, ngrpb
      !$omp parallel do private(ifac, ii, isou, jsou, diipbx, diipby, diipbz, &
      !$omp                     pfac, pip) if(nfabor > thr_n_min)
      do it = 1, nthrdb
        do ifac = iomplb(1,ig,it), iomplb(2,ig,it)

          ii = ifabor(ifac)
          diipbx = diipb(1,ifac)
          diipby = diipb(2,ifac)
          diipbz = diipb(3,ifac)

          ! Terms along U, V, W
          do isou = 1, 3

            pfac = inc*coefaq(isou,ifac)

            ! coefu is a matrix
            do jsou = 1, 3

              pip = vel(jsou,ii)                           &
                  + grdqdm(1,jsou,ii)*diipbx               &
                  + grdqdm(2,jsou,ii)*diipby               &
                  + grdqdm(3,jsou,ii)*diipbz

              pfac = pfac + coefbv(isou,jsou,ifac)*pip

            enddo

            flumab(ifac) = flumab(ifac) + pfac*surfbo(isou,ifac)

          enddo

        enddo
      enddo
    enddo

  endif

! Deallocation
deallocate(grdqdm)
deallocate(qdm, coefaq)

endif

!===============================================================================
! 6.  POUR S'ASSURER DE LA NULLITE DU FLUX DE MASSE AUX LIMITES
!       SYMETRIES PAROIS COUPLEES
!===============================================================================

if (iflmb0.eq.1) then
  ! Force flumab to 0 for velocity
  !$omp parallel do if(nfabor > thr_n_min)
  do ifac = 1, nfabor
    if (isympa(ifac).eq.0) then
      flumab(ifac) = 0.d0
    endif
  enddo
endif

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format('INIMAV APPELE AVEC INIT =',I10)

#else

 1000 format('INIMAV CALLED WITH INIT =',I10)

#endif

!----
! End
!----

return

end subroutine
