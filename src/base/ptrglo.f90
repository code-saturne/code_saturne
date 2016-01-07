!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

! Module for main fortran real array management
! and other generic array resizing procedures

module ptrglo

  !=============================================================================

contains

  !=============================================================================

  ! Resize a real array and synchronize halo

  subroutine resize_sca_real_array ( array )

    use mesh, only: ncel, ncelet

    implicit none

    ! Arguments

    double precision, pointer, dimension(:) :: array

    ! Local variables

    integer iel
    double precision, allocatable, dimension(:) :: buffer

    allocate(buffer(ncel))
    do iel = 1, ncel
      buffer(iel) = array(iel)
    enddo
    deallocate(array)

    allocate(array(ncelet))
    do iel = 1, ncel
      array(iel) = buffer(iel)
    enddo
    deallocate(buffer)

    call synsca (array)
    !==========

  end subroutine resize_sca_real_array

  !=============================================================================

  ! Resize a list of non interleaved real array and synchronize halo

  subroutine resize_n_sca_real_arrays ( size, array )

    use mesh, only: ncel, ncelet

    implicit none

    ! Arguments

    integer size
    double precision, pointer, dimension(:,:) :: array

    ! Local variables

    integer ii, iel
    double precision, allocatable, dimension(:,:) :: buffer

    allocate(buffer(ncel,size))
    do ii = 1, size
      do iel = 1, ncel
        buffer(iel,ii) = array(iel,ii)
      enddo
    enddo
    deallocate(array)

    allocate(array(ncelet,size))
    do ii = 1, size
      do iel = 1, ncel
        array(iel,ii) = buffer(iel,ii)
      enddo
    enddo
    deallocate(buffer)

    do ii = 1, size
      call synsca (array(1,ii))
      !==========
    enddo

  end subroutine resize_n_sca_real_arrays

  !=============================================================================

  ! Resize a vector non-interleaved array and synchronize halo

  subroutine resize_vec_real_array_ni ( array )

    use mesh, only: ncel, ncelet, ndim

    implicit none

    ! Arguments

    double precision, pointer, dimension(:,:) :: array

    ! Local variables

    integer iel, isou
    double precision, allocatable, dimension(:,:) :: buffer

    allocate(buffer(ncel,ndim))
    do isou = 1, ndim
      do iel = 1, ncel
        buffer(iel,isou) = array(iel,isou)
      enddo
    enddo
    deallocate(array)

    allocate(array(ncelet,ndim))
    do isou = 1, ndim
      do iel = 1, ncel
        array(iel,isou) = buffer(iel,isou)
      enddo
    enddo
    deallocate(buffer)

    call synvec (array(1,1), array(1,2), array(1,3))
    !==========

  end subroutine resize_vec_real_array_ni

  !=============================================================================

  ! Resize a vector interleaved array and synchronize halo

  subroutine resize_vec_real_array ( array )

    use mesh, only: ncel, ncelet, ndim

    implicit none

    ! Arguments

    double precision, pointer, dimension(:,:) :: array

    ! Local variables

    integer iel, isou
    double precision, allocatable, dimension(:,:) :: buffer

    allocate(buffer(ndim,ncel))
    do iel = 1, ncel
      do isou = 1, ndim
        buffer(isou,iel) = array(isou,iel)
      enddo
    enddo
    deallocate(array)

    allocate(array(ndim,ncelet))
    do iel = 1, ncel
      do isou = 1, ndim
        array(isou,iel) = buffer(isou,iel)
      enddo
    enddo
    deallocate(buffer)

    call synvin (array)
    !==========

  end subroutine resize_vec_real_array

  !=============================================================================

  ! Resize a tensor non-interleaved array and synchronize halo

  subroutine resize_tens_real_array_ni ( array )

    use mesh, only: ncel, ncelet, ndim

    implicit none

    ! Arguments

    double precision, pointer, dimension(:,:) :: array

    ! Local variables

    integer iel, isou
    double precision, allocatable, dimension(:,:) :: buffer

    allocate(buffer(ncel,ndim*ndim))
    do isou = 1, ndim*ndim
      do iel = 1, ncel
        buffer(iel,isou) = array(iel,isou)
      enddo
    enddo
    deallocate(array)

    allocate(array(ncelet,ndim*ndim))
    do isou = 1, ndim*ndim
      do iel = 1, ncel
        array(iel,isou) = buffer(iel,isou)
      enddo
    enddo
    deallocate(buffer)

    call synten &
    !==========
   ( array(1,1), array(1,2), array(1,3), &
     array(1,4), array(1,5), array(1,6), &
     array(1,7), array(1,8), array(1,9) )

  end subroutine resize_tens_real_array_ni

  !=============================================================================

  ! Resize a tensor interleaved array and synchronize halo

  subroutine resize_tens_real_array ( array )

    use mesh, only: ncel, ncelet, ndim

    implicit none

    ! Arguments

    double precision, pointer, dimension(:,:,:) :: array

    ! Local variables

    integer iel, ii, jj
    double precision, allocatable, dimension(:,:,:) :: buffer

    allocate(buffer(ndim,ndim,ncel))
    do iel = 1, ncel
      do jj = 1, ndim
        do ii = 1, ndim
          buffer(ii,jj,iel) = array(ii,jj,iel)
        enddo
      enddo
    enddo
    deallocate(array)

    allocate(array(ndim,ndim,ncelet))
    do iel = 1, ncel
      do jj = 1, ndim
        do ii = 1, ndim
          array(ii,jj,iel) = buffer(ii,jj,iel)
        enddo
      enddo
    enddo
    deallocate(buffer)

    call syntin (array)
    !==========

  end subroutine resize_tens_real_array

  !=============================================================================

  ! Resize main real array, update pointers and synchronize halo

  subroutine resize_main_real_array ( dt , propce )

    use dimens, only: nproce
    use field
    use mesh  , only: ncel, ncelet
    use numvar
    use parall, only: irangp
    use paramx
    use period, only: iperio

    implicit none

    ! Arguments

    double precision, pointer, dimension(:)   :: dt
    double precision, pointer, dimension(:,:) :: propce

    ! Local variables

    integer iel, iprop

    double precision, allocatable, dimension(:) :: dt0
    double precision, allocatable, dimension(:,:) :: proce0

    ! Buffering array

    allocate(dt0(ncelet))
    allocate(proce0(ncelet,nproce))

    do iel = 1, ncel
      dt0(iel) = dt(iel)
    enddo

    do iprop = 1, nproce
      do iel = 1, ncel
        proce0(iel,iprop) = propce(iel,iprop)
      enddo
    enddo

    ! Reallocate main arrays

    deallocate(dt, propce)

    allocate(dt(ncelet))
    allocate(propce(ncelet,nproce))

    ! Update new main real array  : "real" cells

    do iel = 1, ncel
      dt(iel) = dt0(iel)
    enddo

    do iprop = 1, nproce
      do iel = 1, ncel
        propce(iel,iprop) = proce0(iel,iprop)
      enddo
    enddo

    ! Update new main real array : "ghost" cells
    ! Note that the cs_field API cannot be used at this stage
    ! because the mapping fldtri is not done yet

    if (irangp.ge.0 .or. iperio.eq.1) then

      call synsca(dt)

      do iprop = 1, nproce
        if (iprop.eq.ipproc(irom)) call synsca (propce(1,iprop))
      enddo

    endif

    ! Free buffer

    deallocate(dt0)
    deallocate(proce0)

  return
  end subroutine resize_main_real_array

  !=============================================================================

end module ptrglo
