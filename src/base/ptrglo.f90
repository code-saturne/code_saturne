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
      call synsca (array(:,ii))
      !==========
    enddo

  end subroutine resize_n_sca_real_arrays

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

  ! Resize a symmetric tensor interleaved array and synchronize halo

  subroutine resize_sym_tens_real_array ( array )

    use mesh, only: ncel, ncelet

    implicit none

    ! Arguments

    double precision, pointer, dimension(:,:) :: array

    ! Local variables

    integer iel, isou
    double precision, allocatable, dimension(:,:) :: buffer

    allocate(buffer(6,ncel))
    do iel = 1, ncel
      do isou = 1, 6
        buffer(isou,iel) = array(isou,iel)
      enddo
    enddo
    deallocate(array)

    allocate(array(6,ncelet))
    do iel = 1, ncel
      do isou = 1, 6
        array(isou,iel) = buffer(isou,iel)
      enddo
    enddo
    deallocate(buffer)

    call syntis (array)
    !==========

  end subroutine resize_sym_tens_real_array

  !=============================================================================

end module ptrglo
