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

subroutine ecrlis &
!================

 ( ncelet , ncel   ,                                     &
   dt     , volume )

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE D'ECRITURE DES INFOS LISTING

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! dt   (ncelet)    ! tr ! <-- ! valeur du pas de temps                         !
! volume           ! tr ! <-- ! volume d'un des ncelet elements                !
! (ncelet)         !    !     !                                                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstnum
use cstphy
use albase
use parall
use ppppar
use ppthch
use ppincl
use field
use cs_f_interfaces
use cs_c_bindings

!===============================================================================

implicit none

integer          ncelet, ncel
double precision dt(ncelet), volume(ncelet)

! Local variables

integer          ic, icel, f_id, c_id, f_dim
integer          kval, nfld, f_type
character(len=200) :: chain, chainc, flabel,fname

double precision dervar(9)
double precision varres(9), varnrm(9)

double precision, allocatable, dimension(:) :: w1, w2

double precision, dimension(:), pointer :: field_s_v, field_s_vp
double precision, dimension(:,:), pointer :: field_v_v, field_v_vp

type(solving_info) sinfo

!===============================================================================

! Number of fields
call field_get_n_fields(nfld)

allocate(w1(ncelet), w2(ncelet))

!===============================================================================
! 2. ECRITURE DES CRITERES DE CONVERGENCE
!===============================================================================

write(nfecra,1000)
write(nfecra,1010)
write(nfecra,1011)
write(nfecra,1010)


do f_id = 0, nfld - 1

  call field_get_type(f_id, f_type)
  call field_get_key_int(f_id, keylog, kval)

  ! Is the field of type FIELD_VARIABLE?
  if (iand(f_type, FIELD_VARIABLE).eq.FIELD_VARIABLE) then
    chainc = 'c'
    chain = ' '
    ic = 4
    call field_get_dim (f_id, f_dim)
    call field_get_label(f_id, flabel)
    call field_get_name(f_id, fname)

    chainc(ic:ic+12) = flabel(1:12)
    ic=ic+12
    chain = ' '
    call field_get_key_struct_solving_info(f_id, sinfo)
    write(chain,3000) sinfo%rnsmbr
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+13
    chain = ' '
    write(chain,4000) sinfo%nbivar
    chainc(ic:ic+7) = chain(1:7)
    ic=ic+9
    chain = ' '
    write(chain,3000) sinfo%resvar
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+14
    chain = ' '

    ! Compute the time drift
    !-----------------------

    ! Scalar time drift (except pressure)
    if (f_dim.eq.1.and.(ippmod(icompf).ge.0.or.trim(fname).ne.'pressure')) then
      call field_get_val_s(f_id, field_s_v)
      call field_get_val_prev_s(f_id, field_s_vp)
      dervar(1) = 0.d0
      do icel = 1, ncel
        dervar(1) = dervar(1)                                            &
                  + (field_s_v(icel) - field_s_vp(icel))**2              &
                  *  volume(icel)/dt(icel)
      enddo
      if (irangp.ge.0) call parsom (dervar(1))
      dervar(1) = dervar(1) / voltot

    ! Pressure ime drift (computed in resopv.f90)
    else if (f_dim.eq.1) then
      dervar(1) = sinfo%dervar

    ! Vector or tensor time drift (total time drift)
    else
      call field_get_val_v(f_id, field_v_v)
      call field_get_val_prev_v(f_id, field_v_vp)

      ! Compute the derive
      dervar(1) = 0.d0

      do icel = 1, ncel

        ! Loop over the components
        do c_id = 1, f_dim
          dervar(1) = dervar(1)                                         &
                    + (field_v_v(c_id,icel) - field_v_vp(c_id,icel))**2 &
                    *  volume(icel)/dt(icel)
        enddo
      enddo

      if (irangp.ge.0) call parsom (dervar(1))
      dervar(1) = dervar(1) / voltot

    endif

    write(chain,3000) dervar(1)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+12

    ! L2 time normalized residual
    !----------------------------
    if (f_dim.eq.1) then
      call field_get_val_s(f_id, field_s_v)
      call field_get_val_prev_s(f_id, field_s_vp)
      do icel = 1, ncel
        w1(icel) = volume(icel)*(field_s_v(icel)-field_s_vp(icel))/dt(icel)
        w2(icel) = volume(icel)*field_s_v(icel)
      enddo

      ! L2 error, square root

      varres(1) = sqrt(cs_gres(ncel,volume,w1,w1))
      varnrm(1) = sqrt(cs_gres(ncel,volume,w2,w2))

      if (varnrm(1).gt.0.d0) varres(1) = varres(1)/varnrm(1)
      sinfo%l2residual = varres(1)

    ! Vector or tensor time drift (total L2 time residual)
    else
      call field_get_val_v(f_id, field_v_v)
      call field_get_val_prev_v(f_id, field_v_vp)

      ! Loop over the components
      do c_id = 1, f_dim

        do icel = 1, ncel
          w1(icel) = volume(icel) &
                   * (field_v_v(c_id, icel)-field_v_vp(c_id, icel))/dt(icel)
          w2(icel) = volume(icel)*field_v_v(c_id, icel)
        enddo

        ! L2 error, NO square root

        varres(c_id) = cs_gres(ncel,volume,w1,w1)
        varnrm(c_id) = cs_gres(ncel,volume,w2,w2)

        if (c_id.gt.1) then
          varres(1) = varres(1) + varres(c_id)
          varnrm(1) = varnrm(1) + varnrm(c_id)
        endif

      enddo

      if (varnrm(1).gt.0.d0) varres(1) = varres(1)/varnrm(1)
      sinfo%l2residual = sqrt(varres(1))

    endif

    write(chain,3000) sinfo%l2residual
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+12

    ! Finalize the log of the line
    if (kval.gt.0) write(nfecra,'(a)') chainc(1:ic)

    ! Vector or tensor time drift (by component)
    if (f_dim.gt.1) then
      call field_get_val_v(f_id, field_v_v)
      call field_get_val_prev_v(f_id, field_v_vp)

      ! Loop over the components
      do c_id = 1, f_dim

        ! Compute the time drift
        dervar(c_id) = 0.d0
        do icel = 1, ncel
          dervar(c_id) = dervar(c_id)                                      &
                       + (field_v_v(c_id,icel) - field_v_vp(c_id,icel))**2 &
                       *  volume(icel)/dt(icel)
        enddo
        if (irangp.ge.0) call parsom (dervar(c_id))
        dervar(c_id) = dervar(c_id) / voltot

        chainc = 'c'
        chain = ' '
        ic = 4

        ! Vectors
        if (f_dim.eq.3) then
          if (c_id.eq.1) then
            chain = trim(flabel) // '[X]'
          else if (c_id.eq.2) then
            chain = trim(flabel) // '[Y]'
          else if (c_id.eq.3) then
            chain = trim(flabel) // '[Z]'
          endif
        endif

        ! Symmetric tensors
        if (f_dim.eq.6) then
          if (c_id.eq.1) then
            chain = trim(flabel) // '[XX]'
          else if (c_id.eq.2) then
            chain = trim(flabel) // '[YY]'
          else if (c_id.eq.3) then
            chain = trim(flabel) // '[ZZ]'
          else if (c_id.eq.4) then
            chain = trim(flabel) // '[XY]'
          else if (c_id.eq.5) then
            chain = trim(flabel) // '[YZ]'
          else if (c_id.eq.6) then
            chain = trim(flabel) // '[XZ]'
          endif
        endif

        ! Tensors
        if (f_dim.eq.9) then
          if (c_id.eq.1) then
            chain = trim(flabel) // '[XX]'
          else if (c_id.eq.2) then
            chain = trim(flabel) // '[XY]'
          else if (c_id.eq.3) then
            chain = trim(flabel) // '[XZ]'
          else if (c_id.eq.4) then
            chain = trim(flabel) // '[YX]'
          else if (c_id.eq.5) then
            chain = trim(flabel) // '[YY]'
          else if (c_id.eq.6) then
            chain = trim(flabel) // '[YZ]'
          else if (c_id.eq.7) then
            chain = trim(flabel) // '[ZX]'
          else if (c_id.eq.8) then
            chain = trim(flabel) // '[ZY]'
          else if (c_id.eq.9) then
            chain = trim(flabel) // '[ZZ]'
          endif
        endif

        chainc(ic:ic+12) = chain(1:12)
        ic=ic+12
        chainc(ic:ic+12) = ' '
        ic=ic+13
        chainc(ic:ic+7) = ' '
        ic=ic+9
        chainc(ic:ic+12) = ' '
        ic=ic+14
        chain = ' '
        write(chain,3000) dervar(c_id)
        chainc(ic:ic+12) = chain(1:12)
        ic=ic+12

        ! Print the time drift of the component
        if (kval.gt.0) write(nfecra,'(a)') chainc(1:ic)

      enddo
    endif

    ! Store the time drift and the l2residual
    call field_set_key_struct_solving_info(f_id, sinfo)
  endif

enddo

write(nfecra,1010)
write(nfecra,*) ' '
write(nfecra,*) ' '

deallocate(w1, w2)

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format (/,3X,'** INFORMATIONS SUR LA CONVERGENCE',/,        &
          3X,'   -------------------------------')
 1011 format ('   Variable    Norm 2nd mb.',                      &
        '  Nbiter  Residu norme        Derive Residu en temps')
 1010 format ('---------------------------',                      &
        '----------------------------------------------------')

 3000 format (e12.5)
 4000 format (i7)

#else

 1000 format (/,3X,'** INFORMATION ON CONVERGENCE',/,             &
          3X,'   --------------------------')
 1011 format ('   Variable    Rhs norm    ',                      &
        '  N_iter  Norm. residual      Drift  Time residual')
 1010 format ('---------------------------',                      &
        '--------------------------------------------------')

 3000 format (e12.5)
 4000 format (i7)

#endif
return
end subroutine
