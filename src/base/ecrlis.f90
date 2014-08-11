!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

 ( nvar   , ncelet , ncel   ,                                     &
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
! nvar             ! e  ! <-- ! nombre de variables                            !
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

!===============================================================================

implicit none

integer          nvar, ncelet, ncel
double precision dt(ncelet), volume(ncelet)

! Local variables

integer          ic, icel, ivar, ipp, f_id, f_id_prv, c_id, f_dim
integer          ippf, kval
logical          interleaved
character(len=200) :: chain, chainc

double precision, dimension(:), pointer :: field_s_v, field_s_vp
double precision, dimension(:,:), pointer :: field_v_v, field_v_vp

integer, save :: keypp = -1

!===============================================================================

if (keypp.lt.0) then
  call field_get_key_id("post_id", keypp)
endif

f_id = -1
c_id = 1

do ivar = 1, nvar
  if (ivar.eq.ipr.and.ippmod(icompf).lt.0) cycle

  f_id_prv = f_id
  f_id = ivarfl(ivar)
  if (f_id.eq.f_id_prv) then
    c_id = c_id + 1
  else
    c_id = 1
  endif
  call field_get_key_int(f_id, keypp, ippf)
  if (ippf.le.1) cycle

  call field_get_dim(f_id, f_dim, interleaved)
  if (f_dim.gt.1) then
    call field_get_val_v(f_id, field_v_v)
    call field_get_val_prev_v(f_id, field_v_vp)
  else if (f_dim.eq.1) then
    call field_get_val_s(f_id, field_s_v)
    call field_get_val_prev_s(f_id, field_s_vp)
  endif

  if (f_dim.gt.1) then
    call field_get_key_int(f_id, keylog, kval)
    if (kval.gt.0) then
      call field_get_key_int(f_id, keypp, ipp)
      dervar(ipp) = 0.d0
      do icel = 1, ncel
        dervar(ipp) = dervar(ipp)                                            &
                + (field_v_v(c_id,icel) - field_v_vp(c_id,icel))**2          &
                *  volume(icel)/dt(icel)
      enddo
      if (irangp.ge.0) call parsom (dervar(ipp))
      dervar(ipp) = dervar(ipp) / voltot
    endif
  else if (f_dim.eq.1) then
    call field_get_key_int(f_id, keylog, kval)
    if (kval.gt.0) then
      call field_get_key_int(f_id, keypp, ipp)
      dervar(ipp) = 0.d0
      do icel = 1, ncel
        dervar(ipp) = dervar(ipp)                                            &
                + (field_s_v(icel) - field_s_vp(icel))**2                    &
                *  volume(icel)/dt(icel)
      enddo
      if (irangp.ge.0) call parsom (dervar(ipp))
      dervar(ipp) = dervar(ipp) / voltot
    endif
  endif
enddo

if (ippmod(icompf).lt.0) then
  ipp = ipprtp(ipr)
  if (dervar(ipp).lt.epzero) then
    dervar(ipp) = -1.d0
  endif
  dervar(ipp) = rnsmbr(ipp) / dervar(ipp)
endif

!===============================================================================
! 2. ECRITURE DES CRITERES DE CONVERGENCE
!===============================================================================

write(nfecra,1000)
write(nfecra,1010)
write(nfecra,1011)
write(nfecra,1010)

f_id = -1
c_id = 1

do ivar = 1, nvar

  ipp = ipprtp(ivar)

  f_id_prv = f_id
  f_id = ivarfl(ivar)
  if (f_id.eq.f_id_prv) then
    c_id = c_id + 1
  else
    c_id = 1
  endif

  call field_get_key_int(ivarfl(ivar), keylog, kval)
  if (kval.gt.0) then
    chainc = 'c'
    chain = ' '
    ic = 4
    call field_get_dim (f_id, f_dim, interleaved)
    call field_get_label(f_id, chain)
    if (f_dim.gt.1) then
      if (c_id.eq.1) then
        chain = trim(chain) // 'X'
      else if (c_id.eq.2) then
        chain = trim(chain) // 'Y'
      else if (c_id.eq.3) then
        chain = trim(chain) // 'Z'
      endif
    endif
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+12
    chain = ' '
    write(chain,3000) rnsmbr(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+13
    chain = ' '
    write(chain,4000) nbivar(ipp)
    chainc(ic:ic+7) = chain(1:7)
    ic=ic+9
    chain = ' '
    write(chain,3000) resvar(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+14
    chain = ' '
    write(chain,3000) dervar(ipp)
    chainc(ic:ic+12) = chain(1:12)
    ic=ic+12
    write(nfecra,'(a)') chainc(1:ic)
  endif

enddo

write(nfecra,1010)
write(nfecra,*) ' '
write(nfecra,*) ' '

!--------
! Formats
!--------

#if defined(_CS_LANG_FR)

 1000 format (/,3X,'** INFORMATIONS SUR LA CONVERGENCE',/,        &
          3X,'   -------------------------------')
 1011 format ('   Variable    Norm 2nd mb.',                      &
        '  Nbiter  Residu norme        derive')
 1010 format ('---------------------------',                      &
        '------------------------------------')

 3000 format (e12.5)
 4000 format (i7)

#else

 1000 format (/,3X,'** INFORMATION ON CONVERGENCE',/,             &
          3X,'   --------------------------')
 1011 format ('   Variable    Rhs norm    ',                      &
        '  N_iter  Norm. residual      derive')
 1010 format ('---------------------------',                      &
        '------------------------------------')

 3000 format (e12.5)
 4000 format (i7)

#endif
return
end subroutine
