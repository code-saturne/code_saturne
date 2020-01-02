!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

!> \file clipke.f90
!> \brief clipping of the turbulent kinetic energy and the turbulent
!> dissipation.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     ncelet        number of extended (real + ghost) cells
!> \param[in]     ncel          number of cells
!> \param[in]     iclip         indicator = 0 if viscl0 is used
!>                              otherwise viscl is used.
!______________________________________________________________________________

subroutine clipke &
 ( ncelet, ncel, iclip )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use cstphy
use cstnum
use entsor
use optcal
use parall
use cs_c_bindings
use field

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel
integer          iclip

! Local variables

integer          iclpke, iel, iclpk2, iclpe2
integer          ivar, ii, iwarnk
integer          iclpmn(2), iclpmx(1)
integer          kclipp, clip_e_id, clip_k_id

double precision xepmin,xepm,xe,xkmin,xkm,xk,var,epz2
double precision vmin(2), vmax(2)

double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_var
double precision, dimension(:), pointer :: viscl
double precision, dimension(:), pointer :: cpro_k_clipped
double precision, dimension(:), pointer :: cpro_e_clipped

type(var_cal_opt) :: vcopt

!===============================================================================

call field_get_val_s(icrom, crom)
call field_get_val_s(ivarfl(ik), cvar_k)
call field_get_val_s(ivarfl(iep), cvar_ep)

call field_get_val_s(iviscl, viscl)

call field_get_key_struct_var_cal_opt(ivarfl(ik), vcopt)
iwarnk = vcopt%iwarni

! Initialization to avoid compiler warnings

ivar = 0
iclpmx(1) = 0
! Une petite valeur pour eviter des valeurs exactement nulles.

epz2 = epzero**2

call field_get_key_id("clipping_id", kclipp)

! Postprocess clippings?
call field_get_key_int(ivarfl(ik), kclipp, clip_k_id)
if (clip_k_id.ge.0) then
  call field_get_val_s(clip_k_id, cpro_k_clipped)
endif

call field_get_key_int(ivarfl(iep), kclipp, clip_e_id)
if (clip_e_id.ge.0) then
  call field_get_val_s(clip_e_id, cpro_e_clipped)
endif

!===============================================================================
! ---> Stockage Min et Max pour log
!===============================================================================

do ii = 1, 2

  iclpmn(ii) = 0

  if (ii.eq.1) then
    cvar_var => cvar_k
  elseif(ii.eq.2) then
    cvar_var => cvar_ep
  endif

  vmin(ii) =  grand
  vmax(ii) = -grand
  do iel = 1, ncel
    var = cvar_var(iel)
    vmin(ii) = min(vmin(ii),var)
    vmax(ii) = max(vmax(ii),var)
  enddo

enddo

do iel = 1, ncel
  if (clip_k_id.ge.0) &
    cpro_k_clipped(iel) = 0.d0
  if (clip_e_id.ge.0) &
    cpro_e_clipped(iel) = 0.d0
enddo

!===============================================================================
! ---> Detection des valeurs hors norme "physiques"
!       uniquement pour avertissement
!       ou dans le cas ICLKEP = 1
!===============================================================================

if (iwarnk.ge.2.or.iclkep.eq.1) then

  if (iclip.eq.1) then

    xkm = 1296.d0*sqrt(cmu)/almax**2
    xepm = 46656.d0*cmu/almax**4
    iclpke = 0
    do iel= 1, ncel
      xk = cvar_k(iel)
      xe = cvar_ep(iel)
      xkmin = xkm * (viscl(iel) / crom(iel))**2
      xepmin = xepm * (viscl(iel) / crom(iel))**3
      if (xk.le.xkmin.or.xe.le.xepmin) then
        if(iclkep.eq.1) then
          if (clip_k_id.ge.0) &
            cpro_k_clipped(iel) = xkmin - xk
          cvar_k(iel)  = xkmin
          if (clip_e_id.ge.0) &
            cpro_e_clipped(iel) = xepmin - xe
          cvar_ep(iel) = xepmin
        endif
        iclpke = iclpke + 1
      endif
    enddo

  elseif (iclip.eq.0) then

    xkmin = 1296.d0 * sqrt(cmu)/almax**2 * (viscl0/ro0)**2
    xepmin = 46656.d0 * cmu/almax**4 * (viscl0/ro0)**3
    iclpke = 0
    do iel=1,ncel
      xk = cvar_k(iel)
      xe = cvar_ep(iel)
      if (xk.le.xkmin.or.xe.le.xepmin) then
        if (iclkep.eq.1) then
          cvar_k(iel)  = xkmin
          if (clip_k_id.ge.0) &
            cpro_k_clipped(iel) = xkmin - xk
          cvar_ep(iel) = xepmin
          if (clip_e_id.ge.0) &
            cpro_e_clipped(iel) = xepmin - xe
        endif
        iclpke = iclpke + 1
      endif
    enddo

  else

    write(nfecra,1000) iclip
    call csexit (1)

  endif

  ! ---  Stockage nb de clippings pour log

  if (iclkep.eq.1) then
    iclpmn(1) = iclpke
    iclpmn(2) = iclpke
  endif

  ! ---  Impression eventuelle

  if (iwarnk.ge.2) then
    if (irangp.ge.0) call parcpt (iclpke)
    write(nfecra,1010)iclpke
  endif

endif

!===============================================================================
! ---> Clipping "standard" ICLKEP = 0
!===============================================================================

if (iclkep.eq.0) then

  iclpk2 = 0
  iclpe2 = 0
  do iel = 1, ncel
    xk = cvar_k(iel)
    xe = cvar_ep(iel)
    if (abs(xk).le.epz2) then
      iclpk2 = iclpk2 + 1
      if (clip_k_id.ge.0) &
        cpro_k_clipped(iel) = epz2 - cvar_k(iel)
      cvar_k(iel) = max(cvar_k(iel),epz2)
    elseif(xk.le.0.d0) then
      iclpk2 = iclpk2 + 1
      if (clip_k_id.ge.0) &
        cpro_k_clipped(iel) = -xk
      cvar_k(iel) = -xk
    endif
    if (abs(xe).le.epz2) then
      iclpe2 = iclpe2 + 1
      if (clip_e_id.ge.0) &
        cpro_e_clipped(iel) = epz2 - cvar_ep(iel)
      cvar_ep(iel) = max(cvar_ep(iel), epz2)
    elseif(xe.le.0.d0) then
      iclpe2 = iclpe2 + 1
      if (clip_e_id.ge.0) &
        cpro_e_clipped(iel) = - xe
      cvar_ep(iel) = - xe
    endif
  enddo

  ! Stockage nb de clippings pour log

  iclpmn(1) = iclpk2
  iclpmn(2) = iclpe2

endif

do ii = 1, 2

  if (ii.eq.1) then
    ivar = ik
  elseif(ii.eq.2) then
    ivar = iep
  endif

  call log_iteration_clipping_field(ivarfl(ivar), iclpmn(ii), 0,  &
                                    vmin(ii:ii), vmax(ii:ii), iclpmn(ii),&
                                    iclpmx(1))

enddo

!===============================================================================
! ---> Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET DANS clipke                           ',/,&
'@    =========                                               ',/,&
'@     APPEL DE clipke              AVEC OPTION = ',I10        ,/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut pas etre execute.                       ',/,&
'@                                                            ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                           &
 i10,' VALEURS DU K-EPS AU DELA DES ECHELLES BASEES SUR ALMAX')

#else

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN clipke                                ',/,&
'@    ========                                                ',/,&
'@     CALL OF clipke               WITH OPTION = ',I10        ,/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  The calulation will not be run.                           ',/,&
'@                                                            ',/,&
'@  Contact the support.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1010 format(                                                           &
 i10,' K-EPS VALUES BEYOND THE SCALES BASED ON ALMAX')

#endif

return

end subroutine
