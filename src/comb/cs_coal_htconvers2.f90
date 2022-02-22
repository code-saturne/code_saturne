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

!===============================================================================
! Function:
! --------
!>  \file  cs_coal_htconvers2.f90
!>
!> \brief  - Calculating temperature of particles
!>           Function with enthalpy and concentrations
!>           if imode = 1
!>         - Calculating particles enthalpy
!>           Function with temperature and concentrations
!>           if imode = -1
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     mode          -1: t -> h; 1: h -> t
!> \param[in]     icla          class number
!> \param[in,out] enthal        mass enthalpy (in \f$ j . kg^{-1}) \f$
!> \param[in]     xsolid        mass fraction of components
!> \param[in,out] temper        temperature (in kelvin)
!> \param[in]     t1            coal inlet temperature
!______________________________________________________________________________!

subroutine cs_coal_htconvers2 &
 ( mode , icla , enthal , xsolid , temper , t1)

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use entsor
use cstnum
use cstphy
use ppppar
use ppthch
use coincl
use cpincl
use ppincl

!===============================================================================

implicit none

! Arguments

integer          mode , icla , icha

double precision xsolid(nsolim)
double precision temper , enthal , t1

! Local variables

integer          it , isol , ihflt2

double precision eh1 , eh0 , x2

!===============================================================================
! Conversion mode
ihflt2 = 1

if ( ihflt2.eq.0 ) then

!===============================================================================
! 1. H2 linear function
!===============================================================================

  icha = ichcor(icla)

  if ( mode.eq.-1 ) then

  ! --> Temperature law -> enthalpy (mode = -1)

    enthal = h02ch(icha) + cp2ch(icha)*(temper-trefth)

  elseif ( mode.eq.1 ) then

  ! --> Enthalpy law -> temperature (mode = 1)

    temper =  (enthal-h02ch(icha))/cp2ch(icha) + trefth

  else

    write(nfecra,1000) mode
    call csexit (1)
    !==========

  endif


elseif( ihflt2.ne.0 ) then

!===============================================================================
! 2. H2 tabulated
!===============================================================================

  if ( mode.eq.-1 ) then

    ! --> Temperature law -> enthalpy (mode = -1)

    it = npoc
    if ( temper.ge.thc(it) ) then
      enthal = zero
      do isol = 1, nsolid
        enthal = enthal + xsolid(isol)*ehsoli(isol,it)
      enddo
      go to 11
    endif

    it = 1
    if ( temper.le.thc(it) ) then
      enthal = zero
      do isol = 1, nsolid
        enthal = enthal + xsolid(isol)*ehsoli(isol,it)
      enddo
      go to 11
    endif
    it = 1
 10       continue

    it = it + 1
    if ( temper.le.thc(it) ) then
      eh0 = zero
      eh1 = zero
      do isol = 1, nsolid
        eh0 = eh0 + xsolid(isol)*ehsoli(isol,it-1)
        eh1 = eh1 + xsolid(isol)*ehsoli(isol,it  )
      enddo
      enthal = eh0                                                &
             + (eh1-eh0)*(temper-thc(it-1))                       &
                        /(thc(it)-thc(it-1))
      goto 11
    endif
    goto 10
 11       continue

  elseif ( mode.eq.1 ) then

    ! --> Enthalpy law -> temperature (mode = 1)

    x2 = 0.d0
    do isol = 1, nsolid
      x2 = x2 + xsolid(isol)
    enddo

    if ( x2 .gt. epsicp ) then
      it  = npoc-1
      eh1 = zero
      do isol = 1, nsolid
        eh1 = eh1 + xsolid(isol)*ehsoli(isol,it+1)
      enddo
      if ( enthal.ge.eh1 ) temper = thc(it+1)

      it  = 1
      eh0 = zero
      do isol = 1, nsolid
        eh0 = eh0 + xsolid(isol)*ehsoli(isol,it  )
      enddo
      if ( enthal.le.eh0 ) temper = thc(it)

      do it = 1, npoc-1
        eh0 = zero
        eh1 = zero
        do isol = 1, nsolid
          eh0 = eh0 + xsolid(isol)*ehsoli(isol,it  )
          eh1 = eh1 + xsolid(isol)*ehsoli(isol,it+1)
        enddo
        if ( enthal.ge.eh0 .and. enthal.le.eh1 )                  &
          temper = thc(it)                                        &
            + (enthal-eh0)*(thc(it+1)-thc(it))/(eh1-eh0)

      enddo

    else
      temper = t1
    endif

  else

    write(nfecra,1000) mode
    call csexit (1)

  endif

endif

!--------
! Formats
!--------

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ Error: in cs_coal_htconver2                             ',/,&
'@    ======                                                  ',/,&
'@    Incorrect value of the argument mode                    ',/,&
'@    it must be an integer equal to 1 or -1                  ',/,&
'@    its value here is ', i10                                ,/,&
'@                                                            ',/,&
'@  The calculation can not run.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return
end subroutine
