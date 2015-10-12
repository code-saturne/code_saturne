!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

!> \file precdi.f90
!>
!> \brief Management of the injection of particles formed by precipitation.
!>
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     dt            time step (per cell)
!> \param[in]     iprev         time step indicator for fields
!> \param[out]    val           number of particles to inject (with weight)
!_______________________________________________________________________________

subroutine precdi          &
 (  dt  , iprev , val )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use numvar
use optcal
use cstnum
use parall
use lagran
use lagdim
use mesh
use field
use entsor

!===============================================================================

implicit none
! Arguments

double precision dt(ncelet)
integer          iprev
double precision val

! Local variables

integer          iel , ip , k , npt

double precision, dimension(:,:), pointer :: vela
double precision, dimension(:), pointer :: cscalt
double precision, dimension(:), pointer :: cscal, cscal2

double precision  dintrf(1)

integer                       nbprec_tot, vv
integer                       nbdiss_tot
integer                       value_part, nbprec2
double  precision             d3 , pis6

integer, dimension(:), allocatable :: cell
integer, dimension(:), allocatable :: nbdiss
double  precision, dimension(:), allocatable :: mp
double  precision, dimension(:), allocatable :: mp_diss_t

!===============================================================================

if (iprev.eq.0) then
  call field_get_val_v(ivarfl(iu), vela)
  if (itherm.eq.1.or.itherm.eq.2) then
    call field_get_val_s(ivarfl(isca(iscalt)), cscalt)
  endif
else if (iprev.eq.1) then
  call field_get_val_prev_v(ivarfl(iu), vela)
  if (itherm.eq.1.or.itherm.eq.2) then
    call field_get_val_prev_s(ivarfl(isca(iscalt)), cscalt)
  endif
endif

call field_get_val_s(ivarfl(isca(1)), cscal)
call field_get_val_prev_s(ivarfl(isca(1)), cscal2)

!===============================================================================
! 1. INITIALISATION
!===============================================================================

pis6 = pi / 6.d0

allocate(nbdiss(nbrclas))
allocate(mp(nbrclas))
allocate(mp_diss_t(ncelet))

!number of dissolved particles
nbdiss = 0
nbdiss_tot = 0

!mass of dissolved particles
mp_diss_t = 0.d0
!number of precipated particles
nbprec_tot = 0

nbprec2 = 0

!===============================================================================
! 2. GESTION DES PARTICULES
!===============================================================================
do iel = 1 , ncel
  nbprec2 = nbprec2 + nbprec(iel)
enddo

if (nbprec2 .ge. 1.d6) then
  write(nfecra,1000) nbprec2
  call csexit (1)
endif

allocate(cell(nbprec2))
cell = 0

do iel = 1 , ncel

  !Precipitation (Add particles)
  if (nbprec(iel) .gt. 0) then
    cell(nbprec_tot + 1 : nbprec_tot + nbprec(iel)) = iel
    nbprec_tot = nbprec_tot + nbprec(iel)
  endif

  do k = 1, nbrclas
    mp_diss_t(iel) = mp_diss_t(iel) + mp_diss(iel,k)
  enddo

  !Dissolution (Remove particles)
  if (mp_diss_t(iel) .gt. 0) then
    mp = 0.d0
    do npt  = 1 , nbpart
      do k = 1,  nbrclas
        if (      ipepa(jisor,npt) .eq. iel             &
            .and. eptp(jdp,npt) .eq. ruslag(k,1,idpt)   &
            .and. (mp(k) .lt. mp_diss(iel,k))) then
          ! Removing of particles due to dissolution
          ipepa(jisor,npt) = 0
          mp(k)  = mp(k) + pepa(jrpoi,npt) &
                         * (pi/6.d0 * (eptp(jdp,npt))**3 * rho_preci)
          nbdiss(k) = nbdiss(k) + 1
        endif
      enddo
    enddo
  endif

enddo

do k = 1, nbrclas
  nbdiss_tot = nbdiss_tot + nbdiss(k)
enddo

value_part = nbprec_tot
if (irangp.ge.0) then
  call parcpt(nbprec_tot)
  call parcpt(nbdiss_tot)
endif

npt = nbpart
nbpnew = nbpnew + nbprec_tot
vv = lagr_resize_particle_set(nbpart+nbpnew)

if (nbprec_tot .ge. 1) then

  do ip = npt + 1 , npt + value_part
    !to do: place particle at random location in the cell iel (not always at the cog)
    eptp(jxp,ip) = xyzcen(1,cell(ip - npt))
    eptp(jyp,ip) = xyzcen(2,cell(ip - npt))
    eptp(jzp,ip) = xyzcen(3,cell(ip - npt))
    ipepa(jisor,ip) = cell(ip - npt)
    eptp(juf,ip) =  vela(1, cell(ip - npt))
    eptp(jvf,ip) =  vela(2, cell(ip - npt))
    eptp(jwf,ip) =  vela(3, cell(ip - npt))
    eptp(jup,ip) =  vela(1, cell(ip - npt))
    eptp(jvp,ip) =  vela(2, cell(ip - npt))
    eptp(jwp,ip) =  vela(3, cell(ip - npt))
    eptp(jdp,ip) =  dprec
    d3 = eptp(jdp,ip) * eptp(jdp,ip) * eptp(jdp,ip)
    eptp(jmp,ip) = rho_preci * pis6 * d3
    pepa(jrpoi,ip) = 1.d0
    pepa(jrtsp,ip) = 0.d0

    if ( idepst .eq. 1 ) then

      call zufall(1,dintrf(1))

      pepa(jrinpf,ip) = 5.d0 + 15.d0 * dintrf(1)

      pepa(jryplu,ip) = 1000.d0
      ipepa(jimark,ip) = -1
      ipepa(jdfac,ip)  = 0
      ipepa(jdepo,ip) = 0

    endif

  enddo

endif

val = 0.d0
do ip = npt + 1 , npt + value_part
  val = val + pepa(jrpoi,ip)
enddo

nbpart = nbpart + value_part

deallocate(cell)
deallocate(nbdiss)
deallocate(mp)
deallocate(mp_diss_t)

1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE DE           ',/,&
'@    =========   PRECIPITATION/DISSOLUTION (prec_diss)       ',/,&
'@                                                            ',/,&
'@  Le nombre maximal de particules precipitees               ',/,&
'@    a  ete depasse', I10                                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

end subroutine precdi
