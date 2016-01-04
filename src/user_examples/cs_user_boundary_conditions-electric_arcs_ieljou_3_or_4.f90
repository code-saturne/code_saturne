!-------------------------------------------------------------------------------

!VERS

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

!===============================================================================
! Function:
! ---------
!> \file  cs_user_boundary_conditions-electrics_arcs_ieljou.f90
!> \brief Example of cs_user_boundary_conditions subroutine.f90 for electrics arcs / Joule
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[out]    icodcl        boundary condition code:
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!> \param[in]     itrifb        indirection for boundary faces ordering
!> \param[in,out] itypfb        boundary face types
!> \param[out]    izfppp        boundary face zone number
!> \param[in]     dt            time step (per cell)
!> \param[in,out] rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________

subroutine cs_user_boundary_conditions &
 ( nvar   , nscal  ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   dt     ,                                                       &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ihmpre
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use atincl
use ctincl
use elincl
use cs_fuel_incl
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          icodcl(nfabor,nvarcl)
integer          itrifb(nfabor), itypfb(nfabor)
integer          izfppp(nfabor)

double precision dt(ncelet)
double precision rcodcl(nfabor,nvarcl,3)

! Local variables

!< [loc_var_dec]
integer          ifac  , ii     , iel
integer          i     , ntf    , nb    , id , itrouv
integer          izone
integer          nborne(nbtrmx)
integer          ilelt , nlelt

double precision rnbs2,capaeq
double precision sir(nelemx)   ,sii(nelemx)
double precision sirb(nbtrmx,6),siib(nbtrmx,6)
double precision ur(nbtrmx,6)  ,ui(nbtrmx,6)
double precision sirt(nbtrmx)  ,siit(nbtrmx)
character(len=200) :: chain

integer, allocatable, dimension(:) :: lstelt

double precision, dimension(:), pointer :: cpro_dji, cpro_djr
!< [loc_var_dec]

!===============================================================================

!===============================================================================
! Initialization
!===============================================================================

!< [init]
allocate(lstelt(nfabor))  ! temporary array for boundary faces selection
!< [init]

!===============================================================================
! Assign boundary conditions to boundary faces here

! For each subset:
! - use selection criteria to filter boundary faces of a given subset
! - loop on faces from a subset
!   - set the boundary condition for each face
!===============================================================================

! 1 - Computation of intensity (A/m2) for each electrode
!   -----------------------------------------------------

!< [pre_init]
do i= 1,nbelec
  sir(i) = 0.d0
  sii(i) = 0.d0
enddo

do ntf= 1,nbtrf
  sirt(ntf) = 0.d0
  siit(ntf) = 0.d0
enddo

if (ntcabs.lt.(ntpabs+2)) then
  do ntf = 1,nbtrf
    uroff(ntf) = 0.d0
    uioff(ntf) = 0.d0
  enddo
endif
!< [pre_init]

!     Loop on selected boundary faces

!< [step_1]
do i = 1,nbelec

  chain = ' '
  write(chain,3000) ielecc(i)

  if ( ielect(i).ne. 0 ) then
    call getfbr(chain,nlelt,lstelt)
    !==========

    do ilelt = 1, nlelt

      ifac = lstelt(ilelt)

      iel = ifabor(ifac)

      do id=1,ndimve
        call field_get_val_s(iprpfl(idjr(id)), cpro_djr)
        sir(i) = sir(i)                                           &
             + cpro_djr(iel)*surfbo(id,ifac)
      enddo

      if ( ippmod(ieljou) .eq. 4 ) then
        do id=1,ndimve
          call field_get_val_s(iprpfl(idji(id)), cpro_dji)
          sii(i) = sii(i)                                         &
               + cpro_dji(iel)*surfbo(id,ifac)
        enddo
      endif

    enddo

  endif

enddo
!< [step_1]

! 2 - Definition of Voltage on each termin of transformers
!----------------------------------------------------------------

!  2.1 Computation of Intensity on each termin of transformers

!< [step_2_1]
do i=1,nbelec
  sirb(ielect(i),ielecb(i)) = 0.d0
  if ( ippmod(ieljou) .eq. 4 ) then
    siib(ielect(i),ielecb(i)) = 0.d0
  endif
enddo

do i=1,nbelec
  if ( ielect(i).ne. 0 ) then
    sirb(ielect(i),ielecb(i)) = sirb(ielect(i),ielecb(i))         &
                               +sir(i)
    if ( ippmod(ieljou) .eq. 4 ) then
       siib(ielect(i),ielecb(i)) = siib(ielect(i),ielecb(i))      &
                                  +sii(i)
    endif
  endif
enddo
!< [step_2_1]

!  2.2 RVoltage on each termin

!< [step_2_2]
do ntf=1,nbtrf

!      Primary and Secondary in Triangle

  if (ibrpr(ntf) .eq. 0 .and. ibrsec(ntf) .eq. 0 ) then

    nborne(ntf) = 3

    rnbs2 = 3.d0*rnbs(ntf)*rnbs(ntf)
    ur(ntf,1)=  1.154675d0*tenspr(ntf)/rnbs(ntf)                  &
      + (zr(ntf)*sirb(ntf,1)-zi(ntf)*siib(ntf,1))/rnbs2

    ur(ntf,2)= -0.5773d0*tenspr(ntf)/rnbs(ntf)                    &
      + (zr(ntf)*sirb(ntf,2)-zi(ntf)*siib(ntf,2))/rnbs2
    ur(ntf,3)= -0.5773d0*tenspr(ntf)/rnbs(ntf)                    &
      + (zr(ntf)*sirb(ntf,3)-zi(ntf)*siib(ntf,3))/rnbs2

    ui(ntf,1)=  0.d0                                              &
      + (zi(ntf)*sirb(ntf,1)+zr(ntf)*siib(ntf,1))/rnbs2
    ui(ntf,2)= -1.d0*tenspr(ntf)/rnbs(ntf)                        &
      + (zi(ntf)*sirb(ntf,2)+zr(ntf)*siib(ntf,2))/rnbs2
    ui(ntf,3)=  1.d0*tenspr(ntf)/rnbs(ntf)                        &
      + (zi(ntf)*sirb(ntf,3)+zr(ntf)*siib(ntf,3))/rnbs2

  else

    write(nfecra, *) 'Matrice sur le Transfo a ecrire'
    call csexit(1)

  endif
enddo
!< [step_2_2]

!  2.3 Total intensity for a transformer
!         (zero valued WHEN Offset established)

!< [step_2_3]
do ntf=1,nbtrf
  sirt(ntf) = 0.d0
  if ( ippmod(ieljou) .eq. 4 ) then
    siit(ntf) = 0.d0
  endif
enddo

do i=1,nbelec
  if ( ielect(i).ne. 0 ) then
    sirt(ielect(i)) = sirt(ielect(i)) + sir(i)
    if ( ippmod(ieljou) .eq. 4 ) then
      siit(ielect(i)) = siit(ielect(i)) + sii(i)
    endif
  endif
enddo
!< [step_2_3]

!  2.4 Take in account of Offset

!< [step_2_4]
capaeq = 3.d0

do ntf=1,nbtrf
  uroff(ntf) = uroff(ntf) + sirt(ntf)/capaeq
  if ( ippmod(ieljou) .eq. 4 ) then
    uioff(ntf) = uioff(ntf) + siit(ntf)/capaeq
  endif
enddo

! A reference transformer is assumed to have an Offset zero valued

if ( ntfref .gt. 0 ) then
  uroff(ntfref) = 0.d0
  uioff(ntfref) = 0.d0
endif

do ntf=1,nbtrf
  do nb=1,nborne(ntf)
    ur(ntf,nb) = ur(ntf,nb) + uroff(ntf)
    if ( ippmod(ieljou) .eq. 4 ) then
      ui(ntf,nb) = ui(ntf,nb) + uioff(ntf)
    endif
  enddo
enddo
!< [step_2_4]

! Print of UROFF (real part of offset potential)

write(nfecra,1500)
do ntf=1,nbtrf
  write(nfecra,2000) ntf,uroff(ntf)
enddo
write(nfecra,1501)

!  2.5 Take in account of Boundary Conditions

!     Loop on selected Boundary Faces

!< [step_2_5]
do i=1,nbelec

  CHAIN = ' '
  write(chain,3000) ielecc(i)

  call getfbr(chain,nlelt,lstelt)
  !==========

  do ilelt = 1, nlelt

    ifac = lstelt(ilelt)

    iel = ifabor(ifac)

    itypfb(ifac) = iparoi

!     - Zone number
    izone = i

!      - Allocation of zone number
    izfppp(ifac) = izone

    if ( ielect(i) .ne. 0 ) then
      icodcl(ifac,isca(ipotr))   = 1
      rcodcl(ifac,isca(ipotr),1) = ur(ielect(i),ielecb(i))

      if ( ippmod(ieljou).eq.4  ) then
        icodcl(ifac,isca(ipoti))   = 1
        rcodcl(ifac,isca(ipoti),1) = ui(ielect(i),ielecb(i))
      endif

    else

      ii = ipotr
      icodcl(ifac,isca(ii))   = 3
      rcodcl(ifac,isca(ii),3) = 0.d0

      if ( ippmod(ieljou).eq. 4   ) then
        ii = ipoti
        icodcl(ifac,isca(ii))   = 3
        rcodcl(ifac,isca(ii),3) = 0.d0
      endif
    endif

  enddo

enddo
!< [step_2_5]

! 3 - Test, if not any reference transformer
!      a piece of wall may be at ground.

!< [step_3]
if ( ntfref .eq. 0 ) then

  itrouv = 0
  do ifac = 1, nfabor

    if ( itypfb(ifac) .eq. iparoi ) then

      if (icodcl(ifac,isca(ipotr)) .eq. 1 ) then

        if ( ippmod(ieljou).eq.3 ) then

          if ( abs(rcodcl(ifac,isca(ipotr),1)).lt.1.e-20 ) then
            itrouv = 1
          endif

        else if ( ippmod(ieljou).eq.4 ) then

          if (icodcl(ifac,isca(ipoti)) .eq. 1 ) then

            if (abs(rcodcl(ifac,isca(ipotr),1)).lt.1.e-20         &
           .and.abs(rcodcl(ifac,isca(ipoti),1)).lt.1.e-20 ) then
              itrouv = 1
            endif
          endif
        endif
      endif
    endif
  enddo

  if ( itrouv .eq. 0 ) then
    write(nfecra,1000)
    call csexit (1)
  endif

endif
!< [step_3]

!--------
! Formats
!--------

 1000 format(1X,' ERROR in JOULE : '                          ,/, &
       1X,' ====================   '                          ,/, &
      10X,' Lack of reference : choose a transformer for wich',/, &
      10X,' offset is assumed zero or a face at ground on the',/, &
          ' boundary')
 1500 format(/,2X,' ** INFORMATIONS ON TRANSFOMERS           ',/, &
         2X,'    ---------------------------------------'/    ,/, &
         1X,'      ---------------------------------'         ,/, &
         1X,'      Number of Transfo        UROFF    '        ,/, &
         1X,'      ---------------------------------')
 1501 format(1X,'      ---------------------------------')
 2000 format(10x,i6,12x,e12.5)
 3000 format(i7)

!----
! End
!----

deallocate(lstelt)  ! temporary array for boundary faces selection

return
end subroutine cs_user_boundary_conditions
