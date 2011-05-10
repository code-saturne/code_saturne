!-------------------------------------------------------------------------------

!VERS


!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine uscfpv &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   maxelt , lstelt ,                                              &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     ,                                     &
   ra     )

!===============================================================================
! Purpose:
! -------

!    User subroutine.

!    Set (variable) physical properties for the compressible flow scheme.


! Description
! ===========

! This subroutine replaces the user subroutine 'usphyv' for the
! compressible flow scheme.

! This subroutine is called at the beginning of each time step.

! At the very first time step (not at restart), the only variables that
! have been initialized are those provided:
!   - in the GUI and in the user subroutines 'usini1' and 'uscfx2'; ex.:
!     . the density             (set to ro0)
!     . the molecular viscosity (set to viscl0)
!     . the volumetric molecular viscosity (set to viscv0)
!     . the molecular thermal conductivity (set to visls0(itempk))
!   - in the user subroutines 'usiniv' and 'uscfxi'; ex.:
!     . the unknown variables (null by default)

! This subroutine allows the user to set the cell values for:
!   - the molecular viscosity                            viscl  kg/(m s)
!   - the isobaric specific heat (cp=dh/dT|P)            cp     J/(kg degree)
!   - the molecular thermal conductivity                 lambda W/(m degree)
!   - the molecular diffusivity for user-defined scalars viscls kg/(m s)


! Warnings
! ========

! The density ** must not ** be set here: for the compressible scheme,
! it is one of the unknowns, and it can be initialized as such in the user
! subroutine 'uscfxi' (rtp array).

! The turbulent viscosity ** must not ** be modified here (to modify this
! variable, use the user subroutine 'usvist')

! To set a variable isobaric specific heat, the integer icp must
! have been set to 1: the value for icp is set automatically in the
! subroutine 'uscfth', depending on the thermodynamics laws selected
! by the user.

! To set a variable diffusivity for a given user-defined scalar, the
! variable ivisls(scalar_number) must have been set to 1 in the user
! subroutine 'usini1' or in the GUI (otherwise, a memory problem is
! expected).

! Examples are provided in the present subroutine (but they do not have
! any physical signification).


! Cells identification
! ====================

! Cells may be identified using the 'getcel' subroutine.
! The syntax of this subroutine is described in the 'usclim' subroutine,
! but a more thorough description can be found in the user guide.

! The type of the boundary faces at the previous time step is available
! (except at the first time step, since the arrays itypfb and itrifb have
! not yet been set);


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! maxelt           ! i  ! <-- ! max number of cells and faces (int/boundary)   !
! lstelt(maxelt)   ! ia ! --- ! work array                                     !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and preceding time steps)         !
! propce(ncelet, *)! ra ! <-> ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...3(ncelet    ! ra ! --- ! work arrays                                    !
! ra(*)            ! ra ! --- ! main real work array                           !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use pointe
use numvar
use optcal
use cstphy
use entsor
use period
use ppppar
use ppthch
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          maxelt, lstelt(maxelt)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra
integer          ivart, iel
integer          ipcvis, ipcvsv, ipccp
integer          ipcvsl, ith, iscal, ii, iccfth, imodif
double precision varam, varbm, varcm, vardm
double precision varal, varbl, varcl, vardl
double precision varac, varbc
double precision xrtp

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
!===============================================================================
! 0.  This test allows the user to ensure that the version of this subroutine
!       used is that from his case definition, and not that from the library.
!     However, this subroutine may not be mandatory,
!       thus the default (library reference) version returns immediately.
!===============================================================================

if(1.eq.1) then
  iuscfp = 0
  return
endif


! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. Mandatory initializations
!===============================================================================

! --- Memory initialization

idebia = idbia0
idebra = idbra0

!===============================================================================

! Warning: the examples provided below are physically meaningless.
! =======

! These examples must be adapted by the user. Hence, the default
! (library reference) version stops immediately after each example
! (the 'call csexit(1)' directive must be discarded to use the
! portion of code).

! It is adviced to discard all the examples that are not necessary, so
! as to minimize the risk of error.

! List of examples
! ================

! Ex. 1: molecular viscosity varying with temperature
! Ex. 2: molecular volumetric viscosity varying with temperature
! Ex. 3: isobaric specific heat varying with temperature
! Ex. 4: molecular thermal conductivity varying with temperature
! Ex. 5: molecular diffusivity of user-defined scalars varying with temperature

!===============================================================================


!===============================================================================
! Ex. 1: molecular viscosity varying with temperature
! =====
!    The values of the molecular viscosity are provided as a function of
!    the temperature. All variables are evaluated at the cell centres.
!===============================================================================

! --- Rank of the temperature in the array 'rtp'
!     To refer to the user-defined scalar number 2 instead, for example, use
!     ivart = isca(2)

ivart = isca(itempk)

! --- Rank 'ipcvis' of the molecular dynamic viscosity
!     in the array 'propce' (physical properties at the cell centers)

ipcvis = ipproc(iviscl)

! --- User-defined coefficients for the selected law.
!     The values hereafter are provided as a mere example. They
!     are physically meaningless.

varam = -3.4016d-9
varbm =  6.2332d-7
varcm = -4.5577d-5
vardm =  1.6935d-3

! --- Molecular dynamic viscosity mu at the cell centres, kg/(m s)
!     In this example, mu is provided as a function of the temperature T:
!       mu(T)              =    T  *( T  *( am  * T +  bm  )+ cm  )+ dm
!     that is:
!       propce(iel,ipcvis) =   xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm

do iel = 1, ncel
  xrtp = rtp(iel,ivart)
  propce(iel,ipcvis) =                                          &
       xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm
enddo


! --- Discard the following test so that the code do not stop
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif


!===============================================================================
! Ex. 2: molecular volumetric viscosity varying with temperature
! =====
!    The values of the molecular volumetric viscosity are provided as a function
!    of the temperature. All variables are evaluated at the cell centres.
!===============================================================================

! --- Rank of the temperature in the array 'rtp'
!     To refer to the user-defined scalar number 2 instead, for example, use
!     ivart = isca(2)

ivart = isca(itempk)

! --- Rank 'ipcvsv' of the molecular dynamic viscosity
!     in the array 'propce' (physical properties at the cell centers)

if(iviscv.gt.0) then
  ipcvsv = ipproc(iviscv)
else
  ipcvsv = 0
endif

! --- Stop if the viscosity has not been defined as variable

if(ipcvsv.le.0) then
  write(nfecra,2000) iviscv
  call csexit (1)
endif

! --- User-defined coefficients for the selected law.
!     The values provided hereafter are provided as a mere example. They
!     are physically meaningless.

varam = -3.4016d-9
varbm =  6.2332d-7
varcm = -4.5577d-5
vardm =  1.6935d-3

! --- Molecular dynamic volumetric viscosity kappa at the cell centres, kg/(m s)
!     In this example, kappa is provided as a function of the temperature T:
!       kappa(T)           =    T  *( T  *( am  * T +  bm  )+ cm  )+ dm
!     that is:
!       propce(iel,ipcvsv) =   xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm

do iel = 1, ncel
  xrtp = rtp(iel,ivart)
  propce(iel,ipcvsv) =                                          &
       xrtp*(xrtp*(varam*xrtp+varbm)+varcm)+vardm
enddo

! --- Discard the following test so that the code do not stop
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif


!===============================================================================
! Ex. 3: isobaric specific heat varying with temperature
! =====
!    The values of the isobaric specific heat values are provided as a function
!    of the temperature. All variables are evaluated at the cell centres.
!===============================================================================

! Warning:
! =======
! do not discard the call to the subroutine 'usthht' at the end of this
! example: its purpose is to calculate the isochoric specific heat.
! Indeed, this variable needs to be computed from the isobaric specific heat
! using the thermodynamics laws.

! --- Rank of the temperature in the array 'rtp'
!     To refer to the user-defined scalar number 2 instead, for example, use
!     ivart = isca(2)

ivart = isca(itempk)

! --- Rank 'ipcpp' of the isobaric specific heat
!     in the array 'propce' (physical properties at the cell
!     centers)

if(icp.gt.0) then
  ipccp  = ipproc(icp   )
else
  ipccp  = 0
endif

! --- Stop if the iobaric or iochoric specific heat (cp or cv) has not
!     been defined as variable

if(ipccp.le.0) then
  write(nfecra,1000) icp
  call csexit (1)
endif
if(icv.le.0) then
  write(nfecra,1001) icv
  call csexit (1)
endif

! --- User-defined coefficients for the selected law.
!     The values provided hereafter are provided as a mere example. They
!     are physically meaningless.

varac = 0.00001d0
varbc = 1000.0d0

! --- Isobaric specific heat cp at the cell centres, J/(kg degree)
!     In this example, cp is provided as a function of the temperature T:
!       cp(T)              =      ac * T  + ab
!     that is:
!       propce(iel,ipccp ) =    varac*xrtp+varbc

do iel = 1, ncel
  xrtp = rtp(iel,ivart)
  propce(iel,ipccp ) = varac*xrtp + varbc
enddo

! --- The isochoric specific heat is deduced from the isobaric specific
!     heat using the subroutine 'uscfth'.

iccfth = 432
imodif = 0

call uscfth                                                     &
!==========
 ( nvar   , nscal  ,                                              &
   iccfth , imodif ,                                              &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   propce(1, ipproc(icv) )    , w1     , w2     , w3     )

! --- Discard the following test so that the code do not stop
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif


!===============================================================================
! Ex. 4: molecular thermal conductivity varying with temperature
! =====
!    The values of the molecular thermal conductivity are provided as a function
!    of the temperature. All variables are evaluated at the cell centres.
!===============================================================================

! --- Rank of the temperature in the array 'rtp'
!     To refer to the user-defined scalar number 2 instead, for example, use
!     ivart = isca(2)

ivart = isca(itempk)

! --- Rank 'ipcvsl' of the olecular thermal conductivity
!     in the array 'propce' (physical properties at the cell
!     centers)

if(ivisls(itempk).gt.0) then
  ipcvsl = ipproc(ivisls(itempk))
else
  ipcvsl = 0
endif

! --- Stop if the molecular thermal conductivity has not
!     been defined as variable

if(ipcvsl.le.0) then
  write(nfecra,1010)                                            &
       itempk, itempk, ivisls(itempk)
  call csexit (1)
endif

! --- User-defined coefficients for the selected law.
!     The values provided hereafter are provided as a mere example. They
!     are physically meaningless.

varal = -3.3283d-7
varbl =  3.6021d-5
varcl =  1.2527d-4
vardl =  0.58923d0

! --- Molecular thermal conductivity lambda at the cell centres, W/(m degree)
!     In this example, lambda is provided as a function of the temperature T:
!       lambda(T)          =    T  *( T  *( al  * T +  bl  )+ cl  )+ dl
!     that is:
!       propce(iel,ipcvsl) =   xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl

do iel = 1, ncel
  xrtp = rtp(iel,ivart)
  propce(iel,ipcvsl) =                                          &
       (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)
enddo

! --- Discard the following test so that the code do not stop
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif


!===============================================================================
! Ex. 5: molecular diffusivity of user-defined scalars varying with temperature
! =====
!    The molecular diffusivity can be set for all the user-defined scalars
!    ** except **:
!      - temperature and enthalpy (already dealt with above: for these
!        variables, the 'diffusivity' is the thermal conductivity)
!      - variances of the fluctuations of another scalar variable (the
!        diffusivity is assumed to be equal to that of the associated
!        scalar)
!    The values of the molecular diffusivity are provided as a function
!    of the temperature. All variables are evaluated at the cell centres.
!===============================================================================

! --- Loop on the scalars
do ii = 1, nscaus

! --- Rank of the ii-th scalar in the list of all scalars
  iscal = ii


! --- If the scalar is the temperature, it is marked by ith = 1
!     so that it will be skipped.

  ith = 0
  if (iscal.eq.itempk) ith = 1

! --- If the variable represents the variance of the fluctuations of
!     another scalar variable (iscavr <= 0), it is simply skipped.

  if (ith.eq.0.and.iscavr(iscal).le.0) then

! --- Here, iscal points to any scalar variable except the temperature,
!     the enthalpy and the variance of the fluctuations of another
!     scalar variable.

! --- Rank of the temperature in the array 'rtp'
!     To refer to the user-defined scalar number 2 instead, for example, use
!     ivart = isca(2)

  ivart = isca(itempk)

! --- Rank 'ipcvsl' of the molecular diffusivity of the current scalar iscal
!     in the array 'propce' (physical properties at the cell centers)

    if(ivisls(iscal).gt.0) then
      ipcvsl = ipproc(ivisls(iscal))
    else
      ipcvsl = 0
    endif

! --- Stop if the molecular diffusivity has not been defined as variable

    if(ipcvsl.le.0) then
      write(nfecra,1010) iscal, iscal, ivisls(iscal)
      call csexit (1)
    endif

! --- User-defined coefficients for the selected law.
!     The values provided hereafter are provided as a mere example. They
!     are physically meaningless.

    varal = -3.3283d-7
    varbl =  3.6021d-5
    varcl =  1.2527d-4
    vardl =  0.58923d0

! --- Molecular diffusivity lambda at the cell centres, kg/(m s)
!     In this example, lambda is provided as a function of the temperature T:
!       lambda(T)          =    T  *( T  *( al  * T +  bl  )+ cl  )+ dl
!     that is:
!       propce(iel,ipcvsl) =   xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl

    do iel = 1, ncel
      xrtp = rtp(iel,ivart)
      propce(iel,ipcvsl) =                                        &
           (xrtp*(xrtp*(varal*xrtp+varbl)+varcl)+vardl)
    enddo


  endif
! --- End of the tests on ith and iscavr

enddo
! --- End of the loop on the scalars


! --- Discard the following test so that the code do not stop
if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif


!----
! Formats
!----

 1000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the GUI or in the user subroutine ''usini1'', the',/, &
'@         isobaric specific heat is declared as a property',/,   &
'@         uniform in space: icp = ',i10   ,/,                    &
'@       in the user subroutine ''uscfpv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@',/,                                                            &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usini1'' or ''uscfpv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1001 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the GUI or in the user subroutine ''usini1'', the',/, &
'@         isochoric specific heat is declared as a property',/,  &
'@         uniform in space: icv = ',i10   ,/,                    &
'@       in the user subroutine ''uscfpv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usini1'' or ''uscfpv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 1010 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@     For the scalar ',i10,/,                                    &
'@       in the GUI or in the user subroutine ''usini1'', the',/, &
'@         molecular diffusivity is declared as a property',/,    &
'@         uniform in space: ivisls(',i10   ,') = ',i10   ,/,     &
'@       in the user subroutine ''uscfpv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the GUI input data or the',/, &
'@    user subroutines ''usini1'' or ''uscfpv''.',/,              &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 2000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     The data is inconsistent',/,                               &
'@',/,                                                            &
'@       in the user subroutine ''uscfx2'', the molecular',/,     &
'@         volumetric viscosity is declared as a property',/,     &
'@         uniform in space: iviscv = ',i10,/,                    &
'@       in the user subroutine ''uscfpv'', however, it is',/,    &
'@         assumed to be potentially non uniform in space.',/,    &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Ensure consistency by modifying the user subroutines',/,      &
'@    ''uscfx2'' or ''uscfpv''.',/,                               &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)
 9000 format(                                                     &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/,                                                            &
'@ @@ WARNING:    stop in computation of physical properties',/,  &
'@    =======',/,                                                 &
'@     Call to ''csexit'' from the user subroutine ''uscfpv''.',/,&
'@',/,                                                            &
'@     The subroutine ''csexit'' (run stop) was called from ',/,  &
'@       within the user subroutine ''uscfpv''. The user shall',/,&
'@       ensure that all the default examples provided in the',/, &
'@       reference version of the user subroutine have been',/,   &
'@       discarded. It shall also be checked that there is no',/, &
'@       remaining stopping test at the end of the examples ',/,  &
'@       that have been retained.',/,                             &
'@@',/,                                                           &
'@  The calculation will not be run.',/,                          &
'@',/,                                                            &
'@  Check and modify the user subroutine ''uscfpv''.',/,          &
'@',/,                                                            &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',/)


!----
! End
!----

return
end subroutine
