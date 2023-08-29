!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!> \brief Automatic boundary condition for pulverized coal combution
!>        with Lagrangian module.

!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role
!______________________________________________________________________________!
!> \param[in]     itypfb        boundary face types
!> \param[in]     izfppp        zone number for the boundary face for
!>                                      the specific physic module
!!> \param[in,out] rcodcl        value of the boundary conditions to edge faces
!>
!>                              boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                               -  coefficient (infinite if no exchange)
!>                               -  rcodcl(3) value flux density
!>                               -  (negative if gain) \f$w.m^{-2} \f$ or
!>                               -  roughness in \f$m\f$ if  icodcl=6
!>                                -# for velocity:
!>                                           \f$(\mu+\mu_T)\gradv \vect{u}\f$
!>                                -# for pressure: \f$ \Delta \grad P
!>                                                 \cdot \vect{n} \f$
!>                                -# for scalar:   \f$ C_p \left ( K +
!>                                                 \dfrac{K_T}{\sigma_T} \right)
!>                                                 \grad T \cdot \vect{n} \f$
!______________________________________________________________________________!

subroutine cpltcl &
 ( itypfb , izfppp ,                                              &
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
use dimens, only : nvar
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          izfppp(nfabor)

double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ii, ifac, izone, mode, iel, ige, iok
integer          icha, icke
integer          nbrval
double precision qisqc, viscla, d2s3, uref2, rhomoy, dhy, xiturb
double precision t1
double precision h1    (nozppm)
double precision qimpc (nozppm) , qcalc(nozppm)
double precision coefe (ngazem)
double precision f1mc  (ncharm) , f2mc (ncharm)
double precision, dimension(:), pointer ::  brom
double precision, dimension(:), pointer :: viscl

!===============================================================================
! 0. Initializations
!===============================================================================

call field_get_val_s(ibrom, brom)
call field_get_val_s(iviscl, viscl)

d2s3 = 2.d0/3.d0

!===============================================================================
! 1.  Parallel exchanges for the user data
!===============================================================================
!  In fact this exchange could be avoided by changing uscpcl and by asking
!    the user to give the variables which depend of the area out of the loop
!    on the boundary faces: the variables would be available on all processors.
!  However, it makes the user subroutine a bit more complicated and especially
!    if the user modifies it through, it does not work.
!  We assume that all the provided variables are positive,
!    which allows to use a max for the proceedings know them.
!  If this is not the case, it is more complicated but we can get a max anyway.

if (irangp.ge.0) then
  call parimx(nozapm,iqimp )
  call parimx(nozapm,ientat)
  call parimx(nozapm,ientcp)
  call parrmx(nozapm,qimpat)
  call parrmx(nozapm,timpat)
  nbrval = nozppm*ncharm
  call parrmx(nbrval,qimpcp)
  nbrval = nozppm*ncharm
  call parrmx(nbrval,timpcp)
  nbrval = nozppm*ncharm*nclcpm
  call parrmx(nbrval,distch)
endif

!===============================================================================
! 2.  Correction of the velocities (in norm) for controlling the imposed flow
!       Loop over all inlet faces
!                     =========================
!===============================================================================

! --- Calculated flow
do izone = 1, nozppm
  qcalc(izone) = 0.d0
enddo
do ifac = 1, nfabor
  izone = izfppp(ifac)
  qcalc(izone) = qcalc(izone) - brom(ifac) *             &
                ( rcodcl(ifac,iu,1)*surfbo(1,ifac) +       &
                  rcodcl(ifac,iv,1)*surfbo(2,ifac) +       &
                  rcodcl(ifac,iw,1)*surfbo(3,ifac) )
enddo

if (irangp.ge.0) then
  call parrsm(nozapm,qcalc )
endif

do izone = 1, nozapm
  if (iqimp(izone).eq.0) then
    qimpc(izone) = qcalc(izone)
  endif
enddo

! Rescaling of the velocities (for mass flow)

iok = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  if (iqimp(izone).eq.1) then
    if(abs(qcalc(izone)).lt.epzero) then
      write(nfecra,2001)izone,iqimp(izone),qcalc(izone)
      iok = iok + 1
    endif
  endif
enddo

if (iok.ne.0) then
  call csexit(1)
endif

do ifac = 1, nfabor
  izone = izfppp(ifac)
  if (iqimp(izone).eq.1) then
    qimpc(izone) = qimpat(izone)
    qisqc = qimpc(izone)/qcalc(izone)
    rcodcl(ifac,iu,1) = rcodcl(ifac,iu,1)*qisqc
    rcodcl(ifac,iv,1) = rcodcl(ifac,iv,1)*qisqc
    rcodcl(ifac,iw,1) = rcodcl(ifac,iw,1)*qisqc
  endif
enddo

 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE PHYSIQUES PARTICULIERES              ',/,&
'@    =========                                               ',/,&
'@     COMBUSTION CHARBON PULVERISE COUPLE AU                 ',/,&
'@     TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON :       ',/,&
'@     PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/&
'@                                                            ',/,&
'@  Le debit est impose sur la zone IZONE =     ', I10         ,/,&
'@    puisque                IQIMP(IZONE) =     ', I10         ,/,&
'@  Or, sur cette zone, le produit RHO D S integre est nul :  ',/,&
'@    il vaut                             = ',E14.5            ,/,&
'@    (D est la direction selon laquelle est impose le debit).',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier uscpcl, et en particulier                        ',/,&
'@    - que le vecteur  RCODCL(IFAC,IU,1),             ',/,&
'@                      RCODCL(IFAC,IV,1),             ',/,&
'@                      RCODCL(IFAC,IW,1) qui determine',/,&
'@      la direction de la vitesse est non nul et n''est pas  ',/,&
'@      uniformement perpendiculaire aux face d''entree       ',/,&
'@    - que la surface de l''entree n''est pas nulle (ou que  ',/,&
'@      le nombre de faces de bord dans la zone est non nul)  ',/,&
'@    - que la masse volumique n''est pas nulle               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 3.  Filling the table of the boundary conditions
!===============================================================================

call cs_boundary_conditions_legacy_turbulence(itypfb)

!===============================================================================
! 4.  Filling the boundary conditions table
!     Loop on all input faces
!                     =========================
!     We determine the family and its properties
!     We impose the boundary conditions
!     for the scalars
!===============================================================================

do ii = 1, nzfppp

  izone = ilzppp(ii)

  ! An input ientre must be of type
  ! ientat = 1 (or ientcp = 1 ?)

  if (ientat(izone).eq.1) then

    ! ------ Calculating H1(izone)
    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(io2) = wmole(io2) / (wmole(io2)+xsi*wmole(in2))
    coefe(in2) = 1.d0 - coefe(io2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    t1   = timpat(izone)
    mode = -1
    call cpthp1(mode, h1(izone), coefe, f1mc, f2mc, t1)

  endif

enddo

do ifac = 1, nfabor

  izone = izfppp(ifac)

  if ( itypfb(ifac).eq.ientre ) then

    ! Automatic processing of specific physic scalars

    do icha = 1, ncharb

! ------ CL pour F1M et F2M du charbon ICHA
      rcodcl(ifac,isca(if1m(icha)),1) = zero
      rcodcl(ifac,isca(if2m(icha)),1) = zero

    enddo

! ------ CL pour F3M
    rcodcl(ifac,isca(if3m),1) = zero
! ------ CL pour FP4M
    rcodcl(ifac,isca(if4p2m),1)   = zero
! ------ CL pour HM
    do ige = 1, ngazem
      coefe(ige) = zero
    enddo
    coefe(io2) = wmole(io2) / (wmole(io2)+xsi*wmole(in2))
    coefe(in2) = 1.d0 - coefe(io2)
    do icha = 1, ncharm
      f1mc(icha) = zero
      f2mc(icha) = zero
    enddo
    t1   = timpat(izone)
    mode = -1
    call cpthp1                                                   &
    !==========
    ( mode  , h1(izone) , coefe  , f1mc   , f2mc   ,              &
      t1    )
    rcodcl(ifac,isca(iscalt),1) = h1(izone)

  endif

enddo

!----
! End
!----

return
end subroutine
