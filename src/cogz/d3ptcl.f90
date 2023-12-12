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
! Function :
! --------

!> \file d3ptcl.f90
!> \brief Automatic boundary conditions for 3 PTHEM gas diffusion flame model

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     itypfb        boundary face types
!> \param[out]    izfppp        boundary face zone number
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
!> \param[out]    rcodcl        boundary condition values:
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

subroutine d3ptcl &
 ( itypfb , izfppp , icodcl,                                             &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens, only : nscal, nvar
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          izfppp(nfabor)
integer          icodcl(nfabor,nvar)

double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          igg, ifac, izone, iscal
integer          ii, ifue, ioxy
double precision qimabs, qimpl(1), qisqc
double precision qcalc(nozppm)
double precision coefg(ngazgm)
double precision, dimension(:), pointer ::  brom
double precision, dimension(:), pointer :: viscl

!===============================================================================
! 0.  INITIALISATIONS
!===============================================================================

call field_get_val_s(ibrom, brom)
call field_get_val_s(iviscl, viscl)

do igg = 1, ngazgm
  coefg(igg) = zero
enddo

!===============================================================================
! 1.  ECHANGES EN PARALLELE POUR LES DONNEES UTILISATEUR
!===============================================================================

!  En realite on pourrait eviter cet echange en modifiant usd3pc et en
!    demandant a l'utilisateur de donner les grandeurs dependant de la
!    zone hors de la boucle sur les faces de bord : les grandeurs
!    seraient ainsi disponibles sur tous les processeurs. Cependant,
!    ca rend le sous programme utilisateur un peu plus complique et
!    surtout, si l'utilisateur le modifie de travers, ca ne marche pas.
!  On suppose que toutes les gandeurs fournies sont positives, ce qui
!    permet d'utiliser un max pour que tous les procs les connaissent.
!    Si ce n'est pas le cas, c'est plus complique mais on pourrait
!    s'en tirer avec un max quand meme, sauf qimp qui peut etre negatif.

if (irangp.ge.0) then
  call parmax(tinfue)
  call parmax(tinoxy)
  ! qimabs  = |qimp|, or 0 if no boundary faces on the rank
  do ii = 1, nozapm
    qimabs = abs(qimp(ii))
    qimpl(1) = qimp(ii)
    ! the rank owning the max of qimabs, meaning |qimp|, share
    ! qimp with the others
    call parmxl(1, qimabs, qimpl)
    qimp(ii) = qimpl(1)
  enddo
  call parimx(nozapm,iqimp )
  call parimx(nozapm,ientox)
  call parimx(nozapm,ientfu)
endif

!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES POUR LES SCALAIRES
!
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================

!       ON DETERMINE LA FAMILLE ET SES PROPRIETES
!       ON IMPOSE LES CONDITIONS AUX LIMITES
!===============================================================================

! Fuel and Oxidant inlet are checked
ifue = 0
ioxy = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  if (ientfu(izone).eq.1) then
    ifue = 1
  elseif(ientox(izone).eq.1) then
    ioxy = 1
  endif
enddo
if(irangp.ge.0) then
  call parcmx(ifue)
  call parcmx(ioxy)
endif

! Fuel inlet at TINFUE: HINFUE calculation
if (ifue.eq.1) then
  coefg(1) = 1.d0
  coefg(2) = zero
  coefg(3) = zero
  hinfue = cs_gas_combustion_t_to_h(coefg, tinfue)
endif

! Oxidant inlet at TINOXY: HINOXY calculation
if (ioxy.eq.1) then
  coefg(1) = zero
  coefg(2) = 1.d0
  coefg(3) = zero
  hinoxy = cs_gas_combustion_t_to_h(coefg, tinoxy)
endif


do ifac = 1, nfabor

  izone = izfppp(ifac)

  if (itypfb(ifac).eq.ientre.or.itypfb(ifac).eq.i_convective_inlet) then

    ! Fuel inlet at TINFUE
    if (ientfu(izone).eq.1) then

      ! Neumann for outflow
      if (qimp(izone).lt.0.d0) then

        do iscal = 1, nscal
          icodcl(ifac,isca(iscal)) = 3
          rcodcl(ifac,isca(iscal),3) = 0.d0
        enddo

      ! Dirichlet for inflow
      else
        ! Mean mixture fraction
        rcodcl(ifac,isca(ifm),1)   = 1.d0

        ! Mixture fraction variance
        rcodcl(ifac,isca(ifp2m),1) = 0.d0

        ! Mixture enthalpy
        if (ippmod(icod3p).eq.1) then
          rcodcl(ifac,isca(iscalt),1) = hinfue
        endif

        ! Soot model
        if (isoot.ge.1) then
          rcodcl(ifac,isca(ifsm),1) = 0.d0
          rcodcl(ifac,isca(inpm),1) = 0.d0
        endif

      endif

      ! Density
      brom(ifac) = pther/(cs_physical_constants_r*tinfue/wmolg(1))

    ! Oxydant inlet at TINOXY
    elseif (ientox(izone).eq.1) then

      ! Neumann for outflow
      if (qimp(izone).lt.0.d0) then

        do iscal = 1, nscal
          icodcl(ifac,isca(iscal)) = 3
          rcodcl(ifac,isca(iscal),3) = 0.d0
        enddo

      ! Dirichlet for inflow
      else

        ! Mean mixture fraction
        rcodcl(ifac,isca(ifm),1)   = 0.d0

        ! Mixture fraction variance
        rcodcl(ifac,isca(ifp2m),1) = 0.d0

        ! Mixture enthalpy
        if ( ippmod(icod3p).eq.1 ) then
          rcodcl(ifac,isca(iscalt),1) = hinoxy
        endif

        ! Soot model
        if (isoot.ge.1) then
          rcodcl(ifac,isca(ifsm),1) = 0.d0
          rcodcl(ifac,isca(inpm),1) = 0.d0
        endif

        ! Density
        brom(ifac) = pther/(cs_physical_constants_r*tinoxy/wmolg(2))

      endif

    endif

  endif

enddo

!===============================================================================
! 3.  SI IQIMP = 1 : CORRECTION DES VITESSES (EN NORME) POUR CONTROLER
!                    LES DEBITS IMPOSES
!     SI IQIMP = 0 : CALCUL DE QIMP

!     ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                   =========================
!===============================================================================

! --- Debit calcule

do izone = 1, nozppm
  qcalc(izone) = 0.d0
enddo
do ifac = 1, nfabor
  izone = izfppp(ifac)
  qcalc(izone) = qcalc(izone) - brom(ifac) *             &
     ( rcodcl(ifac,iu,1)*surfbo(1,ifac) +                  &
       rcodcl(ifac,iv,1)*surfbo(2,ifac) +                  &
       rcodcl(ifac,iw,1)*surfbo(3,ifac) )
enddo

if(irangp.ge.0) then
  call parrsm(nozapm,qcalc)
endif

do izone = 1, nozapm
  if (iqimp(izone).eq.0) then
    qimp(izone) = qcalc(izone)
  endif
enddo

! --- Correction des vitesses en norme

do ifac = 1, nfabor
  izone = izfppp(ifac)
  if (iqimp(izone).eq.1) then
    if (abs(qcalc(izone)).gt.epzero) then
      qisqc = qimp(izone)/qcalc(izone)
    else
      qisqc = 0.d0
    endif
    rcodcl(ifac,iu,1) = rcodcl(ifac,iu,1)*qisqc
    rcodcl(ifac,iv,1) = rcodcl(ifac,iv,1)*qisqc
    rcodcl(ifac,iw,1) = rcodcl(ifac,iw,1)*qisqc
  endif
enddo

!===============================================================================
! 4.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES POUR LA TURBULENCE
!===============================================================================

call cs_boundary_conditions_legacy_turbulence(itypfb)

!----
! End
!----

return
end subroutine
