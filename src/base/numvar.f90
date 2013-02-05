!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

!> \file numvar.f90
!> \brief Module for variable numbering

module numvar

  !=============================================================================

  use paramx

  !=============================================================================

  integer, save :: ipr       !< pressure
  integer, save :: iu        !< velocity component \f$ u_x \f$
  integer, save :: iv        !< velocity component \f$ u_y \f$
  integer, save :: iw        !< velocity component \f$ u_z \f$
  integer, save :: ik        !< turbulent kinetic energy \f$ k \f$
  integer, save :: iep       !< turbulent dissipation \f$ \varepsilon \f$
  integer, save :: ir11      !< Reynolds stress component \f$ R_{xx} \f$
  integer, save :: ir22      !< Reynolds stress component \f$ R_{yy} \f$
  integer, save :: ir33      !< Reynolds stress component \f$ R_{zz} \f$
  integer, save :: ir12      !< Reynolds stress component \f$ R_{xy} \f$
  integer, save :: ir23      !< Reynolds stress component \f$ R_{yz} \f$
  integer, save :: ir13      !< Reynolds stress component \f$ R_{zz} \f$
  integer, save :: iphi      !< variable \f$ \phi \f$ of the \f$ \phi-f_b \f$ model
  integer, save :: ifb       !< variable \f$ f_b \f$ of the \f$ \phi-f_b \f$ model
  integer, save :: ial       !< variable \f$ \alpha \f$ of the \f$ Bl-v^2-k \f$ model
  integer, save :: iomg      !< variable \f$ \omega \f$ of the \f$ k-\omega \f$ SST
  integer, save :: inusa     !< variable \f$ \widetilde{\nu}_T \f$ of the Spalart Allmaras
  integer, save :: isca(nscamx) !< isca(i) is the index of the scalar i
  integer, save :: iscapp(nscamx) !< iscapp(i) is the index of the specific physics scalar i
  integer, save :: nscaus    !< number of user scalars
  integer, save :: nscapp    !< number of specific physics scalars
  integer, save :: nscasp    !< number of species scalars
  integer, save :: iuma      !< mesh velocity component \f$ w_x \f$
  integer, save :: ivma      !< mesh velocity component \f$ w_y \f$
  integer, save :: iwma      !< mesh velocity component \f$ w_z \f$

  ! Position des proprietes (physiques ou numeriques)
  !  (dans propce, propfa et propfb)
  !    le numero des proprietes est unique, quelle aue soit la
  !      localisation de ces dernieres (cellule, face, face de bord)
  !    Voir cs_user_boundary_conditions pour quelques exemples

  ! ipproc : pointeurs dans propce
  ! ipprof : pointeurs dans propfa
  ! ipprob : pointeurs dans propfb

  ! irom   : Density at the current time step
  ! iroma  : Density at the previous time step
  ! iviscl : Viscosite moleculaire dynamique en kg/(ms) des phases
  ! ivisct : Viscosite turbulente des phases
  ! ivisla : Viscosite moleculaire dynamique en kg/(ms) des phases au pas
  !          de temps precedent
  ! ivista : Viscosite turbulente des phases au pas de temps precedent
  ! icp    : Chaleur specifique des phases
  ! icpa   : Chaleur specifique des phases au pas de temps precedent
  ! itsnsa : Terme source Navier Stokes des phases au pas de temps precedent
  ! itstua : Terme source des grandeurs turbulentes au pas de temps precedent
  ! itssca : Terme source des scalaires au pas de temps precedent
  ! iestim : Estimateur d'erreur pour Navier-Stokes
  ! ifluma : Flux de masse associe aux variables
  ! ifluaa : Flux de masse explicite (plus vu comme un tableau de travail)
  !          associe aux variables
  ! ismago : constante de Smagorinsky dynamique
  ! icour  : Nombre de Courant des phases
  ! ifour  : Nombre de Fourier des phases
  ! iprtot : Pression totale au centre des cellules Ptot=P*+rho*g.(x-x0)
  !                                                             -  - -
  ! ivisma : Viscosite de maillage en ALE (eventuellement orthotrope)
  ! iustdy : pointer for dilatation source terms
  ! itsrho : pointer for global dilatation source terms
  ! ibetau : pointer for thermal expansion coefficient


  integer, save :: ipproc(npromx), ipprof(npromx), ipprob(npromx), &
                   irom  , iroma , iviscl,                         &
                   ivisct, ivisla, ivista,                         &
                   icp   , icpa  , itsnsa,                         &
                   itstua, itssca(nscamx),                         &
                   iestim(nestmx)         , ifluma(nvarmx),        &
                   ifluaa(nvarmx), ismago, icour ,                 &
                   ifour , iprtot, ivisma(3),                      &
                   iustdy(nscamx), itsrho,                         &
                   ibeta

  ! Position des conditions aux limites
  !  (position dans coefa et coefb des coef (coef. coef.f) relatifs a
  !   une variable donnee)

  ! icoef   : coef numeros 1 (anciens coefa coefb)
  ! icoeff  : coef numeros 2 (anciens coefaf coefbf)
  ! icoefr  : coef number 3 (for the Rij in the momentum eq)
  ! iclrtp  : pointeur dans COEFA et COEFB

  integer, save :: icoef , icoeff , icoefr , iclrtp(nvarmx,3)

  ! Mapping to field structures

  ! ivarfl(i)                  Field id for variable i
  ! iprpfl(i)                  Field id for property i

  integer, save :: ivarfl(nvarmx), iprpfl(npromx)

  !=============================================================================

end module numvar
