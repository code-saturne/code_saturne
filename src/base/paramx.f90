!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

!> \file paramx.f90
!> \brief Module for definition of general parameters

module paramx

  !=============================================================================

  implicit none

  !> \defgroup paramx Module for definition of general parameters

  !> \addtogroup paramx
  !> \{

  !=============================================================================

  !> maximum number of scalars solutions of an
  !> advection equation, apart from the variables of the turbulence model
  !> \f$ (k, \varepsilon, R_{ij}, \omega, \varphi, \overline{f}, \alpha, \nu_t\f$)
  !> , that is to say
  !> the temperature and other scalars (passive or not, user-defined or not)
  integer   nscamx

  !> maximal number of variables = nscamx + 12 (u,v,w,P,Rij,e,alp)
  integer   nvarmx

  parameter(nscamx=200)
  parameter(nvarmx=nscamx+12)

  !> Maximal possible boundary condition types
  integer    ntypmx
  parameter(ntypmx=200)

  !> \anchor iindef
  !> pointer for undefined type face (non-standard case)
  integer   iindef

  !> \anchor ientre
  !> if \ref itypfb=ientre: inlet face.
  !> -  Zero-flux condition for pressure and Dirichlet condition for all
  !> other variables. The value of the Dirichlet must be given in
  !> \ref rcodcl "rcodcl"(ifac,ivar,1) for every value of \c ivar, except for
  !> \c ivar = \ref ipr. The other values of \ref rcodcl and
  !> \ref icodcl are filled automatically.
  integer   ientre

  !> \anchor isolib
  !> if \ref itypfb=isolib: free outlet face (or more precisely free inlet/outlet
  !> with forced pressure)
  !>  - The pressure is always treated with a Dirichlet condition, calculated with the constraint
  !> \f$\displaystyle \frac{\partial }{\partial n}\left(\frac{ \partial P}{\partial \tau}\right)=0\f$.
  !> The pressure is set to \f$P_0\f$ at the first \ref isolib face met.
  !> The pressure calibration is always done on a single face, even if there are
  !> several outlets.
  !>  - if the mass flow is coming in, the velocity is set to zero
  !> and a Dirichlet condition for the scalars and the turbulent quantities is used
  !> (or zero-flux condition if no Dirichlet value has been specified).
  !>  - if the mass flow is going out, zero-flux condition are set for the velocity,
  !> the scalars and the turbulent quantities.
  !>  - Nothing is written in \ref icodcl or \ref rcodcl for the pressure or
  !> the velocity. An optional Dirichlet condition can be specified for the scalars
  !> and turbulent quantities.
  !> \remark A standard \ref isolib outlet face amounts to a Dirichlet
  !> condition (\ref icodcl=1) for the pressure, a free outlet condition
  !> (\ref icodcl=9) for the velocity and a Dirichlet condition
  !> (\ref icodcl=1) if the user has specified a Dirichlet value or a zero-flux
  !> condition (\ref icodcl=3) for the other variables.
  integer   isolib

  !> \anchor isymet
  !> if \ref itypfb=isymet: symmetry face (or wall without friction).
  !> - Nothing to be writen in \ref icodcl and  \ref rcodcl.
  integer   isymet

  !> \anchor iparoi
  !> if \ref itypfb=iparoi: smooth solid wall face, impermeable and with friction.
  integer   iparoi

  !> \anchor iparug
  !> if \ref itypfb=iparug: rough solid wall face, impermeable and with friction.
  integer   iparug

  !> if \ref itypfb=iesicf: imposed inlet/outlet for compressible flow (for example, supersonic inlet).
  !>  - A boundary value has to be given for the following quantities:
  !>         - velocity
  !>         - two of the four thermodynamical properties: density, pressure, total energy, temperature
  !>         - all other variables.
  !>  - Homogeneous Neumann boundary condition for the pressure (seen by the reconstruction
  !> gradients and the diffusion operator).
  !>  - Dirichlet condition for the velocity and the total energy.
  !>  - The boundary convective fluxes of momentum and total energy are computed from a Rusanov scheme
  !> for stability reasons. Note that the pressure boundary value is needed to compute those
  !> two fluxes (seen by the pressure gradient of the momentum equation).
  !>  - If the mass flow is coming in, Dirichlet condition for the scalars and the turbulent quantities
  !> is used (or zero-flux condition if no Dirichlet value has been specified).
  !>  - If the mass flow is going out, zero-flux condition are set for the scalars and the turbulent
  !> quantities.
  integer   iesicf

  !> if \ref itypfb=isspcf: supersonic outlet for compressible flow.
  !>  - Nothing needs to be given. The imposed state at the boundary is the upstream state
  !> (values in boundary cells).
  !>  - Homogeneous Neumann boundary condition for the pressure (seen by the reconstruction
  !> gradients and the diffusion operator).
  !>  - Dirichlet (\ref icodcl=1) for the velocity and the total energy.
  !> (pressure boundary value seen by the pressure gradient of the momentum equation).
  !>  - If the mass flow is coming in, Dirichlet condition for the scalars and the turbulent quantities
  !> is used (or zero-flux condition if no Dirichlet value has been specified).
  !>  - If the mass flow is going out, zero-flux condition are set for the scalars and the turbulent
  !> quantities.
  integer   isspcf

  !> if \ref itypfb=isopcf: mixed outlet for compressible flow with a given pressure.
  !>  - Boundary values are obtained by solving a Riemann problem between an inner (values
  !> at boundary cells center) and an outer state. The given pressure is considered as an
  !> outer value.
  !>  - Homogeneous Neumann boundary condition for the pressure (seen by the reconstruction
  !> gradients and the diffusion operator).
  !>  - Dirichlet (\ref icodcl=1) for the velocity and the total energy.
  !>  - Analytical boundary convective fluxes of momentum and total energy are computed.
  !> Note that the pressure boundary value is needed to compute those two fluxes.
  !> (seen by the pressure gradient of the momentum equation).
  !>  - If the mass flow is coming in, Dirichlet condition for the scalars and the turbulent quantities
  !> is used (or zero-flux condition if no Dirichlet value has been specified).
  !>  - If the mass flow is going out, zero-flux condition are set for the scalars and the turbulent
  !> quantities.
  integer   isopcf

  !> if \ref itypfb=iephcf: mixed inlet for compressible flow with given total pressure
  !>                        and total enthalpy (reservoir boundary conditions).
  !>  - Boundary values are obtained by solving a Riemann problem between an inner (values
  !> at boundary cells center) and an outer state.
  !>  - Homogeneous Neumann boundary condition for the pressure (seen by the reconstruction
  !> gradients and the diffusion operator).
  !>  - Dirichlet (\ref icodcl=1) for velocity and total energy.
  !>  - Analytical boundary convective fluxes of momentum and total energy are computed.
  !> Note that the pressure boundary value is needed to compute those two fluxes
  !> (seen by the pressure gradient of the momentum equation).
  !>   - If the mass flow is coming in, Dirichlet condition for the scalars and the turbulent quantities
  !> is used (or zero-flux condition if no Dirichlet value has been specified).
  !>  - If the mass flow is going out, zero-flux condition are set for the scalars and the turbulent
  !> quantities.
  integer   iephcf

  ! TODO : not available yet.
  integer   ieqhcf

  !> \anchor icscpl
  !> code/code coupling condition
  integer   icscpl

  !> \anchor icscpd
  !> code/code coupling condition with decentered flux
  integer   icscpd

  !> \anchor ifrent
  !> if \ref itypfb=ifrent: free entrance based on Bernoulli equation when
  !> the flow is incoming, standard outlet when outgoing
  integer   ifrent

  !> \anchor ifresf
  !> if \ref itypfb=ifresf: free surface for mobile mesh boundary condition
  !>  - Homogeneous Neumann boundary condition for velocity and total energy (seen by the reconstruction
  !> gradients and the diffusion operator).
  !>  - Dirichlet (\ref icodcl=1) for the pressure.
  !> Alse a boundary condition type for mesh velocity in ALE for modelling
  !> free surface (\f$ \vect{u} \cdot \vect{S} = \vect{w} \cdot \vect{S} \f$).
  integer   ifresf

  !> \anchor i_convective_inlet
  !> if \ref itypfb=i_convective_inlet: inlet face where the total mass flux is
  !>                                    prescribed.
  !> -  Zero-flux condition for pressure and Dirichlet condition for all
  !> other variables. The value of the Dirichlet must be given in
  !> \ref rcodcl "rcodcl"(ifac,ivar,1) for every value of \c ivar, except for
  !> \c ivar = \ref numvar::ipr "ipr". The other values of \ref rcodcl and
  !> \ref icodcl are filled automatically.
  !> The diffusive flux is CANCELLED (therefore the total mass flux is due to
  !> convection only).
  integer   i_convective_inlet


  parameter(iindef=1, ientre=2, isolib=3, isymet=4, iparoi=5,       &
            iparug=6, iesicf=7, isspcf=8, isopcf=9, iephcf=10,      &
            ieqhcf=11, icscpl=12, icscpd=13, ifrent=14, ifresf=15,  &
            i_convective_inlet=16)

  !> maximal number of valuators for Navier-Stokes
  integer    nestmx
  parameter (nestmx=4)

  !> Error estimator for Navier-Stokes.
  !> iest = iespre: prediction, (default name: EsPre).
  !> After the velocity prediction step (yielding \f$\vect{u}^*\f$), the
  !> estimator \f$\eta^{\,pred}_{\,i,k}(\vect{u}^*)\f$, local variable calculated
  !> at every cell \f$ \Omega_i \f$, is created from
  !> \f$\vect{\mathcal R}^{\,pred}(\vect{u}^*)\f$,
  !> which represents the residual of the equation solved during this step:
  !> \f$\vect{u}\f$ and \f$ P \f$:
  !> \f{eqnarray*}{
  !>   \vect{\mathcal R}^{\,pred}(\vect{u}^*)
  !>       & = & \rho^n \dfrac{\vect{u}^*-\vect{u}^n}{\Delta t}
  !>           + \rho^n \vect{u}^n \cdot \gradt (\vect{u}^*)
  !>           - \divv \left((\mu+\mu_t)^n \gradt (\vect{u}^*) \right)
  !>           + \grad(P^n)
  !>   \\  & - & \text{rest of the right-hand member }
  !>            (\vect{u}^n, P^n, \text{other variables}^n)
  !> \f}
  !>  - By definition:
  !> \f$ \eta^{\,pred}_{\,i,k}(\vect{u}^*)= {|\Omega_i|}^{\,(k-2)/2}\ ||\vect{\mathcal R}^{\,pred}(\vect{u}^*)||
  !> _{{IL}^{2}(\Omega_i)} \f$
  !>  - The first family, k=1, suppresses the
  !> volume \f$ |\Omega_i| \f$ which intrinsicly appears  with the norm
  !> \f$ {IL}^{2}(\Omega_i) \f$.
  !>  - The second family, k=2, exactly represents the norm
  !> \f$ {IL}^{2}(\Omega_i) \f$. The size of the cell therefore
  !> appears in its calculation and induces a weighting effect.
  !>  - \f$ \eta^{\,pred}_{\,i,k}(\vect{u}^*)\f$  is ideally equal to zero when the
  !> reconstruction methods are perfect and the associated system is solved exactly.
  integer   iespre

  !> Error estimator for Navier-Stokes.
  !> iest = iesder: drift  (default name: EsDer).
  !> The estimator \f$\eta^{\,der}_{\,i,k}(\vect{u}^{\,n+1})\f$ is based on the
  !> following quantity (intrinsic to the code):
  !> \f{eqnarray*}{
  !> \eta^{\,der}_{\,i,k}(\vect{u}^{\,n+1})
  !>    &=& {|\Omega_i|}^{(k-2)/2}
  !>       || \divs (\text{corrected mass flow after the pressure step})
  !>       - \Gamma||_{{L}^{2}(\Omega_i)}
  !> \\ &=& {|\Omega_i|}^{(1-k)/2}
  !>      | \divs (\text{corrected mass flow after the pressure step})- \Gamma|
  !> \f}
  !>  - Ideally, it is equal to zero when the Poisson equation related to the pressure is
  !> solved exactly.
  integer   iesder

  !> Error estimator for Navier-Stokes.
  !> iest = iescor: correction, (default name: EsCor).
  !> The estimator \f$ \eta^{\,corr}_{\,i,k}(\vect{u}^{\,n+1})\f$ comes directly
  !> from the mass flow calculated with the updated velocity field:
  !> \f{eqnarray*}{
  !> \eta^{\,corr}_{\,i,k}(\vect{u}^{\,n+1})=
  !> |\Omega_i|^{\,\delta_{\,2,k}}\ |div (\rho^n \vect{u}^{n+1}) - \Gamma|
  !> \f}
  !> - The velocities \f$\vect{u}^{n+1}\f$ are taken at the cell centers,
  !> the divergence is calculated after projection on the faces.
  !> \f$ \,\delta_{\,2,k}\f$ represents the Kronecker symbol.
  !> - The first family, k=1, is the absolute raw value of the divergence of the mass flow
  !< minus the mass source term.
  !> The second family, $k=2$, represents a physical property and allows to evaluate
  !> the difference in \f$kg.s^{\,-1}\f$.
  !> - Ideally, it is equal to zero when the Poisson equation is solved exactly and
  !> the projection from the mass flux at the faces to the velocity at the cell
  !> centers is made in a set of  functions with null divergence.
  integer   iescor

  !> Error estimator for Navier-Stokes. iest = iestot: total, (default name: EsTot).
  !> The estimator \f$ \eta^{\,tot}_{\,i,k}(\vect{u}^{\,n+1})\f$, local variable
  !> calculated at every cell \f$\Omega_i\f$, is based on the quantity
  !> \f$\vect{\mathcal R}^{\,tot}(\vect{u}^{\,n+1})\f$, which represents the
  !> residual of the equation using the updated values of
  !> \f$\vect{u}\f$ and \f$P\f$:
  !> \f{eqnarray*}{
  !>   \vect{\mathcal R}^{\,pred}(\vect{u}^*)
  !>       & = & \rho^n \dfrac{\vect{u}^*-\vect{u}^n}{\Delta t}
  !>           + \rho^n \vect{u}^n \cdot \gradt (\vect{u}^*)
  !>           - \divv \left((\mu+\mu_t)^n \gradt (\vect{u}^*) \right)
  !>           + \grad(P^n)
  !>   \\  & - & \text{rest of the right-hand member }
  !>            (\vect{u}^n, P^n, \text{other variables}^n)
  !> \f}
  !> - By definition:
  !> \f$ \eta^{\,tot}_{\,i,k}(\vect{u}^{\,n+1})= {|\Omega_i|}^{\,(k-2)/2}\ ||\vect{\mathcal R}^{\,tot}(\vect{u}^{\,n+1})||
  !> _{{I\hspace{-.25em}L}^{2}(\Omega_i)} \f$
  !> - The mass flux in the convective term is recalculated from \f$\vect{u}^{n+1}\f$
  !> expressed at the cell centers (and not taken from the updated mass flow at the
  !> faces).
  !> - As for the prediction estimator:
  !>   - The first family, k=1, suppresses the
  !> volume \f$ |\Omega_i| \f$ which intrinsicly appears  with the norm
  !> \f$ {IL}^{2}(\Omega_i) \f$.
  !>   - The second family, k=2, exactly represents the norm
  !> \f$ {IL}^{2}(\Omega_i) \f$. The size of the cell therefore
  !> appears in its calculation and induces a weighting effect.
  integer   iestot
  parameter (iespre=1, iesder=2, iescor=3, iestot=4)

  ! conditions aux limites possibles pour la vitesse de maillage en ale

  !> \anchor ibfixe
  !> boundary condition type for mesh velocity in ALE: fixed wall
  integer   ibfixe

  !> \anchor igliss
  !> boundary condition type for mesh velocity in ALE: sliding wall
  integer   igliss

  !> \anchor ivimpo
  !> boundary condition type for mesh velocity in ALE: imposed velocity.
  !> - In the case where all the nodes of a face have a imposed displacement,
  !> it is not necessary to fill the tables with boundary conditions
  !> mesh velocity for this face, they will be erased. In the other case,
  !> the value of the Dirichlet must be given in \ref rcodcl "rcodcl"(ifac,ivar,1)
  !> for every value of \c ivar (\ref iuma, \ref ivma and \ref iwma).
  !> The other boxes of \ref rcodcl and \ref icodcl are completed automatically.
  !> The tangential mesh velocity is taken like a tape speed under the
  !> boundary conditions of wall for the fluid, except if wall fluid velocity
  !> was specified by the user in the interface or \ref cs_user_boundary_conditions
  !> (in which case it is this speed which is considered).
  integer   ivimpo

  parameter(ibfixe=1, igliss=2, ivimpo=3)

  !> maximum number of structures in ALE
  integer    nstrmx
  parameter (nstrmx=200)

  !=============================================================================

  !> \}

  ! Temporary former Fortan subroutine names so that the user soubroutines
  ! do not compile if the deprecated name (now used for the C routine) is used.
  integer cs_user_boundary_conditions, cs_user_parameters,     &
          cs_user_initialization, cs_user_physical_properties, &
          cs_user_extra_operations, ushist, cs_f_user_head_losses

end module paramx
