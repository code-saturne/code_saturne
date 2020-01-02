/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*-----------------------------------------------------------------------------*/

/*!

\page cs_user_mass_source_terms Examples of data settings for mass source terms (cs_user_mass_source_terms.f90)



The subroutine \ref cs_user_mass_source_terms is called at three different stages in the code <tt> (iappel = 1, 2 or 3):</tt>
  - \c iappel = 1: Calculation of the number of cells where a mass source is imposed: \c ncesmp. Called once at the beginning of the calculation.
  - \c iappel = 2: Identification of the cells where a mass source term is imposed: array \c icesmp(ncesmp). Called once at the beginning of the calculation.
  - \c iappel = 3: Calculation of the values of the mass source terms. Called at each time step.


The equation for mass conservation becomes: \f$ \der{\rho}{t}+\divs{(\rho \vect{u})}=\Gamma \f$

The equation for a variable \f$ \varia \f$ becomes:\f$ \frac{\Delta(\rho\varia)}{\Delta t} = ... + \Gamma(\varia^{in} - \varia) \f$  discretized as \f$ \rho \dfrac{\varia^{(n+1)} - \varia^{(n)}}{\Delta t} = ... + \Gamma(\varia^{in} - \varia^{(n+1)}) \f$

\f$ \varia^{in} \f$ is the value of \f$ \varia \f$ associated to the injected mass.

Two options are available:
  - the mass flux is injected with the local value of variable \f$ \varia \f$: \f$ \varia^{in} = \varia^{(n+1)} \f$ (the equation for \f$ \varia \f$ is therefore not modified)
  - the mass flux is injected with a specific value for \f$ \varia \f$: \f$ \varia^{in} \f$ is specified by the user


\section var_user Variables to be specified by the user

 - \c ncesmp: number of cells where a mass source term is imposed

 - \c icetsm(ieltsm): identification of the cells where a mass source
                  term is imposed.
                  For each cell where a mass source term is imposed
                  (\c ielstm in \c [1;ncesmp]), \c icetsm(ieltsm) is the
                  global index number of the corresponding cell
                  (\c icestm(ieltsm) in \c [1;ncel])

 - \c smacel(ieltsm,ipr): value of the injection mass rate gamma (in \f$ kg \cdot m^3 \cdot s^{-1}\f$)
                             in the \c ieltsm cell with mass source term

 - \c itypsm(ieltsm,ivar): type of treatment for variable ivar in the
                       \c ieltsm cell with mass source term.
                     * <tt>itypsm = 0 </tt>--> injection of \c ivar at local value
                     * <tt>itypsm = 1 </tt>--> injection of \c ivar at user
                                      specified value

 - \c smacel(ieltsm,ivar): specified value for variable ivar associated
                       to the injected mass in the \c ieltsm cell with
                       a mass source term
                                  except for <tt> ivar=ipr </tt>

\remark
 - if \c itypsm(ieltsm,ivar)=0, \c smacel(ieltsm,ivar) is not used

 - if \c smacel(ieltsm,ipr)<0, mass is removed from the system,
     therefore \c Code_Saturne automatically considers \f$ \varia^{in}=\varia^{(n+1)}\f$,
     whatever the values of \c itypsm or \c smacel specified by the user

 - if a value \c ivar is not linked for a mass source
     term is imposed, no source term will be taken into account.

 - if a scalar doesn't evolve following the standard equation
     \f$ \dfrac{\partial{(\rho \varia)}}{\partial{dt}} + \divs{(\rho \vect{u} \varia)} = ...\f$
     (alternate convective field for instance), the source term
     set by this routine will not be correct (except in case of
     injection at the local value of the variable). The proper source
     term should be added directly in \ref ustssc.

\remark <b> Idenfication of the cells:</b> \n
The selection of cells where to apply the source term is based on a \ref getcel command. For more info on the syntax of the \ref getcel command, refer to the user manual or to the comments on the similar command \ref getfbr in the routine \ref cs_user_boundary_conditions.

\section loc_var_ms Local variables

\snippet cs_user_mass_source_terms.f90 loc_var

\section init_and_finit Initialization and finalization

The following initialization block needs to be added for the following examples:
\snippet cs_user_mass_source_terms.f90 allocate

At the end of the subroutine, it is recommended to deallocate the work array:
\snippet cs_user_mass_source_terms.f90 deallocate

In theory Fortran 95 deallocates locally-allocated arrays automatically, but deallocating arrays in a symetric manner to their alloacation is good pratice, and avoids using a different logic C and Fortran.

\section one_or_two First or second call

 - First call: <tt> iappel = 1</tt>: \c ncesmp: calculation of the number of cells with mass source term

 - Second call (if \c ncesmp>0): <tt>iappel = 2</tt>: \c icetsm: index number of cells with mass source terms

\warning
  - Do not use \c smacel in this section (it is set on the third call, <tt> iappel = 3 </tt>)
  - Do not use \c icetsm in this section on the first call (\c iappel=1)

This section <tt> (iappel = 1 or 2) </tt> is only accessed at the beginning of a calculation. Should the localization of the mass source terms evolve in time, the user must identify at the beginning all cells that can potentially become a mass source term.

\snippet cs_user_mass_source_terms.f90 one_or_two


\subsection example1_1 Example 1
No mass source term (default)

\snippet cs_user_mass_source_terms.f90 example_1_1

\subsection example1_2 Example 2
Mass source term in the cells that have boundary face of color 3 and the cells with a coordinate X between 2.5 and 5.

In this test in two parts, one must pay attention not to count
the cells twice (a cell with a boundary face of color 3 can
also have a coordinate X between 2.5 and 5).
One should also pay attention that, on the first call, the
array \c icetsm doesn't exist yet. It mustn't be used outside
of tests (\c iappel.eq.2).

\snippet  cs_user_mass_source_terms.f90 example_1_2

\subsection genric_sub Generic subsection: do not modify

  - For \c iappel = 1: specification of \c ncesmp. This block is valid for both examples.

\snippet  cs_user_mass_source_terms.f90 generic_sub


\section three Third call (for ncesmp > 0)

  \c iappel = 3: -  \c itypsm: type of mass source term
                 -  \c smacel : mass source term

\remark If \c itypsm(ieltsm,var) is set to 1, \c smaccel(ieltsm,var) must be set.

\subsection example2_1 Example 1

Simulation of an inlet condition by mass source terms and printing of the total mass rate.

\snippet  cs_user_mass_source_terms.f90 example_2_1

Calculation of the inlet conditions for k and epsilon with standard laws in a circular pipe

\snippet cs_user_mass_source_terms.f90 inlet_cal

\subsection example2_2 Example 2
Simulation of a suction (by a pump for instance) with a total rate of \f$ 80 000 \: kg \cdot s^{-1}\f$. The suction rate is supposed to be uniformly distributed on all the cells selected above.

Calculation of the total volume of the area where the mass source term is imposed (the case of parallel computing is taken into account with the call to \c parsom).

\snippet cs_user_mass_source_terms.f90  calcul_total

The mass suction rate is \f$ \Gamma = - \dfrac{80000}{v_{tot}} \f$ (in \f$ kg \cdot m^{-3} \cdot s^{-1}\f$). It is set below, with a test for cases where \f$ v_{tot} = 0 \f$. The total mass rate is calculated for verification.

\snippet cs_user_mass_source_terms.f90  mass_suction

\subsection end3 End third call

\snippet cs_user_mass_source_terms.f90 end_call_3

\subsection format Formats

\snippet cs_user_mass_source_terms.f90 format

*/
