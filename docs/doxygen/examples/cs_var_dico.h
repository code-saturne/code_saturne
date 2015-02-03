/*============================================================================
 * Code_Saturne documentation page
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2015 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/*!
  \page cs_var_dico Fortran-C naming reference

  \section intro

  This page provides cross-reference tables containing the names of Fortran
  variables and their C counterparts.

  Note: Variables with the same name in Fortran and in C are not listed.

  \section namingcontent Contents

  Cross-reference tables are available for the following categories of
  variables:

   - \subpage mesh
   - \subpage field
   - \subpage local

*/
// _____________________________________________________________________________
/*!

  \page mesh Mesh and mesh quantities variables

  \section mesh_vars Mesh variables

  These variables are defined in the files \ref mesh.f90 and \ref cs_mesh.h.

  - The Fortran variables are global in the code, accessible directly by their
    name.
  - Members of the global C structure \c cs_glob_mesh are accessed as: \n
    <tt>cs_glob_mesh->name</tt>, \n
    \e i.e. adding <tt> cs_glob_mesh-> </tt> in front of the name of the
    variable.

  Fortran code       | C code                                | Description
  ------------------ | ------------------------------------- | ------------
  <tt> \ref ndim     | cs_glob_mesh->dim                     | Space dimension
  \ref ncelet        | cs_glob_mesh->n_cells_with_ghosts     | Total number of cells on the local rank \n (n_cells + n_ghost_cells)
  \ref ncel          | cs_glob_mesh->n_cells                 | Number of cells
  \ref nfac          | cs_glob_mesh->n_i_faces               | Number of interior faces
  \ref nfabor        | cs_glob_mesh->n_b_faces               | Number of boundary faces
  \ref nnod          | cs_glob_mesh->n_vertices              | Number of vertices
                     | cs_glob_mesh->n_b_cells               | Number of boundary cells
  \ref lndfac        | cs_glob_mesh->i_face_vtx_connect_size | Size of the connectivity \n interior faces -> vertices
  \ref lndfbr        | cs_glob_mesh->b_face_vtx_connect_size | Size of the connectivity \n boundary faces -> vertices
  \ref nfml          | cs_glob_mesh->n_families              | Number of families
  \ref ifacel        | cs_glob_mesh->i_face_cells            | Interior faces -> cells connectivity
  \ref ifabor        | cs_glob_mesh->b_face_cells            | Boundary faces -> cells connectivity
  \ref ipnfac        | cs_glob_mesh->i_face_vtx_idx          | Interior faces -> vertices index
  \ref nodfac        | cs_glob_mesh->i_face_vtx_lst          | Interior faces -> vertices connectivity
  \ref ipnfbr        | cs_glob_mesh->b_face_vtx_idx          | Boundary faces -> vertices index
  \ref nodfbr        | cs_glob_mesh->b_face_vtx_lst          | Boundary faces -> vertices connectivity
  \ref ifmfbr        | cs_glob_mesh->b_face_family           | Boundary face family
  \ref ifmcel        | cs_glob_mesh->cell_family             | Cell family
                     | cs_glob_mesh->b_cells                 | Boundary cell list
  \ref xyznod        | cs_glob_mesh->vtx_coord               | Vertex coordinates </tt>

  \section mesh_q_vars Mesh quantities variables

  These variables are defined in the files \ref mesh.f90 and
  \ref cs_mesh_quantities.h.

  - The Fortran variables are global in the code, accessible directly by their
    name.
  - Members of the global C structure \c cs_glob_mesh_quantities are accessed
    as: \n
    <tt>cs_glob_mesh_quantities->name</tt>, \n
    \e i.e. adding <tt> cs_glob_mesh_quantities-> </tt> in front of the name of
    the variable.

  Fortran code       | C code                                  | Description
  ------------------ | --------------------------------------- |-------------
  <tt> \ref isympa   |  cs_glob_mesh_quantities->b_sym_flag    | Symmetry flag for boundary faces
  \ref xyzcen        |  cs_glob_mesh_quantities->cell_cen      | Cell center coordinates
  \ref surfac        |  cs_glob_mesh_quantities->i_face_normal | Surface normal of interior faces
  \ref surfbo        |  cs_glob_mesh_quantities->b_face_normal | Surface normal of border faces
  \ref cdgfac        |  cs_glob_mesh_quantities->i_face_cog    | Center of gravity of interior faces
  \ref cdgfbo        |  cs_glob_mesh_quantities->b_face_cog    | Center of gravity of border faces
  \ref volume        |  cs_glob_mesh_quantities->cell_vol      | Cell volume
  \ref surfan        |  cs_glob_mesh_quantities->i_face_surf   | Surface of interior faces
  \ref surfbn        |  cs_glob_mesh_quantities->b_face_surf   | Surface of boundary faces
  \ref dist          |  cs_glob_mesh_quantities->i_dist        | Distance between the cell center and \n the center of gravity of interior faces
  \ref distb         |  cs_glob_mesh_quantities->b_dist        | Distance between the cell center and \n the center of gravity of border faces
  \ref pond          |  cs_glob_mesh_quantities->weight        | Interior faces weighting factor </tt>


*/
// _____________________________________________________________________________
/*!

  \page field Variables, properties and fields (deprecated \c rtp, \c propce)

  - In Fortran:
    - \ref vars "Variables" are accessed through the deprecated array \c rtp as: \n
      <tt>rtp(iel,index)</tt>.
    - \ref props "Properties" are accessed through the deprecated array \c propce as: \n
      <tt>propce(iel,index)</tt>.
    - The scalar indexes are accessed as \c isca(j), and the scalar values as
      follows: \n
      <tt>rtp(iel,isca(iscalt))</tt>.
    - Now, both variables and properties can be accessed via the
      \ref field.f90 "cs_field" API, as in the following examples: \n
      <tt>call \ref field_get_val_s(ivarfl(ipr), cvar_pr) \n cvar_pr(iel)</tt>, \n\n
      <tt>call \ref field_get_val_v(ivarfl(iu), cvar_vel) \n cvar_vel(isou,iel)</tt>, \n\n
      <tt>call \ref field_get_val_s(iprpfl(icp), cpro_cp) \n cpro_cp(iel)</tt>, \n\n
      where \ref ipr, \ref iu are variable indexes and \ref icp is a property index.
  - In C:
    - Both variables and properties are accessed as: \n
      <tt>CS_F_(name)->val[cell_id]</tt>, \n
      where the \c name is defined in \ref cs_field_pointer.h. \n
    - Indexed variables (such as user scalars) and indexed properties
      are accessed as: \n
      <tt>CS_FI_(name,ii-1)->val[cell_id]</tt>.

  \remark Note that indexes in C begin at 0, while indexes in Fortran begin at 1.
          Thus, Fortran and C loop counters are related in the following as:\n
          <tt>cell_id = iel-1</tt>.

  Cross-reference tables are available for the variables and properties of the
  standard solver and the specific physics features:
  - \ref vars
  - \ref props
  - \ref part
    - \ref atmo
    - \ref comb
    - \ref cfbl
    - \ref elec
    - \ref cogz
    - \ref rayt

  \section vars Variables

  The Fortran variables indexes are defined in the files \ref numvar.f90 (with
  the exception of \c ihm and \c iscal, which are respectively defined in \ref
  ppincl.f90 and \ref optcal.f90) and the C variables names are defined in
  \ref cs_field_pointer.h. \n Note that \c dt is just an \c allocatable array in
  Fortran while it is mapped as a field in C.

  deprecated Fortran code         | Fortran code                                     | C code                       | Description
  ------------------------------- | ------------------------------------------------ | ---------------------------- | ------------
  <tt> dt                         | <tt> dt                                          | CS_F_(dt)->val[cell_id]      | Local time step
  rtp(iel,\ref ipr)               | call field_get_val_s(ivarfl(ipr), cvar_pr)       | CS_F_(p)->val[cell_id]       | Pressure
  <EM>(WARNING: deprecated)</EM> \n rtp(iel,\ref iu) \n rtp(iel,\ref iv) \n rtp(iel,\ref iw) | call field_get_val_v(ivarfl(iu), cvar_vel) \n cvar_vel(1, iel) \n cvar_vel(2, iel) \n cvar_vel(3, iel) | \n CS_F_(u)->val[cell_id][0] \n CS_F_(u)->val[cell_id][1] \n CS_F_(u)->val[cell_id][2] | Velocity
  rtp(iel,\ref ivoidf)            | call field_get_val_s(ivarfl(ivoidf), cvar_voidf) | CS_F_(void_f)->val[cell_id]  | Void fraction for cavitation modelling
  rtp(iel,\ref ik)                | call field_get_val_s(ivarfl(ik  ), cvar_k  )     | CS_F_(k)->val[cell_id]       | Turbulent kinetic energy \f$ k \f$
  rtp(iel,\ref iep)               | call field_get_val_s(ivarfl(iep ), cvar_eps)     | CS_F_(eps)->val[cell_id]     | Turbulent dissipation \f$ \varepsilon \f$
  rtp(iel,\ref ir11)              | call field_get_val_s(ivarfl(ir11), cvar_r11)     | CS_F_(r11)->val[cell_id]     | Reynolds stress component \f$ R_{xx} \f$
  rtp(iel,\ref ir22)              | call field_get_val_s(ivarfl(ir22), cvar_r22)     | CS_F_(r22)->val[cell_id]     | Reynolds stress component \f$ R_{yy} \f$
  rtp(iel,\ref ir33)              | call field_get_val_s(ivarfl(ir33), cvar_r33)     | CS_F_(r33)->val[cell_id]     | Reynolds stress component \f$ R_{zz} \f$
  rtp(iel,\ref ir12)              | call field_get_val_s(ivarfl(ir12), cvar_r12)     | CS_F_(r12)->val[cell_id]     | Reynolds stress component \f$ R_{xy} \f$
  rtp(iel,\ref ir23)              | call field_get_val_s(ivarfl(ir23), cvar_r23)     | CS_F_(r23)->val[cell_id]     | Reynolds stress component \f$ R_{yz} \f$
  rtp(iel,\ref ir13)              | call field_get_val_s(ivarfl(ir13), cvar_r13)     | CS_F_(r13)->val[cell_id]     | Reynolds stress component \f$ R_{xz} \f$
  rtp(iel,\ref iphi)              | call field_get_val_s(ivarfl(iphi), cvar_phi)     | CS_F_(phi)->val[cell_id]     | \f$ \phi \f$ for \f$ \phi-f_b \f$ model
  rtp(iel,\ref ifb)               | call field_get_val_s(ivarfl(ifb ), cvar_fb )     | CS_F_(f_bar)->val[cell_id]   | \f$ f_b \f$ for \f$ \phi-f_b \f$ model
  rtp(iel,\ref ial)               | call field_get_val_s(ivarfl(ial ), cvar_al )     | CS_F_(alpha)->val[cell_id]   | \f$ \alpha \f$ for \f$ Bl-v^2-k \f$ \n or EBRSM model
  rtp(iel,\ref iomg)              | call field_get_val_s(ivarfl(iomg), cvar_omg)     | CS_F_(omg)->val[cell_id]     | \f$ \omega \f$ for \f$ k-\omega \f$ SST model
  rtp(iel,\ref inusa)             | call field_get_val_s(ivarfl(inusa), cvar_nusa)   | CS_F_(nusa)->val[cell_id]    | \f$ \widetilde{\nu}_T \f$ for Spalart-Allmaras
  rtp(iel,\ref iuma) \n rtp(iel,\ref ivma) \n rtp(iel,\ref iwma) | call field_get_val_v(ivarfl(iuma), cvar_mesh_v)  | CS_F_(mesh_u)->val[3*cell_id] \n CS_F_(mesh_u)->val[3*cell_id+1] \n CS_F_(mesh_u)->val[3*cell_id+2] | Mesh velocity
  rtp(iel,\ref isca(\ref ihm))    | call field_get_val_s(ivarfl(isca(ihm)), cvar_hm)       | CS_F_(h)->val[cell_id]       | Enthalpy
  rtp(iel,\ref isca(\ref iscalt)) | call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt) | CS_F_(t)->val[cell_id]       | Temperature </tt>


  \section props Properties

  These properties are defined in the files \ref numvar.f90 and
  \ref cs_field_pointer.h.

  Deprecated Fortran code    | Deprecated Fortran code                                  |C code                       | Description
  -------------------------- | -------------------------------------------------------- |---------------------------- | ------------
  <tt> propce(iel,\ref irom) | call field_get_val_s(iprpfl(irom  ), cpro_rho  )         |CS_F_(rho)->val[cell_id]     | Density at the current time step
  propce(iel,\ref iroma)     | call field_get_val_s(iprpfl(iroma ), cpro_rhoa )         |CS_F_(rho)->val_pre[cell_id] | Density at the previous time step
  propce(iel,\ref iromaa)    | call field_get_val_s(iprpfl(iromaa), cpro_rhoaa)         |<em>not yet implemented</em> | Density at the second previous time
  propce(iel,\ref iviscl)    | call field_get_val_s(iprpfl(iviscl), cpro_viscl)         |CS_F_(mu)->val[cell_id]      | Molecular viscosity
  propce(iel,\ref ivisct)    | call field_get_val_s(iprpfl(ivisct), cpro_visct)         |CS_F_(mu_t)->val[cell_id]    | Turbulent dynamic viscosity
  propce(iel,\ref ivisla)    | call field_get_val_s(iprpfl(iromaa), cpro_romaa)         |CS_F_(mu)->val_pre[cell_id]  | Dynamic molecular viscosity (in kg/(m.s)) \n at the previous time-step
  propce(iel,\ref ivista)    | call field_get_val_s(iprpfl(ivista), cpro_viscta)        |CS_F_(mu_t)->val_pre[cell_id] | Dynamic turbulent viscosity \n at the previous time-step
  propce(iel,\ref icp)       | call field_get_val_s(iprpfl(icp   ), cpro_cp   )         |CS_F_(cp)->val[cell_id]      | Specific heat
  propce(iel,\ref icpa)      | call field_get_val_s(iprpfl(icpa  ), cpro_cpa  )         |CS_F_(cp)->val_pre[cell_id]  | specific heat at the previous time-step
  propce(iel,\ref icrom)     | call field_get_val_s(iprpfl(icrom ), cpro_crom )         |CS_F_(rho)->val[cell_id]     | Density (at cells)
  propce(iel,\ref ibrom)     | call field_get_val_s(iprpfl(ibrom ), bpro_rho  )         |CS_F_(rho_b)->val[cell_id]   | Density (at boundary faces)
  propce(iel,\ref ismago)    | call field_get_val_s(iprpfl(ismago), cpro_smago)         |<em>not yet implemented</em> | Field id of the anisotropic turbulent viscosity
  propce(iel,\ref icour)     | call field_get_val_s(iprpfl(icour ), cpro_cour )         |<em>not yet implemented</em> | Courant number
  propce(iel,\ref ifour)     | call field_get_val_s(iprpfl(ifour ), cpro_four )         |<em>not yet implemented</em> | Fourier number
  propce(iel,\ref iprtot)    | call field_get_val_s(iprpfl(iprtot), cpro_prtot)         |<em>not yet implemented</em> | Total pressure at cell centres
  propce(iel,\ref ivisma(1)) \n propce(iel,\ref ivisma(2)) \n propce(iel,\ref ivisma(3))| <em>not yet implemented</em> | Mesh velocity viscosity for the ALE module
  propce(iel,\ref itsrho)    | call field_get_val_s(iprpfl(itsrho), cpro_tsrho )        |<em>not yet implemented</em> | Global dilatation source terms
  propce(iel,\ref ibeta)     | call field_get_val_s(iprpfl(ibeta ), cpro_beta  )        |<em>not yet implemented</em> | Thermal expansion coefficient
                             | call field_get_val_s(\ref ipri, porosi) \n porosi(iel)   | CS_F_(poro)->val[cell_id]       | Porosity
                             | call field_get_val_v(\ref iprf, porosf) \n porosf(ii,iel)| CS_F_(t_poro)->val[cell_id][ii] | Tensorial porosity
                             | call field_get_val_v(\ref ifrbr, forbr) \n forbr(ii,iel) | <em>not yet implemented</em>    | Field id of the stresses at boundary
                             | call field_get_val_s(\ref iylbr, yplbr) \n yplbr(iel)    | <em>not yet implemented</em>    | Field id of \f$y^+\f$ at boundary
                             | call field_get_val_v(\ref idten, dttens) \n dttens(ii,iel) | <em>not yet implemented</em>    | Field id for the dttens tensor</tt>


  \section part Specific physics

  \subsection atmo Atmospheric

  Defined in \ref optcal.f90, \ref atincl.f90, \ref atvarp.f90 and
  \ref cs_field_pointer.h.

  Fortran code                             | C code                                 | Description
  ---------------------------------------- | -------------------------------------- | ------------
  <tt> rtp(iel,\ref isca(\ref iscalt))     | CS_F_(pot_t)->val[cell_id]             | Potential temperature
  rtp(iel,\ref isca(\ref itotwt))          | CS_F_(totwt)->val[cell_id]             | Total water content
  rtp(iel,\ref isca(\ref intdrp))          | CS_F_(ntdrp)->val[cell_id]             | Total number of droplets
  rtp(iel,\ref isca(\ref isca_chem(iesp))) | CS_FI_(chemistry,iesp-1)->val[cell_id] | Chemistry species (indexed) </tt>


  \subsection comb Coal combustion

  Defined in \ref ppincl.f90, \ref ppcpfu.f90 and \ref cs_field_pointer.h.

  Fortran code                            | C code                           | Description
  --------------------------------------- | -------------------------------- | ------------
  <tt> rtp(iel,\ref isca(\ref inp(iesp))) | CS_FI_(np,iesp-1)->val[cell_id]  | Particles per kg for coal class
  rtp(iel,\ref isca(\ref ixch(iesp)))     | CS_FI_(xch,iesp-1)->val[cell_id] | Reactive coal mass fraction for coal class
  rtp(iel,\ref isca(\ref ixck(iesp)))     | CS_FI_(xck,iesp-1)->val[cell_id] | Coke mass fraction for coal class
  rtp(iel,\ref isca(\ref ixwt(iesp)))     | CS_FI_(xwt,iesp-1)->val[cell_id] | Water mass fraction for coal class
  rtp(iel,\ref isca(\ref ih2(iesp)))      | CS_FI_(h2,iesp-1)->val[cell_id]  | Mass enthalpy for coal class (permeatic case)
  rtp(iel,\ref isca(\ref if1m(iesp)))     | CS_FI_(f1m,iesp-1)->val[cell_id] | Mean value light volatiles for coal class
  rtp(iel,\ref isca(\ref if2m(iesp)))     | CS_FI_(f2m,iesp-1)->val[cell_id] | Mean value heavy volatiles for coal class
  rtp(iel,\ref isca(\ref if4m))           | CS_F_(f4m)->val[cell_id]         | -
  rtp(iel,\ref isca(\ref if5m))           | CS_F_(f5m)->val[cell_id]         | -
  rtp(iel,\ref isca(\ref if6m))           | CS_F_(f6m)->val[cell_id]         | -
  rtp(iel,\ref isca(\ref if7m))           | CS_F_(f7m)->val[cell_id]         | -
  rtp(iel,\ref isca(\ref if8m))           | CS_F_(f8m)->val[cell_id]         | -
  rtp(iel,\ref isca(\ref if9m))           | CS_F_(f9m)->val[cell_id]         | -
  rtp(iel,\ref isca(\ref ifvp2m))         | CS_F_(fvp2m)->val[cell_id]       | -
  rtp(iel,\ref isca(\ref iyco2))          | CS_F_(yco2)->val[cell_id]        | CO2 fraction
  rtp(iel,\ref isca(\ref iyhcn))          | CS_F_(yhcn)->val[cell_id]        | HCN fraction
  rtp(iel,\ref isca(\ref iyno))           | CS_F_(yno)->val[cell_id]         | NO fraction
  rtp(iel,\ref isca(\ref iynh3))          | CS_F_(ynh3)->val[cell_id]        | NH3 enthalpy
  rtp(iel,\ref isca(\ref ihox))           | CS_F_(hox)->val[cell_id]         | Ox enthalpy </tt>


  \subsection cfbl Compressible

  Defined in \ref ppincl.f90 and \ref cs_field_pointer.h.

  Fortran code                         | C code                        | Description
  ------------------------------------ | ----------------------------- | ------------
  <tt> rtp(iel,\ref isca(\ref ienerg)) | CS_F_(energy)->val[cell_id]   | Total energy
  rtp(iel,\ref isca(\ref itempk))      | CS_F_(t_kelvin)->val[cell_id] | Temperature, in Kelvin </tt>


  \subsection elec Electric arcs

  Defined in \ref ppincl.f90 and \ref cs_field_pointer.h.

  Fortran code                          | C code                             | Description
  ------------------------------------- | ---------------------------------- | ------------
  <tt> rtp(iel,\ref isca(\ref ipotr))   | CS_F_(potr)->val[cell_id]          | Electric potential, real part
  rtp(iel,\ref isca(\ref ipoti))        | CS_F_(poti)->val[cell_id]          | Electric potential, imaginary part
  rtp(iel,\ref isca(\ref ipotva(1))) \n rtp(iel,isca(ipotva(2))) \n rtp(iel,isca(ipotva(3))) | CS_F_(potva)->val[cell_id][0] \n CS_F_(potva)->val[cell_id][1] \n CS_F_(potva)->val[cell_id][2] | Vector potential
  rtp(iel,\ref isca(\ref iycoel(iesp))) | CS_FI_(ycoel,iesp-1)->val[cell_id] | Constituent mass fraction </tt>


  \subsection cogz Gas combustion

  Defined in \ref ppincl.f90 and \ref cs_field_pointer.h.

  Fortran code                      | C code                     | Description
  --------------------------------- | -------------------------- | ------------
  <tt> rtp(iel,\ref isca(\ref ifm)) | CS_F_(fm)->val[cell_id]    | Mixture fraction
  rtp(iel,\ref isca(\ref ifp2m))    | CS_F_(fp2m)->val[cell_id]  | Mixture fraction variance
  rtp(iel,\ref isca(\ref ifsm))     | CS_F_(fsm)->val[cell_id]   | Soot mass fraction
  rtp(iel,\ref isca(\ref inpm))     | CS_F_(npm)->val[cell_id]   | Soot precursor number
  rtp(iel,\ref isca(\ref iygfm))    | CS_F_(ygfm)->val[cell_id]  | Fresh gas fraction
  rtp(iel,\ref isca(\ref iyfm))     | CS_F_(yfm)->val[cell_id]   | Mass fraction
  rtp(iel,\ref isca(\ref iyfp2m))   | CS_F_(yfp2m)->val[cell_id] | Mass fraction variance
  rtp(iel,\ref isca(\ref icoyfp))   | CS_F_(coyfp)->val[cell_id] | Mass fraction covariance </tt>


  \subsection rayt Radiative transfer

  Defined in \ref radiat.f90 and \ref cs_field_pointer.h.

  Fortran code                 | C code                               | Description
  ---------------------------- | ------------------------------------ | ------------
  <tt> propce(iel,\ref ilumin) | CS_F_(rad_lumin)->val[cell_id]       | Radiative luminance
  propce(iel,\ref iqx) \n propce(iel,\ref iqy) \n propce(iel,\ref iqz) | CS_F_(rad_q)->val[cell_id][0] \n CS_F_(rad_q)->val[cell_id][1] \n CS_F_(rad_q)->val[cell_id][2] | Radiative flux
  propce(iel,\ref itsre(iesp)) | CS_FI_(rad_ets,iesp-1)->val[cell_id] | Radiative flux explicit source term
  propce(iel,\ref itsri(iesp)) | CS_FI_(rad_its,iesp-1)->val[cell_id] | Radiative flux implicit source term
  propce(iel,\ref iabs(iesp))  | CS_FI_(rad_abs,iesp-1)->val[cell_id] | Radiative absorption
  propce(iel,\ref iemi(iesp))  | CS_FI_(rad_emi,iesp-1)->val[cell_id] | Radiative emission
  propce(iel,\ref icak(iesp))  | CS_FI_(rad_cak,iesp-1)->val[cell_id] | Radiative absorption coefficient
  call field_get_val_s(\ref itparo,btparo) \n btparo(iel) | CS_F_(tparo)->val[cell_id]      | Wall temperature
  call field_get_val_s(\ref iqinci,bqinci) \n bqinci(iel) | CS_F_(qinci)->val[cell_id]      | Radiative incident radiative flux density
  call field_get_val_s(\ref ixlam,bxlam) \n bxlam(iel)    | CS_F_(xlam)->val[cell_id]       | Wall thermal conductivity
  call field_get_val_s(\ref iepa,bepa) \n bepa(iel)       | CS_F_(epa)->val[cell_id]        | Wall thickness
  call field_get_val_s(\ref ieps,beps) \n beps(iel)       | CS_F_(emissivity)->val[cell_id] | Wall emissivity
  call field_get_val_s(\ref ifnet,bfnet) \n bfnet(iel)    | CS_F_(fnet)->val[cell_id]       | Boundary radiative flux
  call field_get_val_s(\ref ifconv,bfconv) \n bfconv(iel) | CS_F_(fconv)->val[cell_id]      | Boundary radiative convective flux
  call field_get_val_s(\ref ihconv,bhconv) \n bhconv(iel) | CS_F_(hconv)->val[cell_id]      | Radiative exchange coefficient </tt>

*/
// _____________________________________________________________________________
/*!

  \page local Common local variables

  The following table provides a non-exhaustive list of local variables which
  are used in the code in a recurring manner.

  Fortran code  | C code     | Description
  ------------- | ---------- | ------------
  <tt> iel      | cell_id    | Cell index
  ifac          | face_id    | Face index
  ig            | g_id       | Interior face number of associated groups
  it            | t_id       | Interior face number of threads
  idimtr        | tr_dim     | Indicator for tensor perodicity of rotation
  flumas        | i_massflux | Mass flux at interior faces
  flumab        | b_massflux | Mass flux at boundary faces
  viscf         | i_visc     | \f$ \mu_\fij \dfrac{S_\fij}{\ipf \jpf} \f$ \n  at interior faces for the r.h.s.
  viscb         | b_visc     | \f$ \mu_\fib \dfrac{S_\fib}{\ipf \centf} \f$ \n  at border faces for the r.h.s.
  smbrp         | rhs        | Right hand side \f$ \vect{Rhs} \f$ </tt>

  \section conv Local naming convention for fields (Fotran and C)

  Rules have been stablished for local names denoting fields, depending on their nature. The convention, applying both in Fortran and in C, is as follows:

  - The first letter of the name indicates the location at which the field values are defined:
    - \b c for values at the cell centers.
    - \b i for values at the interior faces.
    - \b b for values at the boundary faces.
  - The next three letters indicate if the field is a variable (at the current time step or the previous time step) or a property:
    - \b var for variables at the current time step.
    - \b vara for variables at the previous time step.
    - \b pro for properties.
  - An underscore \b _ follows.
  - Finally, the <b> short name </b> of the variable/property is specified. This short name is built from the variable/property Fortran index, removing the \c i at the beginning of the word.

  The following examples ilustrate this convention:

  \c cvar_pr: Values of the variable pressure field defined at the cell centers, at the current time step. \n
  \c cvara_pr: Values of the variable pressure field defined at the cell centers, at the previous time step. \n
  \c cpro_cp: Values of the property specific heat defined field at the cell centers. \n

*/
