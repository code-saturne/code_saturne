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

!> \file entsor.f90
!> \brief Module for input/output

module entsor

  !=============================================================================

  use paramx

  implicit none

  !=============================================================================

  !> \defgroup entsor Module for input/output

  !> \addtogroup entsor
  !> \{

  !> standard output
  integer, save :: nfecra
  !> printing level
  integer, save :: iwarni(nvarmx)

  !> continuation file vortex method (necessarily ASCII , specific structure)
  integer, save :: impmvo
  !> continuation file vortex method (necessarily ASCII , specific structure)
  integer, save :: impvvo
  !> data file vortex method.
  integer, save :: impdvo
  !> name of file, see usvort module.
  character*13, save :: ficdat
  !> save frequency
  integer, save :: ntsuit

  !> in post post processing output nth variable if ichrvr(n)=1
  integer, save :: ichrvr(nvppmx)

  !> \defgroup userfile Additional user files

  !> \addtogroup userfile
  !> \{

  !> thermophysical input file for a specific physics
  !> when Janaf is used (name of the file)
  character*32, save :: ficfpp
  !> logical unit for Janaf file
  integer, save :: impfpp
  !> perform Janaf (=1) or not (=0)
  integer, save :: indjon

  !> Input files for the atmospheric specific physics
  !> (name of the meteo profile file)
  character*32, save :: ficmet
  !> logical unit of the meteo profile file
  integer, save :: impmet

  !> \}

  !> \defgroup history History user files

  !> \addtogroup history
  !> \{

  !> history files
  !> maximum number of history user files
  integer    nushmx
  parameter(nushmx=16)
  !> directory of history files
  character*80, save :: emphis
  !> prefix of history files
  character*80, save :: prehis
  !> logical unit for specific usert history files
  integer, save :: impush(nushmx)
  !> names of specific usert history files
  character*13, save :: ficush(nushmx)
  !> sytock file and mobile structure varibles output unit
  integer, save :: impsth(2)

  !> maximum number of probes
  !> see associated format in \ref ecrhis
  integer    ncaptm
  parameter(ncaptm=100)

  !> time plot format (1: .dat, 2: .csv, 3: both)
  integer, save :: tplfmt
  !> number of probes (<=  ncaptm)
  integer, save :: ncapt
  !> nthist : output frequency : > 0 ou -1 (never) or uninitialized (-999)
  integer, save :: nthist
  !> frhist : output frequency in seconds
  double precision, save :: frhist
  !> nthsav : save period (> 0) or -1. Files are open and closed
  integer, save :: nthsav
  !> ihisvr : number of probes  and index by variable (-999 = uninitialized)
  integer, save :: ihisvr(nvppmx,ncaptm+1)
  !> write indicator (O or 1) for history of internal mobile structures
  integer, save :: ihistr
  !> probes corresponding element
  integer, save :: nodcap(ncaptm)
  !> row of process containing nodcap (in parallel processing)
  integer, save :: ndrcap(ncaptm)
  !> xyzcap : required position for a probe
  double precision, save :: xyzcap(3,ncaptm)
  !> tplflw : time plot flush wall-time interval (none if <= 0)
  double precision, save :: tplflw

  !> \}

  !> \defgroup lagrange Lagrange files

  !> \addtogroup lagrange
  !> \{

  !> name of Lagrange listing
  character*6, save :: ficlal

  !> logical unit of Lagrange listing
  integer, save :: implal
  !> output period of Lagrange listing
  integer, save :: ntlal

  !> Lagrange file
  integer, save :: impla1
  !> Lagrange file
  integer, save :: impla2
  !> Lagrange file
  integer, save :: impla3
  !> Lagrange file
  integer, save :: impla4
  !> Lagrange file
  integer, save :: impla5(15)

  !> \}

  !> \addtogroup userfile
  !> \{

  ! --- Fichiers utilisateurs

  !> maximal number of user files
  integer    nusrmx
  parameter(nusrmx=10)

  !> unit numbers for potential user specified files useful if and only if the user needs files (therefore always useful, by security)
  integer, save ::      impusr(nusrmx)

  !> name of the potential user specified files. In the case of a non-parallel
  !> calculation, the suffix applied the file name is a two digit number:
  !> from \f$ \texttt{usrf01} \f$ to \f$ \texttt{usrf10} \f$ .
  !> In the case of a parallel-running calculation, the four digit processor index-number is
  !> added to the suffix. For instance, for a calculation running on two
  !> processors: from \f$ \texttt{usrf01.n\_0001} \f$ to  \f$ \texttt{usrf10.n\_0001} \f$ and
  !> from \f$ \texttt{usrf01.n\_0002} \f$ to \f$ \texttt{usrf10.n\_0002} \f$ . The opening,
  !> closing, format and location of these files must be managed by the user.
  !> useful if and only if the user needs files (therefore always useful, by security)
  character*13, save :: ficusr(nusrmx)

  !> \}

  !> \defgroup listing Output listing

  !> \addtogroup listing
  !> \{

  !> name of the variables (unknowns, physical properties ...): used in the
  !> execution listing, in the post-processing files, etc.
  !> If not initialised,  the code chooses the manes by default.
  !> It is recommended not to define variable names of more than 16
  !> characters, to get a clear execution listing (some advanced writing
  !> levels take into account only the first 16 characters).
  !> always useful}
  character*80, save :: nomvar(nvppmx)

  !> locator pointer vor variables output
  integer, save :: ipprtp(nvarmx)
  !> locator pointer vor variables output
  integer, save :: ipppro(npromx)
  !> locator pointer vor variables output
  integer, save :: ippdt
  !> locator pointer vor variables output
  integer, save :: ipptx
  !> locator pointer vor variables output
  integer, save :: ippty
  !> locator pointer vor variables output
  integer, save :: ipptz
  !> locator pointer vor variables output
  integer, save :: ipp2ra(nvppmx)

  !> for every quantity (variable, physical or numerical property ...),
  !> indicator concerning the writing in the execution report file
  !> default value (-999) is automatically converted into 1 if the concerned
  !> quantity is one of the main variables (pressure, velocity, turbulence,
  !> scalar), the density, the time step if \c idtvar > 0 or the turbulent
  !> viscosity. Otherwise converted into 0.
  !> = 1: writing in the execution listing.
  !> = 0: no writing.
  integer, save :: ilisvr(nvppmx)
  !> index of a variable if ipp corresponds to a determined variable
  !> (p,u,k...)
  !> 0 si ipp corresponds to an annex variable (cp, mut...) or to nothing.
  integer, save :: itrsvr(nvppmx)
  !> writing period in the execution report file. \c ntlist= -1: no writing
  !> \c ntlist\> 0: period (every \c ntlist time step). The value of
  !> \c ntlist must be adapted according to the number of iterations
  !> carried out in the calculation. Keeping \c ntlist to 1 will indeed provide
  !> a maximum volume of information, but if the number of time steps
  !> is too large the execution report file might become too big and unusable
  !> (problems with disk space, memory problems while opening the file with a
  !> text editor, problems finding the desired information in the file, ...).
  integer, save :: ntlist

  !> \defgroup convergence Convergence information

  !> \addtogroup convergence
  !> \{

  !> number of iterations
  integer, save :: nbivar(nvppmx)
  !> right-hand-side norm
  double precision, save :: rnsmbr(nvppmx)
  !> normed residual
  double precision, save :: resvar(nvppmx)
  !> norm of drift in time
  double precision, save :: dervar(nvppmx)

  !> \}
  !> \}

  !> \defgroup other_output Boundary post-processing

  !> \addtogroup other_output
  !> \{

  !> indicates the data to post-process on the boundary mesh (the boundary mesh must
  !> have been activated with \c ichrbo=1. Its value is
  !> the product of the following integers, depending on the variables
  !> that should be post-processed:
  !> \c ipstyp: \f$  $y^+$  \f$ at the boundary
  !> \c ipstcl: value of the variables at the
  !> boundary (using the boundary conditions but without reconstruction)
  !> \c ipstft}: thermal flux at the boundary
  !> ( \f$ $W\,m^{-2}$ \f$ ), if a thermal scalar has been defined (\c iscalt)
  !> For instance, with \c ipstdv=ipstyp*ipstcl,
  !> \f$ $y^+$  \f$ and the variables will be post-processed at the boundaries.
  !> With \c ipstdv=1, none of these data are post-processed at the boundaries.
  !> always useful if \c ichrbo=1
  integer, save :: ipstdv(6)

  !> post-processed property: Efforts (1: all; 2: tangent; 4: normal)
  integer    ipstfo
  !> post-processed property: yplus
  integer    ipstyp
  !> post-processed property: Tplus
  integer    ipsttp
  !> post-processed property: thermal flux rebuilt
  integer    ipstft
  !> post-processed property: boundary temperature
  integer    ipsttb
  !> post-processed property: Nusselt
  integer    ipstnu
  parameter (ipstfo=1, ipstyp=2, ipsttp= 3, ipstft=4, ipsttb=5, ipstnu=6)

  !> margin in seconds on the remaining CPU time which is necessary to allow
  !> the calculation to stop automatically and write all the required results
  !> (for the machines having a queue manager).
  !> = -1: calculated automatically,
  !> = 0: margin defined by the user.
  !> Always useful, but the default value should not be changed.
  double precision, save :: tmarus
  !> \}
  !> \}

  !=============================================================================

end module entsor
