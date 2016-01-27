!-------------------------------------------------------------------------------

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

!> \file radiat.f90
!> Module for Radiation

module radiat

  !===========================================================================

  use ppppar

  implicit none

  !> \defgroup radiat Module for Radiative transfer

  !> \addtogroup radiat
  !> \{

  !===========================================================================

  !> Activation of the radiative transfer module:
  !>  - 0: not activated
  !>  - 1: DOM
  !>  - 2: P1
  integer, save :: iirayo

  !> Phase which radiats (Bulk by default, but may be coal class or fuel
  !> drolplets phase)

  integer, save :: nrphas

  !> Verbosity level in the listing concerning the calculation of
  !> the wall temperatures:
  !>  - 0: no display
  !>  - 1: standard
  !>  - 2: complete
  !> useful if and only if the radiation module is activated
  integer, save :: iimpar

  !> Verbosity level in the listing concerning the solution of
  !> the radiative transfer equation:
  !>  - 0: no display
  !>  - 1: standard
  !>  - 2: complete
  !> Useful if and only if the radiation module is activated
  integer, save :: iimlum

  !> When gas or coal combustion is activated, \ref imodak indicates whether the
  !> absorption coefficient shall be calculated ``automatically'' (=1) or read from
  !> the data file (=0)
  !> Useful if the radiation module is activated
  integer, save :: imodak

  !>ADF model:
  !> - 0 no ADF model
  !> - 1 ADF model with 8 intervals of wave length
  !> - 2 ADF model with 50 intervals of wave length
  integer, save :: imoadf

  !>FSCK model:
  !> - 0 no FSCK model
  !> - 1 FSCK model activated
  integer, save :: imfsck

  !--> pointeur dans le macrotableau propce :

  !                       ITSRE --> Terme source explicite
  !                       ITSRI --> Terme source implicite
  !                       IQX,IQY,IQZ --> Composantes du vecteur densite de flux radiatif
  !                       IABSO --> part d'absorption dans le terme source explicite
  !                       IEMI --> part d'emission dans le terme source explicite
  !                       ICAK --> coefficient d'absorption
  !                       ILUMIN --> POINTEUR QUI PERMET DE REPERER L INTEGRALE DE LA
  !                                  LUMINANCE DANS LA TABLEAU propce

  integer, save ::  itsre(1+nclcpm) , itsri(1+nclcpm) ,                      &
                    iqx   ,   iqy   , iqz   ,                                &
                    iabso(1+nclcpm) , iemi(1+nclcpm)  , icak(1+nclcpm)  ,    &
                    ilumin

  !--> field ids for specific  boundary fields
  !                       IQINCI --> densite de flux incident radiatif
  !                       IXLAM  --> conductivite thermique de la paroi
  !                       IEPA   --> epaisseur de la paroi
  !                       IEPS   --> emissivite de la paroi
  !                       IFNET  --> Flux Net radiatif
  !                       IFCONV --> Flux Convectif
  !                       IHCONV --> Coef d'echange fluide
  !                       IQINSP --> densite de flux incident radiatif spectral

  integer, save ::  iqinci = -1
  integer, save ::  ixlam  = -1
  integer, save ::  iepa   = -1
  integer, save ::  ieps   = -1
  integer, save ::  ifnet  = -1
  integer, save ::  ifconv = -1
  integer, save ::  ihconv = -1
  integer, save ::  iqinsp = -1

  !--> XNP1MX : pour le modele P-1,
  !     pourcentage de cellules pour lesquelles on admet que l'epaisseur
  !     optique depasse l'unite bien que ce ne soit pas souhaitable
  !> With the P-1 model (\ref iirayo =2), \ref xnp1mx is the percentage of cells of
  !> the calculation domain for which it is acceptable that the optical
  !> thickness is lower than unity (more precisely, where \f$ KL \f$ is lower than
  !> 1, where \f$ K \f$ is the absorption coefficient of the medium and \f$ L \f$ is a
  !> characteristic length of the domain), although it is not to be desired
  !> Useful if and only if the radiation module is activated with the P-1 method
  double precision, save ::  xnp1mx

  !--> ISTPP1 : pour le modele P-1,
  !     indicateur d'arret mis a 1 dans ppcabs si le pourcentage de cellules
  !     pour lesquelles l'epaisseur optique depasse l'unite est superieur a
  !     XNP1MX  (on s'arrete a la fin du pas de temps)
  integer, save ::           istpp1

  !> Indicates the method used to calculate the radiative source term:
  !>  - 0: semi-analytic calculation (compulsory with transparent media)
  !>  - 1: conservative calculation
  !>  - 2: semi-analytic calculation corrected in order to be globally conservative
  !> Useful if and only if the radiation module is activated
  !> \remark If the medium is transparent, the choice has no effect on the calculation
  integer, save ::           idiver

  !> Index of the quadrature and number of directions for a single octant
  !> - Quadrature Sn (n(n+2) directions)
  !>
  !>   - 1: S4 (24 directions)
  !>   - 2: S6 (48 directions)
  !>   - 3: S8 (80 directions)
  !>
  !> - Quadrature Tn (8n^2 directions)
  !>
  !>   - 4: T2 (32 directions)
  !>   - 5: T4 (128 directions)
  !>   - 6: Tn (8*ndirec^2 directions)
  integer, save :: i_quadrature

  !> Parameter assiociated to the Tn
  integer, save :: ndirec

  !> For the Tn quadrature, \ref ndirec squared
  integer, save :: ndirs

  !--> directions of angular values of the quadrature sx, sy, sz
  !    and weight of the solid angle associated

  double precision, dimension(:), allocatable :: sx, sy, sz, angsol

  !> Indicates whether the radiation variables should be initialized (=0) or read
  !> from a restart file (=1)
  !> Useful if and only if the radiation module is activated (in this case, a
  !> restart file rayamo must be available)
  integer, save :: isuird

  !> Period of the radiation module.
  !> The radiation module is called every \ref nfreqr time steps (more precisely,
  !> every time \ref optcal::ntcabs "ntcabs" is a multiple of \ref nfreqr).
  !> Also, in order to have proper initialization of the variables, whatever
  !> the value of \ref nfreqr, the radiation module is called at
  !> the first time step of a calculation (restart or not).
  !> Useful if and only if the radiation module is activated}
  integer, save ::           nfreqr

  !> Spectral radiation models (ADF and FSCK)
  !> Number of ETRs to solve
  integer, save ::           nwsgg
  !> Weights of the Gaussian quadrature
  double precision, dimension(:), allocatable :: wq

  !--> Informations sur les zones frontieres

  ! NBZRDM Nombre max. de  zones frontieres
  ! NOZRDM Numero max. des zones frontieres

  integer    nbzrdm
  parameter (nbzrdm=2000)
  integer    nozrdm
  parameter (nozrdm=2000)

  ! NZFRAD Nombre de zones de bord (sur le proc courant)
  ! ILZRAY Liste des numeros de zone de bord (du proc courant)
  ! NOZARM Numero de zone de bord atteint max
  !   exemple zones 1 4 2 : NZFRAD=3,NOZARM=4

  integer, save ::           nozarm, nzfrad, ilzrad(nbzrdm)

  !--> Types de condition pour les temperatures de paroi :
  !       ITPIMP Profil de temperature imposee
  !       IPGRNO Parois grises ou noires
  !       IPREFL Parois reflechissante
  !       IFGRNO Flux de conduction impose dans la paroi
  !                   ET paroi non reflechissante (EPS non nul)
  !       IFREFL Flux de conduction impose dans la paroi
  !                   ET paroi reflechissante     (EPS = 0)
  !       ITPT1D Resolution de l'equation de la chaleur (module tp1d)

  integer   itpimp   , ipgrno   , iprefl   , ifgrno   , ifrefl   , itpt1d
  parameter(itpimp=1 , ipgrno=21, iprefl=22, ifgrno=31, ifrefl=32, itpt1d=4)

  !> \}

  !=============================================================================

contains

  !=============================================================================

  ! Allocate arrays

  subroutine init_quadrature(ndirs)

    ! Arguments

    integer, intent(in) :: ndirs

    ! Local variables

    integer :: err = 0

    if (.not.allocated(sx)) then
      allocate(sx(ndirs), stat=err)
    endif

    if (.not.allocated(sy)) then
      allocate(sy(ndirs), stat=err)
    endif

    if (.not.allocated(sz)) then
      allocate(sz(ndirs), stat=err)
    endif

    if (.not.allocated(angsol)) then
      allocate(angsol(ndirs), stat=err)
    endif

    if (err /= 0) then
      write (*, *) "Error allocating array."
      call csexit(err)
    endif

    return

  end subroutine init_quadrature

  !=============================================================================

  ! Free related arrays

  subroutine finalize_quadrature

    if (allocated(sx)) then
      deallocate(sx)
    endif

    if (allocated(sy)) then
      deallocate(sy)
    endif

    if (allocated(sz)) then
      deallocate(sz)
    endif

    if (allocated(angsol)) then
      deallocate(angsol)
    endif

  end subroutine finalize_quadrature

  !=============================================================================

end module radiat
