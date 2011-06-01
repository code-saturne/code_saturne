!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

! Module for pointer variables

module pointe

  !=============================================================================

  use paramx

  !=============================================================================

  !... Auxiliaires

  ! Array ! Dimension               ! Description

  ! cocg   ! ncelet*9                ! stockage pour gradient
  ! cocgb  ! ncelbr*9                ! stockage pour gradient bord
  ! coci   ! ncelet*9                ! stockage pour gradient si init. par mc
  ! cocib  ! ncelbr*9                ! stockage pour gradient bord si init. par
  ! dispar ! ncelet                  ! distance a la face de type 5 (phase 1) la
  !                                    plus proche
  ! yplpar ! ncelet                  ! yplus associe (LES only)
  ! forbr  ! nfabor*3                ! efforts aux bords (si posttraite)
  ! yplbr  ! nfabor                  ! yplus bord (si post-traite)
  ! uetbor ! nfabor                  ! uetbor  bord (si LES+VanDriest)
  !        !                         ! ou si le modele de depot est actif
  ! idfstr ! nfabor                  ! tableau d'indirection pour les structures
  !                                    mobiles EN ALE

  double precision, allocatable, dimension(:,:,:) :: cocg, cocgb, coci, cocib
  double precision, allocatable, dimension(:,:) :: forbr
  double precision, allocatable, dimension(:) :: dispar, yplpar
  double precision, allocatable, dimension(:) :: yplbr, uetbor

  ! dudxy  ! (ncelet-ncel,3,3)   ! sauvegarde du gradient de la
  !        !                     ! vitesse en cas de rotation
  ! wdudxy ! (ncelet-ncel,3,3)   ! tableau de travail lie a dudxyz
  ! drdxy  ! (ncelet-ncel,6,3)   ! sauvegarde du gradient de rij
  !        !                     ! en cas de rotation
  ! wdrdxy ! (ncelet-ncel,6,3)   ! tableau de travail lie a drdxyz

  double precision, allocatable, dimension(:,:,:) :: dudxy, wdudxy
  double precision, allocatable, dimension(:,:,:) :: drdxy, wdrdxy

  ! itypfb ! nfabor                  ! type des faces de bord
  ! itrifb ! nfabor                  ! indirection pour tri faces de bord
  ! izfppp ! nfabor                  ! pour reperage des zones frontieres
  !        !                         ! associees aux faces de bord (phys. part.)
  ! izfrad ! nfabor                  ! pour reperage des zones frontieres
  !        !                         ! associees aux faces de bord (radiat.)
  ! isympa ! nfabor                  ! zero pour annuler le flux de masse
  !        !                         !   (symetries et parois avec cl couplees)
  !        !                         ! un sinon
  ! ifapat ! ncelet                  ! numero de face de bord 5 la plus proche
  ! Peut etre serait il plus approprie de le verser dans pointe

  integer, allocatable, dimension(:) :: itypfb, itrifb, izfppp, izfrad
  integer, allocatable, dimension(:) :: isympa, ifapat

  integer, allocatable, dimension(:) :: idfstr

  ! s2kw   ! ncelet                  ! stockage de 2 Sij.Sij en k-omega
  ! divukw ! ncelet                  ! stockage de divu en k-omega (en meme
  !                                    temps que s2kw)

  double precision, allocatable, dimension(:) :: s2kw , divukw

  !... Parametres du module thermique 1D
  ! nfpt1d !                         ! nb faces de bord avec module thermique 1D
  ! inppt1 ! nfpt1d                  ! nombre de mailles dans la paroi
  ! ieppt1 ! nfpt1d                  ! epaisseur de la paroi
  ! irgpt1 ! nfpt1d                  ! raison du maillage
  ! iifpt1 ! nfpt1d                  ! numero de la face
  ! itppt1 ! nfpt1d                  ! temperature de paroi
  ! iiclt1 ! nfpt1d                  ! type de condition limite
  ! itept1 ! nfpt1d                  ! temperature exterieure
  ! ihept1 ! nfpt1d                  ! coefficient d'echange exterieur
  ! ifept1 ! nfpt1d                  ! flux thermique exterieur
  ! ixlmt1 ! nfpt1d                  ! diffusivite thermique
  ! ircpt1 ! nfpt1d                  ! rho*Cp
  ! idtpt1 ! nfpt1d                  ! pas de temps

  integer, save :: nfpt1d , nmxt1d , inppt1 , iifpt1 , iiclt1 ,    &
                   ieppt1 , irgpt1 , itppt1 ,                      &
                   itept1 , ihept1 , ifept1 ,                      &
                   ixlmt1 , ircpt1 , idtpt1

  !... Auxiliaires accessibles directement dans ia, ra
  !     tous ces tableaux sont (dimension)

  ! Pointeur Dimension       Description
  ! ncepdc !                         ! nombre de cellules avec pdc
  ! iicepd ! ncepdc                  ! numero des cellules avec pertedecharge
  ! ickupd ! (ncepdc,6)              ! valeur des coeff de pdc
  ! ncetsm !                         ! nombre de cellules avec tsm
  ! iicesm ! ncetsm                  ! numero des cellules avec tsmasse
  ! iitpsm ! ncetsm                  ! type de tsm
  ! ismace ! ncetsm                  ! valeur de tsm

  integer, save :: ncepdc,                 &
                   iicepd, ickupd, ncetsm, &
                   iicesm, iitpsm, ismace

contains

  !=============================================================================

  ! Initialize auxilliary arrays

  subroutine init_aux_arrays &

( ncelet , ncel   , ncelbr , nfac  , nfabor , &
  iverif )

    use paramx
    use parall
    use period
    use optcal
    use entsor
    use ppincl
    use lagran
    use radiat
    use albase
    use ihmpre

    implicit none

    ! Arguments

    integer, intent(in) :: ncelet, ncel, ncelbr, nfac, nfabor
    integer, intent(in) :: iverif

    ! Boundary-face related arrays

    allocate(itrifb(nfabor), itypfb(nfabor))
    if (ippmod(iphpar).ge.1 .or. iihmpr.eq.1) then
      allocate(izfppp(nfabor))
    endif
    if (iirayo.gt.0) then
      allocate(izfrad(nfabor))
    endif

    ! Symmetry faces

    allocate(isympa(nfabor))

    ! ALE array for structure definition

    if (iale.eq.1) then
      allocate(idfstr(nfabor))
    endif

    ! Gradient calculation

    allocate(cocg(ncelet,3,3), cocgb(ncelbr,3,3))
    if (imrgra.eq.1 .or. iverif.eq.1) then
      allocate(coci(ncelet,3,3), cocib(ncelbr,3,3))
    endif

    ! Wall-distance calculation

    if (ineedy.eq.1 .and. abs(icdpar).eq.1) then
      allocate(dispar(ncelet))
      if (     (itytur.eq.4 .and. idries.eq.1) &
          .or. (iilagr.ge.1 .and. iroule.eq.2) ) then
        allocate(yplpar(ncelet))
      endif
    endif
    if (ineedy.eq.1 .and. abs(icdpar).eq.2) then
      allocate(ifapat(ncelet))
    endif

    ! Periodicity (work arrays for rotation)

    if (iperot.gt.0) then
      allocate(dudxy(ncelet-ncel,3,3), wdudxy(ncelet-ncel,3,3))
      if (itytur.eq.3) then
        allocate(drdxy(ncelet-ncel,6,3), wdrdxy(ncelet-ncel,6,3))
      endif
    endif

    ! Forces on boundary faces

    if (ineedf.eq.1) then
      allocate(forbr(3,nfabor))
    endif

    ! Friction velocity on boundary faces

    if (     (itytur.eq.4 .and. idries.eq.1) &
        .or. (iilagr.ge.1 .and. idepst.gt.0) ) then
      allocate(uetbor(nfabor))
    endif

    ! Non-dimensional distance to the wall (for post-processing)

    if (mod(ipstdv,ipstyp).eq.0) then
      allocate(yplbr(nfabor))
    endif

    ! Temporary storage arrays for k-omega model

    if (iturb.eq.60) then
      allocate(s2kw(ncelet))
      allocate(divukw(ncelet))
    endif

    return

  end subroutine init_aux_arrays

  !=============================================================================

  ! Free auxilliary arrays

  subroutine finalize_aux_arrays

    deallocate(itrifb, itypfb)
    if (allocated(izfppp)) deallocate(izfppp)
    if (allocated(izfrad)) deallocate(izfrad)
    if (allocated(idfstr)) deallocate(idfstr)
    deallocate(cocg, cocgb)
    if (allocated(coci)) deallocate(coci, cocib)
    if (allocated(dispar)) deallocate(dispar)
    if (allocated(yplpar)) deallocate(yplpar)
    if (allocated(ifapat)) deallocate(ifapat)
    if (allocated(dudxy)) deallocate(dudxy, wdudxy)
    if (allocated(drdxy)) deallocate(drdxy, wdrdxy)
    if (allocated(forbr)) deallocate(forbr)
    if (allocated(uetbor)) deallocate(uetbor)
    if (allocated(yplbr)) deallocate(yplbr)
    if (allocated(s2kw)) deallocate(s2kw, divukw)

    return

  end subroutine finalize_aux_arrays

  !=============================================================================

end module pointe
