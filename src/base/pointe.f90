!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2012 EDF S.A.
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

! Module for pointer variables

module pointe

  !=============================================================================

  use paramx

  !=============================================================================

  !... Auxiliaires

  ! Array  ! Dimension               ! Description

  ! dispar ! ncelet                  ! distance a la face de type 5 (phase 1) la
  !                                    plus proche
  ! yplpar ! ncelet                  ! yplus associe (LES only)
  ! forbr  ! nfabor*3                ! efforts aux bords (si posttraite)
  ! yplbr  ! nfabor                  ! yplus bord (si post-traite)
  ! uetbor ! nfabor                  ! uetbor  bord (si LES+VanDriest)
  !        !                         ! ou si le modele de depot est actif
  ! idfstr ! nfabor                  ! tableau d'indirection pour les structures
  !                                    mobiles EN ALE

  double precision, allocatable, dimension(:,:) :: forbr
  double precision, allocatable, dimension(:) :: dispar, yplpar
  double precision, allocatable, dimension(:) :: yplbr, uetbor

  ! Specific arrays for the coupled case

  ! coefau ! (3,nfabor)     ! explicit Boundary conditions for the velocity
  ! coefbu ! (3,3,nfabor)   ! implicit Boundary conditions for the velocity
  ! claale ! (3,nfabor)     ! explicit Boundary conditions for the mesh velocity
  ! clbale ! (3,3,nfabor)   ! implicit Boundary conditions for the mesh velocity
  ! cfaale ! (3,nfabor)     ! explicit Boundary conditions for the mesh velocity
  ! cfbale ! (3,3,nfabor)   ! implicit Boundary conditions for the mesh velocity
  ! coefaut! (3,nfabor)     ! explicit Boundary conditions for the turbulent heat flux
  ! coefbut! (3,3,nfabor)   ! implicit Boundary conditions for the turbulent heat flux
  ! cofarut! (3,nfabor)     ! explicit Boundary conditions for the turbulent heat flux
  !                         !  in the temperature transport equation
  ! cofbrut! (3,3,nfabor)   ! implicit Boundary conditions for the turbulent heat flux
  !                         !  in the temperature transport equation

  double precision, dimension(:,:), allocatable :: coefau, cofafu
  double precision, dimension(:,:,:), allocatable :: coefbu, cofbfu

  double precision, dimension(:,:), allocatable :: cfaale, claale
  double precision, dimension(:,:,:), allocatable :: cfbale, clbale

  double precision, dimension(:,:), allocatable :: coefaut, cofafut, cofarut
  double precision, dimension(:,:,:), allocatable :: coefbut, cofbfut, cofbrut

  ! dudxy  ! (ncelet-ncel,3,3)   ! sauvegarde du gradient de la
  !        !                     ! vitesse en cas de rotation
  ! wdudxy ! (ncelet-ncel,3,3)   ! tableau de travail lie a dudxyz
  ! drdxy  ! (3,6,ncelet-ncel)   ! sauvegarde du gradient de rij
  !        !                     ! en cas de rotation
  ! wdrdxy ! (3,6,ncelet-ncel)   ! tableau de travail lie a drdxyz

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
  ! straio ! ncelet                  ! strain rate tensor at the
  !                                    previous time step

  double precision, allocatable, dimension(:) :: s2kw , divukw
  double precision, allocatable, dimension(:,:) :: straio

  !... Parametres du module thermique 1D

  ! nfpt1d !                         ! nb faces de bord avec module thermique 1D
  ! izft1d ! nfabor                  ! zone de t1d
  ! nppt1d ! nfpt1d                  ! nombre de mailles dans la paroi
  ! eppt1d ! nfpt1d                  ! epaisseur de la paroi
  ! rgpt1d ! nfpt1d                  ! raison du maillage
  ! ifpt1d ! nfpt1d                  ! numero de la face
  ! tppt1d ! nfpt1d                  ! temperature de paroi
  ! iclt1d ! nfpt1d                  ! type de condition limite
  ! tept1d ! nfpt1d                  ! temperature exterieure
  ! hept1d ! nfpt1d                  ! coefficient d'echange exterieur
  ! fept1d ! nfpt1d                  ! flux thermique exterieur
  ! xlmbt1 ! nfpt1d                  ! diffusivite thermique
  ! rcpt1d ! nfpt1d                  ! rho*Cp
  ! dtpt1d ! nfpt1d                  ! pas de temps

  integer, save :: nfpt1d, nmxt1d
  integer, allocatable, dimension(:) :: izft1d, nppt1d , ifpt1d , iclt1d
  double precision, allocatable, dimension(:) :: eppt1d , rgpt1d , tppt1d
  double precision, allocatable, dimension(:) :: tept1d , hept1d , fept1d
  double precision, allocatable, dimension(:) :: xlmbt1 , rcpt1d , dtpt1d

  !... Auxiliaires

  ! ncepdc !                         ! nombre de cellules avec pdc
  ! icepdc ! ncepdc                  ! numero des cellules avec pertedecharge
  ! izcpdc ! ncelet                  ! zone de pdc
  ! ckupdc ! (ncepdc,6)              ! valeur des coeff de pdc

  integer, save :: ncepdc
  integer, allocatable, dimension(:) :: icepdc, izcpdc
  double precision, allocatable, dimension(:,:) :: ckupdc

  ! ncetsm !                         ! nombre de cellules avec tsm
  ! icetsm ! ncetsm*nvar             ! numero des cellules avec tsmasse
  ! itypsm ! ncetsm                  ! type de tsm
  ! izctsm ! ncelet                  ! zone de tsm

  integer, save :: ncetsm
  integer, allocatable, dimension(:) :: icetsm, izctsm
  integer, allocatable, dimension(:,:) :: itypsm
  double precision, allocatable, dimension(:,:) :: smacel

  ! porosi ! ncelet                  ! value of the porosity
  double precision, allocatable, dimension(:) :: porosi

  ! visten ! ncelet                  ! symmetric tensor cell visco
  double precision, allocatable, dimension(:,:) :: visten

contains

  !=============================================================================

  ! Initialize auxilliary arrays

  subroutine init_aux_arrays &

( ncelet , ncel   , ncelbr , nfac  , nfabor )

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

    ! Local variables
    integer                iok, ivar

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
      if (ivelco.eq.1) then
        allocate(cfaale(3,nfabor), claale(3,nfabor))
        allocate(cfbale(3,3,nfabor), clbale(3,3,nfabor))
      endif
    endif

    ! Boundary condition for the velocity when components are coupled

    if (ivelco.eq.1) then
      allocate(coefau(3,nfabor),cofafu(3,nfabor))
      allocate(coefbu(3,3,nfabor),cofbfu(3,3,nfabor))
    endif

    ! Boundary condition for the turbulent heat flux when transport equations are used
    if (ityturt.eq.3) then
      allocate(coefaut(3,nfabor),cofafut(3,nfabor),cofarut(3,nfabor))
      allocate(coefbut(3,3,nfabor),cofbfut(3,3,nfabor),cofbrut(3,3,nfabor))
    endif

    ! Porosity array when needed

    if (iporos.eq.1) then
      allocate(porosi(ncelet))
    endif

    ! Symmetric cell diffusivity when needed
    iok = 0
    do ivar = 1, nvarmx
      if (idften(ivar).eq.6) iok = iok + 1
    enddo

    if (iok.ge.1) then
      allocate(visten(6,ncelet))
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
        allocate(drdxy(3,6,ncelet-ncel), wdrdxy(3,6,ncelet-ncel))
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

   ! Strain rate tensor at the previous time step
   ! if rotation curvature correction
   if (irccor.eq.1) then
      if (idtvar.ge.0) then
        allocate(straio(ncelet,6))
      endif
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
    if (allocated(izcpdc)) deallocate(izcpdc)
    if (allocated(izctsm)) deallocate(izctsm)
    if (allocated(izft1d)) deallocate(izft1d)
    if (allocated(coefau)) deallocate(coefau, cofafu, coefbu, cofbfu)
    if (allocated(coefaut)) deallocate(coefaut, cofafut, cofarut, &
                                       coefbut, cofbfut, cofbrut)
    if (allocated(porosi)) deallocate(porosi)
    if (allocated(visten)) deallocate(visten)
    if (allocated(cfaale)) deallocate(cfaale, cfbale, claale, clbale)
    if (allocated(dispar)) deallocate(dispar)
    if (allocated(yplpar)) deallocate(yplpar)
    if (allocated(ifapat)) deallocate(ifapat)
    if (allocated(dudxy)) deallocate(dudxy, wdudxy)
    if (allocated(drdxy)) deallocate(drdxy, wdrdxy)
    if (allocated(forbr)) deallocate(forbr)
    if (allocated(uetbor)) deallocate(uetbor)
    if (allocated(yplbr)) deallocate(yplbr)
    if (allocated(s2kw)) deallocate(s2kw, divukw)
    if (allocated(straio))  deallocate(straio)

    return

  end subroutine finalize_aux_arrays

  !=============================================================================

  subroutine init_kpdc

    allocate(icepdc(ncepdc))
    allocate(ckupdc(ncepdc,6))

  end subroutine init_kpdc

  !=============================================================================

  subroutine finalize_kpdc

    deallocate(icepdc)
    deallocate(ckupdc)

  end subroutine finalize_kpdc

  !=============================================================================

  subroutine init_tsma ( nvar )

    implicit none

    integer :: nvar

    allocate(icetsm(ncetsm))
    allocate(itypsm(ncetsm,nvar))
    allocate(smacel(ncetsm,nvar))

  end subroutine init_tsma

  !=============================================================================

  subroutine finalize_tsma

    deallocate(icetsm)
    deallocate(itypsm)
    deallocate(smacel)

  end subroutine finalize_tsma

  !=============================================================================

  subroutine init_pt1d

    use cstnum

    allocate(nppt1d(nfpt1d), ifpt1d(nfpt1d), iclt1d(nfpt1d))
    allocate(eppt1d(nfpt1d), rgpt1d(nfpt1d), tppt1d(nfpt1d))
    allocate(tept1d(nfpt1d), hept1d(nfpt1d), fept1d(nfpt1d))
    allocate(xlmbt1(nfpt1d), rcpt1d(nfpt1d), dtpt1d(nfpt1d))

    !---> INITIALISATION DES TABLEAUX
    !     a des valeurs sortant en erreur dans vert1d
    !     sauf pour les variables de conditions aux limites
    !     qui sont initialisees a des valeurs standard
    !     (flux nul, Timpose=0, coef d'echange infini)

    ifpt1d(:) = -999
    nppt1d(:) = -999
    iclt1d(:) = 3
    eppt1d(:) = -999.d0
    rgpt1d(:) = -999.d0
    tppt1d(:) = 0.d0
    tept1d(:) = 0.d0
    hept1d(:) = rinfin
    fept1d(:) = 0.d0
    xlmbt1(:) = -999.d0
    rcpt1d(:) = -999.d0
    dtpt1d(:) = -999.d0

  end subroutine init_pt1d

  !=============================================================================

  subroutine finalize_pt1d

    deallocate(nppt1d, ifpt1d, iclt1d)
    deallocate(eppt1d, rgpt1d, tppt1d)
    deallocate(tept1d, hept1d, fept1d)
    deallocate(xlmbt1, rcpt1d, dtpt1d)

  end subroutine finalize_pt1d

end module pointe
