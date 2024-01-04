!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

!> \file scalai.f90
!> \brief Resolution of source term convection diffusion equations
!> for scalars in a time step.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nscal         total number of scalars
!> \param[in]     iterns        Navier-Stokes iteration number
!> \param[in]     dt            time step (per cell)
!______________________________________________________________________________

subroutine scalai &
 ( nscal  ,                                              &
   iterns , dt     )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use entsor
use optcal
use cstphy
use cstnum
use pointe
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use atchem
use sshaerosol, only : iaerosol, CS_ATMO_AEROSOL_OFF
use field
use cs_c_bindings
use cs_cf_bindings
use cfpoin, only: hgn_relax_eq_st
use cs_nz_condensation

!===============================================================================

implicit none

! Arguments

integer          nscal, iterns
double precision dt(ncelet)

! Local variables

integer          iscal, ivar, iel, isou, ifac
integer          ii, iisc, itspdv, icalc, iappel, ifcvsl
integer          ispecf, scal_id, f_id, f_dim

double precision fmb, hint, pimp, hext, cmax, cmid, cmin
double precision, allocatable, dimension(:) :: viscf, viscb

double precision, dimension(:), pointer :: cvar_var, cvara_var
double precision, dimension(:), pointer :: cvar_fm, cpro_viscls
double precision, dimension(:,:), pointer :: cvar_vav, cvara_vav
double precision, dimension(:), pointer :: coefap, coefbp, cofafp, cofbfp
double precision, dimension(:), pointer :: coefa_fm, coefb_fm

integer :: keyvar

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

procedure() :: kinrates, ppinv2, cs_coal_masstransfer
procedure() :: set_dirichlet_scalar, max_mid_min_progvar, elflux
procedure() :: compute_gaseous_chemistry, elreca

interface

  subroutine cs_cf_energy                            &
     (field_id)                                      &
    bind(C, name='cs_cf_energy')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: field_id
  end subroutine cs_cf_energy

   subroutine solve_equation_scalar                   &
     (f_id, ncesmp,                                   &
     iterns, itspdv,                                  &
     icetsm, itypsm,                                  &
     smacel,                                          &
     viscf,  viscb)                                   &
    bind(C, name='cs_f_solve_equation_scalar')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: f_id, ncesmp, iterns, itspdv
    integer(c_int), dimension(*) :: icetsm, itypsm
    real(kind=c_double), dimension(*) :: smacel
    real(kind=c_double), dimension(*), intent(inout) :: viscf, viscb
  end subroutine solve_equation_scalar

  subroutine solve_equation_vector                   &
    (id,     ncesmp,                                 &
     iterns, icetsm,                                 &
     itypsm, smacel,                                 &
     viscf, viscb)                                   &
    bind(C, name='cs_f_solve_equation_vector')
    use, intrinsic :: iso_c_binding
    implicit none
    integer(c_int), value :: id, ncesmp, iterns
    integer(c_int), dimension(*) :: icetsm, itypsm
    real(kind=c_double), dimension(*) :: smacel
    real(kind=c_double), dimension(*), intent(inout) :: viscf, viscb
  end subroutine solve_equation_vector

end interface

!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_key_id("scalar_id", keyvar)

! Allocate temporary arrays for the species resolution
allocate(viscf(nfac), viscb(nfabor))

ipass = ipass + 1

! Atmospheric chemistry => all chemical fields are not buoyant
if (      ichemistry.ge.1 .and. iaerosol.eq.CS_ATMO_AEROSOL_OFF &
    .and. nespg.gt.0 .and. iterns.eq.-1) then
  ! Computation of kinetics rates
  call kinrates()
endif

!===============================================================================
! 2. Handle model or specific physics scalars.
!===============================================================================

if (nscapp.gt.0) then

  ! Initialize variable values from BC's.

  if (ippmod(iphpar).ge.1) then

    call ppinv2

    ! TODO: check that the following loop may be removed, as it is
    !       handled in tridim.
    !
    !       For safety, current values are copied to previous values.
    !       This should not be required, as it is handled before,
    !       but may in fact used in cpflux and a few others.

    if (ipass.eq.1.and.isuite.eq.0) then
      if (nscapp.gt.0) then
        do ii = 1, nscapp
          iscal = iscapp(ii)
          ivar  = isca(iscal)
          call field_get_dim(ivarfl(ivar), f_dim)
          if (f_dim.eq.1) then
            call field_get_val_s(ivarfl(ivar), cvar_var)
            call field_get_val_prev_s(ivarfl(ivar), cvara_var)
            do iel = 1, ncelet
              cvara_var(iel) = cvar_var(iel)
            enddo
          else
            call field_get_val_v(ivarfl(ivar), cvar_vav)
            call field_get_val_prev_v(ivarfl(ivar), cvara_vav)
            do iel = 1, ncelet
              do isou = 1, 3
                cvara_vav(isou,iel)= cvar_vav(isou,iel)
              enddo
            enddo
          endif
        enddo
      endif
    endif

  endif

! ---> Calculs TS relatifs a la physique du charbon
!      GMDEV1, GMDEV2, GMHET, GMDCH
  if (ippmod(iccoal).ne.-1) then

    call cs_coal_masstransfer(ncelet, ncel, volume)

  endif

!    ATTENTION : POUR LE CLIPPING AVEC ICLP = 1, IL FAUT

!                QUE LES SCALAIRES SOIENT RESOLUS AVANT
!                LEUR VARIANCE ASSOCIEE

! ---> Boucle sur les scalaires physiques particulieres.
!      On peut imaginer a la place des resolutions couplees.
!      Ici, on ne donne qu'un exemple.

  do ii = 1, nscapp

    iscal = iscapp(ii)
    ivar  = isca(iscal)

!     Schema compressible sans choc :
! ---> Traitement special pour la masse volumique,
!                     la temperature et l'energie
!     L'indicateur ISPECF sera non nul si l'on ne doit pas resoudre
!       le scalaire plus bas avec solve_equation_scalar.

    ispecf = 0

    if (ippmod(icompf).ge.0 .and. iterns.eq.-1) then

      if ( iscal.eq.itempk ) then
        ispecf = 1
      elseif ( iscal.eq.ienerg ) then
        ispecf = 2
      endif

! ---> Masse volumique : deja resolue
! ---> Temperature     : n'est pas une variable transportee
! ---> Enthalpie       : a resoudre

      if (ispecf.eq.2) then
         call cs_cf_energy(ivarfl(isca(iscal)))
      endif

    endif

!     Pour le compressible, on ne resout pas celles deja resolues ou
!       qui ne doivent pas l'etre
!     Pour les autres physiques, on resout les scalaires dans l'ordre
!       (pas de tests a faire)

    if (ispecf.eq.0) then

! ---> Variances et scalaires

!   iscavr dira scalaire/variance
!   itspdv dira si calcul des termes de prod et dissip supplementaires
!         ou non (1 : oui, 0 : non)

!         si iscavr = 0
!           scalaire
!           itspdv = 0
!         sinon
!           variance
!           si iscavr > 0 et iscavr < nscal+1
!             itspdv = 1
!           sinon
!             pour le moment, on s'arrete
!             a terme, les combustionnistes pourront donner leur propre
!               grandeur scalaire associee a la variance et
!               eventuellement reconstruite en dehors
!               de solve_equation_scalar.
!               il faudra alors peut etre declarer des tableaux
!               de travail suppl
!           fin si
!         fin si

      iisc = iscal
      if (iscavr(iisc).eq.0) then
        itspdv = 0
      elseif (iscavr(iisc).gt.0.and.iscavr(iisc).le.nscal) then
        itspdv = 1
      else
        write(nfecra,9000)iisc,iisc,nscal,iscavr(iisc)
        call csexit (1)
      endif

      ! Specific treatment BC for gas combustion: steady laminar flamelet
      if (ippmod(islfm).ge.0) then
        call field_get_coefa_s(ivarfl(isca(ifm)), coefa_fm)
        call field_get_coefb_s(ivarfl(isca(ifm)), coefb_fm)
        call field_get_val_s(ivarfl(isca(ifm)), cvar_fm)
        call field_get_coefa_s(ivarfl(ivar), coefap)
        call field_get_coefb_s(ivarfl(ivar), coefbp)
        call field_get_coefaf_s(ivarfl(ivar), cofafp)
        call field_get_coefbf_s(ivarfl(ivar), cofbfp)

        call field_get_key_int (ivarfl(ivar), kivisl, ifcvsl)
        if (ifcvsl.ge.0) then
          call field_get_val_s(ifcvsl, cpro_viscls)
        endif
        ! Make sure that mixture fraction is solved
        if (mode_fp2m.eq.1) then
          if (ivar.eq.isca(ifsqm)) then
            do ifac = 1, nfabor
              iel = ifabor(ifac)
              fmb = coefa_fm(ifac) + coefb_fm(ifac)*cvar_fm(iel)
              hint = (cpro_viscls(iel))/distb(ifac)
              pimp = fmb**2.d0
              hext = rinfin
              call set_dirichlet_scalar                               &
                ( coefap(ifac), cofafp(ifac),                         &
                  coefbp(ifac), cofbfp(ifac),                         &
                  pimp        , hint        , hext )
            enddo
          endif
        endif

        if (ippmod(islfm).ge.2) then
          if (ivar.eq.isca(ipvm)) then
            do ifac = 1, nfabor
              iel = ifabor(ifac)
              fmb = coefa_fm(ifac) + coefb_fm(ifac)*cvar_fm(iel)
              hint = (cpro_viscls(iel))/distb(ifac)

              call max_mid_min_progvar(fmb,cmax,cmid,cmin)
              pimp = cmid
              hext = rinfin
              call set_dirichlet_scalar                               &
                ( coefap(ifac), cofafp(ifac),                         &
                  coefbp(ifac), cofbfp(ifac),                         &
                  pimp        , hint        , hext )
            enddo
          endif
        endif
      endif

      call field_get_dim(ivarfl(isca(iscal)), f_dim)

      if (f_dim.eq.1) then

        call solve_equation_scalar                          &
             (ivarfl(isca(iisc)), ncetsm, iterns,           &
             itspdv, icetsm, itypsm, smacel,                &
             viscf, viscb)

     else

       call solve_equation_vector                           &
            (ivarfl(isca(iisc)), ncetsm, iterns, icetsm,    &
            itypsm, smacel, viscf, viscb)

      endif

! ---> Versions Electriques
!             Effet Joule
!             Arc Electrique
!             Conduction ionique

!     On calcule ici j, E, j.E reels et imagimaires

      if ((ippmod(ieljou).ge.1 .or. ippmod(ielarc).ge.1)  &
        .and. iterns.eq.-1) then

!     On utilise le  fait que les scalaires sont dans l'ordre
!       H, PotR, [PotI], [A] pour faire le calcul de j, E et j.E
!       apres  la determination de PotR [et PotI].

        icalc = 0
!            On peut y aller apres PotR si on est en arc
!                                         ou en Joule sans PotI

        if (ippmod(ielarc).ge.1.or.ippmod(ieljou).eq.1             &
             .or.ippmod(ieljou).eq.3) then
          call field_get_id('elec_pot_r', f_id)
          call field_get_key_int(f_id, keyvar, scal_id)
          if (iscal.eq.scal_id) then
            icalc = 1
          endif
        endif

        !  On y va apres PotI si on est en Joule avec PotI
        if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
          call field_get_id('elec_pot_i', f_id)
          call field_get_key_int(f_id, keyvar, scal_id)
          if (iscal.eq.scal_id) then
            icalc = 1
          endif
        endif

        if (icalc.eq.1) then

          ! Compute j, E and j.E
          iappel = 1

          call elflux(iappel)

          ! Readjust electric variables j, j.E (and Pot, E)

          if (ielcor .eq.1  .and. ntcabs .gt. 1) then
            call elreca(dt)
          endif

        endif

      endif

    endif

  enddo  ! End of loop on specific physical models
endif

! Electric arcs:
! computation of magnetic field B and Laplace effect jxB
if (ippmod(ielarc).ge.1 .and. iterns.eq.-1) then
  ! On utilise le  fait que les scalaires sont dans l'ordre
  !   H, PotR, [PotI], [A] pour faire le calcul de A, B, jxB
  !   apres la determination et le recalage de j
  iappel = 2
  call elflux(iappel)
endif

! Compressible homogeneous two-phase model:
! return to equilibrium source term step for volume, mass, energy fractions
if (ippmod(icompf).eq.1.and.hgn_relax_eq_st.ge.0) then
  call cs_cf_hgn_source_terms_step
endif

!===============================================================================
! 3. TRAITEMENT DES SCALAIRES UTILISATEURS STANDARD
!     On voit que meme s'ils sont numerotes en premier, on peut les
!       traiter en dernier si on veut.
!     On peut imaginer aussi de by-passer cette phase si le modele a
!       physique particuliere le demande.
!===============================================================================

if (nscaus.gt.0) then

! ---> Boucle sur les scalaires utilisateur.

  do ii = 1, nscaus

    iscal = ii
    ivar  = isca(iscal)

! ---> Variances et scalaires

!   iscavr dira scalaire/variance
!   itspdv dira si calcul des termes de prod et dissip supplementaires
!         ou non (1 : oui, 0 : non)

!         si iscavr = 0
!           scalaire
!           itspdv = 0
!         sinon
!           variance
!           si iscavr > 0 et iscavr < nscal+1
!             itspdv = 1
!           sinon
!             pour le moment, on s'arrete
!             a terme, les combustionnistes pourront donner leur propre
!               grandeur scalaire associee a la variance et
!               eventuellement reconstruite en dehors de
!               solve_equation_scalar il faudra alors
!               peut etre declarer des tableaux de travail suppl
!           fin si
!         fin si

    iisc = iscal
    if (iscavr(iisc).eq.0) then
      itspdv = 0
    elseif (iscavr(iisc).gt.0.and.iscavr(iisc).le.nscal) then
      itspdv = 1
    else
      write(nfecra,9000)iisc,iisc,nscal,iscavr(iisc)
      call csexit(1)
    endif

    call field_get_dim(ivarfl(isca(iscal)), f_dim)

    if (f_dim.eq.1) then

      call solve_equation_scalar                            &
          (ivarfl(isca(iisc)), ncetsm, iterns,              &
           itspdv, icetsm, itypsm, smacel,                  &
           viscf, viscb)
    else

      call solve_equation_vector                            &
           (ivarfl(isca(iisc)), ncetsm, iterns, icetsm,     &
            itypsm, smacel, viscf, viscb)

    endif

  enddo ! End of loop on user-defined scalars.

endif

! Atmospheric gaseous chemistry
! Resolution of chemical evolution of species
if (ichemistry.ge.1 .and. iaerosol.eq.CS_ATMO_AEROSOL_OFF .and. nespg.gt.0 .and. iterns.eq.-1) then
  call compute_gaseous_chemistry(dt)
endif

! Atmospheric gas + aerosol chemistry
if (ichemistry.ge.1 .and. iaerosol.ne.CS_ATMO_AEROSOL_OFF .and. iterns.eq.-1) then
  call cs_atmo_aerosol_time_advance()
endif

! Free memory
deallocate(viscf, viscb)

!===============================================================================
! 4.  FORMATS
!===============================================================================

 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE SOLVING SCALARS EQUATIONS          ',/,&
'@    ========                                                ',/,&
'@    SCALAR NUMBER ',I10                                      ,/,&
'@    iscavr(',I10   ,') MUST BE A POSITIVE OR NULL INTEGER   ',/,&
'@      AND LOWER OR EQUAL THAN NSCAL = ', I10                 ,/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculaton will not be run.                           ',/,&
'@                                                            ',/,&
'@  If iscavr(I) is null, the scalar I is not a variance      ',/,&
'@  If iscavr(I) is positive, the scalar I is a variance:     ',/,&
'@    it is the variance of the fluctuations of the scalar J  ',/,&
'@    whose number is iscavr(I)                               ',/,&
'@                                                            ',/,&
'@  Check parameters.                                         ',/,&
'@  Contact support.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----
return

end subroutine
