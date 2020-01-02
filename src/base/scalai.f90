!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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
!> \param[in]     nvar          total number of variables
!> \param[in]     nscal         total number of scalars
!> \param[in]     iterns        Navier-Stokes iteration number
!> \param[in]     dt            time step (per cell)
!______________________________________________________________________________

subroutine scalai &
 ( nvar   , nscal  ,                                              &
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
use siream
use field
use cs_cf_bindings
use cfpoin, only: hgn_relax_eq_st

!===============================================================================

implicit none

! Arguments

integer          nvar, nscal, iterns
double precision dt(ncelet)

! Local variables

integer          iscal, ivar, iel, isou
integer          ii, iisc, itspdv, icalc, iappel
integer          ispecf, scal_id, f_id, f_dim

double precision, allocatable, dimension(:) :: dtr
double precision, allocatable, dimension(:) :: viscf, viscb

double precision, dimension(:), pointer :: cvar_var, cvara_var
double precision, dimension(:,:), pointer :: cvar_vav, cvara_vav
integer :: keyvar

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1. Initialization
!===============================================================================

call field_get_key_id("scalar_id", keyvar)

! Allocate temporary arrays for the species resolution
allocate(dtr(ncelet))
allocate(viscf(nfac), viscb(nfabor))

ipass = ipass + 1

! Atmospheric chemistry

if (ichemistry.ge.1 .and. nscal.gt.0) then
  ! Computation of kinetics rates
  call kinrates()
  !==========
endif

!===============================================================================
! 2. Handle model or specific physics scalars.
!===============================================================================

if (nscapp.gt.0) then

  ! Initialize variable values from BC's.

  if (ippmod(iphpar).ge.1) then

    call ppinv2(nvar, nscal, dt)
    !==========

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

    call cs_coal_masstransfer &
    !=======================
   ( ncelet , ncel   ,        &
     volume )

  endif

! ---> Calculs TS relatifs a la physique du fuel
!      GAMEVA, GAMHTF

  if (ippmod(icfuel).ne.-1) then

    call cs_fuel_masstransfer(ncelet, ncel)

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

    ! ---> Pas de temps (avec facteur multiplicatif eventuel)

    if (cdtvar(ivar).ne.1.d0) then
      do iel = 1, ncel
        dtr(iel) = dt(iel)*cdtvar(ivar)
      enddo
    else
      do iel = 1, ncel
        dtr(iel) = dt(iel)
      enddo
    endif

!     Schema compressible sans choc :
! ---> Traitement special pour la masse volumique,
!                     la temperature et l'energie
!     L'indicateur ISPECF sera non nul si l'on ne doit pas resoudre
!       le scalaire plus bas avec covofi.

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

        call cfener &
        !==========
 ( nvar   , nscal  ,                                              &
   ncepdc , ncetsm ,                                              &
   iscal  ,                                                       &
   icepdc , icetsm , itypsm ,                                     &
   dt     ,                                                       &
   ckupdc , smacel )

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
!               eventuellement reconstruite en dehors de covofi.
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

      call field_get_dim(ivarfl(isca(iscal)), f_dim)

      if (f_dim.eq.1) then

! ---> Appel a covofi pour la resolution

        call covofi                                                 &
        !==========
   ( nvar   , nscal  ,                                              &
     ncepdc , ncetsm , nfbpcd , ncmast ,                            &
     iterns , iisc   , itspdv ,                                     &
     icepdc , icetsm , ifbpcd , ltmast ,                            &
     itypsm , itypcd , itypst ,                                     &
     dtr    , tslagr ,                                              &
     ckupdc , smacel , spcond , svcond , flxmst ,                   &
     viscf  , viscb  )

      else

! ---> Appel a covofv pour la resolution

        call covofv                                                 &
        !==========
   ( nvar   , nscal  ,                                              &
     ncepdc , ncetsm ,                                              &
     iterns , iisc   ,                                              &
     icepdc , icetsm ,                                              &
     itypsm ,                                                       &
     dtr    ,                                                       &
     ckupdc , smacel ,                                              &
     viscf  , viscb  )

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
!     On y va apres PotI si on est en Joule avec PotI
        if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
          call field_get_id('elec_pot_i', f_id)
          call field_get_key_int(f_id, keyvar, scal_id)
          if (iscal.eq.scal_id) then
            icalc = 1
          endif
        endif

        if (icalc.eq.1) then

!     Calcul de j, E et j.E
          iappel = 1

          call elflux(iappel)
          !==========


!     Recalage des variables electriques j, j.E (et Pot, E)

          if ( ielcor .eq.1  .and. ntcabs .gt. 1 ) then

            call elreca(dt)
            !==========

          endif

        endif

      endif

    endif


! ---> Fin de la Boucle sur les scalaires physiques particulieres.
  enddo
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

! ---> Pas de temps (avec facteur multiplicatif eventuel)

    if (cdtvar(ivar).ne.1.d0) then
      do iel = 1, ncel
        dtr(iel) = dt(iel)*cdtvar(ivar)
      enddo
    else
      do iel = 1, ncel
        dtr(iel) = dt(iel)
      enddo
    endif


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
!               eventuellement reconstruite en dehors de covofi.
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

    call field_get_dim(ivarfl(isca(iscal)), f_dim)

    if (f_dim.eq.1) then

! ---> Appel a covofi pour la resolution

      call covofi                                                   &
      !==========
   ( nvar   , nscal  ,                                              &
     ncepdc , ncetsm , nfbpcd , ncmast ,                            &
     iterns , iisc   , itspdv ,                                     &
     icepdc , icetsm , ifbpcd , ltmast ,                            &
     itypsm , itypcd , itypst ,                                     &
     dtr    , tslagr ,                                              &
     ckupdc , smacel , spcond , svcond , flxmst ,                   &
     viscf  , viscb  )

    else

! ---> Appel a covofv pour la resolution

        call covofv                                                 &
        !==========
   ( nvar   , nscal  ,                                              &
     ncepdc , ncetsm ,                                              &
     iterns , iisc   ,                                              &
     icepdc , icetsm ,                                              &
     itypsm ,                                                       &
     dtr    ,                                                       &
     ckupdc , smacel ,                                              &
     viscf  , viscb  )

      endif


! ---> Fin de la Boucle sur les scalaires utilisateurs.
  enddo

endif

! Atmospheric gaseous chemistry
! Resolution of chemical evolution of species
if (ichemistry.ge.1 .and. nscal.gt.0 .and. iterns.eq.-1) then
  call compute_gaseous_chemistry(dt)
endif

! Atmospheric aerosol chemistry
if (iaerosol.eq.1 .and. nscal.gt.0 .and. iterns.eq.-1) then
  call compute_siream(dt)
endif

! Free memory
deallocate(dtr)
deallocate(viscf, viscb)

!===============================================================================
! 4.  FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 9000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA RESOLUTION DES SCALAIRES         ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE NUMERO ',I10                                    ,/,&
'@    iscavr(',I10   ,') DOIT ETRE UN ENTIER                  ',/,&
'@      POSITIF OU NUL ET                                     ',/,&
'@      INFERIEUR OU EGAL A NSCAL = ',I10                      ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Si iscavr(I) est nul, le scalaire I n est pas une variance',/,&
'@  Si iscavr(I) est positif, le scalaire I est une variance :',/,&
'@    il s agit de la variance des fluctuations du scalaire J ',/,&
'@    dont le numero est iscavr(I)                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@  Contacter l''assistance.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

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

#endif

!----
! End
!----
return

end subroutine
