!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

!===============================================================================
! Function :
! ----------

!> \file schtmp.f90
!>
!> \brief Management of the mass flux, the viscosity, the density, the specific
!> heat  and the tsnsa array in case of a theta-scheme.
!>
!> Please refer to the
!> <a href="../../theory.pdf#massflux"><b>mass flux</b></a> section
!> of the theory guide for more informations.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Arguments                                                                    !
!______________________________________________________________________________!
!  mode           name               role                                      !
!______________________________________________________________________________!
!> \param[in]     nscal              total number of scalars
!> \param[in]     iappel             call number (before of after phyvar
!______________________________________________________________________________

subroutine schtmp &
 ( nscal  , iappel )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use entsor
use cstphy
use ppincl
use pointe
use parall
use period
use mesh
use field

!===============================================================================

implicit none

! Arguments

integer          nscal  , iappel

! Local variables

integer          iel    , ifac   , iscal
integer          iflmas , iflmab
integer          f_id
integer          key_t_ext_id, icpext
integer          iviext
integer          iroext

double precision flux   , theta  , aa, bb, viscos, varcp

double precision, dimension(:), pointer :: i_mass_flux, b_mass_flux
double precision, dimension(:), pointer :: i_mass_flux_prev, b_mass_flux_prev
double precision, dimension(:), pointer :: brom, broma, bromaa
double precision, dimension(:), pointer :: crom, croma, cromaa
double precision, dimension(:), pointer :: cpro_viscl, cpro_visct
double precision, dimension(:), pointer :: cpro_cp, cpro_visls
double precision, dimension(:), pointer :: cproa_cp, cproa_visls, cproa_visct
double precision, dimension(:), pointer :: cpro_rho_mass, bpro_rho_mass

!===============================================================================

!===============================================================================
! 0. Initialisation
!===============================================================================

call field_get_key_id("time_extrapolated", key_t_ext_id)

call field_get_id_try("density_mass", f_id)
if (f_id.ge.0) &
  call field_get_val_s(f_id, cpro_rho_mass)

call field_get_id_try("boundary_density_mass", f_id)
if (f_id.ge.0) &
  call field_get_val_s(f_id, bpro_rho_mass)

!===============================================================================
! 1. At the really beginning of the time step
!===============================================================================

if (iappel.eq.1) then

  ! --- Store the previous mass flux (n-1->n) in *_mass_flux_prev
  ! Note: if there is no previous values, nothing is done
  ! for explicit schemes (istmpf=0) a specific treatment is done
  ! because previous value is used as a work array...
  if (istmpf.ne.0 .or. staggered.eq.1) then
    call field_get_key_int(ivarfl(iu), kimasf, iflmas)
    call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
    if (itpcol.eq.0) then
      call field_current_to_previous(iflmas)
      call field_current_to_previous(iflmab)
    else
      call field_get_val_s(iflmas, i_mass_flux)
      call field_get_val_s(iflmab, b_mass_flux)
      call field_get_val_prev_s(iflmas, i_mass_flux_prev)
      call field_get_val_prev_s(iflmab, b_mass_flux_prev)
      do ifac = 1 , nfac
        flux                   =                          i_mass_flux(ifac)
        i_mass_flux(ifac)      = 2.d0*i_mass_flux(ifac) - i_mass_flux_prev(ifac)
        i_mass_flux_prev(ifac) = flux
      enddo

      do ifac = 1 , nfabor
        flux                   =                          b_mass_flux(ifac)
        b_mass_flux(ifac)      = 2.d0*b_mass_flux(ifac) - b_mass_flux_prev(ifac)
        b_mass_flux_prev(ifac) = flux
      enddo
    endif
  endif

  ! If required, the density at time step n-1 is updated
  ! Note that for VOF and dilatable algorithms, density at time step n-2
  ! is also updated
  ! Note that at the begining of the calculation, previous values have been
  ! initialized by inivar
  if (irovar.gt.0) then
    if (ischtp.eq.2.and.itpcol.eq.1) then
      call field_get_val_s(icrom, crom)
      call field_get_val_prev_s(icrom, croma)
      call field_get_val_prev2_s(icrom, cromaa)

      do iel = 1, ncelet
        cpro_rho_mass(iel) = 3.d0*crom(iel) - 3.d0*croma(iel)     &
                            + cromaa(iel)
      enddo

      call field_get_val_s(ibrom, brom)
      call field_get_val_prev_s(ibrom, broma)
      call field_get_val_prev2_s(ibrom, bromaa)
      do ifac = 1, nfabor
        bpro_rho_mass(ifac) = 3.d0*brom(ifac) - 3.d0*broma(ifac)  &
                             + bromaa(ifac)
      enddo
    endif

    call field_current_to_previous(icrom)
    call field_current_to_previous(ibrom)

    if (((idilat.gt.1.or.ivofmt.gt.0).and.itpcol.eq.0)            &
        .or.ippmod(icompf).eq.3) then
      ! Save the latest density field seen by the mass equation
      ! it will be updated in the correction step.
      call field_get_val_s(icrom, crom)
      do iel = 1, ncelet
        cpro_rho_mass(iel) = crom(iel)
      enddo

      call field_get_val_s(ibrom, brom)
      do ifac = 1, nfabor
        bpro_rho_mass(ifac) = brom(ifac)
      enddo

    endif
  endif

  call field_get_key_int(iviscl, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_current_to_previous(iviscl)
  endif
  call field_get_key_int(ivisct, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_current_to_previous(ivisct)
  endif
  if (icp.ge.0) then
    call field_get_key_int(icp, key_t_ext_id, icpext)
    if (icpext.gt.0) then
      call field_current_to_previous(icp)
    endif
  endif

  ! Remarque : si on faisait cette operation pour tous les
  ! scalaires, on la ferait plusieurs fois pour les scalaires
  ! ayant une variance
  if (nscal.ge.1) then
    do iscal = 1, nscal
      ! Diffusivity
      call field_get_key_int (ivarfl(isca(iscal)), kivisl, f_id)
      if (f_id.ge.0.and.iscavr(iscal).le.0) then
        call field_get_key_int(f_id, key_t_ext_id, iviext)
        if (iviext.gt.0) then
          call field_current_to_previous(f_id)
        endif
      endif
      ! Densisty
      call field_get_key_int (ivarfl(isca(iscal)), kromsl, f_id)
      if (f_id.ge.0.and.iscavr(iscal).le.0) then
        call field_get_key_int(f_id, key_t_ext_id, iroext)
        if (iroext.gt.0) then
          call field_current_to_previous(f_id)
        endif
      endif
    enddo
  endif

  return


!===============================================================================
! 2.  JUSTE APRES PHYVAR (ET DONC AVANT NAVSTO)
!===============================================================================


elseif (iappel.eq.2) then

! 2.1 MISE A JOUR DES VALEURS ANCIENNES
! =====================================

!     On passe ici dans le cas ou l'on suspecte que la valeur portee
!       par les variables "anciennes" n'est pas satisfaisante pour
!       extrapoler.
!     On passe ici au premier pas de temps et lorsque le fichier suite
!       ne comportait pas grandeur requise.

  if (initro.ne.1) then
    initro = 1
    call field_current_to_previous(icrom)
    call field_current_to_previous(ibrom)
  endif
  if (initvi.ne.1) then
    initvi = 1
    call field_current_to_previous(iviscl)
    call field_current_to_previous(ivisct)
  endif
  if (initcp.ne.1) then
    initcp = 1
    if (icp.gt.0) then
      call field_current_to_previous(icp)
    endif
  endif

  ! Remarque : si on faisant cette operation pour tous les
  ! scalaires, on la ferait plusieurs fois pour les scalaires
  ! ayant une variance
  if (nscal.ge.1) then
    do iscal = 1, nscal
      ! Diffusivity
      if (initvs(iscal).ne.1) then
        initvs(iscal) = 1
        call field_get_key_int (ivarfl(isca(iscal)), kivisl, f_id)
        if (f_id.ge.0.and.iscavr(iscal).le.0) then
          call field_current_to_previous(f_id)
        endif
      endif
    enddo
  endif


! 2.2 EXTRAPOLATION DES NOUVELLES VALEURS
! =======================================

! --- Extrapolation de la viscosite dans le cas d'un
!     theta schema
!     A partir de Fn-1 et Fn on calcule Fn+theta
!     On conserve les nouvelles valeurs dans l'ancien tableau pour
!     retablir en fin de pas de temps

   call field_get_key_int(iviscl, key_t_ext_id, iviext)
   if (iviext.gt.0) then
    call field_get_val_s(iviscl, cpro_viscl)
    call field_get_val_prev_s(iviscl, cproa_visls)
    theta  = thetvi
    do iel = 1, ncel
      viscos = cpro_viscl(iel)
      cpro_viscl(iel) = (1.d0+theta) * cpro_viscl(iel)    &
           -       theta  * cproa_visls(iel)
      cproa_visls(iel) = viscos
    enddo
  endif

  call field_get_key_int(ivisct, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_get_val_s(ivisct, cpro_visct)
    call field_get_val_prev_s(ivisct, cproa_visct)
    theta  = thetvi
    do iel = 1, ncel
      viscos = cpro_visct(iel)
      cpro_visct(iel) = (1.d0+theta) * cpro_visct(iel)    &
           -       theta  * cproa_visct(iel)
      cproa_visct(iel) = viscos
    enddo
  endif
  if (icp.gt.0) then
    call field_get_key_int(icp, key_t_ext_id, icpext)
    if (icpext.gt.0) then
      call field_get_val_s(icp, cpro_cp)
      call field_get_val_prev_s(icp, cproa_cp)
      theta  = thetcp
      do iel = 1, ncel
        varcp  = cpro_cp(iel)
        cpro_cp(iel) = (1.d0+theta) * cpro_cp(iel)      &
                      -      theta  * cproa_cp(iel)
        cproa_cp(iel) = varcp
      enddo
    endif
  endif

  !     Remarque : si on faisant cette operation pour tous les
  !       scalaires, on la ferait plusieurs fois pour les scalaires
  !       ayant une variance ET CE SERAIT FAUX
  if (nscal.ge.1) then
    do iscal = 1, nscal
      ! Diffusivity
      call field_get_key_int (ivarfl(isca(iscal)), kivisl, f_id)
      if (f_id.ge.0.and.iscavr(iscal).le.0) then
        call field_get_key_int(f_id, key_t_ext_id, iviext)
        if (iviext.gt.0) then
          theta  = thetvs(iscal)
          call field_get_val_s(f_id, cpro_visls)
          call field_get_val_prev_s(f_id, cproa_visls)
          do iel = 1, ncel
            viscos = cpro_visls(iel)
            cpro_visls(iel) = (1.d0+theta)*cpro_visls(iel) &
                                    -theta *cproa_visls(iel)
            cproa_visls(iel) = viscos
          enddo
        endif
      endif

    enddo
  endif

  return


!===============================================================================
! 3.  JUSTE APRES NAVSTO, DANS LES BOUCLES U/P ET ALE
!===============================================================================

elseif (iappel.eq.3) then

!     On traite ici le flux de masse uniquement
!        On suppose qu'il n'y en a qu'un seul.

!     si istmpf = 1 : standard : on ne fait rien
!     si istmpf = 2 : ordre 2 : on ne fait rien
!     si istmpf = 0 : explicite : on remet F(n) dans
!       i_mass_flux sauf a la derniere iteration (un traitement
!       complementaire sera fait en iappel=4)

!     Dans le cas ou on itere sur navsto, on passe ici
!     - a toutes les iterations si ISTMPF.NE.0
!     - a toutes les iterations sauf la derniere si ISTMPF.EQ.0

  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, i_mass_flux)
  call field_get_val_s(iflmab, b_mass_flux)
  call field_get_val_prev_s(iflmas, i_mass_flux_prev)
  call field_get_val_prev_s(iflmab, b_mass_flux_prev)

  if (istmpf.eq.2.and.itpcol.eq.1) then
    theta  = 0.5d0
    aa = 1.d0/(2.d0-theta)
    bb = (1.d0-theta)/(2.d0-theta)
    do ifac = 1 , nfac
      i_mass_flux(ifac) = aa * i_mass_flux(ifac) + bb * i_mass_flux_prev(ifac)
    enddo
    do ifac = 1 , nfabor
      b_mass_flux(ifac) = aa * b_mass_flux(ifac) + bb * b_mass_flux_prev(ifac)
    enddo
  else if (istmpf.eq.0) then
    do ifac = 1 , nfac
      i_mass_flux(ifac) = i_mass_flux_prev(ifac)
    enddo
    do ifac = 1 , nfabor
      b_mass_flux(ifac) = b_mass_flux_prev(ifac)
    enddo
  endif

  return

!===============================================================================
! 3.  JUSTE APRES NAVSTO (ET STRDEP), HORS LES BOUCLES U/P ET ALE
!===============================================================================

elseif (iappel.eq.4) then

!     On traite ici le flux de masse uniquement
!        On suppose qu'il n'y en a qu'un seul.

!     Si istmpf = 1 : standard : on ne fait rien
!     Si istmpf = 2 : ordre 2 : on ne fait rien
!     Si istmpf = 0 : explicite
!       On sauvegarde F_(n+1) dans i_mass_flux_prev, mais on continue
!       les calculs avec F_(n) mis dans i_mass_flux
! TODO it would be simpler to use the theta scheme of istmpf=2 with theta=0!

!     On retablira au dernier appel de schtmp pour istmpf = 0

!     Dans le cas ou on itere sur navsto, on passe ici
!       - a toutes les iterations            si istmpf.ne.0
!       - uniquement a la derniere iteration si istmpf.eq.0
!         (ce faisant, a partir de la deuxieme sous-iteration,
!          le calcul sera fait avec F(n+1) et plus F(n), mais on
!          suppose que l'utilisateur a choisi de faire des sous-iter
!          aussi pour impliciter le flux de masse)

  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, i_mass_flux)
  call field_get_val_s(iflmab, b_mass_flux)
  call field_get_val_prev_s(iflmas, i_mass_flux_prev)
  call field_get_val_prev_s(iflmab, b_mass_flux_prev)

  ! Already checked
  if (istmpf.eq.0) then
    do ifac = 1 , nfac
      flux = i_mass_flux(ifac)
      i_mass_flux(ifac) = i_mass_flux_prev(ifac)
      i_mass_flux_prev(ifac) = flux
    enddo
    do ifac = 1 , nfabor
      flux = b_mass_flux(ifac)
      b_mass_flux(ifac) = b_mass_flux_prev(ifac)
      b_mass_flux_prev(ifac) = flux
    enddo
  endif

  return

!===============================================================================
! 3.  JUSTE APRES SCALAI
!===============================================================================

elseif (iappel.eq.5) then

! 3.1 RETABLISSEMENT POUR LE FLUX DE MASSE
! ========================================

!     On corrige les manipulations sur le flux de masse faites dans
!       l'appel precedent afin d'etre pret pour le pas de temps suivant.

!     Si istmpf = 1 : standard : on ne fait rien
!     Si istmpf = 2 : ordre 2 : on ne fait rien
!     Si istmpf = 0 : explicite
!       on remet F_(n+1) (stocke dans i_mass_flux_prev) dans i_mass_flux
!       de sorte que les deux flux contiennent la meme chose

  if (istmpf.eq.0) then
    call field_get_key_int(ivarfl(iu), kimasf, iflmas)
    call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
    call field_get_val_s(iflmas, i_mass_flux)
    call field_get_val_s(iflmab, b_mass_flux)
    call field_get_val_prev_s(iflmas, i_mass_flux_prev)
    call field_get_val_prev_s(iflmab, b_mass_flux_prev)
    do ifac = 1 , nfac
      i_mass_flux(ifac) = i_mass_flux_prev(ifac)
    enddo
    do ifac = 1 , nfabor
      b_mass_flux(ifac) = b_mass_flux_prev(ifac)
    enddo
  endif

! 3.1 RETABLISSEMENT POUR LES PROPRIETES PHYSIQUES
! ================================================

  call field_get_key_int(iviscl, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_get_val_s(iviscl, cpro_viscl)
    call field_get_val_prev_s(iviscl, cproa_visls)
    do iel = 1, ncel
      cpro_viscl(iel) = cproa_visls(iel)
    enddo
  endif

  call field_get_key_int(ivisct, key_t_ext_id, iviext)
  if (iviext.gt.0) then
    call field_get_val_s(ivisct, cpro_visct)
    call field_get_val_prev_s(ivisct, cproa_visct)
    do iel = 1, ncel
      cpro_visct(iel) = cproa_visct(iel)
    enddo
  endif

  if (icp.gt.0) then
    call field_get_key_int(icp, key_t_ext_id, icpext)
    if (icpext.gt.0) then
      call field_get_val_s(icp, cpro_cp)
      call field_get_val_prev_s(icp, cproa_cp)
      do iel = 1, ncel
        cpro_cp(iel) = cproa_cp(iel)
      enddo
    endif
  endif

  ! Remarque : si on faisant cette operation pour tous les
  ! scalaires, on la ferait plusieurs fois pour les scalaires
  ! ayant une variance
  if (nscal.ge.1) then
    do iscal = 1, nscal
      ! Diffusivity
      call field_get_key_int (ivarfl(isca(iscal)), kivisl, f_id)
      if (f_id.ge.0.and.iscavr(iscal).le.0) then
        call field_get_key_int(f_id, key_t_ext_id, iviext)
        if (iviext.gt.0) then
          call field_get_val_s(f_id, cpro_visls)
          call field_get_val_prev_s(f_id, cproa_visls)
          do iel = 1, ncel
            cpro_visls(iel) = cproa_visls(iel)
          enddo
        endif
      endif
    enddo
  endif

  return

endif

!===============================================================================

!----
! End
!----

end subroutine
