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
integer          ifcvsl
double precision flux   , theta  , aa, bb, viscos, xmasvo, varcp
double precision, dimension(:), pointer :: i_mass_flux, b_mass_flux
double precision, dimension(:), pointer :: i_mass_flux_prev, b_mass_flux_prev
double precision, dimension(:), pointer :: brom, crom, broma, croma
double precision, dimension(:), pointer :: cpro_viscl, cpro_visct
double precision, dimension(:), pointer :: cpro_cp, cpro_visls
double precision, dimension(:), pointer :: cproa_cp, cproa_visls, cproa_visct
double precision, dimension(:), pointer :: cpro_romls
double precision, dimension(:), pointer :: cproa_romls

!===============================================================================

!===============================================================================
! 0. INITIALISATION
!===============================================================================

!===============================================================================
! 1.  AU TOUT DEBUT DE LA BOUCLE EN TEMPS
!===============================================================================

if (iappel.eq.1) then

! --- Application du schema en temps sur le flux de masse
!     Soit F le flux de masse
!       - si istmpf = 2 (schema non standard, theta>0 en fait = 0.5)
!         i_mass_flux      contient F_(n-2+theta) et on y met F(n-1+theta)
!         i_mass_flux_prev contient F_(n-1+theta) et
!                                    on y met une extrapolation en n+theta
!       - si istmpf = 1 (schema non standard, theta=0)
!         i_mass_flux      contient deja F_n et
!         i_mass_flux_prev n'est pas utilise : on ne fait rien
!       - sinon : istmpf = 0 (schema standard, theta= -999)
!         i_mass_flux      et i_mass_flux_prev contiennent tous les deux F(n)
!                                            : on ne fait rien


!     Ordre 2 en temps pour le flux (theta = 0.5) a entrer dans navsto
!       Flux convectif = 2F(n-1+theta)-F(n-2+theta)
!       Flux conserve  =  F(n-1+theta)
!     Au premier pas de temps, l'ancien a ete initialise dans inivar (a 0)
!       en suite de calcul, les deux ont ete relus.

  if (istmpf.eq.2) then
    call field_get_key_int(ivarfl(iu), kimasf, iflmas)
    call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
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

!     Les valeurs courantes ecrasent les valeurs anterieures
!       en cas d'extrapolation en temps (theta > 0)
!       Pour RHO, on le fait en double si ICALHY = 1 (et sur NCELET)
!     Au debut du calcul les flux nouveau et ancien ont ete initialises inivar
  if (iroext.gt.0) then
    call field_current_to_previous(icrom)
    call field_current_to_previous(ibrom)
  endif
  if (iviext.gt.0) then
    call field_get_val_s(iviscl, cpro_viscl)
    call field_get_val_s(ivisct, cpro_visct)
    call field_get_val_prev_s(iviscl, cproa_visls)
    call field_get_val_prev_s(ivisct, cproa_visct)
    do iel = 1, ncel
      cproa_visls(iel) = cpro_viscl(iel)
      cproa_visct(iel) = cpro_visct(iel)
    enddo
  endif
  if (icpext.gt.0) then
    if (icp.ge.0) then
      call field_get_val_s(icp, cpro_cp)
      call field_get_val_prev_s(icp, cproa_cp)
      do iel = 1, ncel
        cproa_cp(iel) = cpro_cp(iel)
      enddo
    endif
  endif

  ! Remarque : si on faisait cette operation pour tous les
  ! scalaires, on la ferait plusieurs fois pour les scalaires
  ! ayant une variance
  if (nscal.ge.1) then
    do iscal = 1, nscal
      ! Diffusivity
      call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
      if (ifcvsl.ge.0.and.iscavr(iscal).le.0) then
        if (ivsext(iscal).gt.0) then
          call field_get_val_s(ifcvsl, cpro_visls)
          call field_get_val_prev_s(ifcvsl, cproa_visls)
          do iel = 1, ncel
            cproa_visls(iel) = cpro_visls(iel)
          enddo
        endif
      endif
      ! Densisty
      call field_get_key_int (ivarfl(isca(iscal)), kromsl, ifcvsl)
      if (ifcvsl.ge.0.and.iscavr(iscal).le.0) then
        if (iroext.gt.0) then
          call field_get_val_s(ifcvsl, cpro_romls)
          call field_get_val_prev_s(ifcvsl, cproa_romls)
          do iel = 1, ncel
            cproa_romls(iel) = cpro_romls(iel)
          enddo
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
    if (iroext.gt.0) then
      call field_current_to_previous(icrom)
      call field_current_to_previous(ibrom)
    endif
  endif
  if (initvi.ne.1) then
    initvi = 1
    if (iviext.gt.0) then
      call field_get_val_s(iviscl, cpro_viscl)
      call field_get_val_s(ivisct, cpro_visct)
      call field_get_val_prev_s(iviscl, cproa_visls)
      call field_get_val_prev_s(ivisct, cproa_visct)
      do iel = 1, ncel
        cproa_visls(iel) = cpro_viscl(iel)
        cproa_visct(iel) = cpro_visct(iel)
      enddo
    endif
  endif
  if (initcp.ne.1) then
    initcp = 1
    if (icpext.gt.0) then
      if (icp   .gt.0) then
        call field_get_val_s(icp, cpro_cp)
        call field_get_val_prev_s(icp, cproa_cp)
        do iel = 1, ncel
          cproa_cp(iel) = cpro_cp(iel)
        enddo
      endif
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
        call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
        if (ifcvsl.ge.0.and.iscavr(iscal).le.0) then
          if (ivsext(iscal).gt.0) then
            call field_get_val_s(ifcvsl, cpro_visls)
            call field_get_val_prev_s(ifcvsl, cproa_visls)
            do iel = 1, ncel
              cproa_visls(iel) = cpro_visls(iel)
            enddo
          endif
        endif
      endif
    enddo
  endif


! 2.2 EXTRAPOLATION DES NOUVELLES VALEURS
! =======================================

! --- Extrapolation de la viscosite et de la masse volumique dans le cas d'un
!     theta schema
!     A partir de Fn-1 et Fn on calcule Fn+theta
!     On conserve les nouvelles valeurs dans l'ancien tableau pour
!     retablir en fin de pas de temps

!     Le calcul pour Rho est fait sur NCELET afin d'economiser un echange.

  if (iroext.gt.0) then
    call field_get_val_s(icrom, crom)
    call field_get_val_prev_s(icrom, croma)
    theta  = thetro
    do iel = 1, ncelet
      xmasvo = crom(iel)
      crom(iel) = (1.d0+theta)*crom(iel) - theta*croma(iel)
      croma(iel) = xmasvo
    enddo
    call field_get_val_s(ibrom, brom)
    call field_get_val_prev_s(ibrom, broma)
    do ifac = 1, nfabor
      xmasvo = brom(ifac)
      brom(ifac) = (1.d0+theta)*brom(ifac) - theta*broma(ifac)
      broma(ifac) = xmasvo
    enddo
  endif
  if (iviext.gt.0) then
    call field_get_val_s(iviscl, cpro_viscl)
    call field_get_val_s(ivisct, cpro_visct)
    call field_get_val_prev_s(iviscl, cproa_visls)
    call field_get_val_prev_s(ivisct, cproa_visct)
    theta  = thetvi
    do iel = 1, ncel
      viscos = cpro_viscl(iel)
      cpro_viscl(iel) = (1.d0+theta) * cpro_viscl(iel)    &
           -       theta  * cproa_visls(iel)
      cproa_visls(iel) = viscos
      viscos = cpro_visct(iel)
      cpro_visct(iel) = (1.d0+theta) * cpro_visct(iel)    &
           -       theta  * cproa_visct(iel)
      cproa_visct(iel) = viscos
    enddo
  endif
  if (icpext.gt.0) then
    if (icp.ge.0) then
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
      call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
      if (ifcvsl.ge.0.and.iscavr(iscal).le.0) then
        if (ivsext(iscal).gt.0) then
          theta  = thetvs(iscal)
          call field_get_val_s(ifcvsl, cpro_visls)
          call field_get_val_prev_s(ifcvsl, cproa_visls)
          do iel = 1, ncel
            viscos = cpro_visls(iel)
            cpro_visls(iel) = (1.d0+theta)*cpro_visls(iel) &
                                    -theta *cproa_visls(iel)
            cproa_visls(iel) = viscos
          enddo
        endif
      endif
      ! Density
      call field_get_key_int (ivarfl(isca(iscal)), kromsl, ifcvsl)
      if (ifcvsl.ge.0.and.iscavr(iscal).le.0) then
        if (iroext.gt.0) then
          theta = thetvs(iscal)
          call field_get_val_s(ifcvsl, cpro_romls)
          call field_get_val_prev_s(ifcvsl, cproa_romls)
          do iel = 1, ncel
            viscos = cpro_romls(iel)
            cpro_romls(iel) = (1.d0+theta)*cpro_romls(iel) &
                                    -theta *cproa_romls(iel)
            cproa_romls(iel) = viscos
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
!     si istmpf = 2 : ordre 2 (thetfl > 0 : = 0.5)
!       on calcule f(n+theta) par interpolation a partir
!       de F_(n-1+theta) et f(n+1) et on le met dans i_mass_flux
!     si istmpf = 0 : explicite (thetfl = 0) : on remet F(n) dans
!       i_mass_flux sauf a la derniere iteration (un traitement
!       complementaire sera fait en iappel=4)

!     Dans le cas ou on itere sur navsto, on passe ici
!     - a toutes les iterations si ISTMPF.NE.0
!     - a toutes les iterations sauf la derniere si ISTMPF.EQ.0

  call field_get_key_int(ivarfl(iu), kimasf, iflmas)
  call field_get_key_int(ivarfl(iu), kbmasf, iflmab)
  call field_get_val_s(iflmas, i_mass_flux)
  call field_get_val_s(iflmab, b_mass_flux)

  if (istmpf.ne.1) then
    call field_get_val_prev_s(iflmas, i_mass_flux_prev)
    call field_get_val_prev_s(iflmab, b_mass_flux_prev)
  endif

  if (istmpf.eq.2) then
    theta  = thetfl
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
!     Si istmpf = 2 : ordre 2 (thetfl > 0 : = 0.5)
!       on calcule F(n+theta) par interpolation a partir
!       de F_(n-1+theta) et F(n+1) et on le met dans i_mass_flux)
!     Si istmpf = 0 : explicite (thetfl = 0)
!       On sauvegarde F_(n+1) dans i_mas_flux_prev,mais on continue
!       les calculs avec F_(n) mis dans i_mass_flux

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
!     Si istmpf = 2 : ordre 2 (thetfl > 0 : = 0.5) : on ne fait rien
!     Si istmpf = 0 : explicite (thetfl = 0)
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

!     Le calcul pour Rho est fait sur NCELET afin d'economiser un echange.

  if (iroext.gt.0) then
    call field_get_val_s(icrom, crom)
    call field_get_val_prev_s(icrom, croma)
    do iel = 1, ncelet
      crom(iel) = croma(iel)
    enddo
    call field_get_val_s(ibrom, brom)
    call field_get_val_prev_s(ibrom, broma)
    do ifac = 1, nfabor
      brom(ifac) = broma(ifac)
    enddo
  endif
  if (iviext.gt.0) then
    call field_get_val_s(iviscl, cpro_viscl)
    call field_get_val_s(ivisct, cpro_visct)
    call field_get_val_prev_s(iviscl, cproa_visls)
    call field_get_val_prev_s(ivisct, cproa_visct)
    do iel = 1, ncel
      cpro_viscl(iel) = cproa_visls(iel)
      cpro_visct(iel) = cproa_visct(iel)
    enddo
  endif
  if (icpext.gt.0) then
    if (icp.ge.0) then
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
      call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
      if (ifcvsl.ge.0.and.iscavr(iscal).le.0) then
        if (ivsext(iscal).gt.0) then
          call field_get_val_s(ifcvsl, cpro_visls)
          call field_get_val_prev_s(ifcvsl, cproa_visls)
          do iel = 1, ncel
            cpro_visls(iel) = cproa_visls(iel)
          enddo
        endif
      endif
      ! Density
      call field_get_key_int (ivarfl(isca(iscal)), kivisl, ifcvsl)
      if (ifcvsl.ge.0.and.iscavr(iscal).le.0) then
        if (iroext.gt.0) then
          call field_get_val_s(ifcvsl, cpro_romls)
          call field_get_val_prev_s(ifcvsl, cproa_romls)
          do iel = 1, ncel
            cpro_romls(iel) = cproa_romls(iel)
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
