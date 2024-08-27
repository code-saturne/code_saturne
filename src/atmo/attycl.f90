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
!> \file attycl.f90
!> \brief Automatic boundary conditions for atmospheric module
!>   (based on meteo file)

!> \brief Automatically compute the boundary conditions from the meteo file
!>      or from the imbrication profiles
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   itypfb          boundary face types
!> \param[out]  icodcl          face boundary condition code
!>                               - 1 Dirichlet
!>                               - 2 Radiative outlet
!>                               - 3 Neumann
!>                               - 4 sliding and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 5 smooth wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 6 rough wall and
!>                                 \f$ \vect{u} \cdot \vect{n} = 0 \f$
!>                               - 9 free inlet/outlet
!>                                 (input mass flux blocked to 0)
!>                               - 13 Dirichlet for the advection operator and
!>                                    Neumann for the diffusion operator
!> \param[out]  rcodcl          Boundary conditions value
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!
!-------------------------------------------------------------------------------
subroutine attycl ( itypfb, icodcl, rcodcl )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only: nvar
use entsor
use parall
use ppppar
use ppthch
use ppincl
use mesh
use field
use atincl
use atchem
use atimbr
use sshaerosol
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

procedure() :: mscrss

integer          itypfb(nfabor)
integer          izfppp(nfabor)
integer          icodcl(nfabor,nvar)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, iel, ilelt
integer          ii, nbrsol, nelts
integer          jsp, isc, ivar
integer          fid_axz
double precision d2s3, zent, vs, xuent, xvent, xwent
double precision vel_dir(3), shear_dir(3)
double precision xkent, xeent, tpent, qvent,ncent
double precision xcent
double precision viscla, uref2, rhomoy, dhy, xiturb
double precision rscp, pp, dum

integer, dimension(:), pointer :: elt_ids

double precision, dimension(:), pointer :: brom, coefap, viscl
double precision, dimension(:,:), pointer :: cpro_met_vel
double precision, dimension(:), pointer :: cpro_met_potemp
double precision, dimension(:), pointer :: cpro_met_qv, cpro_met_nc
double precision, dimension(:), pointer :: cpro_met_k, cpro_met_eps
double precision, dimension(:), pointer :: cpro_met_p
double precision, dimension(:), pointer :: cpro_met_rho
double precision, dimension(:), pointer :: cpro_met_axz
double precision, pointer, dimension(:)   :: bvar_temp_sol
double precision, pointer, dimension(:)   :: bvar_tempp
double precision, pointer, dimension(:)   :: bvar_total_water

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

d2s3 = 2.d0/3.d0

xuent = 0.d0
xvent = 0.d0
xkent = 0.d0
xeent = 0.d0
tpent = 0.d0

rscp = rair/cp0

call field_get_val_s(ibrom, brom)
call field_get_val_s(iviscl, viscl)

if (imeteo.ge.2) then
  call field_get_val_s_by_name('meteo_pot_temperature', cpro_met_potemp)
  call field_get_val_v_by_name('meteo_velocity', cpro_met_vel)
  call field_get_val_s_by_name('meteo_tke', cpro_met_k)
  call field_get_val_s_by_name('meteo_eps', cpro_met_eps)
  call field_get_val_s_by_name('meteo_pressure', cpro_met_p)
  call field_get_val_s_by_name('meteo_density', cpro_met_rho)

  if (ippmod(iatmos).eq.2) then
    call field_get_val_s_by_name('meteo_humidity', cpro_met_qv)
    call field_get_val_s_by_name('meteo_drop_nb', cpro_met_nc)
  endif

endif
call field_get_id_try('meteo_shear_anisotropy', fid_axz)
if (fid_axz.ne.-1) then
  call field_get_val_s(fid_axz, cpro_met_axz)
endif

! Soil atmosphere boundary conditions
!------------------------------------
if (iatsoil.ge.1) then
  call field_get_val_s_by_name("soil_temperature", bvar_temp_sol)
  call field_get_val_s_by_name("soil_pot_temperature", bvar_tempp)
  call field_get_val_s_by_name("soil_total_water", bvar_total_water)
  call atmo_get_soil_zone(nelts, nbrsol, elt_ids)

  do ilelt = 1, nelts

    ifac = elt_ids(ilelt) + 1 ! C > Fortran

    ! Rough wall if no specified
    ! Note: roughness and thermal roughness are computed in solmoy
    if (itypfb(ifac).eq.0) itypfb(ifac) = iparug

    if (iscalt.ne.-1) then
      ! If not yet specified
      if (rcodcl(ifac,isca(iscalt),1).gt.rinfin*0.5d0)  then
        ! Dirichlet with wall function Expressed directly in term of
        ! potential temperature
        icodcl(ifac,isca(iscalt))   = -6
        rcodcl(ifac,isca(iscalt),1) = bvar_tempp(ilelt)
      endif
    endif
    if (ippmod(iatmos).eq.2) then
      ! If not yet specified
      if (rcodcl(ifac,isca(iymw),1).gt.rinfin*0.5d0)  then
        icodcl(ifac, isca(iymw)) = 6
        rcodcl(ifac, isca(iymw),1) = bvar_total_water(ilelt)
      endif
    endif

  enddo
endif

!===============================================================================
! 2.  SI IPROFM = 1 : CHOIX ENTREE/SORTIE SUIVANT LE PROFIL METEO SI
!                       ITYPFB N'A PAS ETE MODIFIE
!                     VARIABLES TIREES DU PROFIL METEO SI
!                       RCODCL(IFAC,IVAR,1) N'A PAS ETE MODIFIE

!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
! ==============================================================================
!
if (imbrication_flag) then

  call summon_cressman(ttcabs)

  if (cressman_u) then
    call mscrss(id_u,2,rcodcl(1,iu,1))
    if (imbrication_verbose) then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,ubord=", cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                               &
             cdgfbo(3,ifac),                                               &
             rcodcl(ifac, iu, 1)
      enddo
    endif
  endif

  if (cressman_v) then
    call mscrss(id_v,2,rcodcl(1,iv,1))
    if(imbrication_verbose)then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,vbord=", cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                               &
             cdgfbo(3,ifac),                                               &
             rcodcl(ifac, iv, 1)
      enddo
    endif
  endif

  if (cressman_tke) then
    call mscrss(id_tke,2,rcodcl(1,ik,1))
    if(imbrication_verbose)then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,tkebord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                                &
             cdgfbo(3,ifac),                                                &
             rcodcl(ifac, ik, 1)
      enddo
    endif
  endif

  if (cressman_eps) then
    call mscrss(id_eps,2,rcodcl(1,iep,1))
    if(imbrication_verbose)then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,epsbord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                                &
             cdgfbo(3,ifac),                                                &
             rcodcl(ifac,iep,1)
      enddo
    endif
  endif

  if (cressman_theta .and. ippmod(iatmos).ge.1) then
    call mscrss(id_theta, 2, rcodcl(1,isca(iscalt),1))
    if (imbrication_verbose) then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,thetabord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                                  &
             cdgfbo(3,ifac),                                                  &
             rcodcl(ifac,isca(iscalt),1)
      enddo
    endif
  endif

  if (cressman_qw .and. ippmod(iatmos).ge.2) then
    call mscrss(id_qw, 2, rcodcl(1,isca(iymw),1))
    if (imbrication_verbose) then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,qwbord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                               &
             cdgfbo(3,ifac),                                               &
             rcodcl(ifac,isca(iymw),1)
      enddo
    endif
  endif

  if (cressman_nc .and. ippmod(iatmos).ge.2) then
    call mscrss(id_nc, 2, rcodcl(1,isca(intdrp),1))
    if (imbrication_verbose) then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,ncbord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                               &
             cdgfbo(3,ifac),                                               &
             rcodcl(ifac,isca(intdrp),1)
      enddo
    endif
  endif

endif ! imbrication_flag

! ==============================================================================
! Note: the interpolated field are already in rcodcl(ifac,*,1)
! .. and will replace the values from the standard meteo profile
! ==============================================================================

call field_get_coefa_s(ivarfl(ipr), coefap)

do ifac = 1, nfabor

  iel = ifabor(ifac)

  ! If meteo profile is on, we take the value and store it in rcolcl if not
  ! already modified
  ! It will be used for inlet or backflows
  if (imeteo.ge.1) then
    zent = cdgfbo(3,ifac)

    ! If specified by the user or by code-code coupling
    if (rcodcl(ifac,iu,1).lt.rinfin*0.5d0) then
      xuent = rcodcl(ifac,iu,1)
    else if (imeteo.eq.1) then
      call intprf &
      (nbmetd, nbmetm,                                               &
        zdmet, tmmet, umet , zent  , ttcabs, xuent )
    else
      xuent = cpro_met_vel(1, iel)
    endif

    xwent = 0.d0
    if (rcodcl(ifac,iw,1).lt.rinfin*0.5d0) then
      xwent = rcodcl(ifac,iw,1)
    endif
    if (rcodcl(ifac,iv,1).lt.rinfin*0.5d0) then
      xvent = rcodcl(ifac,iv,1)
    else if (imeteo.eq.1) then
      call intprf &
      (nbmetd, nbmetm,                                               &
       zdmet, tmmet, vmet , zent  , ttcabs, xvent )
    else
      xvent = cpro_met_vel(2, iel)
      xwent = cpro_met_vel(3, iel)
    endif

    if (ik.ge.1) then
      if (rcodcl(ifac,ik,1).lt.rinfin*0.5d0) then
        xkent = rcodcl(ifac,ik,1)
      else if (imeteo.eq.1) then
        call intprf &
        (nbmetd, nbmetm,                                               &
        zdmet, tmmet, ekmet, zent  , ttcabs, xkent )
      else
        xkent = cpro_met_k(iel)
      endif
    else if (imeteo.eq.1) then
      call intprf &
      (nbmetd, nbmetm,                                               &
       zdmet, tmmet, ekmet, zent  , ttcabs, xkent )
    else
      xkent = cpro_met_k(iel)
    endif

    if (iep.ge.1) then
      if (rcodcl(ifac,iep,1).lt.rinfin*0.5d0) then
        xeent = rcodcl(ifac,iep,1)
      else if (imeteo.eq.1) then
        call intprf &
        (nbmetd, nbmetm,                                               &
        zdmet, tmmet, epmet, zent  , ttcabs, xeent )
      else
        xeent = cpro_met_eps(iel)
      endif
    else if (imeteo.eq.1) then
      call intprf &
      (nbmetd, nbmetm,                                               &
       zdmet, tmmet, epmet, zent  , ttcabs, xeent )
    else
      xeent = cpro_met_eps(iel)
    endif

    if (ippmod(iatmos).ge.1) then
      if (rcodcl(ifac,isca(iscalt),1).lt.rinfin*0.5d0) then
        tpent = rcodcl(ifac,isca(iscalt),1)
      else if (imeteo.eq.1) then
        call intprf &
          (nbmett, nbmetm,                                               &
          ztmet, tmmet, tpmet, zent  , ttcabs, tpent )
      else
        tpent = cpro_met_potemp(iel)
      endif
    endif

    vs = xuent*surfbo(1,ifac) + xvent*surfbo(2,ifac)

    ! Velocity direction, will be normalized afterwards
    vel_dir(1) = xuent
    vel_dir(2) = xvent
    vel_dir(3) = xwent
    shear_dir(1) = 0.d0
    shear_dir(2) = 0.d0
    if (fid_axz.eq.-1) then
      shear_dir(3) = -sqrt(cmu) ! Rxz/k
    else
      shear_dir(3) = cpro_met_axz(iel) ! Rxz/k
    endif

    ! On met a jour le type de face de bord s'il n'a pas ete specifie
    !   par l'utilisateur.
    ! Pour une entree, on remplit la condition de Dirichlet si elle n'a pas
    ! ete  specifiee par utilisateur.

    if (iautom(ifac).ge.1) then
      if (vs.gt.epzero) then
        itypfb(ifac) = isolib
      else
        if (itypfb(ifac).eq.0) itypfb(ifac) = ientre
      endif
    endif

    if (itypfb(ifac).eq.ientre.or.itypfb(ifac).eq.i_convective_inlet) then

      if (rcodcl(ifac,iu,1).gt.rinfin*0.5d0) rcodcl(ifac,iu,1) = xuent
      if (rcodcl(ifac,iv,1).gt.rinfin*0.5d0) rcodcl(ifac,iv,1) = xvent
      if (rcodcl(ifac,iw,1).gt.rinfin*0.5d0) rcodcl(ifac,iw,1) = xwent

      call turbulence_bc_set_uninit_inlet_k_eps(ifac, xkent, xeent, &
                                                vel_dir, shear_dir, rcodcl)

      if (iscalt.ne.-1) then

        if (rcodcl(ifac,isca(iscalt),1).gt.rinfin*0.5d0) &
          rcodcl(ifac,isca(iscalt),1) = tpent

        !  Humid Atmosphere
        if ( ippmod(iatmos).eq.2 ) then
          if (rcodcl(ifac,isca(iymw),1).gt.rinfin*0.5d0)  then
            if (imeteo.eq.1) then
              call intprf &
                (nbmett, nbmetm, ztmet, tmmet, qvmet, zent, ttcabs, qvent )
            else
              qvent = cpro_met_qv(iel)
            endif
            rcodcl(ifac,isca(iymw),1) = qvent
          endif

          if (rcodcl(ifac,isca(intdrp),1).gt.rinfin*0.5d0)  then
            if (imeteo.eq.1) then
              call intprf &
                (nbmett, nbmetm, ztmet, tmmet, ncmet, zent, ttcabs, ncent )
            else
              ncent = cpro_met_nc(iel)
            endif
            rcodcl(ifac,isca(intdrp),1) = ncent
          endif
        endif

      endif

    endif

    if (iatmst.gt.0) then

      if (iautom(ifac).ge.1) then

        ! Dirichlet on the pressure: expressed in solved pressure directly
        ! (not in total pressure), that is why -1 is used
        ! (transformed as 1 in cs_boundary_conditions_type).
        icodcl(ifac,ipr) = -1
        rcodcl(ifac,ipr,1) = coefap(ifac)

        ! Dirichlet on turbulent variables
        call turbulence_bc_set_uninit_inlet_k_eps(ifac, xkent, xeent, &
                                                  vel_dir, shear_dir, rcodcl)

        if (iautom(ifac).eq.1) then

          ! Homogeneous Neumann on the velocity
          icodcl(ifac,iu) = 3
          icodcl(ifac,iv) = 3
          icodcl(ifac,iw) = 3
          rcodcl(ifac,iu,3) = 0.d0
          rcodcl(ifac,iv,3) = 0.d0
          rcodcl(ifac,iw,3) = 0.d0

        else if (iautom(ifac).eq.2) then

          ! Dirichlet on the velocity
          icodcl(ifac,iu) = 1
          icodcl(ifac,iv) = 1
          icodcl(ifac,iw) = 1
          rcodcl(ifac,iu,1) = xuent
          rcodcl(ifac,iv,1) = xvent
          rcodcl(ifac,iw,1) = xwent

        endif

      endif

    endif

  endif

  ! Conversion Temperature to potential temperature for Dirichlet and
  ! wall boundary conditions
  !
  ! if icodcl < 0 it is directly expressed in term of potential temperature
  ! so no need of conversion.
  if (iscalt.ne.-1) then

    ivar = isca(iscalt)
    if (icodcl(ifac,ivar).eq.-1) then
      icodcl(ifac,ivar) = 1
    else if (icodcl(ifac,ivar).eq.-5) then
      icodcl(ifac,ivar) = 5
    else if (icodcl(ifac,ivar).eq.-6) then
      icodcl(ifac,ivar) = 6
    else if ((icodcl(ifac,ivar).eq.1 &
      .or.icodcl(ifac,ivar).eq.5     &
      .or.icodcl(ifac,ivar).eq.6)    &
      .and.rcodcl(ifac,ivar,1).lt.rinfin*0.5d0)  then

      zent = cdgfbo(3,ifac)
      if (imeteo.eq.0) then
        call atmstd(zent, pp, dum, dum)
      else if (imeteo.eq.1) then
        ! Pressure profile from meteo file:
        call intprf(nbmett, nbmetm, ztmet , tmmet , phmet , zent, ttcabs, pp)
      else
        pp = cpro_met_p(iel) - cpro_met_rho(iel) * gz * (xyzcen(3, iel) - cdgfbo(3,ifac))
      endif

      ! Convert from temperature in Kelvin to potential temperature
      rcodcl(ifac,ivar,1) = rcodcl(ifac,ivar,1) * (ps/pp)**rscp

    endif

  endif
enddo

! Atmospheric gaseous chemistry
if (ichemistry.ge.1) then

  do ifac = 1, nfabor

    if (itypfb(ifac).eq.ientre) then

      zent = cdgfbo(3,ifac)

      ! For species present in the concentration profiles chemistry file,
      ! profiles are used here as boundary conditions if boundary conditions have
      ! not been treated earlier (eg, in cs_user_boundary_conditions)
      do ii = 1, nespgi
        if (rcodcl(ifac,isca(isca_chem(idespgi(ii))),1).gt.0.5d0*rinfin) then
          call intprf &
            (nbchmz, nbchim,                                               &
            zproc, tchem, espnum(1+(ii-1)*nbchim*nbchmz), zent  , ttcabs, xcent )
          ! The first nespg user scalars are supposed to be chemical species
          rcodcl(ifac,isca(isca_chem(idespgi(ii))),1) = xcent
        endif
      enddo

      ! For other species zero Dirichlet conditions are imposed,
      ! unless they have already been treated earlier (eg, in cs_user_boundary_conditions)
      do ii =1 , nespg
        if (rcodcl(ifac,isca(isca_chem(ii)),1).gt.0.5d0*rinfin) then
          rcodcl(ifac,isca(isca_chem(ii)),1) = 0.0d0
        endif
      enddo

    endif

  enddo

endif

! Atmospheric aerosol chemistry
if (iaerosol.ne.CS_ATMO_AEROSOL_OFF) then

  do ifac = 1, nfabor

    if (itypfb(ifac).eq.ientre) then

      do jsp = 1, nlayer_aer*n_aer+n_aer
        isc = isca_chem(nespg + jsp)
        if (rcodcl(ifac,isca(isc),1).gt.0.5d0*rinfin) &
            rcodcl(ifac,isca(isc),1) = dlconc0(jsp)
      enddo

      ! For other species zero dirichlet conditions are imposed,
      ! unless they have already been treated earlier (eg, in usatcl)
      do ii = 1, nlayer_aer*n_aer+n_aer
        isc = isca_chem(nespg + ii)
        if (rcodcl(ifac,isca(isc),1).gt.0.5d0*rinfin) &
            rcodcl(ifac,isca(isc),1) = 0.0d0
      enddo

      ! For gaseous species which have not been treated earlier
      ! (for example species not present in the third gaseous scheme,
      ! which can be treated in usatcl of with the file chemistry)
      ! zero dirichlet conditions are imposed
      do ii = 1, nespg
        isc = isca_chem(ii)
        if (rcodcl(ifac,isca(isc),1).gt.0.5d0*rinfin) &
          rcodcl(ifac,isca(isc),1) = 0.0d0
      enddo

    endif

  enddo

endif

!--------
! Formats
!--------


!----
! End
!----

return
end subroutine attycl
