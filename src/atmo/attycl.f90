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
!> \param[in]   izfppp          boundary face zone number for atmospheric module
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
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradv \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!
!-------------------------------------------------------------------------------
subroutine attycl ( itypfb, izfppp, icodcl, rcodcl )

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
use siream
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          izfppp(nfabor)
integer          icodcl(nfabor,nvar)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, iel, izone
integer          ii
integer jsp, isc
double precision d2s3, zent, vs, xuent, xvent
double precision xkent, xeent, tpent, qvent,ncent
double precision xcent
double precision viscla, uref2, rhomoy, dhy, xiturb
double precision, dimension(:), pointer :: brom, coefap, viscl

! arrays for cressman interpolation
double precision , dimension(:),allocatable :: u_bord
double precision , dimension(:),allocatable :: v_bord
double precision , dimension(:),allocatable :: tke_bord
double precision , dimension(:),allocatable :: eps_bord
double precision , dimension(:),allocatable :: theta_bord
double precision , dimension(:),allocatable :: qw_bord
double precision , dimension(:),allocatable :: nc_bord

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

d2s3 = 2.d0/3.d0

xuent = 0.d0
xvent = 0.d0
xkent = 0.d0
xeent = 0.d0
tpent = 0.d0

call field_get_val_s(ibrom, brom)
call field_get_val_s(iviscl, viscl)

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
    allocate(u_bord(nfabor))
    call mscrss(id_u,2,u_bord)
    if (imbrication_verbose) then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,ubord=", cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                               &
             cdgfbo(3,ifac),                                               &
             u_bord(ifac)
      enddo
    endif
  endif

  if (cressman_v) then
    allocate(v_bord(nfabor))
    call mscrss(id_v,2,v_bord)
    if(imbrication_verbose)then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,vbord=", cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                               &
             cdgfbo(3,ifac),                                               &
             v_bord(ifac)
      enddo
    endif
  endif

  if (cressman_tke) then
    allocate(tke_bord(nfabor))
    call mscrss(id_tke,2,tke_bord)
    if(imbrication_verbose)then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,tkebord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                                &
             cdgfbo(3,ifac),                                                &
             tke_bord(ifac)
      enddo
    endif
  endif

  if (cressman_eps) then
    allocate(eps_bord(nfabor))
    call mscrss(id_eps,2,eps_bord)
    if(imbrication_verbose)then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,epsbord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                                &
             cdgfbo(3,ifac),                                                &
             eps_bord(ifac)
      enddo
    endif
  endif

  if (cressman_theta .and. ippmod(iatmos).ge.1) then
    allocate(theta_bord(nfabor))
    call mscrss(id_theta, 2, theta_bord)
    if (imbrication_verbose) then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,thetabord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                                  &
             cdgfbo(3,ifac),                                                  &
             theta_bord(ifac)
      enddo
    endif
  endif

  if (cressman_qw .and. ippmod(iatmos).ge.2) then
    allocate(qw_bord(nfabor))
    call mscrss(id_qw, 2, qw_bord)
    if (imbrication_verbose) then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,qwbord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                               &
             cdgfbo(3,ifac),                                               &
             qw_bord(ifac)
      enddo
    endif
  endif

  if (cressman_nc .and. ippmod(iatmos).ge.2) then
    allocate(nc_bord(nfabor))
    call mscrss(id_nc, 2, nc_bord)
    if (imbrication_verbose) then
      do ifac = 1, max(nfabor,100), 1
        write(nfecra,*)"attycl::xbord,ybord,zbord,ncbord=",cdgfbo(1,ifac), &
             cdgfbo(2,ifac),                                               &
             cdgfbo(3,ifac),                                               &
             nc_bord(ifac)
      enddo
    endif
  endif

endif ! imbrication_flag

! ==============================================================================
! the interpolated field u_bord,v_bord,  .. will replace the values from the
! standard meteo profile
! ==============================================================================

call field_get_coefa_s(ivarfl(ipr), coefap)

do ifac = 1, nfabor

  izone = izfppp(ifac)

  if (iprofm(izone).eq.1.and.imeteo.eq.1) then

!     On recupere les valeurs du profil et on met a jour RCODCL s'il n'a pas
!       ete modifie. Il servira si la face est une face d'entree ou si c'est une
!       face de sortie (si le flux est rentrant).
    zent = cdgfbo(3,ifac)

    if (imbrication_flag .and.cressman_u) then
      xuent = u_bord(ifac)
    else
      call intprf                                                    &
      !==========
      (nbmetd, nbmetm,                                               &
       zdmet, tmmet, umet , zent  , ttcabs, xuent )
    endif

    if (imbrication_flag .and.cressman_v) then
      xvent = v_bord(ifac)
    else
      call intprf                                                    &
      !==========
      (nbmetd, nbmetm,                                               &
       zdmet, tmmet, vmet , zent  , ttcabs, xvent )
    endif

    if (imbrication_flag .and.cressman_tke) then
      xkent = tke_bord(ifac)
    else
      call intprf                                                    &
      !==========
      (nbmetd, nbmetm,                                               &
       zdmet, tmmet, ekmet, zent  , ttcabs, xkent )
    endif

    if (imbrication_flag .and.cressman_eps) then
      xeent = eps_bord(ifac)
    else
      call intprf                                                    &
      !==========
      (nbmetd, nbmetm,                                               &
       zdmet, tmmet, epmet, zent  , ttcabs, xeent )
    endif

    if(imbrication_flag .and.cressman_theta                          &
       .and. ippmod(iatmos).ge.1 ) then
        tpent = theta_bord(ifac)
    else
      call intprf                                                    &
      !==========
      (nbmett, nbmetm,                                               &
       ztmet, tmmet, tpmet, zent  , ttcabs, tpent )
    endif

    vs = xuent*surfbo(1,ifac) + xvent*surfbo(2,ifac)

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
      if (rcodcl(ifac,iw,1).gt.rinfin*0.5d0) rcodcl(ifac,iw,1) = 0.d0

      call turbulence_bc_set_uninit_inlet_k_eps(ifac, xkent, xeent, rcodcl)

      if (iscalt.ne.-1) then

        if (rcodcl(ifac,isca(iscalt),1).gt.rinfin*0.5d0) &
             rcodcl(ifac,isca(iscalt),1) = tpent

          !  Humid Atmosphere
          if ( ippmod(iatmos).eq.2 ) then
            if (rcodcl(ifac,isca(iymw),1).gt.rinfin*0.5d0)  then
              if (imbrication_flag .and. cressman_qw)then
                qvent = qw_bord(ifac)
              else
                call intprf &
                !==========
                (nbmett, nbmetm, ztmet, tmmet, qvmet, zent, ttcabs, qvent )
              endif
              rcodcl(ifac,isca(iymw),1) = qvent
            endif

            if (rcodcl(ifac,isca(intdrp),1).gt.rinfin*0.5d0)  then
              if (imbrication_flag .and. cressman_nc)then
                ncent = nc_bord(ifac)
              else
                call intprf &
                !==========
                (nbmett, nbmetm, ztmet, tmmet, ncmet, zent, ttcabs, ncent )
              endif
              rcodcl(ifac,isca(intdrp),1) = ncent
            endif
          endif

      endif

    endif

    if (iatmst.gt.0) then

      if (iautom(ifac).ge.1) then

        ! Dirichlet on the pressure: expressed in solved pressure directly
        ! (not in total pressure), that is why -1 is used (transformed as 1 in typecl.f90).
        icodcl(ifac,ipr) = -1
        rcodcl(ifac,ipr,1) = coefap(ifac)

        ! Dirichlet on turbulent variables
        call turbulence_bc_set_uninit_inlet_k_eps(ifac, xkent, xeent, rcodcl)

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
          rcodcl(ifac,iw,1) = 0.d0

        endif

      endif

    endif

  else

    if (itypfb(ifac).eq.ientre.or.itypfb(ifac).eq.i_convective_inlet) then

      if (icalke(izone).eq.1) then
        uref2 =   rcodcl(ifac,iu,1)**2                       &
                + rcodcl(ifac,iv,1)**2                       &
                + rcodcl(ifac,iw,1)**2
        uref2 = max(uref2,epzero)
        rhomoy = brom(ifac)
        iel    = ifabor(ifac)
        viscla = viscl(iel)
        dhy    = dh(izone)
        call turbulence_bc_inlet_hyd_diam(ifac,                            &
                                          uref2, dhy, rhomoy, viscla,      &
                                          rcodcl)
      else if (icalke(izone).eq.2) then
        uref2 =   rcodcl(ifac,iu,1)**2                       &
                + rcodcl(ifac,iv,1)**2                       &
                + rcodcl(ifac,iw,1)**2
        uref2 = max(uref2,epzero)
        dhy    = dh(izone)
        xiturb = xintur(izone)
        call turbulence_bc_inlet_turb_intensity(ifac,                      &
                                                uref2, xiturb, dhy,        &
                                                rcodcl)
      endif

    endif

  endif

enddo

! Atmospheric gaseous chemistry
if (ifilechemistry.ge.1) then

 do ifac = 1, nfabor

  if (itypfb(ifac).eq.ientre) then

   izone = izfppp(ifac)

   if (iprofc(izone).eq.1) then

    zent = cdgfbo(3,ifac)

    ! For species present in the concentration profiles file,
    ! profiles are used here as boundary conditions if boundary conditions have
    ! not been treated earier (eg, in usatcl)
    do ii = 1, nespgi
      if (rcodcl(ifac,isca(isca_chem(idespgi(ii))),1).gt.0.5d0*rinfin) then
        call intprf                                                    &
        !==========
        (nbchmz, nbchim,                                               &
        zproc, tchem, espnum(1+(ii-1)*nbchim*nbchmz), zent  , ttcabs, xcent )
        ! The first nespg user scalars are supposed to be chemical species
        rcodcl(ifac,isca(isca_chem(idespgi(ii))),1) = xcent
      endif
    enddo

   endif

   ! For other species zero dirichlet conditions are imposed,
   ! unless they have already been treated earlier (eg, in usatcl)
   do ii =1 , nespg
    if (rcodcl(ifac,isca(isca_chem(ii)),1).gt.0.5d0*rinfin) then
      rcodcl(ifac,isca(isca_chem(ii)),1) = 0.0d0
    endif
   enddo

  endif

 enddo

endif

! Atmospheric aerosol chemistry
if (iaerosol.eq.1) then

  do ifac = 1, nfabor

    if (itypfb(ifac).eq.ientre) then

     izone = izfppp(ifac)

      if (iprofa(izone).eq.1) then
        do jsp = 1, nesp_aer*nbin_aer+nbin_aer
          isc = (isca_chem(1) - 1)  + nespg_siream + jsp
          if (rcodcl(ifac,isca(isc),1).gt.0.5d0*rinfin) &
              rcodcl(ifac,isca(isc),1) = dlconc0(jsp)
        enddo
      endif

      ! For other species zero dirichlet conditions are imposed,
      ! unless they have already been treated earlier (eg, in usatcl)
      do ii = 1, nesp_aer*nbin_aer+nbin_aer
        isc = (isca_chem(1) - 1)  + nespg_siream + ii
        if (rcodcl(ifac,isca(isc),1).gt.0.5d0*rinfin) &
            rcodcl(ifac,isca(isc),1) = 0.0d0
      enddo

       ! For gaseous species which have not been treated earlier
       ! (for example species not present in the third gaseous scheme,
       ! which can be treated in usatcl of with the file chemistry)
       ! zero dirichlet conditions are imposed
       do ii = 1, nespg_siream
        isc = (isca_chem(1) - 1)  + ii
        if (rcodcl(ifac,isca(isc),1).gt.0.5d0*rinfin) &
            rcodcl(ifac,isca(isc),1) = 0.0d0
       enddo

    endif

  enddo

endif

!----
! FORMATS
!----

! ---------------------------------
! clean up the 'imbrication'
! ---------------------------------
if (imbrication_flag)then
  if(cressman_u) then
    deallocate(u_bord)
  endif
  if(cressman_v) then
    deallocate(v_bord)
  endif
  if(cressman_tke) then
    deallocate(tke_bord)
  endif
  if(cressman_eps) then
    deallocate(eps_bord)
  endif
  if(cressman_theta .and. ippmod(iatmos).ge.1 ) then
    deallocate(theta_bord)
  endif
  if(cressman_qw .and. ippmod(iatmos).ge.2 ) then
    deallocate(qw_bord)
  endif
  if(cressman_nc .and. ippmod(iatmos).ge.2 ) then
    deallocate(nc_bord)
  endif
endif

!----
! FIN
!----

return
end subroutine attycl
