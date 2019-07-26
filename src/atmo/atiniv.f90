!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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
!> \file atiniv.f90
!> \brief Initialisation of calculation variables for the atmospheric module,
!> it is the counterpart of usiniv.f90.
!>
!> Initialise for example the meteorological field for each cell of
!> the domain by interpolation of the data from the meteo file
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]   nvar        total number of variables
!> \param[in]   nscal       total number of scalars
!> \param[in]   dt          time step value
!-------------------------------------------------------------------------------

subroutine atiniv  ( nvar, nscal, dt )

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use ppppar
use ppthch
use ppincl
use atincl
use mesh
use atchem
use siream
use field

!===============================================================================

implicit none

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

integer          imode, iel
double precision d2s3
double precision zent,xuent,xvent,xkent,xeent,tpent,qvent,ncent

integer k,ii, isc
double precision xcent

double precision, dimension(:,:), pointer :: vel
double precision, dimension(:), pointer :: cvar_k, cvar_ep, cvar_phi
double precision, dimension(:), pointer :: cvar_fb, cvar_omg, cvar_nusa
double precision, dimension(:,:), pointer :: cvar_rij
double precision, dimension(:), pointer :: cvar_r11, cvar_r22, cvar_r33
double precision, dimension(:), pointer :: cvar_r12, cvar_r13, cvar_r23
double precision, dimension(:), pointer :: cvar_despgi, cvar_sc
double precision, dimension(:), pointer :: cvar_scalt, cvar_totwt, cvar_ntdrp

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

d2s3 = 2.d0/3.d0

if (itytur.eq.2) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (itytur.eq.3) then
  if (irijco.eq.1) then
    call field_get_val_v(ivarfl(irij), cvar_rij)
  else
    call field_get_val_s(ivarfl(ir11), cvar_r11)
    call field_get_val_s(ivarfl(ir22), cvar_r22)
    call field_get_val_s(ivarfl(ir33), cvar_r33)
    call field_get_val_s(ivarfl(ir12), cvar_r12)
    call field_get_val_s(ivarfl(ir23), cvar_r23)
    call field_get_val_s(ivarfl(ir13), cvar_r13)
  endif
  call field_get_val_s(ivarfl(iep), cvar_ep)
elseif (iturb.eq.50) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iep), cvar_ep)
  call field_get_val_s(ivarfl(iphi), cvar_phi)
  call field_get_val_s(ivarfl(ifb), cvar_fb)
elseif (iturb.eq.60) then
  call field_get_val_s(ivarfl(ik), cvar_k)
  call field_get_val_s(ivarfl(iomg), cvar_omg)
elseif (iturb.eq.70) then
  call field_get_val_s(ivarfl(inusa), cvar_nusa)
endif

!===============================================================================
! 2. READING THE METEO PROFILE FILE (IF IMETEO = 1 DEFAULT OPTION):
!===============================================================================

if (imeteo.gt.0) then

  imode = 1
  call atlecm &
  !==========
  ( imode )

endif

if (iatra1.gt.0) then

  imode = 1
  call usatdv &
  !==========
  ( imode )

endif

! Atmospheric gaseous chemistry
if (ifilechemistry.ge.1) then

  ! Second reading of chemical profiles file
  imode = 1
  call atlecc                                                     &
  !==========
  ( imode)

  ! Computation of the conversion factor matrix used for
  ! the reaction rates jaccobian matrix
  do ii = 1, nespg
    do k = 1, nespg
      conv_factor_jac((chempoint(k)-1)*nespg+chempoint(ii)) = dmmk(ii)/dmmk(k)
    enddo
  enddo

  ! Volume initilization with profiles for species present
  ! in the chemical profiles file
  if (isuite.ne.0.and.init_at_chem.eq.1) then
    do k = 1, nespgi
      call field_get_val_s(ivarfl(isca(isca_chem(idespgi(k)))), cvar_despgi)

      do iel = 1, ncel
        zent = xyzcen(3,iel)
        call intprf                                                         &
        !==========
        (nbchmz, nbchim,                                                    &
         zproc, tchem, espnum(1+(k-1)*nbchim*nbchmz), zent  , ttcabs, xcent )
        ! The first nespg user scalars are supposed to be chemical species
        cvar_despgi(iel) = xcent
      enddo

    enddo
  endif
endif

! Atmospheric aerosol chemistry
if (iaerosol.eq.1) then

  ! Reading intial concentrations and numbers
  call atleca()

  ! Initialization
  if (init_at_chem.eq.1) then
    do ii = 1, nesp_aer*nbin_aer + nbin_aer
      isc = (isca_chem(1) - 1) + nespg_siream + ii
      call field_get_val_s(ivarfl(isca(isc)), cvar_sc)
      do iel = 1, ncel
        cvar_sc(iel) = dlconc0(ii)
      enddo
    enddo
  endif

endif

! Check simulation times used by atmo
! radiative transfer or chemistry models
if (     (iatra1.eq.1.or.ifilechemistry.ge.1)              &
    .and.(syear.eq.-999.or.squant.eq.-999.or.shour.eq.-999 &
          .or.smin.eq.-999.or.ssec.eq.-999)) then
  if (iatra1.eq.1) write(nfecra,1000)
  if (ifilechemistry.ge.1) write(nfecra,1001)
  call csexit (1)
endif

! Check radiative module latitude / longitude
if (iatra1.eq.1 .and. (xlat.ge.rinfin*0.5 .or. xlon.ge.rinfin*0.5)) then
  write(nfecra,1002)
  call csexit (1)
endif

! Check latitude / longitude from meteo file
if (imeteo.gt.0) then
  if (maxval(xmet).ge.rinfin*0.5 .or. maxval(ymet).ge.rinfin*0.5) then
    write(nfecra,1003)
    call csexit (1)
  endif
endif

! Check latitude / longitude from chemistry file
if (ifilechemistry.ge.1) then
  if (maxval(xchem).ge.rinfin*0.5 .or. maxval(ychem).ge.rinfin*0.5) then
    write(nfecra,1004)
    call csexit (1)
  endif
endif


!===============================================================================
! 3. Dry atmosphere: default initialization of potential temperature
!===============================================================================

! Only if the simulation is not a restart from another one
if (isuite.eq.0) then

  if (initmeteo.eq.1) then

    if (ippmod(iatmos).eq.1) then
      call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
    else if (ippmod(iatmos).eq.2) then
      call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)
      call field_get_val_s(ivarfl(isca(itotwt)), cvar_totwt)
      call field_get_val_s(ivarfl(isca(intdrp)), cvar_ntdrp)
    endif

    if (imeteo.eq.0) then

      if (ippmod(iatmos).eq.1) then
        ! The thermal scalar is potential temperature
        do iel = 1, ncel
          cvar_scalt(iel) = t0
        enddo
      endif

      if (ippmod(iatmos).eq.2) then
        ! The thermal scalar is liquid potential temperature
        do iel = 1, ncel
          cvar_scalt(iel) = t0
          cvar_totwt(iel) = 0.d0
          cvar_ntdrp(iel) = 0.d0
        enddo
      endif

    ! Only if meteo file is present:
    else

      do iel = 1, ncel

        zent = xyzcen(3,iel)

        call intprf                                                   &
        !==========
       (nbmetd, nbmetm,                                               &
        zdmet, tmmet, umet , zent  , ttcabs, xuent )

        call intprf                                                   &
        !==========
       (nbmetd, nbmetm,                                               &
        zdmet, tmmet, vmet , zent  , ttcabs, xvent )

        call intprf                                                   &
        !==========
       (nbmetd, nbmetm,                                               &
        zdmet, tmmet, ekmet, zent  , ttcabs, xkent )

        call intprf                                                   &
        !==========
       (nbmetd, nbmetm,                                               &
        zdmet, tmmet, epmet, zent  , ttcabs, xeent )

        vel(1,iel) = xuent
        vel(2,iel) = xvent
        vel(3,iel) = 0.d0

    !     ITYTUR est un indicateur qui vaut ITURB/10
        if    (itytur.eq.2) then

          cvar_k(iel)  = xkent
          cvar_ep(iel) = xeent

        elseif (itytur.eq.3) then

          if (irijco.eq.1) then
            cvar_rij(1,iel) = d2s3*xkent
            cvar_rij(2,iel) = d2s3*xkent
            cvar_rij(3,iel) = d2s3*xkent
            cvar_rij(4,iel) = 0.d0
            cvar_rij(5,iel) = 0.d0
            cvar_rij(6,iel) = 0.d0
          else
            cvar_r11(iel) = d2s3*xkent
            cvar_r22(iel) = d2s3*xkent
            cvar_r33(iel) = d2s3*xkent
            cvar_r12(iel) = 0.d0
            cvar_r13(iel) = 0.d0
            cvar_r23(iel) = 0.d0
          endif
          cvar_ep(iel)  = xeent

        elseif (iturb.eq.50) then

          cvar_k(iel)   = xkent
          cvar_ep(iel)  = xeent
          cvar_phi(iel) = d2s3
          cvar_fb(iel)  = 0.d0

        elseif (iturb.eq.60) then

          cvar_k(iel)   = xkent
          cvar_omg(iel) = xeent/cmu/xkent

        elseif (iturb.eq.70) then

          cvar_nusa(iel) = cmu*xkent**2/xeent

        endif


        if (ippmod(iatmos).eq.1) then
          ! The thermal scalar is potential temperature
            call intprf                                                 &
            !==========
         (nbmett, nbmetm,                                               &
          ztmet, tmmet, tpmet, zent  , ttcabs, tpent )

            cvar_scalt(iel) = tpent
        endif

        if (ippmod(iatmos).eq.2) then
          ! The thermal scalar is liquid potential temperature
            call intprf                                                 &
            !==========
         (nbmett, nbmetm,                                               &
          ztmet, tmmet, tpmet, zent  , ttcabs, tpent )
            cvar_scalt(iel) = tpent

            call intprf                                                 &
            !==========
         (nbmett, nbmetm,                                               &
          ztmet, tmmet, qvmet, zent  , ttcabs, qvent )
            cvar_totwt(iel) = qvent

            call intprf                                                 &
            !==========
         (nbmett, nbmetm,                                               &
          ztmet, tmmet, ncmet, zent  , ttcabs, ncent )
            cvar_ntdrp(iel) = ncent
        endif

      enddo

    endif
  endif

endif

!===============================================================================
! 4. USER  OPTIONS
!===============================================================================

call cs_user_f_initialization &
!==========================
( nvar   , nscal  ,                                            &
  dt     )

!----
! FORMATS
!----


#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@    MODELE DE RAYONNEMENT (IATRA1) DEMANDE                  ',/,&
'@                                                            ',/,&
'@    Le temps de la simulation est mal defini                ',/,&
'@    Revoir les variables syear, squant, shour, smin, ssec   ',/,&
'@                                                            ',/,&
'@    Par priorite decroissante ces variables peuvent        ',/,&
'@    etre definies dans cs_user_parameters.f90 ou le fichier ',/,&
'@    meteo ou le fichier chimie eventuel                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@    MODULE DE CHIMIE (ICHEMISTRY) DEMANDE                   ',/,&
'@                                                            ',/,&
'@    Le temps de la simulation est mal defini                ',/,&
'@    Revoir les variables syear, squant, shour, smin, ssec   ',/,&
'@                                                            ',/,&
'@    Par priorite decroissante ces variables peuvent         ',/,&
'@    etre definies dans cs_user_parameters.f90 ou le fichier ',/,&
'@    meteo ou le fichier chimie                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@    MODELE DE RAYONNEMENT (IATRA1) DEMANDE                  ',/,&
'@                                                            ',/,&
'@    Les coordonnees xlat, xlon du domaine sont mal definies ',/,&
'@                                                            ',/,&
'@    Voir cs_user_parameters.f90                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@                                                            ',/,&
'@    Les coordonnees xmet, ymet du profile meteo sont mal    ',/,&
'@    definies.                                               ',/,&
'@                                                            ',/,&
'@    Voir cs_user_parameters.f90                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (ATMOSPHERIQUE) DEMANDEE          ',/,&
'@    MODULE DE CHIMIE (ICHEMISTRY) DEMANDE                   ',/,&
'@                                                            ',/,&
'@    Les coordonnees xchem, ychem pour les profiles de       ',/,&
'@    concentration (modèle de chimie) sont mal definies.     ',/,&
'@                                                            ',/,&
'@    Voir cs_user_parameters.f90                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                RADITIVE MODEL (IATRA1)                     ',/,&
'@                                                            ',/,&
'@    The simulation time is wrong                            ',/,&
'@    Check variables syear, squant, shour, smin, ssec        ',/,&
'@                                                            ',/,&
'@    By decreasing priority these variablse can be defined   ',/,&
'@    in cs_user_parameters or the meteo file                 ',/,&
'@    or the chemistry file                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@    The simulation time is wrong                            ',/,&
'@    Check variables syear, squant, shour, smin, ssec        ',/,&
'@                                                            ',/,&
'@    By decreasing priority these variablse can be defined   ',/,&
'@    in cs_user_parameters or the meteo file                 ',/,&
'@    or the chemistry file                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                ATMOSPHERIC  MODULE                         ',/,&
'@                RADITIVE MODEL (IATRA1)                     ',/,&
'@                                                            ',/,&
'@    Wrong xlat and xlon coordinates.                        ',/,&
'@                                                            ',/,&
'@    See cs_user_parameters.f90                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1003 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC MODULE                                    ',/,&
'@                                                            ',/,&
'@    Wrong coordinates xmet, ymet for the meteo profile.     ',/,&
'@                                                            ',/,&
'@    See cs_user_parameters.f90                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1004 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@    Wrong xchem, ychem coordinates for the concentration    ',/,&
'@    profiles (chemistry model).                             ',/,&
'@                                                            ',/,&
'@    See cs_user_parameters.f90                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----

return

end subroutine atiniv
