!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine atiniv &
!================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL
!    POUR LA PHYSIQUE PARTICULIERE : ECOULEMENTS ATMOSPHERIQUES
!    PENDANT DE USINIV.F

! Cette routine est appelee en debut de calcul (suite ou non)
!     avant le debut de la boucle en temps

! Elle permet d'INITIALISER ou de MODIFIER (pour les calculs suite)
!     les variables de calcul,
!     les valeurs du pas de temps


! On dispose ici de ROM et VISCL initialises par RO0 et VISCL0
!     ou relues d'un fichier suite
! On ne dispose des variables VISCLS, CP (quand elles sont
!     definies) que si elles ont pu etre relues dans un fichier
!     suite de calcul

! LA MODIFICATION DES PROPRIETES PHYSIQUES (ROM, VISCL, VISCLS, CP)
!     SE FERA EN STANDARD DANS LE SOUS PROGRAMME PPPHYV
!     ET PAS ICI

! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! dt(ncelet)       ! tr ! <-- ! valeur du pas de temps                         !
! rtp              ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules                                    !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

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

double precision dt(ncelet), rtp(ncelet,nflown:nvar), propce(ncelet,*)

! Local variables

integer          imode, iel
double precision d2s3
double precision zent,xuent,xvent,xkent,xeent,tpent,qvent,ncent

integer k,ii, isc
double precision xcent

double precision, dimension(:,:), pointer :: vel

!===============================================================================

! Map field arrays
call field_get_val_v(ivarfl(iu), vel)

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================

d2s3 = 2.d0/3.d0

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
  if (init_at_chem.eq.1) then
    do iel = 1, ncel

      zent = xyzcen(3,iel)

      do k = 1, nespgi
        call intprf                                                         &
        !==========
        (nbchmz, nbchim,                                                    &
         zproc, tchem, espnum(1+(k-1)*nbchim*nbchmz), zent  , ttcabs, xcent )
        ! The first nespg user scalars are supposed to be chemical species
        rtp(iel,isca(idespgi(k))) = xcent
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
    do iel = 1, ncel
      do ii = 1, nesp_aer*nbin_aer + nbin_aer
        isc = (isca_chem(1) - 1) + nespg_siream + ii
        rtp(iel,isca(isc)) = dlconc0(ii)
      enddo
    enddo
  endif

endif

! Verifications
if ((iatra1.eq.1.or.ichemistry.ge.1).and.(syear.eq.-999.or.squant.eq.-999.or.shour.eq.-999&
.or.smin.eq.-999.or.ssec.eq.-999)) then
  if (iatra1.eq.1) write(nfecra,1000)
  if (ichemistry.ge.1) write(nfecra,1001)
  call csexit (1)
endif

if ((iatra1.eq.1.or.ichemistry.ge.1).and.(xlat.ge.rinfin*0.5.or.xlon.ge.rinfin*0.5)) then
  if (iatra1.eq.1) write(nfecra,1002)
  if (ichemistry.ge.1) write(nfecra,1003)
  call csexit (1)
endif

!===============================================================================
! 3. Dry atmosphere: default initialization of potential temperature
!===============================================================================

! Only if the simulation is not a restart from another one
if (isuite.eq.0) then

  if (initmeteo.eq.1) then
    if (imeteo.eq.0) then

      if (ippmod(iatmos).eq.1) then
        ! The thermal scalar is potential temperature
        do iel = 1, ncel
          rtp(iel,isca(iscalt)) = t0
        enddo
      endif

      if (ippmod(iatmos).eq.2) then
        ! The thermal scalar is liquid potential temperature
        do iel = 1, ncel
          rtp(iel,isca(iscalt)) = t0
          rtp(iel,isca(itotwt)) = 0.d0
          rtp(iel,isca(intdrp)) = 0.d0
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

          rtp(iel,ik)  = xkent
          rtp(iel,iep) = xeent

        elseif (itytur.eq.3) then

          rtp(iel,ir11) = d2s3*xkent
          rtp(iel,ir22) = d2s3*xkent
          rtp(iel,ir33) = d2s3*xkent
          rtp(iel,ir12) = 0.d0
          rtp(iel,ir13) = 0.d0
          rtp(iel,ir23) = 0.d0
          rtp(iel,iep)  = xeent

        elseif (iturb.eq.50) then

          rtp(iel,ik)   = xkent
          rtp(iel,iep)  = xeent
          rtp(iel,iphi) = d2s3
          rtp(iel,ifb)  = 0.d0

        elseif (iturb.eq.60) then

          rtp(iel,ik)   = xkent
          rtp(iel,iomg) = xeent/cmu/xkent

        elseif (iturb.eq.70) then

          rtp(iel,inusa) = cmu*xkent**2/xeent

        endif


        if (ippmod(iatmos).eq.1) then
          ! The thermal scalar is potential temperature
            call intprf                                                 &
            !==========
         (nbmett, nbmetm,                                               &
          ztmet, tmmet, tpmet, zent  , ttcabs, tpent )

            rtp(iel,isca(iscalt)) = tpent
        endif

        if (ippmod(iatmos).eq.2) then
          ! The thermal scalar is liquid potential temperature
            call intprf                                                 &
            !==========
         (nbmett, nbmetm,                                               &
          ztmet, tmmet, tpmet, zent  , ttcabs, tpent )
            rtp(iel,isca(iscalt)) = tpent

            call intprf                                                 &
            !==========
         (nbmett, nbmetm,                                               &
          ztmet, tmmet, qvmet, zent  , ttcabs, qvent )
            rtp(iel,isca(itotwt)) = qvent

            call intprf                                                 &
            !==========
         (nbmett, nbmetm,                                               &
          ztmet, tmmet, ncmet, zent  , ttcabs, ncent )
            rtp(iel,isca(intdrp)) = ncent
        endif

      enddo

    endif
  endif

endif

!===============================================================================
! 4. USER  OPTIONS
!===============================================================================

call cs_user_initialization &
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
'@    Les coordonnees xlat et xlon du domaine sont mal definies',/,&
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
'@    MODULE DE CHIMIE (ICHEMISTRY) DEMANDE                   ',/,&
'@                                                            ',/,&
'@    Les coordonnees xlat et xlon du domaine sont mal definies',/,&
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
'@    Wrong xlat and xlon coordinates                         ',/,&
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
'@      ATMOSPHERIC CHEMISTRY                                 ',/,&
'@                                                            ',/,&
'@    Wrong xlat and xlon coordinates                         ',/,&
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
