!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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
   dt     , rtp    , propce , propfb )

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

! Les proprietes physiaues sont accessibles dans le tableau
!     PROPCE (prop au centre), PROPFB (prop aux faces de bord)
!     Ainsi,
!      PROPCE(IEL,IPPROC(IROM  )) designe ROM   (IEL)
!      PROPCE(IEL,IPPROC(IVISCL)) designe VISCL (IEL)
!      PROPCE(IEL,IPPROC(ICP   )) designe CP    (IEL)
!      PROPCE(IEL,IPPROC(IVISLS(ISCAL))) designe VISLS (IEL ,ISCAL)

!      PROPFB(IFAC,IPPROB(IROM  )) designe ROMB  (IFAC)

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
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
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
use dimens, only: ndimfb
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

!===============================================================================

implicit none

integer          nvar   , nscal


double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfb(ndimfb,*)

! Local variables

integer          imode, iel
double precision d2s3
double precision zent,xuent,xvent,xkent,xeent,tpent,qvent,ncent

!===============================================================================
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

        rtp(iel,iu) = xuent
        rtp(iel,iv) = xvent
        rtp(iel,iw) = 0.d0

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
  dt     , rtp    , propce , propfb )

!----
! FORMATS
!----

!----
! FIN
!----

return

end subroutine atiniv
