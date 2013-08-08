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

subroutine ppiniv &
!================

 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb , coefa  , coefb  )

!===============================================================================
! FONCTION :
! --------

! INITIALISATION DES VARIABLES DE CALCUL POUR
!      LA PHYSIQUE PARTICULIERE

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

! Les proprietes physiques sont accessibles dans le tableau
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
! coefa coefb      ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
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
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh

!===============================================================================

implicit none

integer          nvar   , nscal


double precision dt(ncelet), rtp(ncelet,*), propce(ncelet,*)
double precision propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables



!===============================================================================

!===============================================================================
! 1.  INITIALISATION VARIABLES LOCALES
!===============================================================================


!===============================================================================
! 2. AIGUILLAGE VERS LE MODELE ADEQUAT
!===============================================================================


! ---> Combustion gaz
!      Flamme de diffusion : chimie 3 points

 if ( ippmod(icod3p).ge.0 ) then
  call d3pini                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb )
  endif

! ---> Combustion gaz
!      Flamme de premelange : modele EBU

 if ( ippmod(icoebu).ge.0 ) then
  call ebuini                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb , coefa  , coefb  )
endif

! ---> Combustion gaz
!      Flamme de premelange : modele LWC

 if ( ippmod(icolwc).ge.0 ) then
  call lwcini                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb , coefa  , coefb  )
endif

! ---> Combustion charbon pulverise

if ( ippmod(iccoal).ge.0 ) then
  call cs_coal_varini                                             &
  !==================
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb , coefa  , coefb  )
endif


! ---> Combustion charbon pulverise couples Lagrangien

if (ippmod(icpl3c).ge.0) then
  call cplini(nvar, nscal, dt, rtp, coefa, coefb)
  !==========
endif

! ---> Combustion fuel

if  (ippmod(icfuel).ge.0 ) then
  call cs_fuel_varini                                             &
  !==================
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb , coefa  , coefb  )
endif

! ---> Compressible

if ( ippmod(icompf).ge.0 ) then
  call cfiniv                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb , coefa  , coefb  )
endif

! ---> Version electrique
!      Effet Joule
!      Conduction ionique

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then

  call eliniv                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb , coefa  , coefb  )

endif

! ---> Ecoulements atmospheriques

if ( ippmod(iatmos).ge.0 ) then

  call atiniv                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb )

endif

! ---> Ecoulements aerorefrigerants

if ( ippmod(iaeros).ge.0 ) then

  call ctiniv                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   dt     , rtp    , propce , propfb , coefa  , coefb  )

endif

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
