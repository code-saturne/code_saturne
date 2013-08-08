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

subroutine ppphyv &
!================

 ( nvar   , nscal  ,                                              &
   ibrom  ,                                                       &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  )

!===============================================================================
! FONCTION :
! --------

! ROUTINE PHYSIQUE PARTICULIERE : REMPLISSAGE DES VARIABLES PHYSIQUES
! ATTENTION :
! =========





! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique ppvist devra etre creee)


!  Il FAUT AVOIR PRECISE ICP = 1
!     ==================
!    si on souhaite imposer une chaleur specifique
!    CP variable (sinon: ecrasement memoire).


!  Il FAUT AVOIR PRECISE IVISLS(Numero de scalaire) = 1
!     ==================
!     i on souhaite une diffusivite VISCLS variable
!     pour le scalaire considere (sinon: ecrasement memoire).




! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usipsu :
!             . la masse volumique (initialisee a RO0)
!             . la viscosite       (initialisee a VISCL0)
!      - dans usppiv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)


! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!        !    !     !                                                !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
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
use entsor
use pointe
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          nvar   , nscal

integer          ibrom

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)

! Local variables


!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! --- Initialisation memoire



!===============================================================================
! 2. AIGUILLAGE VERS LE MODELE ADEQUAT
!===============================================================================

! ---> Flamme de diffusion chimie 3 points

  if ( ippmod(icod3p).ge.0 ) then

    call d3pphy                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  )

  endif

! ---> Flamme de diffusion chimie equilibre

!        IF ( IPPMOD(ICODEQ).GE.0 )


! ---> Flamme de premelange : Modele EBU

  if ( ippmod(icoebu).ge.0 ) then

    call ebuphy                                                   &
    !==========
 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  )
  endif

! ---> Flamme de premelange : Modele BML

!        IF ( IPPMOD(ICOBML).GE.0 )
!     &     CALL BMLPHY

! ---> Flamme de premelange : Modele LWC

  if ( ippmod(icolwc).ge.0 ) then

     call lwcphy                                                  &
     !==========
 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa  , propce , propfb ,                    &
   coefa  , coefb  )
  endif

! ---> Flamme charbon pulverise

   if ( ippmod(iccoal).ge.0 ) then

     call cs_coal_physprop                                        &
     !====================
 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb   )

   endif


! ---> Flamme charbon pulverise couplee Transport Lagrangien
!      des particules de charbon

  if ( ippmod(icpl3c).ge.0 ) then

     call cplphy                                                  &
     !==========
 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  )

   endif

! ---> Flamme fuel

  if ( ippmod(icfuel).ge.0 ) then

     call cs_fuel_physprop                                        &
     !====================
 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  )

   endif

! ---> Compressible

  if ( ippmod(icompf).ge.0 ) then

     call cfphyv                                                  &
     !==========
 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  )

   endif

! ---> Physique particuliere : Versions electriques
!          Effet Joule
!          Arc electrique
!          Conduction ionique

if ( ippmod(ieljou).ge.1 .or.                                     &
     ippmod(ielarc).ge.1 .or.                                     &
     ippmod(ielion).ge.1       ) then

!     En Joule, on impose a l'utilisateur de programmer ses lois
!        sur les proprietes (masse volumique , ...)
!        Des exemples physiques sont fournis dans uselph.
!     En arc electrique, on lit un fichier de donnees et on interpole.

  call elphyv                                                     &
  !==========
 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  )

endif

! ---> Aerorefrigerants

if ( ippmod(iaeros).ge.0 ) then

   call ctphyv                                                    &
   !==========
 ( nvar   , nscal  ,                                              &
   ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   , propce , propfb ,                   &
   coefa  , coefb  )

endif

! ---> Atmospheric Flows

if ( ippmod(iatmos).ge.1 ) then

   call atphyv                                                    &
   !==========
 ( ibrom  , izfppp ,                                              &
   dt     , rtp    , rtpa   ,                                     &
   propce , propfb ,                                              &
   coefa  , coefb  )

endif

!----
! FIN
!----

return
end subroutine
