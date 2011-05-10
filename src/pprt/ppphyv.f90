!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine ppphyv &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  ,                                                       &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )

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
!    dans usini1 si on souhaite imposer une chaleur specifique
!    CP variable (sinon: ecrasement memoire).


!  Il FAUT AVOIR PRECISE IVISLS(Numero de scalaire) = 1
!     ==================
!     dans usini1 si on souhaite une diffusivite VISCLS variable
!     pour le scalaire considere (sinon: ecrasement memoire).




! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usini1 :
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
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! ibrom            ! te ! <-- ! indicateur de remplissage de romb              !
!        !    !     !                                                !
! ia(*)            ! ia ! --- ! main integer work array                        !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp, rtpa        ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current and previous time steps)          !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfa(nfac, *)  ! ra ! <-- ! physical properties at interior face centers   !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
! coefa, coefb     ! ra ! <-- ! boundary conditions                            !
!  (nfabor, *)     !    !     !                                                !
! w1...8(ncelet    ! tr ! --- ! tableau de travail                             !
! ra(*)            ! ra ! --- ! main real work array                           !
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

integer          idbia0 , idbra0
integer          nvar   , nscal

integer          ibrom
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision w1(ncelet),w2(ncelet),w3(ncelet),w4(ncelet)
double precision w5(ncelet),w6(ncelet),w7(ncelet),w8(ncelet)
double precision ra(*)

! Local variables

integer          idebia, idebra , ifinia , ifinra
integer          if3max, iw9    , iw10

!===============================================================================

!===============================================================================
! 1. INITIALISATIONS
!===============================================================================

! --- Initialisation memoire

idebia = idbia0
idebra = idbra0


!===============================================================================
! 2. AIGUILLAGE VERS LE MODELE ADEQUAT
!===============================================================================

! ---> Flamme de diffusion chimie 3 points

  if ( ippmod(icod3p).ge.0 ) then

    call d3pphy                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )

  endif

! ---> Flamme de diffusion chimie equilibre

!        IF ( IPPMOD(ICODEQ).GE.0 )


! ---> Flamme de premelange : Modele EBU

  if ( ippmod(icoebu).ge.0 ) then

    call ebuphy                                                   &
    !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )
  endif

! ---> Flamme de premelange : Modele BML

!        IF ( IPPMOD(ICOBML).GE.0 )
!     &     CALL BMLPHY

! ---> Flamme de premelange : Modele LWC

  if ( ippmod(icolwc).ge.0 ) then

     call lwcphy                                                  &
     !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa  , propce , propfa , propfb ,           &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3    , w4     ,                             &
   w5     , w6     , w7    , w8     ,                             &
   ra     )
  endif

! ---> Flamme charbon pulverise

  if ( ippmod(icp3pl).ge.0 ) then

    ifinia = idebia

    if3max = idebra
    iw9    = if3max + ncelet
    iw10   = iw9    + ncelet
    ifinra = iw10   + ncelet
    call rasize('pppphy',ifinra)

     call cpphyv                                                  &
     !==========
 ( ifinia , ifinra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   ra(if3max),                                                    &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra(iw9), ra(iw10),                                             &
   ra     )

   endif

! ---> Flamme charbon pulverise couplee Transport Lagrangien
!      des particules de charbon

  if ( ippmod(icpl3c).ge.0 ) then

     call cplphy                                                  &
     !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )

   endif

! ---> Flamme fuel

  if ( ippmod(icfuel).ge.0 ) then

     call fuphyv                                                  &
     !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )

   endif

! ---> Compressible

  if ( ippmod(icompf).ge.0 ) then

     call cfphyv                                                  &
     !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     ,                                     &
   ra     )

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
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )

endif

! ---> Aerorefrigerants

if ( ippmod(iaeros).ge.0 ) then

   call ctphyv                                                    &
   !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )

endif

! ---> Atmospheric Flows

if ( ippmod(iatmos).ge.1 ) then

   call atphyv                                                    &
   !==========
 ( idebia , idebra ,                                              &
   nvar   , nscal  ,                                              &
   ibrom  , ia(iizfpp) ,                                          &
   ia     ,                                                       &
   dt     , rtp    , rtpa   ,                                     &
   propce , propfa , propfb ,                                     &
   coefa  , coefb  ,                                              &
   w1     , w2     , w3     , w4     ,                            &
   w5     , w6     , w7     , w8     ,                            &
   ra     )

endif

!----
! FIN
!----

return
end subroutine
