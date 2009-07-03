!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

subroutine cplpro &
!================

 ( ipropp , ipppst )

!===============================================================================
!  FONCTION  :
!  ---------


!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!      INIT DES POSITIONS DES VARIABLES D'ETAT
!        (DANS VECTEURS PROPCE, PROPFA, PROPFB)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ipropp           ! e  ! ->  ! numero de la derniere case utlisee             !
!                  !    !     ! dans ipproc, ipprob, ipprof                    !
! ipppst           ! e  ! <-- ! pointeur indiquant le rang de la               !
!                  !    !     !  derniere grandeur definie aux                 !
!                  !    !     !  cellules (rtp,propce...) pour le              !
!                  !    !     !  post traitement                               !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "dimens.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "cstnum.h"
include "ppppar.h"
include "ppthch.h"
include "coincl.h"
include "cpincl.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer       ipropp, ipppst

! VARIABLES LOCALES

integer       iprop, ige

!===============================================================================


! ---> Definition des pointeurs relatifs aux variables d'etat

iprop = ipropp

!    Phase continue
iprop   = iprop + 1
itemp1  = iprop
do ige = 1, (ngaze-2*ncharb)
! ---- Cf. definition de NGAZE dans cplecd.F
  iprop     = iprop + 1
  iym1(ige) = iprop
enddo
iprop = iprop + 1
immel = iprop

! ---- Nb de variables algebriques (ou d'etat)
!         propre a la physique particuliere NSALPP
!         total NSALTO

nsalpp = iprop - ipropp
nsalto = iprop


! ---> Positionnement dans le tableau PROPCE

iprop         = nproce

!    Phase continue (melange gazeux)
iprop           = iprop + 1
ipproc(itemp1)  = iprop
ipppst          = ipppst + 1
ipppro(iprop)   = ipppst

do ige = 1, (ngaze-2*ncharb)
! ---- Cf. definition de NGAZE dans cplecd.F
  iprop             = iprop + 1
  ipproc(iym1(ige)) = iprop
  ipppst            = ipppst + 1
  ipppro(iprop)     = ipppst
enddo

iprop         = iprop + 1
ipproc(immel) = iprop
 ipppst       = ipppst + 1
ipppro(iprop) = ipppst

nproce = iprop


! ---> Positionnement dans le tableau PROPFB
!      Au centre des faces de bord

iprop = nprofb
nprofb = iprop

! ---> Positionnement dans le tableau PROPFA
!      Au centre des faces internes (flux de masse)

iprop = nprofa
nprofa = iprop

return
end
