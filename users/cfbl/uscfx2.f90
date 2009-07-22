!-------------------------------------------------------------------------------

!VERS


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

                  subroutine uscfx2
!================


!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES OPTIONS DES VARIABLES POUR LE COMPRESSIBLE SANS CHOC
!   EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1


!    CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE
!    =============================================


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
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
include "ppincl.h"

!===============================================================================

integer          iphas

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if(1.eq.1) then
  write(nfecra,9000)
  call csexit (1)
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@                      MODULE COMPRESSIBLE                   ',/,&
'@                                                            ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uscfx2 DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@     Ce sous-programme utilisateur permet de definir les    ',/,&
'@       options generales. Il est indispensable.             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================
! 1. PROPRIETES PHYSIQUES
!===============================================================================

! --> Pour chaque phase

iphas = 1

! --> Conductivite thermique laminaire

!       Conductivite thermique constante : IVISLS = 0
!       Conductivite thermique variable  : IVISLS = 1

ivisls(itempk(iphas)) = 0

!       Conductivite thermique de reference :

!       VISLS0 = LAMBDA0  (conductivite thermique en W/(m K))


!       ATTENTION ! IL FAUT QUE VISLS0 SOIT STRICTEMENT POSITIF
!         (donner une valeur meme si la conductivite est variable)

visls0(itempk(iphas)) = 3.d-2

!       Si la conductivite thermique est variable, il faut donner
!       sa loi de variation dans uscfpv.F


! --> Viscosite en volume

!       Viscosite en volume de reference :

!       VISCV0 = KAPPA0  (viscosite en volume en kg/(m s))
!       IVISCV = 0 : uniforme en espace et constant en temps
!              = 1 : variable en espace et  en temps

iviscv(iphas) = 0
viscv0(iphas) = 0.d0

!       Si la Viscosite en volume est variable, il faut donner
!       sa loi de variation dans uscfpv.F


!----
! FIN
!----

return
end
