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

!> \file cscini.f90
!> \brief Initialization of main variables for code_saturne / code_saturne
!> coupling.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!> \param[in]     nvar          total number of variables
!______________________________________________________________________________

subroutine cscini &
 ( nvar  )

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
use period
use albase
use cplsat
use turbomachinery

!===============================================================================

implicit none

! Arguments

integer          nvar

! Local variables

integer          numcpl
integer          ialemx , nvcpmx, ifcpmx

!===============================================================================


do numcpl = 1, nbrcpl

  ! L'interpolation face/face doit être définie pour tous les couplages
  ! de manière identique.

  call mxicpl(numcpl, ifaccp, ifcpmx)
  !==========

  ifaccp = ifcpmx

  ! On vérifie si l'une des instances est en résolution en repère relatif

  call mxicpl(numcpl, icorio, icormx(numcpl))
  !==========

  ! De la même manière, si l'on a une approche ALE sur l'un des
  ! maillages, on doit mettre à jour la localisation.

  call mxicpl(numcpl, iale  , ialemx)
  !==========

  ! Si on est en turbomachine avec maillages glissant, on doit aussi
  ! mettre à jour la localisation

  if (ialemx.eq.1.or.iturbo.eq.2) then
    imajcp(numcpl) = 1
  else
    imajcp(numcpl) = 0
  endif

  ! Détermination du nombre de variables couplées entre les deux
  ! instances du couplage NUMCPL. Toutes les variables d'une instance
  ! sont couplées, SAUF dans le cas de l'ALE où la vitesse de maillage
  ! ne sera pas couplée.
  ! Il faudrait faire quelque en revanche pour les physiques particulières.

  if (iale.eq.0) then
    nvarcp(numcpl) = nvar
  else
    nvarcp(numcpl) = nvar - 3
  endif

  ! Nombre total de variable envoyées: max des variables de chaque
  ! exécutable

  call mxicpl(numcpl, nvarcp(numcpl), nvcpmx)
  !==========

  nvarto(numcpl) = nvcpmx

  ! Cohérence des modèles de turbulence entre chaque instance de CS ;
  ! pour l'instant, on ne traite que les cas de couplage entre
  ! modeles RANS et laminaires, sauf pour le modele v2f (dans ce cas
  ! il n'y a que du couplage mono-modele)

  call tbicpl(numcpl, 1, 1, iturb, iturcp(numcpl))
  !==========

  if (iturb.eq.50.and.iturcp(numcpl).ne.50) then
    write(nfecra,1000) numcpl
    call csexit(1)
    !==========
  elseif (iturb.eq.51.and.iturcp(numcpl).ne.51) then
    write(nfecra,1002) numcpl
    call csexit(1)
    !==========
  elseif (itytur.eq.4.and.                               &
       iturcp(numcpl)/10.ne.4) then
    write(nfecra,1001) numcpl
    call csexit(1)
    !==========
  endif

  ! Cohérence des referentiels de resolution

  if (icorio.ne.icormx(numcpl)) then
    write(nfecra,1100) numcpl
    call csexit(1)
    !==========
  endif

enddo

!--------
! FORMAT
!--------

 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES MODELES DE TURBULENCE POUR LE COUPLAGE ' ,I10        ,/,&
'@    SONT DIFFERENTS ALORS QUE L UN DES MODELES EST LE       ',/,&
'@    V2F PHI_FBAR. CE CAS DE FIGURE N''EST PAS PRIS          ',/,&
'@    EN COMPTE POUR LE MOMENT.                               ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usipph (cs_user_parameters.f90)                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES MODELES DE TURBULENCE POUR LE COUPLAGE ' ,I10        ,/,&
'@    SONT DIFFERENTS ALORS QUE L UN DES MODELES EST LE       ',/,&
'@    V2F BL-V2/K. CE CAS DE FIGURE N''EST PAS PRIS           ',/,&
'@    EN COMPTE POUR LE MOMENT.                               ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usipph (cs_user_parameters.f90)                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE COUPLAGE ', I10, ' EST UN COUPLAGE RANS/LES.         ',/,&
'@    CE CAS DE FIGURE N''EST PAS PRIS EN COMPTE POUR         ',/,&
'@    LE MOMENT.                                              ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usipph (cs_user_parameters.f90)                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1100 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES REFERENTIEL DE RESOLUTION POUR LE COUPLAGE ' ,I10    ,/,&
'@    SONT DIFFERENTS. CE CAS DE FIGURE N''EST PAS PRIS       ',/,&
'@    EN COMPTE.                                              ',/,&
'@    UTILISER PLUTOT UN MODELE TURBOMACHINE.                 ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usipph (cs_user_parameters.f90) ou definir un    ',/,&
'@    rotor de turbomachine (cs_user_turbomachinery.f90)      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine
