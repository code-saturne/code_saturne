!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2024 EDF S.A.
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

procedure() :: mxicpl, tbicpl, csexit

!===============================================================================


do numcpl = 1, nbrcpl

  ! Face to face interpolation should be defined for all couplings
  ! in the same manner.

  call mxicpl(numcpl, ifaccp, ifcpmx)

  ifaccp = ifcpmx

  ! Check if one of the instances is solved in a relative reference frame.

  call mxicpl(numcpl, icorio, icormx(numcpl))

  ! In the same manner, if one of the meshes is defored (ALE) or moves
  ! (such as with turbomachinery), the location (mesh mapping)
  ! should be updated.

  call mxicpl(numcpl, iale, ialemx)

  if (ialemx.eq.1.or.iturbo.eq.2) then
    imajcp(numcpl) = 1
  else
    imajcp(numcpl) = 0
  endif

  ! Determine the number of coupled variables betwwen instances of
  ! coupling numcpl. All variables of a givne instance are coupled,
  ! EXCEPT for ALE, where the mesh velocity is not coupled.
  ! Something should be done for specific physical models.

  if (iale.eq.0) then
    nvarcp(numcpl) = nvar
  else
    nvarcp(numcpl) = nvar - 3
  endif

  ! Total number of sent variables. max number of variables for each
  ! domain.

  call mxicpl(numcpl, nvarcp(numcpl), nvcpmx)

  nvarto(numcpl) = nvcpmx

  ! Coherencey between turbulence models between each code_saturne instance ;
  ! Currently, only couplings between RANS and laminar domains is handled,
  ! except for v2f (where both domains must use the same model).

  call tbicpl(numcpl, 1, 1, iturb, iturcp(numcpl))

  if (iturb.eq.50 .and. iturcp(numcpl).ne.50) then
    write(nfecra,1000) numcpl
    call csexit(1)
  elseif (iturb.eq.51 .and. iturcp(numcpl).ne.51) then
    write(nfecra,1002) numcpl
    call csexit(1)
  elseif (itytur.eq.4 .and. iturcp(numcpl)/10.ne.4) then
    write(nfecra,1001) numcpl
    call csexit(1)
  endif

  ! Coherency between reference frames.

  if (icorio.ne.icormx(numcpl)) then
    write(nfecra,1100) numcpl
    call csexit(1)
  endif

enddo

!--------
! Format
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
