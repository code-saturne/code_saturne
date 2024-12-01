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

subroutine cs_coal_param

!===============================================================================
!  FONCTION  :
!  ---------

!       INIT DES OPTIONS DES VARIABLES TRANSPORTEES
!                     ET DES VARIABLES ALGEBRIQUES
!               COMBUSTION CHARBON PULVERISE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use ppcpfu
use cs_coal_incl
use field
use cs_c_bindings

!===============================================================================

implicit none

integer          ii , jj , iok
integer          icha , isc, krvarfl, kscacp

double precision wmolme, turb_schmidt

type(var_cal_opt) :: vcopt

!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! --> Nature des scalaires transportes

call field_get_key_id("is_temperature", kscacp)

do isc = 1, nscapp
  call field_set_key_int(ivarfl(isca(iscapp(isc))), kscacp, 0)
enddo

! ---- On resout en enthalpie avec un CP constant (Cf. cpvarp)

itherm = 2
call field_set_key_int(ivarfl(isca(iscalt)), kscacp, 0)

call field_get_key_id("variance_dissipation", krvarfl)

! --> Donnees physiques ou numeriques propres aux scalaires CP

do isc = 1, nscapp

  jj = iscapp(isc)

  if (jj .ne. iscalt .and. iscavr(jj).le.0) then

    ! En combustion on considere que la viscosite turbulente domine
    ! ON S'INTERDIT DONC LE CALCUL DES FLAMMES LAMINAIRES AVEC Le =/= 1

    call field_set_key_double(ivarfl(isca(jj)), kvisl0, viscl0)

  endif

! ------ Schmidt ou Prandtl turbulent

  turb_schmidt = 0.7d0
  call field_set_key_double(ivarfl(isca(jj)), ksigmas, turb_schmidt)

!   - Interface code_saturne:
!     ======================

  ii = isca(iscapp(isc))

  call field_get_key_struct_var_cal_opt(ivarfl(ii), vcopt)

  if (vcopt%isstpc == -999) then

    vcopt%blencv = 0.d0
    vcopt%ischcv = 1
    vcopt%isstpc = 0
    vcopt%ircflu = 0

    call field_set_key_struct_var_cal_opt(ivarfl(ii), vcopt)

  endif

enddo

!===============================================================================
! 2. INFORMATIONS COMPLEMENTAIRES
!===============================================================================

! ---- Calcul de RO0 a partir de T0 et P0
!        (loi des gaz parfaits applliquee a l'air)

!    On initialise RO0 avec l'oxydant 1 qui est sense etre
!    l'oxydant majoritaire

wmolme = ( wmole(io2) *oxyo2(1) +wmole(in2) *oxyn2(1)             &
          +wmole(ih2o)*oxyh2o(1)+wmole(ico2)*oxyco2(1))           &
        /(oxyo2(1)+oxyn2(1)+oxyh2o(1)+oxyco2(1))

ro0 = p0*wmolme / (cs_physical_constants_r*t0)

! ---- Initialisation pour la masse volumique du coke

do icha = 1, ncharb
  rhock(icha) = rho0ch(icha)
enddo

! ---> Masse volumique variable et viscosite constante (pour les suites)
irovar = 1
ivivar = 0

!===============================================================================
! 4. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
!===============================================================================

iok = 0
call cs_coal_verify (iok)

if (iok.gt.0) then
  write(nfecra,9999)iok
  call csexit (1)
else
  write(nfecra,9998)
endif

!--------
! Formats
!--------

 9998 format(                                                     &
'',                                                             /,&
' No error detected during coal data verification.',            /)
 9999 format(                                                     &
'@',/,&
'@',                                                           /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                           /,&
'@ @@ ATTENTION : stop in setup/data input.',                  /,&
'@    =========',                                              /,&
'@  The coal model parameters are incoherent or incomplete.',  /,&
'@',                                                           /,&
'@  The computation will not be run (', i10,' errors).',       /,&
'@',                                                           /,&
'@  Check messages above for detailed information.',           /,&
'@',                                                           /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                           /)

!----
! End
!----

return
end subroutine
