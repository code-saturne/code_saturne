!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2022 EDF S.A.
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

subroutine coini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!         INIT DES OPTIONS DES VARIABLES POUR
!              POUR LA COMBUSTION
!        FLAMME DE DIFFUSION ET DE PREMELANGE
!   EN COMPLEMENT DE CE QUI A DEJA ETE FAIT DANS USINI1

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
use radiat
use field
use cs_c_bindings

!===============================================================================

implicit none

! Local variables

integer          jj, iok
integer          isc
double precision wmolme, turb_schmidt
!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! 1.1 Definition des scamin et des scamax des variables transportees
! ==================================================================

! --> Flamme de diffusion : chimie 3 points

if (ippmod(icod3p).ge.0) then

  ! ---- Variance du taux de melange
  !        Type de clipping superieur pour la variance
  !        0 pas de clipping, 1 clipping var max de fm, 2 clipping a SCAMAX
  !iclvfl(ifp2m) = 1
  iclvfl(ifp2m) = 2
  !        scamin(ifp2m) = 0.d0
  !        scamax(ifp2m) = 0.25D0

endif

! --> Flamme de diffusion : Steady laminar flamelet

if (ippmod(islfm).ge.0) then

  ! Manière de calculer la variance de fraction de mélange
  ! mode_fp2m 0: Variance transport equation(VTE)
  ! mode_fp2m 1: 2nd moment of mixutre fraction transport equation (STE)

  ! ---- Variance du taux de melange
  !        Type de clipping superieur pour la variance
  !        0 pas de clipping, 1 clipping var max de fm, 2 clipping a SCAMAX
  if (mode_fp2m .eq. 0) then
    iclvfl(ifp2m) = 1
    !        scamin(ifp2m) = 0.d0
    !        scamax(ifp2m) = 0.25D0
  endif
endif


! --> Flamme de premelange : modele LWC

if (ippmod(icolwc).ge.0) then
  iclvfl(ifp2m) = 0
  iclvfl(iyfp2m) = 0
endif

! 1.4 Donnees physiques ou numeriques propres aux scalaires COMBUSTION
! ====================================================================

do isc = 1, nscapp

  jj = iscapp(isc)

  if (iscavr(jj).le.0) then

! ---- En combustion on considere que la viscosite turbulente domine
!      ON S'INTERDIT DONC LE CALCUL DES FLAMMES LAMINAIRES AVEC Le =/= 1

    call field_set_key_double(ivarfl(isca(jj)), kvisl0, viscl0)

  endif

  ! Schmidt ou Prandtl turbulent

  turb_schmidt = 0.7d0
  call field_set_key_double(ivarfl(isca(jj)), ksigmas, turb_schmidt)

  ! Coeff dissipation des fluctuations

  rvarfl(jj) = 0.8d0

enddo

!===============================================================================
! 2. INFORMATIONS COMPLEMENTAIRES
!===============================================================================

! --> Calcul de RO0 a partir de T0 et P0

if (ippmod(icod3p).ne.-1 .or.                                   &
    ippmod(icoebu).ne.-1 .or.                                   &
    ippmod(icolwc).ne.-1) then
  wmolme = wmolg(2)
  ro0 = pther*wmolme / (cs_physical_constants_r*t0)
  roref = ro0
else if (ippmod(islfm).ne.-1) then
  ro0 = flamelet_library(FLAMELET_RHO, 1, 1, 1, 1)
  roref = ro0
endif

! On met les constantes a -GRAND pour obliger l'utilisateur a les definir
!  (ou les laisser) dans cs_user_combustion.
! --> Constante modele EBU par defaut
cebu   =-grand

! --> Constantes modele LWC par defaut
vref  =-grand
lref  =-grand
ta    =-grand
tstar =-grand

! --> Coefficient de relaxation de la masse volumique
!      RHO(n+1) = SRROM * RHO(n) + (1-SRROM) * RHO(n+1)
srrom =-grand

! --> Viscosite laminaire associee au scalaire enthalpie
!       DIFTL0 (diffusivite dynamique en kg/(m s))
diftl0      =-grand

! --> Reference temperature for fuel and oxydant (K)
tinfue = -grand
tinoxy = -grand

! --> Diffusion 3 points, tableaux HH HF THF
!           (generes dans d3pphy)

nmaxf = 9
nmaxh = 9

! ---> Masse volumique variable
irovar = 1

!===============================================================================
! 3. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

! GUI
if (ippmod(icoebu).ge.0) then
  call uicpi1(srrom, diftl0)
  cebu   = 2.5d0
else if (ippmod(icod3p).ge.0) then
  call uicpi1(srrom, diftl0)
  call uicpi2(tinoxy, tinfue)
! else if (ippmod(icolwc).ge.0) then
  !TODO no GUI yet
endif

! User subroutines
call cs_user_combustion

!===============================================================================
! 4. VERIFICATION DES DONNERS FOURNIES PAR L'UTLISATEUR
!===============================================================================

iok = 0
if (ippmod(icoebu).ge.0) then

  call ebuver(iok)
  if (iok.gt.0) then
    write(nfecra,9999)iok
    call csexit(1)
  else
    write(nfecra,9998)
  endif

else if (ippmod(icod3p).ge.0) then
  call d3pver(iok)
  if (iok.gt.0) then
    write(nfecra,9991)iok
    call csexit(1)
  else
    write(nfecra,9990)
  endif

else if (ippmod(islfm).ge.0) then
  call cs_steady_laminar_flamelet_verify(iok)
  if (iok.gt.0) then
    write(nfecra,9991)iok
    call csexit(1)
  else
    write(nfecra,9990)
  endif

else if (ippmod(icolwc).ge.0) then
  call lwcver(iok)
  if (iok.gt.0) then
    write(nfecra,9993)iok
    call csexit(1)
  else
    write(nfecra,9992)
  endif

endif

 9998 format(                                                     &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                                                    (usebu1).',/)
 9999 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier usebu1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9990 format(                                                     &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                                        (cs_user_combustion).',/)
 9991 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier cs_user_combustion.'                              ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9992 format(                                                     &
'                                                             ',/,&
' Pas d erreur detectee lors de la verification des donnees   ',/,&
'                                                    (uslwc1).',/)
 9993 format(                                                     &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier uslwc1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


return
end subroutine
