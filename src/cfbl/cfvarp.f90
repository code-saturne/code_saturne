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

subroutine cfvarp
!================


!===============================================================================
!  FONCTION  :
!  ---------

!              INIT DES POSITIONS DES VARIABLES
!            POUR LE COMPRESSIBLE SANS CHOC SELON
! REMPLISSAGE DES PARAMETRES (DEJA DEFINIS) POUR LES SCALAIRES PP

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
use ppincl

!===============================================================================

implicit none

! Local variables

integer          ii, iphas, iprop, iok, iccfth, imodif
integer          iit(1)
double precision dblpre(1)

!===============================================================================
!===============================================================================
! 1. DEFINITION DES POINTEURS
!===============================================================================


if ( ippmod(icompf).ge.0 ) then

  iprop =0

  do iphas = 1, nphas

! ---- Masse volumique
    iprop = iprop + 1
    irho(iphas) = iscapp(iprop)
!     Alias pour les C.L.
    irun(iphas) = irho(iphas)

! ---- Energie totale
    iprop = iprop + 1
    ienerg(iphas) = iscapp(iprop)
!     Alias pour les C.L.
    irunh(iphas) = ienerg(iphas)

! ---- Temperature (post)
    iprop = iprop + 1
    itempk(iphas) = iscapp(iprop)

! ---- Phase associee
    iphsca(irho  (iphas)) = iphas
    iphsca(ienerg(iphas)) = iphas
    iphsca(itempk(iphas)) = iphas

! ---- Viscosite dynamique de reference relative au scalaire IRHO
    ivisls(irho  (iphas)) = 0
    visls0(irho  (iphas)) = epzero

! ---- Viscosite dynamique de reference relative au scalaire ITEMPK
    ivisls(itempk(iphas)) = 0
    visls0(itempk(iphas)) = epzero

! ---- Initialisation par defaut de la viscosite en volume (cste)
    iviscv(iphas) = 0
    viscv0(iphas) = 0.d0


!===============================================================================
! 2. OPTIONS DE CALCUL
!===============================================================================

! --> Cv constant ou variable (par defaut : constant)
    icv(iphas) = 0
    cv0(iphas) = 0.d0

    iccfth = -1
    imodif = 0
    ii     = 1
    iit(1)    = 1
    dblpre(1) = 0.d0
    call uscfth                                                   &
    !==========
 ( ii , ii ,                                                      &
   ii , ii , ii ,                                                 &
   iccfth , imodif  , iphas   ,                                   &
   ii , ii , ii , ii ,                                            &
   iit , iit , iit ,                                              &
   dblpre , dblpre , dblpre , dblpre , dblpre , dblpre ,          &
   dblpre , dblpre ,                                              &
   dblpre , dblpre , dblpre , dblpre ,                            &
   dblpre , dblpre , dblpre )

! --> Utilisation d'un flux de masse specifique pour la vitesse

!     ATTENTION   PAS ENCORE IMPLEMENTE
!========   LAISSER IFLMAU(IPHAS) = 0

    iflmau(iphas) = 0

  enddo


!===============================================================================
! 3. ON REDONNE LA MAIN A L'UTILISATEUR
!===============================================================================

  call uscfx2
  !==========


!===============================================================================
! 4. TRAITEMENT ET VERIFICATION DES DONNEES FOURNIES PAR L'UTILISATEUR
!===============================================================================

! ---- Viscosite dynamique de reference relative au scalaire IENERG
  do iphas = 1, nphas
    if(ivisls(itempk(iphas)).gt.0 .or. icv(iphas).gt.0) then
      ivisls(ienerg(iphas)) = 1
    else
      ivisls(ienerg(iphas)) = 0
    endif

    visls0(ienerg(iphas)) = epzero

  enddo

  iok = 0

  do iphas = 1, nphas

    if(visls0(itempk(iphas)).le.0.d0) then
      write(nfecra,1000) iphas, visls0(itempk(iphas))
      iok = 1
    endif

    if(viscv0(iphas).lt.0.d0) then
      write(nfecra,2000) iphas, viscv0(iphas)
      iok = 1
    endif

  enddo

  if(iok.gt.0) call csexit (1)

endif

!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========   MODULE COMPRESSIBLE                         ',/,&
'@                                                            ',/,&
'@    LA CONDUCTIVITE THERMIQUE (PHASE ',I6,') DOIT ETRE      ',/,&
'@    UN REEL POSITIF STRICTEMENT                             ',/,&
'@    ELLE A POUR VALEUR ',E12.4                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========   MODULE COMPRESSIBLE                         ',/,&
'@                                                            ',/,&
'@    LA CONDUCTIVITE THERMIQUE (PHASE ',I6,') DOIT ETRE      ',/,&
'@                                                            ',/,&
'@    LA VISCOSITE EN VOLUME (PHASE ',I6,') DOIT ETRE         ',/,&
'@    UN REEL POSITIF                                         ',/,&
'@    ELLE A POUR VALEUR ',E12.4                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx2.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine

