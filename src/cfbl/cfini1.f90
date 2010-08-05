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

subroutine cfini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!         INIT DES OPTIONS DES VARIABLES POUR
!              LE COMPRESSIBLE SANS CHOC
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

implicit none

!===============================================================================
! Common blocks
!===============================================================================

include "paramx.f90"
include "dimens.f90"
include "numvar.f90"
include "optcal.f90"
include "cstphy.f90"
include "entsor.f90"
include "cstnum.f90"
include "ppppar.f90"
include "ppthch.f90"
include "ppincl.f90"

!===============================================================================

! Local variables

integer          ipp , ii
integer          iphas, iok

!===============================================================================
!===============================================================================
! 0. VERIFICATION ISCALT, ISCSTH
!===============================================================================
!     L'utilisateur ne doit pas y avoir touche.

do iphas = 1, nphas
  if(iscalt(iphas).ne.-1) then
    write(nfecra,1000)iphas,iscalt(iphas)
    call csexit (1)
    !==========
  endif
enddo
do ii = 1, nscapp
  if(iscsth(iscapp(ii)).ne.-10) then
    write(nfecra,1001)ii,iscapp(ii),iscapp(ii),iscsth(iscapp(ii))
    call csexit (1)
    !==========
  endif
enddo
!===============================================================================
! 1. VARIABLES TRANSPORTEES
!===============================================================================

! 1.1 Definition des scamin et des scamax des variables transportees
! ==================================================================


do iphas = 1, nphas

  if(  (abs(scamin(irho  (iphas))+grand).gt.epzero).or.           &
       (abs(scamin(ienerg(iphas))+grand).gt.epzero).or.           &
       (abs(scamin(itempk(iphas))+grand).gt.epzero).or.           &
       (abs(scamax(irho  (iphas))-grand).gt.epzero).or.           &
       (abs(scamax(ienerg(iphas))-grand).gt.epzero).or.           &
       (abs(scamax(itempk(iphas))-grand).gt.epzero) ) then
    write(nfecra,2000)                                            &
         scamin(irho  (iphas)),scamax(irho  (iphas)),             &
         scamin(ienerg(iphas)),scamax(ienerg(iphas)),             &
         scamin(itempk(iphas)),scamax(itempk(iphas))
    call csexit (1)
  endif
!        SCAMIN(IRHO  (IPHAS))   = -GRAND
!        SCAMAX(IRHO  (IPHAS))   =  GRAND
!        SCAMIN(IENERG(IPHAS))   = -GRAND
!        SCAMAX(IENERG(IPHAS))   =  GRAND
!        SCAMIN(ITEMPK(IPHAS))   = -GRAND
!        SCAMAX(ITEMPK(IPHAS))   =  GRAND

enddo

! 1.2 Nature des scalaires transportes
! ====================================

! ---- Type de scalaire (0 passif, 1 temperature en K
!                                 -1 temperature en C
!                                  2 enthalpie en J
!                                  3 energie totale en J)
!      La distinction -1/1 sert pour le rayonnement

do iphas = 1, nphas

  iscsth(irho  (iphas)) = 0
  iscsth(ienerg(iphas)) = 3
  iscsth(itempk(iphas)) = 0

  iscalt(iphas) = ienerg(iphas)

enddo


!         - Schema convectif % schema 2ieme ordre
!           = 0 : upwind
!           = 1 : second ordre
do ii = 1, nvarmx
  blencv(ii) = 0.d0
enddo

!         Upwind necessaire pour le schema utilise


! 1.3 Variable courante : nom, sortie chrono, suivi listing, sortie hist
! ======================================================================

!     Comme pour les autres variables,
!       si l'on n'affecte pas les tableaux suivants,
!       les valeurs par defaut seront utilisees

!     NOMVAR( ) = nom de la variable
!     ICHRVR( ) = sortie chono (oui 1/non 0)
!     ILISVR( ) = suivi listing (oui 1/non 0)
!     IHISVR( ) = sortie historique (nombre de sondes et numeros)
!     si IHISVR(.,1)  = -1 sortie sur toutes les sondes

!     NB : Les 8 premiers caracteres du noms seront repris dans le
!          listing 'developpeur'

! ======================================================================

do iphas = 1, nphas

  ipp = ipprtp(isca(irho  (iphas)))
  NOMVAR(IPP)  = 'Rho'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

  ipp = ipprtp(isca(ienerg(iphas)))
  NOMVAR(IPP)  = 'EnergieT'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

  ipp = ipprtp(isca(itempk(iphas)))
  NOMVAR(IPP)  = 'Temp K'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1

enddo


!===============================================================================
! 2. PARAMETRES GLOBAUX
!===============================================================================

! --- Couplage vitesse/pression (0 : algorithme classique,
!                                1 : couplage instationnaire)
!     Uniquement en monophasique et en incompressible

if( ipucou.ne.0 ) then
  write(nfecra,3000) ipucou
  call csexit (1)
endif


! --- Estimateurs pour Navier-Stokes

!     Interdits en compressible

do iphas = 1, nphas
  if( (iescal(iespre,iphas).ne.0) .or.                            &
      (iescal(iesder,iphas).ne.0) .or.                            &
      (iescal(iescor,iphas).ne.0) .or.                            &
      (iescal(iestot,iphas).ne.0) ) then
    write(nfecra,4000)
    call csexit (1)
  endif
!       IESCAL(IESPRE,IPHAS) = 0
!       IESCAL(IESDER,IPHAS) = 0
!       IESCAL(IESCOR,IPHAS) = 0
!       IESCAL(IESTOT,IPHAS) = 0
enddo


!===============================================================================
! 3. OPTIONS DE CALCUL PAR DEFAUT
!===============================================================================

!     Pour chaque phase

do iphas = 1, nphas

! --> Conditions aux limites prenant en compte l'equilibre hydrostatique
!     (oui = 1 , non = 0)

  icfgrp(iphas) = 1


! ---> Masse volumique variable et viscosite constante (pour les suites)
  irovar(iphas) = 1
  ivivar(iphas) = 0

enddo


!===============================================================================
! 4. ON REDONNE LA MAIN A L'UTLISATEUR
!===============================================================================

call uscfx1
!==========

!===============================================================================
! 5. OPTIONS DE CALCUL OBLIGATOIRES
!     qui pourront etre remontees au dessus de uscfx1
!     selon les developpements
!===============================================================================

!     Pour chaque phase

do iphas = 1, nphas

  idiff(isca(irho(iphas))) = 1

! --> Implicitation du terme de convection de l'equation de masse
!     (oui = 1 , non = 0)
!     On choisit 0 ; c'est la seule option qui a ete testee. Elle
!       facilite le codage pour le respect du flux de masse au bord.

  iconv(isca(irho(iphas))) = 0

! --> Prise en compte de la pression predite pour resoudre Navier-Stokes
!     (oui = 1 , non = 0)

  igrdpp(iphas) = 0

! --> Prediction de pression par une equation d'evolution

!     ATTENTION   PAS ENCORE IMPLEMENTE
!========   LAISSER IPPRED(IPHAS) = 0

  ippred(iphas) = 0


enddo

!===============================================================================
! 6. VERIFICATIONS
!===============================================================================

iok = 0
do iphas = 1, nphas
  if(icfgrp(iphas).ne.0.and.icfgrp(iphas).ne.1) then
    WRITE(NFECRA,5000)IPHAS,'ICFGRP',ICFGRP(IPHAS)
    iok = 1
  endif
enddo

if (iok.ne.0) then
  call csexit (1)
endif

!--------
! FORMATS
!--------

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@  La valeur de ISCALT est renseignee automatiquement.       ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas la renseigner dans usini1, or  ',/,&
'@    elle a ete affectee comme suit :                        ',/,&
'@    ISCALT(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@  Les valeurs de ISCSTH sont renseignees automatiquement.   ',/,&
'@                                                            ',/,&
'@  L''utilisateur ne doit pas les renseigner dans usini1, or ',/,&
'@    pour le scalaire ',I10   ,' correspondant au scalaire   ',/,&
'@    physique particuliere ',I10   ,' on a                   ',/,&
'@    ISCSTH(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@  Les bornes des variables rho, energie ou temperature      ',/,&
'@    ont ete modifiees dans usini1 :                         ',/,&
'@                                                            ',/,&
'@                      SCAMIN        SCAMAX                  ',/,&
'@  rho         ',2E14.5                                       ,/,&
'@  energie     ',2E14.5                                       ,/,&
'@  temperature ',2E14.5                                       ,/,&
'@                                                            ',/,&
'@  Les bornes de ces variables ne doivent pas etre modifiees ',/,&
'@  dans usini1. On peut modifier les bornes des variables    ',/,&
'@  rho et energie dans uscfx1, mais ce n''est pas conseille. ',/,&
'@  Il est preferable de gerer les depassements éventuels     ',/,&
'@  au moyen du sous programme uscfth (arret du calcul en fin ',/,&
'@  de pas de temps en cas de depassement).                   ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@  L''option IPUCOU = ',I10                                   ,/,&
'@    n''est pas compatible avec le module compressible       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Imposer IPUCOU = 0 dans usini1.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@  Les estimateurs ne sont pas compatibles avec le module    ',/,&
'@    compressible                                            ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Imposer IESCAL(.,.) = 0 dans usini1.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHYSIQUE PARTICULIERE (COMPRESSIBLE) DEMANDEE           ',/,&
'@                                                            ',/,&
'@    PHASE ',I10                                              ,/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier uscfx1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

return
end subroutine
