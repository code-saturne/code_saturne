!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine elveri &
!================

 ( iok    )

!===============================================================================
!  FONCTION  :
!  ---------

! VERIFICATION DES PARAMETRES DE CALCUL POUR LE MODULE ELECTRIQUE
!     APRES INTERVENTION UTILISATEUR (COMMONS)

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
use elincl
use radiat
use field

!===============================================================================

implicit none

! Arguments

integer          iok

! Local variables

integer          imono , ntf, ifcvsl

!===============================================================================

! --> Verifications tardives ...
if( ippmod(ielarc).ne.-1.and.ippmod(ielarc).ne.2) then
  write(nfecra,1000)ippmod(ielarc)
  iok = iok + 1
endif
if( ippmod(ieljou).ne.-1.and.ippmod(ieljou).ne.1.and.             &
                             ippmod(ieljou).ne.2.and.             &
                             ippmod(ieljou).ne.3.and.             &
                             ippmod(ieljou).ne.4     ) then
  write(nfecra,1001)ippmod(ieljou)
  iok = iok + 1
endif
if( ippmod(ielion).ne.-1) then
  write(nfecra,1002)ippmod(ielion)
  iok = iok + 1
endif

!===============================================================================
! 2. OPTIONS DU CALCUL : TABLEAUX DE ppincl.h
!                                    elincl.h
!===============================================================================

! --> Coefficient de relaxation de la masse volumique
if( srrom.lt.0.d0 .or. srrom.ge.1.d0) then
  WRITE(NFECRA,2000)'SRROM ', SRROM
  iok = iok + 1
endif

! --> Recalage des variables electriques
if( ielcor.ne.0 .and. ielcor.ne.1) then
  WRITE(NFECRA,2001)'IELCOR', IELCOR
  iok = iok + 1
endif

! --> Si on recale
if( ielcor.eq.1) then

!   -- en arc
  if(ippmod(ielarc).ge.1) then
!        Intensite de courant imposee
    if( couimp.lt.0.d0 ) then
      WRITE(NFECRA,2002)'COUIMP', COUIMP
      iok = iok + 1
    endif
!        Difference de potentiel initiale imposee
    if( dpot  .lt.0.d0 ) then
      WRITE(NFECRA,2002)'DPOT  ', DPOT
      iok = iok + 1
    endif
  endif

!   -- en Joule
  if( ippmod(ieljou).ge.1) then
!         Puissance dissipee imposee
    if( puisim.lt.0.d0 ) then
      WRITE(NFECRA,2002)'PUISIM', PUISIM
      iok = iok + 1
    endif
!         Coefficient multiplicatif correcteur initial
    if( coejou.lt.0.d0 ) then
      WRITE(NFECRA,2002)'COEJOU', COEJOU
      iok = iok + 1
    endif
!         IF( ABS(COEJOU-1.D0).GT.EPZERO ) THEN
!           WRITE(NFECRA,2003)'COEJOU', COEJOU
!           IOK = IOK + 1
!         ENDIF
!        Difference de potentiel initiale imposee
    if( dpot  .lt.0.d0 ) then
      WRITE(NFECRA,2002)'DPOT  ', DPOT
      iok = iok + 1
    endif
  endif

endif

! --> Rho et viscosite variables
if(irovar.ne.1) then
  write(nfecra,2010)irovar
  iok = iok + 1
endif
if(ivivar.ne.1) then
  write(nfecra,2011)ivivar
  iok = iok + 1
endif

! --> Cp et visls (iscalt) variables
if(icp.le.0) then
  write(nfecra,2012)'       icp',    icp,               &
                    '       icp','la chaleur massique     ',   &
                    '       icp'
  iok = iok + 1
endif

call field_get_key_int (ivarfl(isca(iscalt)), kivisl, ifcvsl)
if (ifcvsl.lt.0) then
  write(nfecra,2013) ifcvsl
  iok = iok + 1
endif

! --> Conductivites electriques variables
!     Pour le potentiel reel
call field_get_key_int (ivarfl(isca(ipotr)), kivisl, ifcvsl)
if (ifcvsl.lt.0) then
  write(nfecra,2014)'reel', ifcvsl, 'ipotr'
  iok = iok + 1
endif
!     Pour le potentiel imaginaire (Joule)
if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4 ) then
  call field_get_key_int (ivarfl(isca(ipoti)), kivisl, ifcvsl)
  if (ifcvsl.lt.0) then
    write(nfecra,2014)'imaginaire', ifcvsl, 'ipoti'
    iok = iok + 1
  endif
endif

! --> verif option sur les transfos

if(ippmod(ieljou).eq.3 .or. ippmod(ieljou).eq.4 ) then

  imono = 1
  do ntf=1,nbtrf

    if ( ibrpr(ntf) .ne. 2 .or.                                   &
         ibrsec(ntf).ne. 2       ) then
      imono = 0
    endif
  enddo

  if ( imono .eq. 1 .and. ippmod(ieljou).ne.3 ) then
    write(nfecra,2020)
    iok = iok + 1
  endif

  if ( imono .eq. 0 .and. ippmod(ieljou).ne.4 ) then
    write(nfecra,2021)
    iok = iok + 1
  endif
endif

! ---> verif rayonnement/Arc electrique


if ( ippmod(ielarc).ge.1 .and. ixkabe.eq.2 .and. iirayo.gt.0 ) then
  write(nfecra,2022)
  iok = iok + 1
endif

if ( ippmod(ielarc).ge.1 .and. ixkabe.ne.1 .and. iirayo.gt.0 ) then
  write(nfecra,2023)
  iok = iok + 1
endif



!===============================================================================
! 5. FORMATS VERIFICATION
!===============================================================================

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    IPPMOD(IELARC) NE PEUT PRENDRE QUE LES VALEURS -1 ET 2  ',/,&
'@                                                            ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  IPPMOD(IELARC) = -1 : module arc electrique desactive     ',/,&
'@                 =  2 : module arc electrique 3D            ',/,&
'@                  avec equation sur le potentiel vecteur A  ',/,&
'@                                                            ',/,&
'@                  Aucune autre valeur n''est permise. En    ',/,&
'@                    particulier, la version avec calcul de B',/,&
'@                    par le theoreme d''Ampere n''est pas    ',/,&
'@                    disponible.                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usppmo.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@  IPPMOD(IELJOU) NE PEUT PRENDRE QUE LES VALEURS -1, 1 ET 2 ',/,&
'@                                                            ',/,&
'@  IL VAUT ICI ',I10                                          ,/,&
'@                                                            ',/,&
'@  IPPMOD(IELJOU) = -1 : module Joule desactive              ',/,&
'@                 =  1 : module Joule avec potentiel reel    ',/,&
'@                 =  2 : module Joule avec potentiel complexe',/,&
'@                                                            ',/,&
'@                  Aucune autre valeur n''est permise.       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usppmo.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@  IPPMOD(IELION) NE PEUT PRENDRE QUE LA VALEUR -1           ',/,&
'@                                                            ',/,&
'@  IL VAUT ICI ',I10                                          ,/,&
'@                                                            ',/,&
'@  IPPMOD(IELION) = -1 : module mobilite ionique desactive   ',/,&
'@                                                            ',/,&
'@                  Aucune autre valeur n''est permise.       ',/,&
'@                  Le module n''est pas disponible.          ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usppmo.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    ',A6,' DOIT ETRE UN REEL SUPERIEUR OU EGAL A ZERO ET    ',/,&
'@                             INFERIEUR STRICTEMENT A 1.     ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier useli1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER COMPRIS EGAL A 0 OU 1        ',/,&
'@                                                            ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier useli1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    ',A6,' DOIT ETRE UN REEL STRICTEMENT POSITIF            ',/,&
'@                                                            ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier useli1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
!2003 format(
!    &'@                                                            ',/,
!    &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!    &'@                                                            ',/,
!    &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
!    &'@    =========                                               ',/,
!    &'@                MODULE ELECTRIQUE                           ',/,
!    &'@                                                            ',/,
!    &'@    ',A6,' DOIT ETRE UN REEL EGAL A 1                       ',/,
!    &'@                                                            ',/,
!    &'@    IL VAUT ICI ',E14.5                                      ,/,
!    &'@                                                            ',/,
!    &'@  Le calcul ne peut etre execute.                           ',/,
!    &'@                                                            ',/,
!    &'@  Verifier useli1.                                          ',/,
!    &'@                                                            ',/,
!    &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!    &'@                                                            ',/)

 2010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    IROVAR DOIT ETRE UN ENTIER EGAL A 1                     ',/,&
'@                                                            ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  IROVAR = 1 indique que la masse volumique       est       ',/,&
'@    variable. Elle est donnee par une loi dans uselph ou    ',/,&
'@    par fichier de donnees.                                 ',/,&
'@  IROVAR ne doit pas etre modifie par l''utilisateur        ',/,&
'@    lorsque le module electrique est active.                ',/,&
'@                                                            ',/,&
'@  Verifier useli1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    IVIVAR DOIT ETRE UN ENTIER EGAL A 1                     ',/,&
'@                                                            ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  IVIVAR = 1 indique que la viscosite moleculaire est       ',/,&
'@    variable. Elle est donnee par une loi dans uselph ou    ',/,&
'@    par fichier de donnees.                                 ',/,&
'@  IVIVAR ne doit pas etre modifie par l''utilisateur        ',/,&
'@    lorsque le module electrique est active.                ',/,&
'@                                                            ',/,&
'@  Verifier useli1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2012 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    ',A13,      ' DOIT ETRE UN ENTIER STRICTEMENT POSITIF   ',/,&
'@                                                            ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  ',A13,      ' > 0 indique que ',A26                        ,/,&
'@    est variable. Elle est donnee par une loi dans uselph ou',/,&
'@    par fichier de donnees.                                 ',/,&
'@  ',A13,      ' ne doit pas etre modifie par                ',/,&
'@    l''utilisateur lorsque le module electrique est active. ',/,&
'@                                                            ',/,&
'@  Verifier useli1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2013 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@  Le numero de champ associe a la propriete lambda/Cp       ',/,&
'@    doit etre positif ou nul                                ',/,&
'@                                                            ',/,&
'@    Il vaut ici ', i10                                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  field_set_key_int(ivarfl(isca(iscalt)), kivisl, ...)      ',/,&
'@    ne doit pas etre appele par l''utilisateur lorsque      ',/,&
'@    le module electrique est active.                        ',/,&
'@                                                            ',/,&
'@  Verifier useli1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2014 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@  Le numero de champ associe a la conductivite electrique   ',/,&
'@    pour le potentiel ', a, ' doit etre positif ou nul      ',/,&
'@                                                            ',/,&
'@    Il vaut ici ', i10                                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  field_set_key_int(ivarfl(isca(',a,')), kivisl, ...)      ',/,&
'@    ne doit pas etre appele par l''utilisateur lorsque      ',/,&
'@    le module electrique est active.                        ',/,&
'@                                                            ',/,&
'@  Verifier useli1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    Tout vos branchements sont mono et vous avez choisi     ',/,&
'@    un Type d''ecoulement avec transport du potentiel       ',/,&
'@    imaginaire                                              ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  Verifier usppmo et votre fichier dp_elec                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2021 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    Certains de vos branchements ne sont pas mono et        ',/,&
'@    vous avez choisis un Type d''ecoulement sans            ',/,&
'@    transport du potentiel imaginaire                       ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  Verifier usppmo et votre fichier dp_elec                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2022 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    VOUS NE POUVEZ PAS ACTIVER LE MODULE RAYONNEMENT        ',/,&
'@    ET INTRODUIRE UN TERME SOURCE RADIATIF DANS LE          ',/,&
'@    FICHIER DP_ELEC                                         ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  Verifier votre fichier dp_elec ou desactiver le module    ',/,&
'@  rayonnement                                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2023 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                MODULE ELECTRIQUE                           ',/,&
'@                                                            ',/,&
'@    POUR ACTIVER LE MODULE RAYONNEMENT AVEC LE MODULE       ',/,&
'@    ARC ELECTRIQUE IL FAUT DONNER LE COEFFICIENT            ',/,&
'@    D''ABSORPTION DANS LE FICHIER DP_ELEC                   ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  Verifier votre fichier dp_elec ou desactiver le module    ',/,&
'@  rayonnement                                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!===============================================================================
! 6. SORTIE
!===============================================================================

return
end subroutine
