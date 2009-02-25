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

                  subroutine usray1
!================



!===============================================================================
!  FONCTION  :
!  ---------


!   ============================================================
!   ROUTINE UTILISATEUR : INITIALISATION DES COMMONS RAYONNEMENT
!   ============================================================



!  CALCUL DES FLUX ET DU TERME SOURCE RADIATIFS

!  METHODE DES VOLUMES FINIS.


! 1/ DONNEES DES LUMINANCES ENTRANTES AUX LIMITES DU DOMAINE
!       (C    .L : REFLEXION ET EMISSION ISOTROPE)

!                              ->  ->           ->
! 2/ CALCUL DE LA LUMINANCE L( X , S ) AU POINT X

!                                   D L
!    PAR RESOLUTION DE L'EQUATION : --- = -TK.L +TS
!                                   D S
!                       ->                o
!    OU ENCORE : DIV (L.S ) = -TK.L + TK.L

!                                 ->   /    ->  ->  ->
! 3/ CALCUL DES DENSITES DE FLUX  Q = /  L( X , S ).S DOMEGA
!                                    /4.PI

!                                      /    ->  ->
!        ET DE L'ABSORPTION       SA= /  L( X , S ).  DOMEGA
!                                    /4.PI

!    PAR INTEGRATION DES LUMINANCES SUR LES ANGLES SOLIDES.

!    N . B : CA SERT A CALCULER LE TAUX D'ECHAUFFEMENT
!    -----
!                                      /    ->  ->  ->  ->
! 4/ CALCUL DU FLUX INCIDENT QINCID = /  L( X , S ).S . N DOMEGA
!                                    /->->
!       ->                          / S.N >0
!       N NORMALE FLUIDE VERS PAROI




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
include "entsor.h"
include "optcal.h"
include "cstphy.h"
include "parall.h"
include "period.h"
include "radiat.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

integer          iphas
character*2      num


!===============================================================================



! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN



!===============================================================================
! 1. UTILISATION DU MODULE DE TRANSFERTS RADIATIFS


!     DANS LE CAS DES PHYSIQUES PARTICULIERES (COMBUSTION/CHARBON/
!                                                            ELEC/FIOUL)

!         LES RENSEIGNEMENTS NE DOIVENT PAS ETRE FOURNIS
!                            ==============
!         (ils sont completes dans le fichiers de donnees)

!     DANS LES AUTRES CAS
!         LES RENSEIGNEMENTS SUIVANTS SONT OBLIGATOIRES



!===============================================================================

!--> Par defaut, il n'y a qu'une seule phase traitee.

iphas = 1


if( ippmod(iphpar).le.1 ) then

!-->  IRAYON = 0 : PAS DE TRANSFERTS RADIATIFS (PAR DEFAUT)
!            = 1 : TRANSFERTS RADIATIFS, METHODE DES ORDONNEES DISCRETES
!            = 2 : TRANSFERTS RADIATIFS, APPROXIMATION P-1

  irayon(iphas) = 1

endif

!===============================================================================
! 2. PARAMETRES DU MODULE DE TRANSFERTS RADIATIFS
!===============================================================================

!-->  INDICATEUR SUITE DE CALCUL (LECTURE DU FICHIER SUITE DE RAYONNEMENT)
!     (0      : PAS DE LECTURE D'UN FICHIER SUITE DE RAYONNEMENT
!      1      : RELECTURE D'UN FICHIER SUITE DE RAYONNEMENT
!      ISUITE : RELECTURE D'UN FICHIER SUITE DE RAYONNEMENT SI LE CALCUL FLUIDE EST
!               AUSSI UNE SUITE )

isuird = isuite

!-->  FREQUENCE DE PASSAGE DANS LE MODULE DE RAYONNEMENT

nfreqr = 1

!-->  NOMBRE DE DIRECTIONS : 32 OU 128 (UTILE UNIQUEMENT SI IRAYON=1)

ndirec = 32

!-->  INITIALISATION DU MODE DE CALCUL DU TERME SOURCE RADIATIF EXPLICITE
!     IDIVER = 0 => CALCUL SEMI-ANALYTIQUE
!     IDIVER = 1 => CALCUL CONSERVATIF
!     IDIVER = 2 => CALCUL SEMI-ANALYTIQUE CORRIGE POUR ETRE CONSERVATIF
!     (EN RAYONNEMENT TRANSPARENT, LE CHOIX EST SANS INFLUENCE)

idiver = 2

!--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS TEMPERATURE DE PAROI

iimpar = 1

!--> NIVEAU D'AFFICHAGE (0,1,2) DES RENSEIGNEMENTS SOLVEUR

iimlum = 0

!===============================================================================
! 3. RENSEIGNEMENTS POUR LE POST-PROCESSING
!===============================================================================

! ATTENTION : les grandeurs physiques liees au milieu semi-transparent,
!             controlees par le mot-cle IRAYVP, respectent la frequence
!             de sortie chronologique imposée via NTCHR (dans USINI1.F),
!             alors que les grandeurs liees a la frontiere ne sont
!             disponibles que pour le dernier pas de temps.
!             En cas de besoin de visualisation de l'evolution en temps
!             de grandeurs surfaciques contacter l'equipe de
!             developpement.

do iphas = 1, nphas
  WRITE(NUM,'(I1)') IPHAS

!===============================================================================
! 3.1 VARIABLE DU MILIEU SEMI-TRANSPARENT
!===============================================================================

!=======================================================================
!    * IL FAUT METTRE LA VALEUR DE IRAYVP A 1 POUR LA VISUALISATION *
!=======================================================================

!--> TERME SOURCE RADIATIF (ANALYTIQUE/CONSERVATIF/SEMI-ANALYTIQUE)

  NBRVAP(ITSRAY,IPHAS) = 'Srad_'//NUM
  irayvp(itsray,iphas) = -1

!--> VECTEUR DENSITE DE FLUX RADIATIF

  NBRVAP(IQRAYP,IPHAS) = 'Qrad_'//NUM
  irayvp(iqrayp,iphas) = -1

!--> PART DE L'ABSORPTION DANS LE TERME SOURCE RADIATIF

  NBRVAP( IABSP,IPHAS) = 'Absorp_'//NUM
  irayvp( iabsp,iphas) = -1

!--> PART DE L'EMISSION DANS LE TERME SOURCE RADIATIF

  NBRVAP( IEMIP,IPHAS) = 'Emiss_'//NUM
  irayvp( iemip,iphas) = -1

!--> COEFFICIENT D'ABSORPTION DU MILIEU SEMI-TRANSPARENT

  NBRVAP(  ICAKP,IPHAS) = 'CoefAb_'//NUM
  irayvp(  icakp,iphas) = -1

!===============================================================================
! 3.2 VARIABLE SUR LES FRONTIERES DE TYPE PAROI DU DOMAINE DE CALCUL
!===============================================================================

!=======================================================================
!    * IL FAUT METTRE LA VALEUR DE IRAYVF A 1 POUR LA VISUALISATION *
!=======================================================================

!--> TEMPERATURE DES FACES FRONTIERES DE PAROI

  NBRVAF(ITPARP,IPHAS) = 'Temp_paroi_'//NUM
  irayvf(itparp,iphas) = -1

!--> FLUX INCIDENT RADIATIF RECU PAR LES FACES FRONTIERES DE PAROI

  NBRVAF(IQINCP,IPHAS) = 'Flux_incident_'//NUM
  irayvf(iqincp,iphas) = -1

!--> CONDUCTIVITE THERMIQUES DES FACES FRONTIERES DE PAROIS

  NBRVAF(IXLAMP,IPHAS) = 'Conductivite_th_'//NUM
  irayvf(ixlamp,iphas) = -1

!--> EPAISSEUR DES FACES FRONTIERES DE PAROIS

  NBRVAF( IEPAP,IPHAS) = 'Epaisseur_'//NUM
  irayvf( iepap,iphas) = -1

!--> EMISSIVITE DES FACES FRONTIERES DE PAROIS

  NBRVAF( IEPSP,IPHAS) = 'Emissivite_'//NUM
  irayvf( iepsp,iphas) = -1

!--> FLUX NET RADIATIF AUX FACES FRONTIERES DE PAROIS

  NBRVAF(IFNETP,IPHAS) = 'Flux_net_'//NUM
  irayvf(ifnetp,iphas) = -1

!--> FLUX CONVECTIF AUX FACES FRONTIERES DE PAROIS

  NBRVAF(IFCONP,IPHAS) = 'Flux_convectif_'//NUM
  irayvf(ifconp,iphas) = -1

!--> COEFFICIENT D'ECHANGE CONVECTIF AUX FACES FRONTIERES DE PAROIS

  NBRVAF(IHCONP,IPHAS) = 'Coef_ech_convectif_'//NUM
  irayvf(ihconp,iphas) = -1

!===============================================================================

enddo

return

end
