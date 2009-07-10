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

!  FONCTION  :
!  ---------

! ROUTINES UTILISATEUR POUR ENTREE DES PARAMETRES DE CALCUL (COMMONS)
!  ON PASSE ICI POUR TOUS LES CALCULS

! SI L'ON UTILISE L'IHM DE CODE_SATURNE, LE FICHIER PRESENT N'EST
!   PAS INDISPENSABLE (IL COMPLETE OU PREND LE PAS SUR LES
!   PARAMETRES ENTRES PAR L'INTERFACE)

! ON TROUVERA PLUSIEURS ROUTINES DANS LE FICHIER PRESENT, CHACUNE
!   DESTINEE A RENSEIGNER CERTAINS PARAMETRES PARTICULIERS

! POUR MODIFIER LA VALEUR PAR DEFAUT DE PARAMETRES N'APPARAISSANT
!   PAS DANS LES EXEMPLES FOURNIS ICI, INTERVENIR CI-DESSOUS DANS
!   - USIPSU   POUR LES OPTIONS NUMERIQUES ET PHYSIQUES
!   - USIPES   POUR LES OPTIONS RELATIVES AUX ENTREES-SORTIES

! PAR CONVENTION "PHYSIQUE PARTICULIERE" RENVOIE AUX MODULES
!   SUIVANTS UNIQUEMENT :
!           CHARBON PULVERISE, COMBUSTION GAZ, ELECTRIQUE


!-------------------------------------------------------------------------------


!===============================================================================


subroutine usipph &
!================

 ( nphmax, nphas , iihmpu, nfecra , iturb , icp , iverif )


!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE UTILISATEUR POUR ENTREE
!   DES PARAMETRES DEPENDANT DU NOMBRE DE PHASES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nphmax           ! e  ! <-- ! nombre de phases maximum                       !
! nphas            ! e  ! <-- ! nombre de phases actives                       !
! iihmpu           ! e  ! <-- ! indique si un fichier de donnees de            !
!                  !    !     ! l'ihm est utilise (1: oui, 0: non)             !
! nfecra           ! e  ! <-- ! numero d'unite std pour impression             !
! iturb(nphmax)    ! te ! <-- ! modele de turbulence                           !
! icp  (nphmax)    ! te ! <-- ! indicateur de cp uniforme ou non               !
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
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


!     AUCUN COMMON NE DOIT APPARAITRE ICI


!===============================================================================

! Arguments

integer nphmax, nphas, iihmpu, nfecra
integer iturb(nphmax), icp(nphmax)
integer iverif

! VARIABLES LOCALES

integer iphas

!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
!       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
!       LA BIBLIOTHEQUE)
!===============================================================================

if (iverif.eq.0) then
  if(iihmpu.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usipph DOIT ETRE COMPLETE',/,&
'@       DANS LE FICHIER usini1.F                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================


!     ON NE RENSEIGNERA DANS CE SOUS-PROGRAMME QUE LES PARAMETRES

!       qui y apparaissent deja, A L'EXCLUSION DE tout autre.
!                                ================


!     SI L'ON NE DISPOSE PAS DE L'INTERFACE DE CODE_SATURNE :

!       on renseignera TOUS les parametres qui apparaissent dans
!                      ====                          ce sous-programme.


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :

!       on trouvera dans les sous-programmes utilisateur
!       des exemples commentes sur le modele de la section presente.

!       l'utilisateur pourra, si necessaire, les decommenter et les
!       adapter a ses besoins.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!===============================================================================

! --- Turbulence (pour chaque phase)
!       0...Laminaire
!      10...Longueur de melange
!      20...k-epsilon
!      21...k-epsilon a production lineaire
!      30...Rij-epsilon standard (LRR)
!      31...Rij-epsilon SSG
!      40...LES (Smagorinsky)
!      41...LES (Dynamique)
!      42...LES (WALE)
!      50...v2f (phi-model)
!      60...k-omega SST
!  Pour 10, contacter l'equipe de developpement avant utilisation

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphas = 1
iturb(iphas) = 20

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- Chaleur massique variable (ICP=1) ou non (ICP=0)
!       pour chaque phase IPHAS

!     A renseigner uniquement si les physiques particulieres (charbon,
!       combustion, electrique) NE SONT PAS activees.

!     Pour ces physiques particulieres, il NE FAUT PAS modifier ICP ici
!       et les choix suivants sont imposes
!          charbon et combustion : CP constant ;
!          electrique            : CP variable.

!     Attention : completer usphyv avec la loi donnant Cp
!     =========   si et seulement si on a choisi CP variable ici
!                                                   (avec ICP(IPHAS)=1)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphas = 1
icp(iphas) = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!----
! FORMATS
!----



return
end


!===============================================================================


subroutine usinsc &
!================

 ( iihmpu, nfecra , nscaus , iverif )


!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE UTILISATEUR POUR ENTREE
!   DU NOMBRE DE SCALAIRES UTILISATEUR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! iihmpu           ! e  ! <-- ! indique si un fichier de donnees de            !
!                  !    !     ! l'ihm est utilise (1: oui, 0: non)             !
! nfecra           ! e  ! <-- ! numero d'unite std pour impression             !
! nscaus           ! e  ! <-- ! nombre de scalaires utilisateur                !
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
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


!     AUCUN COMMON NE DOIT APPARAITRE ICI


!===============================================================================

! Arguments

integer iihmpu, nfecra
integer nscaus
integer iverif

! VARIABLES LOCALES


!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
!       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
!       LA BIBLIOTHEQUE)
!===============================================================================

if (iverif.eq.0) then
  if(iihmpu.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usinsc DOIT ETRE COMPLETE',/,&
'@       DANS LE FICHIER usini1.F                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================


!     ON NE RENSEIGNERA DANS CE SOUS-PROGRAMME QUE LES PARAMETRES

!       qui y apparaissent deja, A L'EXCLUSION DE tout autre.
!                                ================


!     SI L'ON NE DISPOSE PAS DE L'INTERFACE DE CODE_SATURNE :

!       on renseignera TOUS les parametres qui apparaissent dans
!                      ====                          ce sous-programme.


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :

!       on trouvera dans les sous-programmes utilisateur
!       des exemples commentes sur le modele de la section presente.

!       l'utilisateur pourra, si necessaire, les decommenter et les
!       adapter a ses besoins.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!===============================================================================

! --- Nombre de scalaires UTILISATEUR (thermiques ou non et quelle que
!       soit leur phase porteuse). Il s'agit des scalaires qui viennent
!       en supplement des scalaires "de base" suivants (naturellement
!       inclus dans la modelisation) :
!        - pression
!        - grandeurs turbulentes
!        - NSCAPP scalaires introduits par l'eventuel modele de
!          combustion, charbon ou electrique mis en oeuvre

!     Ainsi, pour un calcul sans modele de combustion, de charbon
!       pulverise ni de modele electrique, les scalaires utilisateurs
!       pourront etre par exemple :
!        - la temperature ou l'enthalpie,
!        - des fractions massiques de scalaires transportes
!        - la moyenne du carre des fluctuations d'un autre scalaire
!          utilisateur

!     Le nombre de scalaires maximal est donne par NSCAMX dans paramx.h :
!       il s'agit de la valeur maximale admissible pour NSCAUS + NSCAPP.


!     Imposer NSCAUS = 0 s'il n'y a pas de scalaire utilisateur.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

nscaus = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!----
! FORMATS
!----



return
end


!===============================================================================


subroutine usipsc &
!================

 ( nscmax, nscaus, iihmpu, nfecra, iscavr, ivisls , iverif )


!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE UTILISATEUR POUR ENTREE
!   DE PARAMETRES DEPENDANT DU NOMBRE DE SCALAIRES UTILISATEUR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nscmax           ! e  ! <-- ! nombre de scalaires maximal                    !
! nscaus           ! e  ! <-- ! nombre de scalaires utilisateur                !
! iihmpu           ! e  ! <-- ! indique si un fichier de donnees de            !
!                  !    !     ! l'ihm est utilise (1: oui, 0: non)             !
! nfecra           ! e  ! <-- ! numero d'unite std pour impression             !
! iscavr(nscmax    ! te ! <-- ! numero du scalaire correspondant               !
!                  !    !     ! pour les scalaires variance                    !
! ivisls(nscmax    ! te ! <-- ! diffusivite des scalaires uniforme             !
!                  !    !     ! ou non                                         !
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
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


!     AUCUN COMMON NE DOIT APPARAITRE ICI


!===============================================================================

! Arguments

integer nscmax, nscaus, iihmpu, nfecra
integer iscavr(nscmax), ivisls(nscmax)
integer iverif

! VARIABLES LOCALES

integer iutile, iscal

!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
!       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
!       LA BIBLIOTHEQUE)
!===============================================================================

if (iverif.eq.0) then
  if(iihmpu.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usipsc DOIT ETRE COMPLETE',/,&
'@       DANS LE FICHIER usini1.F                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================


!     ON NE RENSEIGNERA DANS CE SOUS-PROGRAMME QUE LES PARAMETRES

!       qui y apparaissent deja, A L'EXCLUSION DE tout autre.
!                                ================


!     SI L'ON NE DISPOSE PAS DE L'INTERFACE DE CODE_SATURNE :

!       on renseignera TOUS les parametres qui apparaissent dans
!                      ====                          ce sous-programme.


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :

!       on trouvera dans les sous-programmes utilisateur
!       des exemples commentes sur le modele de la section presente.

!       l'utilisateur pourra, si necessaire, les decommenter et les
!       adapter a ses besoins.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!===============================================================================

! --- Moyenne du carre des fluctuations d'un scalaire UTILISATEUR :
!     Si on souhaite qu'un scalaire utilisateur j represente la moyenne
!      du carre des fluctuations du scalaire utilisateur k, on indique
!      ISCAVR(j) = k.
!     Les valeurs prises par ISCAVR sont donc naturellement superieures
!      ou egales a 1 et inferieures ou egales au nombre total de scalaires.
!      Ainsi, si on impose ISCAVR(j) = k, il faut
!            0<j<NSCAUS+1, 0<k<NSCAUS+1 et j different de k.

!     Par exemple pour que le scalaire utilisateur 3 soit la moyenne
!      du carre des fluctuations du scalaire utilisateur 2, on impose :
!                         ISCAVR(3) = 2
!      avec NSCAUS au moins egal a 3.

!     Ne pas intervenir si l'on ne souhaite pas inclure explicitement a
!       la simulation la moyenne du carre des fluctuations d'un scalaire
!       utilisateur.

!     Pour les scalaires non utilisateur relatifs a des physiques
!       particulieres, (charbon, combustion, electrique : voir usppmo)
!       implicitement definis selon le modele, les informations sont
!       donnees automatiquement par ailleurs : on ne modifie pas ISCAVR.


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     Le test sur IUTILE permet de desactiver les instructions (qui
!       ne sont fournies qu'a titre d'exemple a adapter)

iutile = 0
if(iutile.eq.1) then
  iscavr(3) = 2
endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN




! --- Diffusivite variable (IVISLS=1) ou constante (IVISLS=0) pour
!       chaque scalaire UTILISATEUR     hormis ceux qui
!       representent la moyenne du carre des fluctuations d'un autre.

!     Pour les scalaires utilisateur ISCAL qui representent la moyenne du
!       carre des fluctuations d'un autre scalaire utilisateur,
!                                on ne renseigne pas IVISLS(ISCAL) ici.
!       C'est l'objet du test sur ISCAVR(ISCAL) dans l'exemple ci-dessous.
!       En effet, la diffusivite de la moyenne du carre des fluctuations
!       d'un scalaire est supposee avoir le meme comportement que la
!       diffusivite de ce scalaire.

!     Pour les scalaires non utilisateur relatifs a des physiques
!       particulieres, (charbon, combustion, electrique : voir usppmo)
!       implicitement definis selon le modele,
!       les informations sont donnees automatiquement par ailleurs :
!                                          on ne modifie pas IVISLS ici.

!     Attention : completer usphyv avec la loi donnant la diffusivite
!     =========   si et seulement si on a choisi IVISLS=1 ici




! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

do iscal = 1, nscaus

!     Pour les scalaires utilisateur qui ne representent pas
!       la moyenne du carre des fluctuations d'un autre scalaire
  if(iscavr(iscal).le.0) then

    ivisls(iscal) = 0

  endif

enddo

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!----
! FORMATS
!----



return
end


!===============================================================================


subroutine usipgl &
!================

 ( nphmax, nesmax,                                                &
   iespre, iesder, iescor, iestot,                                &
   nphas , iihmpu, nfecra,                                        &
   idtvar, ipucou, iphydr, ialgce , iescal , iverif )


!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE UTILISATEUR POUR ENTREE
!   DE PARAMETRES GLOBAUX

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nphmax           ! e  ! <-- ! nombre de phases maximal                       !
! nesmax           ! e  ! <-- ! nombre d'estimateurs d'erreur                  !
!                  !    !     !   maximal pour une phase                       !
! iespre           ! e  ! <-- ! numero de l'estimateur d'erreur                !
!                  !    !     !   prediction                                   !
! iesder           ! e  ! <-- ! numero de l'estimateur d'erreur                !
!                  !    !     !   derive                                       !
! iescor           ! e  ! <-- ! numero de l'estimateur d'erreur                !
!                  !    !     !   correction                                   !
! iestot           ! e  ! <-- ! numero de l'estimateur d'erreur                !
!                  !    !     !   total                                        !
! nphas            ! e  ! <-- ! nombre de phases actives                       !
! iihmpu           ! e  ! <-- ! indique si un fichier de donnees de            !
!                  !    !     ! l'ihm est utilise (1: oui, 0: non)             !
! nfecra           ! e  ! <-- ! numero d'unite std pour impression             !
! idtvar           ! e  ! --> ! indicateur pas de temps variable               !
! ipucou           ! e  ! --> ! indicateur couplage u-p renforce               !
! iphydr           ! e  ! --> ! indicateur de prise en compte de               !
!                  !    !     ! l'equilibre entre le gradient de               !
!                  !    !     ! pression et les termes de gravite et           !
!                  !    !     ! de perte de charge                             !
! ialgce           ! e  ! <-- ! indicateur de methode de calcul des            !
!                  !    !     !  centres de gravite des cellules               !
! iescal           ! te ! <-- ! indicateur d'activation des                    !
!(nesmax,nphmax    !    !     ! estimateurs d'erreur pour navier               !
!                  !    !     ! stokes                                         !
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
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


!     AUCUN COMMON NE DOIT APPARAITRE ICI


!===============================================================================

! Arguments

integer nphmax, nesmax
integer iespre, iesder, iescor, iestot
integer nphas , iihmpu, nfecra
integer idtvar, ipucou, iphydr
integer iescal(nesmax,nphmax)
integer iverif

! VARIABLES LOCALES

integer iphas, ialgce

!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
!       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
!       LA BIBLIOTHEQUE)
!===============================================================================

if (iverif.eq.0) then
  if(iihmpu.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usipgl DOIT ETRE COMPLETE',/,&
'@       DANS LE FICHIER usini1.F                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================


!     ON NE RENSEIGNERA DANS CE SOUS-PROGRAMME QUE LES PARAMETRES

!       qui y apparaissent deja, A L'EXCLUSION DE tout autre.
!                                ================


!     SI L'ON NE DISPOSE PAS DE L'INTERFACE DE CODE_SATURNE :

!       on renseignera TOUS les parametres qui apparaissent dans
!                      ====                          ce sous-programme.


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :

!       on trouvera dans les sous-programmes utilisateur
!       des exemples commentes sur le modele de la section presente.

!       l'utilisateur pourra, si necessaire, les decommenter et les
!       adapter a ses besoins.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!===============================================================================

! --- Pas de temps (0 : uniforme et constant
!                   1 : variable en temps et uniforme en espace
!                   2 : variable en espace et en temps
!                  -1 : algorithme stationnaire)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

idtvar = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- Couplage vitesse/pression (0 : algorithme classique,
!                                1 : couplage instationnaire)
!     Uniquement en monophasique

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

ipucou = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- Prise en compte de la pression hydrostatique
!                               (0 : algorithme usuel
!                                1 : prise en compte explicite)
!     Uniquement en monophasique

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphydr = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- Estimateurs pour Navier-Stokes (champ de vitesse non fige)
!     On conseille de realiser une suite de calcul sur quelques pas de
!       temps en activant les plus parlants ci dessous
!        (=2 pour activer, =0 pour desactiver).

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphas = 1
!       div(rho u) -Gamma
iescal(iescor,iphas) = 0
!       precision de la resolution de la quantite de mouvement
iescal(iestot,iphas) = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


!----
! FORMATS
!----



return
end


!===============================================================================


subroutine usipsu &
!================

 ( nmodpp , iverif )


!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE UTILISATEUR POUR ENTREE
!   DE PARAMETRES UTILISATEUR SUPPLEMENTAIRES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nmodpp           ! e  ! <-- ! nombre de modeles phys.part. actives           !
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
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
include "cstnum.h"
include "dimens.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "vector.h"
include "parall.h"
include "period.h"
include "ihmpre.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"
include "coincl.h"
include "cpincl.h"
include "elincl.h"

!===============================================================================

! Arguments

integer nmodpp
integer iverif

! VARIABLES LOCALES

integer iphas, iutile, ii, jj, imom

!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
!       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
!       LA BIBLIOTHEQUE)
!===============================================================================

if (iverif.eq.0) then
  if(iihmpr.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usipsu DOIT ETRE COMPLETE',/,&
'@       DANS LE FICHIER usini1.F                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================


!     CE SOUS-PROGRAMME PERMET DE RENSEIGNER LES PARAMETRES

!       qui n'apparaissent pas deja dans les autres sous-programmes
!       du present fichier.


!     IL EST POSSIBLE D'AJOUTER OU DE RETRANCHER DES PARAMETRES


!     LE NUMERO DES PROPRIETES PHYSIQUES ET DES  VARIABLES EST CONNU


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :

!       on trouvera dans les sous-programmes utilisateur
!       des exemples commentes sur le modele de la section presente.

!       l'utilisateur pourra, si necessaire, les decommenter et les
!       adapter a ses besoins.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!===============================================================================


! Options du calcul (optcal.h)
! ============================

! --- Suite de calcul : ISUITE ( = 1) ou non (0)
!     Avec relecture du fichier suite auxiliaire ILEAUX ( = 1) ou non (0)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

isuite = 0
ileaux = 1

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Duree
!       NTMABS = numero absolu du dernier pas de temps desire
!         si on a deja fait 10 pas de temps
!         et qu'on veut en faire 10 autres,
!           il faut imposer NTMABS a 10 + 10 = 20

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

ntmabs = 10

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Pas de temps de reference
!     L'exemple donne ici est probablement inadapte a votre cas

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

dtref  = 0.01d0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Pas de temps maximal DTMAX
!     Imposer une valeur a partir de grandeurs caracteristiques du cas
!      sinon, le code prendra par defaut un multiple de dtref
!     Ex. avec
!        Ld : longueur "dynamique" ie par exemple longueur du domaine
!        Ud : Vitesse debitante caracteristique
!        Lt : longueur thermique ie par exemple la hauteur du domaine
!                                                    (selon la gravite)
!        Delta_rho/rho : ecart relatif de masse volumique
!        g : acceleration de la pesanteur


!     DTMAX = Min(Ld/Ud,Sqrt(Lt/(g Delta_rho/rho)))








! --- Temperature ou enthalpie



!   Lorsque des physiques particulieres sont activees
!                                       (charbon, combustion, electrique)
!     on NE renseigne PAS cette section : on NE modifie ni ISCALT ni ISCSTH
!                                   (le test IF(NMODPP.EQ.0) sert a cela).


!   Par contre, si ces physiques particulieres ne sont PAS activees :

!     Si un scalaire UTILISATEUR represente la temperature ou l'enthalpie
!       (de la phase IPHAS) :
!          on donne le numero de ce scalaire dans ISCALT(IPHAS) et
!          on indique  ISCSTH(ISCALT(IPHAS)) = 1 si c'est la temperature
!                  ou  ISCSTH(ISCALT(IPHAS)) = 2 si c'est l'enthalpie

!     Si aucun scalaire ne represente la temperature ou l'enthalpie (de
!       la phase IPHAS)
!          on indique ISCALT(IPHAS) = -1
!          et on ne renseigne pas ISCSTH(ISCALT(IPHAS))


!     Pour le module rayonnement mis en oeuvre en dehors des physiques
!      charbon, combustion, electrique :
!      si l'on a choisi de resoudre en temperature (c'est-a-dire si
!      ISCSTH(ISCALT(IPHAS)) = 1), la temperature du fluide est alors
!      supposee en KELVIN (attention aux conditions aux limites et a
!      l'expression des proprietes physiques dependantes de la
!      temperature).
!      Neanmoins, bien que ce soit deconseille, si l'on souhaite
!      que le solveur fluide travaille avec une temperature en degres
!      Celsius, on doit imposer ISCSTH(ISCALT(IPHAS)) = -1.
!      Ce choix est source d'erreurs pour l'utilisateur. En effet, les
!      conditions aux limites pour la temperature du fluide seront alors
!      en degres Celsius, alors que les conditions aux limites imposees
!      pour le rayonnement dans usray2 devront, pour leur part, etre
!      obligatoirement en Kelvin.



!    Si les physiques particulieres ne sont pas activees
!       (charbon, combustion, electrique : voir usppmo)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

if(nmodpp.eq.0) then

  iphas = 1
!         Numero du scalaire representant la temperature ou l'enthalpie
!                                               ou -1 s'il n'y en a pas
!         Lorsque le choix est fait par l'interface de Code_Saturne,
!           le scalaire representant la temperature ou l'enthalpie
!           est obligatoirement le premier
  iscalt(iphas) = -1

!         S'il y a une variable temperature ou enthalpie
  if(iscalt(iphas).gt.0) then
!           on indique si c'est la temperature (=1) ou l'enthalpie (=2)
    iscsth(iscalt(iphas)) = 1
  endif

endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Calcul (suite) a champ de vitesse fige (1 oui, 0 non)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iccvfg = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Methode des vortex pour les conditions d'entree en L.E.S.
!              (0 : non activee,  1 : activee)
!     La methode des vortex ne concerne que les modeles L.E.S.
!       et n'est valable qu'avec une seule phase
!     Pour utiliser la methode des vortex, veuillez renseigner
!       le sous-programme usvort.F

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphas = 1
if (itytur(iphas).eq.4) then
  ivrtex = 0
endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Schema convectif

!       BLENCV = 0 pour upwind (ordre 1 en espace, "stable mais diffusif")
!              = 1 pour centre/second order (ordre 2 en espace)
!           on peut utiliser des valeurs reelles intermediaires
!       Ici on choisit
!         pour la vitesse de la phase 1 et les scalaires utilisateurs :
!            un schema centre upwind a 100% de centre (BLENCV=1)
!         pour les autres variables
!            la valeur par defaut du code (upwind en standard, centre en LES)

!     En particulier, pour les scalaires utilisateur
!       si l'on suspecte un niveau de diffusion numerique trop grand
!         sur une variable IVAR representant un scalaire utilisateur
!         ISCAL (avec IVAR=ISCA(ISCAL)), il peut etre utile d'imposer
!         BLENCV(IVAR) = 1.0D0 pour utiliser un schema d'ordre 2 en
!         espace pour la convection. Pour la temperature ou l'enthalpie,
!         en particulier, on pourra donc choisir dans ce cas :
!          BLENCV(ISCA(ISCALT(IPHAS))) = 1.0D0

!       Pour les scalaires non utilisateur relatifs a des physiques
!         particulieres, (charbon, combustion, electrique : voir usppmo)
!         implicitement definis selon le modele,
!         les informations sont donnees automatiquement par ailleurs :
!                                         on ne modifie pas BLENCV ici.


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphas = 1

blencv(iu(iphas)) = 1.0d0
blencv(iv(iphas)) = 1.0d0
blencv(iw(iphas)) = 1.0d0
if(nscaus.ge.1) then
  do ii = 1, nscaus
    blencv(isca(ii)) = 1.0d0
  enddo
endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- Parametres  multigrille algebrique

!     IMGR = 0 : PAS DE MULTIGRILLE
!     IMGR = 1 : MULTIGRILLE ALGEBRIQUE

!     Uniquement disponible pour la pression

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphas = 1
imgr(ipr(iphas)) = 1

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


!=========================================================================

! --- Stabilisation en regime turbulent

!     Pour les cas difficiles, une stabilisation peut etre obtenue en
!     ne reconstruisant pas les flux de convection et de diffusion
!     pour les variables du modele de turbulence, soit
!       en k-epsilon : IF (ITYTUR(IPHAS).EQ.2) THEN
!          IRCFLU(IK(IPHAS))   = 0 et IRCFLU(IEP(IPHAS))  = 0
!       en Rij-epsilon : IF (ITYTUR(IPHAS).EQ.3) THEN
!          IRCFLU(IR11(IPHAS)) = 0,   IRCFLU(IR22(IPHAS)) = 0,
!          IRCFLU(IR33(IPHAS)) = 0,
!          IRCFLU(IR12(IPHAS)) = 0,   IRCFLU(IR23(IPHAS)) = 0,
!          IRCFLU(IR23(IPHAS)) = 0,
!                                  et IRCFLU(IEP(IPHAS))  = 0
!     (noter que la variable ITYTUR vaut ITURB/10)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     Le test sur IUTILE permet de desactiver les instructions (qui
!       ne sont fournies qu'a titre d'exemple a adapter)

iutile = 0
if(iutile.eq.1) then

  iphas = 1
  if (iturb(iphas).eq.20) then
    ircflu(ik(iphas))   = 0
    ircflu(iep(iphas))  = 0
  endif

endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! Constantes physiques (cstphy.h)
! ===============================

! --- gravite (g en m/s2, avec le signe dans le repere de calcul)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

gx = 0.d0
gy = 0.d0
gz = 0.d0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Proprietes  de reference du fluide (pour chaque phase)

!       RO0        : masse volumique en kg/m3
!       VISCL0     : viscosite dynamique en kg/(m s)
!       CP0        : chaleur specifique en J/(degres kg)
!       T0         : temperature de reference en Kelvin
!       P0         : pression totale de reference en Pascal
!                    le calcul est par contre effectue sur
!                    une pression reduite P*=Ptot-ro0*g.(x-xref)
!                    (sauf en compressible)
!       XYZP0(3,.) : coordonnees du point de reference pour
!                    la pression totale (la ou elle vaut P0)

!     En general, il n'est pas necessaire de fournir un point
!       de reference XYZ0. S'il y a des sorties, le code
!       prendra le centre de la face de sortie de reference.
!       Par contre, si on prevoit de fixer explicitement des
!       conditions de Dirichlet pour la pression, mieux vaut
!       indiquer a quelle reference les valeurs se rapportent
!       (pour une meilleure resolution de la pression reduite)


!     Les autres proprietes seront fournies par defaut dans tous
!       les cas    .

!     Neanmoins, on peut noter que :

!       Dans les cas standard (ni combustion gaz, ni charbon, ni electrique,
!                              ni compressible) :
!       ---------------------
!         RO0, VISCL0 et CP0
!             sont utiles et representent soit les proprietes du fluide
!             si elles sont constantes, soit de simples valeurs moyennes
!             pour l'initialisation si les proprietes sont variables et
!             donnees dans usphyv.
!         T0  est inutile
!         P0  est utile mais n'est pas utilise dans une loi d'etat. P0
!             est une valeur de reference pour le solveur incompressible
!             qui servira a caler la pression a la sortie (eventuelle)
!             du domaine. On peut la prendre nulle ou egale a une valeur
!             physique en Pascal.

!       En modules electriques :
!       ----------------------
!         RO0, VISCL0 et CP0
!             sont utiles mais representent de simples valeurs moyennes
!             initiales ; la masse volumique, la viscosite dynamique
!             moleculaire et la chaleur massique sont donnees
!             obligatoirement dans PROPCE (qu'elles soient physiquement
!             variables ou non) : voir uselph pour le module effet Joule
!             et le fichier de donnees dp_ELE en arc electrique.
!         T0  est utile et doit etre en Kelvin (> 0) mais represente
!             une simple valeur d'initialisation.
!         P0  est utile mais n'est pas utilise dans une loi d'etat. P0
!             est une valeur de reference pour le solveur incompressible
!             qui servira a caler la pression a la sortie (eventuelle)
!             du domaine. On peut la prendre nulle ou egale a une valeur
!             physique en Pascal.

!       En combustion gaz :
!       -----------------
!         RO0 est inutile (il est recalcule automatiquement par la
!             loi des gaz parfaits a partir de T0 et P0).
!         VISCL0 est indispensable : c'est la viscosite dynamique
!             moleculaire, supposee constante pour le fluide
!         CP0 est indispensable : c'est la chaleur massique,
!             supposee constante (modelisation des termes sources faisant
!             intervenir un Nusselt local dans
!             le module lagrangien, valeur de reference permettant
!             de calculer un couple (temperature, coefficient d'echange)
!             en rayonnement)
!         T0  est indispensable et doit etre en Kelvin (> 0)
!         P0  est indispensable et doit etre en Pascal (> 0)

!       En charbon pulverise :
!       --------------------
!         RO0 est inutile (il est recalcule automatiquement par la
!             loi des gaz parfaits a partir de T0 et P0).
!         VISCL0 est indispensable : c'est la viscosite dynamique
!             moleculaire, supposee constante pour le fluide (son
!             effet est a priori faible devant les effets turbulents)
!         CP0 est indispensable : c'est la chaleur massique,
!             supposee constante (modelisation des termes sources faisant
!             intervenir un Nusselt local dans le module charbon ou
!             le module lagrangien, valeur de reference permettant
!             de calculer un couple (temperature, coefficient d'echange)
!             en rayonnement)
!         T0  est indispensable et doit etre en Kelvin (> 0)
!         P0  est indispensable et doit etre en Pascal (> 0)

!       En compressible :
!       --------------------
!         RO0 est inutile, stricto sensu ; neanmoins, comme l'experience
!             montre que les utilisateurs se servent souvent de cette
!             variable, on demande de lui affecter ici une valeur
!             strictement positive (par exemple, une valeur initiale)
!         VISCL0 est utile et represente la viscosite dynamique
!             moleculaire, lorsqu'elle est constante ou une valeur qui
!             servira lors des initialisations (ou dans les conditions
!             d'entree de la turbulence, selon le choix de l'utilisateur)
!         CP0 est indispensable : c'est la chaleur massique,
!             supposee constante dans la thermodynamique disponible par
!             defaut
!         T0  est indispensable et doit etre en Kelvin (> 0)
!         P0  est indispensable et doit etre en Pascal (> 0)
!             Avec la loi thermodynamique disponible par defaut,
!             T0 et P0 servent a l'initialisation de la masse volumique.
!         XYZP0 est inutile car la variable pression represente directement
!             la presion totale


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphas = 1

ro0   (iphas) = 0.235d0
viscl0(iphas) = 0.84d-6
cp0   (iphas) = 1219.d0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

iphas = 1

t0    (iphas) = 1000.d0 + 273.15d0
p0    (iphas) = 1.013d5
!     On ne specifie XYZ0 que si on compte explicitement fixer
!       des conditions de Dirichlet sur la pression
!      XYZP0(1,IPHAS) = 0.D0
!      XYZP0(2,IPHAS) = 0.D0
!      XYZP0(3,IPHAS) = 0.D0


! --- IROVAR IVIVAR : masse volumique et viscosite constantes ou non

!     Lorsqu'un modele physique particuliere a ete active
!                  (charbon, combustion, modules electriques, compressible)
!       on ne renseigne PAS les variables IROVAR et IVIVAR ici
!                                  elles sont definies automatiquement.
!       Pour le compressible, cependant, IVIVAR peut être modifie dans
!                                  le sous-programme utilisateur uscfx1.

!     Lorsqu'aucun modele physique particuliere n'a ete active,
!         (i.e. charbon, combustion, electrique, compressible : voir usppmo)
!       il faut indiquer si la masse volumique et la viscosite moleculaire
!         sont constantes (IROVAR=0, IVIVAR=0)
!           ou variables  (IROVAR=1, IVIVAR=1)

!       si elles sont variables , il faut alors definir la loi dans usphyv.
!       si elles sont constantes, elles prennent les valeurs RO0 et VISCL0

!       a titre d'exemple, on suppose ci-dessous qu'elles sont constantes

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

if(nmodpp.eq.0) then
  iphas = 1
  irovar(iphas) = 0
  ivivar(iphas) = 0
endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- Valeurs min SCAMIN et max SCAMAX admissibles pour
!        chaque scalaire UTILISATEUR :

!      Les resultats sont clipppes a la fin de chaque pas de temps

!      Si SCAMIN > SCAMAX, on ne clippe pas.

!      Pour un scalaire JJ representant la moyenne du carre des
!         fluctuations d'un autre, on pourra s'abstenir de renseigner
!         ces valeurs (un clipping par defaut est mis en place).
!         C'est l'objet du test sur ISCAVR(JJ) dans l'exemple ci-dessous.


!      Pour les scalaires non utilisateur relatifs a des physiques
!         particulieres, (charbon, combustion, electrique : voir usppmo)
!         implicitement definis selon le modele, les informations sont
!         donnees automatiquement par ailleurs : on ne renseigne pas
!         SCAMIN SCAMAX.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     S'il y a des scalaires utilisateur
if(nscaus.gt.0) then

!       On boucle sur les scalaires utilisateurs :
  do jj = 1, nscaus
!         Pour les scalaires qui ne sont pas des variances
    if(iscavr(jj).le.0) then
!           On definit la borne min et la borne max
      scamin(jj) =-grand
      scamax(jj) =+grand
    endif
  enddo

endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Diffusivite de reference VISLS0 en kg/(m s) pour
!        chaque scalaire UTILISATEUR   hormis  pour ceux qui
!        representent la moyenne du carre des fluctuations d'un autre.

!     Pour les scalaires non utilisateur relatifs a des physiques
!       particulieres, (charbon, combustion, electrique : voir usppmo)
!       implicitement definis selon le modele,
!       les informations sont donnees automatiquement par ailleurs :
!                                        on ne modifie pas VISLS0 ici.

!     Pour les scalaires utilisateur JJ qui representent la moyenne du
!       carre des fluctuations d'un autre scalaire utilisateur,
!                                   on ne renseigne pas VISLS0(JJ) ici.
!       C'est l'objet du test sur ISCAVR(JJ) dans l'exemple ci-dessous.
!       En effet la diffusivite de la moyenne du carre des fluctuations
!       d'un scalaire est supposee identique à la diffusivite de ce
!       scalaire.

!     Lorsqu'on n'a pas active de physique particuliere
!       (charbon, combustion, electrique) et si un scalaire utilisateur
!        represente la temperature ou l'enthalpie,
!        on definit ici pour la temperature ou l'enthalpie :
!                                     VISLS0(ISCALT(IPHAS)) = Lambda/Cp

!     Ici, a titre d'exemple, on affecte a la VISCL0 la viscosite de
!       la phase porteuse, ce qui convient a des traceurs passifs qui
!       suivent le fluide.


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     S'il y a des scalaires utilisateur
if(nscaus.gt.0) then

!       On boucle sur les scalaires utilisateurs :
  do jj = 1, nscaus
!         Pour les scalaires qui ne sont pas des variances
    if(iscavr(jj).le.0) then
!           On definit la diffusivite
      visls0(jj) = viscl0(iphsca(jj))
    endif
  enddo

endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- Vitesse de reference pour l'initialisation de la turbulence (m2/s)
!      (utile seulement en turbulence)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphas = 1
uref(iphas)    = 1.d0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Echelle de longueur de reference en metres pour
!      initialisation de epsilon (et clipping particulier
!      de la turbulence, mais ce n'est pas l'option par defaut)
!      Donner une valeur de l'ordre de la plus grande
!      dimension du domaine physique dans lequel l'ecoulement
!      peut se developper
!      (utile seulement en turbulence)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iphas = 1
almax(iphas) = -grand

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- Definition des moments
!     (au maximum NBMOMX moments, correlations d'ordre maximum NDGMOX)

!     On calcule des moyennes temporelles du type <f1*f2*f3*...*fn>
!     Les fi sont des variables definies aux cellules (tableaux RTP et PROPCE)

!        IDFMOM(i,IMOM) repere la variable fi du moment IMOM
!          si IDFMOM > 0 c'est une variable resolue (RTP)
!          si IDFMOM < 0 c'est une variable auxiliaire (PROPCE)
!        IMOOLD(IMOM) donne en cas de suite le numero, dans l'ancien calcul
!                du moment a utiliser pour initialiser le moment IMOM du
!                nouveau calcul (par defaut IMOOLD(IMOM)=IMOM).
!                La valeur -1 indique qu'on doit reinitialiser le moment IMOM
!        NTDMOM(IMOM) donne le pas de temps de debut du calcul du moment

!     On donne ci dessous l'exemple du calcul des moments <u> et <rho u v>
!       le moment <u> est relu dans le fichier suite si on est en suite,
!         le moment <rho u v> est reinitialise a zero
!       le moment <u> est calcule a partir du pas de temps 1000
!         le moment <rho u v> est calcule a partir du pas de temps 10000


!     Le test sur IUTILE permet de desactiver les instructions (qui
!       ne sont fournies qu'a titre d'exemple a adapter)

iutile = 0
if(iutile.eq.1) then

!     Premier moment : <u>
  imom  = 1
  iphas = 1
  idfmom(1,imom) =  iu(iphas)
  ntdmom(imom)   =  1000
!     Second moment : <rho u v>
  imom  = 2
  iphas = 1
  idfmom(1,imom) = -irom(iphas)
  idfmom(2,imom) =  iu(iphas)
  idfmom(3,imom) =  iv(iphas)
  imoold(imom)   = -1
  ntdmom(imom)   =  10000

endif

!----
! FORMATS
!----



return
end


!===============================================================================


subroutine usipes &
!================

 ( nmodpp , iverif )


!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE UTILISATEUR POUR ENTREE
!   DE PARAMETRES UTILISATEUR SUPPLEMENTAIRES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! nmodpp           ! e  ! <-- ! nombre de modeles phys.part. actives           !
! iverif           ! e  ! <-- ! indicateur des tests elementaires              !
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
include "cstnum.h"
include "dimens.h"
include "numvar.h"
include "optcal.h"
include "cstphy.h"
include "entsor.h"
include "vector.h"
include "parall.h"
include "period.h"
include "ihmpre.h"
include "ppppar.h"
include "ppthch.h"
include "ppincl.h"

!===============================================================================

! Arguments

integer nmodpp
integer iverif

! VARIABLES LOCALES

integer iphas, ipp, imom

!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
!       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
!       LA BIBLIOTHEQUE)
!===============================================================================

if (iverif.eq.0) then
  if(iihmpr.eq.1) then
    return
  else
    write(nfecra,9000)
    call csexit (1)
  endif
endif

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR usipes DOIT ETRE COMPLETE',/,&
'@       DANS LE FICHIER usini1.F                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================


!     CE SOUS-PROGRAMME PERMET DE RENSEIGNER LES PARAMETRES

!       qui n'apparaissent pas deja dans les autres sous-programmes
!       du present fichier.


!     IL EST POSSIBLE D'AJOUTER OU DE RETRANCHER DES PARAMETRES


!     LE NUMERO DES PROPRIETES PHYSIQUES ET DES  VARIABLES EST CONNU


! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     SI L'ON    DISPOSE     DE L'INTERFACE DE CODE_SATURNE :

!       on trouvera dans les sous-programmes utilisateur
!       des exemples commentes sur le modele de la section presente.

!       l'utilisateur pourra, si necessaire, les decommenter et les
!       adapter a ses besoins.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!===============================================================================

!===============================================================================
! 1. ENTREE-SORTIE (optcal.h)
!===============================================================================

! --- ecriture du fichier suite auxiliaire IECAUX = 1 oui, 0 non

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

iecaux = 1

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- pas des sorties listing

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

ntlist = 1

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

! --- post traitement

!     ICHRVL : post traitement du domaine fluide (oui 1/non 0)
!     ICHRBO : post traitement du bord du domaine (oui 1/non 0)
!     ICHRSY : Post traitement des zones couplees avec Syrthes
!              (1 oui, 0 non)
!     ICHRMD : indique si les maillages ecrits seront :
!               0 : fixes,
!               1 : deformables a topologie constante,
!               2 : modifiables (pourront etre completement redefinis en
!                   cours de calcul via le sous-programme USMPST).
!              10 : comme INDMOD = 0, avec champ de déplacement
!              11 : comme INDMOD = 1, avec champ de déplacement
!              12 : comme INDMOD = 2, avec champ de déplacement

!     FMTCHR : format de sortie, parmi
!              'EnSight Gold', 'MED_fichier', ou 'CGNS'
!     OPTCHR : options associees au format de sortie, separees
!              par des virgules, parmi
!              'text'              (format texte, pour EnSight)
!              'binary'            (format binaire, choix par defaut)
!              'big_endian'        (force les sorties EnSight binaires
!                                   en mode 'big-endian')
!              'discard_polygons'  (ignore faces de type polygone)
!              'discard_polyhedra' (ignore cellules de type polyedre)
!              'divide_polygons'   (decoupe faces de type polygone)
!              'divide_polyhedra'  (decoupe cellules de type polyedre)
!              'split_tensors'     (ecrit les tenseurs en tant que
!                                   scalaires separes)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

ichrvl = 1
ichrbo = 0
ichrsy = 0

ichrmd = 0

FMTCHR = 'EnSight Gold'
OPTCHR = 'binary'

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- pas des chrono (-1 : une seule sortie en fin de calcul)
!                    (valeur strictement positive : periode)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

ntchr = -1

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- pas des sorties historiques

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

nthist = 1

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- Nombre de sondes et positions (limite a NCAPTM=100)

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

ncapt  = 4
xyzcap(1,1) = 0.30d0
xyzcap(2,1) = 0.15d0
xyzcap(3,1) = 0.01d0

xyzcap(1,2) = 0.30d0
xyzcap(2,2) = 0.00d0
xyzcap(3,2) = 0.01d0

xyzcap(1,3) = 0.30d0
xyzcap(2,3) =-0.08d0
xyzcap(3,3) = 0.01d0

xyzcap(1,4) = 0.60d0
xyzcap(2,4) =-0.05d0
xyzcap(3,4) = 0.01d0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


! --- variable courante

!     Comme pour les autres variables,
!       si l'on n'affecte pas les tableaux suivants,
!       les valeurs par defaut seront utilisees

!     NOMVAR( ) = nom de la variable
!     ICHRVR( ) = sortie chono (oui 1/non 0)
!     ILISVR( ) = suivi listing (oui 1/non 0)
!     IHISVR( ) = sortie historique (nombre de sondes et numeros)
!     si IHISVR(.,1)  = -1 sortie sur toutes les sondes

!     NB : Seuls les 8 premiers caracteres du nom seront repris dans le
!          listing le plus detaille



! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     Variables dynamiques courantes

!     Exemples pour la phase 1
iphas = 1

!     variable pression
ipp = ipprtp(ipr   (iphas))
NOMVAR(IPP)   = 'Pression'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     variable v1x
ipp = ipprtp(iu    (iphas))
NOMVAR(IPP)   = 'VitesseX'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     variable v1y
ipp = ipprtp(iv    (iphas))
NOMVAR(IPP)   = 'VitesseY'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     variable v1z
ipp = ipprtp(iw    (iphas))
NOMVAR(IPP)   = 'VitesseZ'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

if(itytur(iphas).eq.2) then

!     energie turbulente
  ipp = ipprtp(ik    (iphas))
  NOMVAR(IPP)   = 'EnerTurb'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     dissipation turbulente
  ipp = ipprtp(iep   (iphas))
  NOMVAR(IPP)   = 'Dissip'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif(itytur(iphas).eq.3) then

!     tensions de Reynolds
  ipp = ipprtp(ir11  (iphas))
  NOMVAR(IPP)   = 'Tens.R11'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     tensions de Reynolds
  ipp = ipprtp(ir22  (iphas))
  NOMVAR(IPP)   = 'Tens.R22'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     tensions de Reynolds
  ipp = ipprtp(ir33  (iphas))
  NOMVAR(IPP)   = 'Tens.R33'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     tensions de Reynolds
  ipp = ipprtp(ir12  (iphas))
  NOMVAR(IPP)   = 'Tens.R12'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     tensions de Reynolds
  ipp = ipprtp(ir13  (iphas))
  NOMVAR(IPP)   = 'Tens.R13'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     tensions de Reynolds
  ipp = ipprtp(ir23  (iphas))
  NOMVAR(IPP)   = 'Tens.R23'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     dissipation turbulente
  ipp = ipprtp(iep   (iphas))
  NOMVAR(IPP)   = 'Dissip'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif(iturb(iphas).eq.50) then

!     energie turbulente
  ipp = ipprtp(ik    (iphas))
  NOMVAR(IPP)   = 'EnerTurb'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     dissipation turbulente
  ipp = ipprtp(iep   (iphas))
  NOMVAR(IPP)   = 'Dissip'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     phi
  ipp = ipprtp(iphi  (iphas))
  NOMVAR(IPP)   = 'Phi'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!     f_barre
  ipp = ipprtp(ifb   (iphas))
  NOMVAR(IPP)   = 'f_barre'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

elseif(iturb(iphas).eq.60) then

!-->  energie turbulente
  ipp = ipprtp(ik    (iphas))
  NOMVAR(IPP)   = 'EnerTurb'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

!-->  omega
  ipp = ipprtp(iomg  (iphas))
  NOMVAR(IPP)   = 'Omega'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1

endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!     Variables scalaires utilisateur

!     On peut modifier ici les tableaux relatifs aux scalaires utilisateurs
!       mais les scalaires reserves pour les physiques particulieres sont
!       geres automatiquement. D'ou les tests sur NSCAUS qui assurent que
!       les scalaires vises sont effectivement des scalaires utilisateurs.
!       Par physique particuliere, on entend
!        uniquement celles definies par des modules specifiques
!        du code tels que charbon, combustion, electrique : voir usppmo)
!       Pour les scalaires non utilisateur relatifs a des physiques
!         particulieres, (charbon, combustion, electrique : voir usppmo)
!         implicitement definis selon le modele, les informations sont
!         donnees automatiquement par ailleurs : on ne modifie pas BLENCV.

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

if(isca(1).gt.0.and.nscaus.ge.1) then
 ipp = ipprtp(isca  (1))
  NOMVAR(IPP)  = 'scal 1'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

if(isca(2).gt.0.and.nscaus.ge.2) then
  ipp = ipprtp(isca  (2))
  NOMVAR(IPP)  = 'scal 2'
  ichrvr(ipp)  = 1
  ilisvr(ipp)  = 1
  ihisvr(ipp,1)= -1
endif

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


!     Autres variables

iphas = 1

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     variable masse volumique (sortie en post uniquement si variable
!                               ou si physique particuliere)
ipp = ipppro(ipproc(irom  (iphas)))
NOMVAR(IPP)   = 'masse vol'
ichrvr(ipp)   = max(irovar(iphas),nmodpp)
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     chaleur specifique
if(icp   (iphas).gt.0) then
  ipp = ipppro(ipproc(icp   (iphas)))
  NOMVAR(IPP)   = 'chal. spec.'
  ichrvr(ipp)   = 0
  ilisvr(ipp)   = 0
  ihisvr(ipp,1) = 0
endif

!     viscosite laminaire
ipp = ipppro(ipproc(iviscl(iphas)))
NOMVAR(IPP)   = 'visc. laminaire'
ichrvr(ipp)   = 0
ilisvr(ipp)   = 0
ihisvr(ipp,1) = 0

!     viscosite turbulente
ipp = ipppro(ipproc(ivisct(iphas)))
NOMVAR(IPP)   = 'visc. turb1'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     nombre de Courant
ipp = ipppro(ipproc(icour(iphas)))
NOMVAR(IPP)   = 'Nb Courant'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

!     nombre de Fourier
ipp = ipppro(ipproc(ifour(iphas)))
NOMVAR(IPP)   = 'Nb Fourier'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 0
ihisvr(ipp,1) = -1

!     variable CSMAGO dans le cas des modeles L.E.S. dynamiques
!     (carre de la "constante" de Smagorinsky)
if(ismago(iphas).gt.0) then
  ipp = ipppro(ipproc(ismago(iphas)))
  NOMVAR(IPP)   = 'Csdyn2'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
endif

!     moyennes temporelles (exemple pour le moment 1)
if(nbmomt.gt.0) then
  imom = 1
  ipp = ipppro(ipproc(icmome(imom)))
  NOMVAR(IPP) = 'MoyTps01'
  ichrvr(ipp) = 1
  ilisvr(ipp) = 1
  ihisvr(ipp,1) = -1
endif

!     pression totale (non definie en compressible)
if (ippmod(icompf).lt.0) then
  ipp = ipppro(ipproc(iprtot(iphas)))
  NOMVAR(IPP)   = 'Pression totale'
  ichrvr(ipp)   = 1
  ilisvr(ipp)   = 1
  ihisvr(ipp,1) = -1
endif

!     pas de temps local
ipp = ippdt
NOMVAR(IPP)   = 'pdt local'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

!     temps caracteristique du couplage
!        instationnaire vitesse/pression
ipp = ipptx
NOMVAR(IPP)   = 'Tx'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

ipp = ippty
NOMVAR(IPP)   = 'Ty'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

ipp = ipptz
NOMVAR(IPP)   = 'Tz'
ichrvr(ipp)   = 1
ilisvr(ipp)   = 1
ihisvr(ipp,1) = -1

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


!----
! FORMATS
!----



return
end


!===============================================================================


subroutine ustbtr &
!================

 ( ncel   , ncelet , nfac   , nfabor , nnod   ,                   &
   longia , longra ,                                              &
   nideve , nituse , nrdeve , nrtuse )

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE UTILISATEUR POUR DEFINIR LES DIMENSIONS
!   DES MACROS TABLEAUX IA ET RA
!   DES TABLEAUX UTILISATEUR ITUSER ET RTUSER
!   DES TABLEAUX DEVELOPPEUR IDEVEL ET RDEVEL

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! ncel             ! e  ! <-- ! nombre de cellules hors halo                   !
! ncelet           ! e  ! <-- ! nombre de cellules halo inclus                 !
! nfac             ! e  ! <-- ! nombre de faces internes                       !
! nfabor           ! e  ! <-- ! nombre de faces de bord                        !
! nnod             ! e  ! <-- ! nombre de noeuds                               !
! longia           ! e  ! --> ! dimension du tableau ia                        !
! longra           ! e  ! --> ! dimension du tableau ra                        !
! nideve           ! e  ! --> ! dimension du tableau idevel                    !
! nituse           ! e  ! --> ! dimension du tableau ituser                    !
! nrdeve           ! e  ! --> ! dimension du tableau rdevel                    !
! nrtuse           ! e  ! --> ! dimension du tableau rtuser                    !
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

!===============================================================================

! Arguments

integer          ncel  , ncelet, nfac  , nfabor, nnod
integer          longia, longra
integer          nideve, nituse, nrdeve, nrtuse

!===============================================================================

! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
!===============================================================================

if(1.eq.1) return

!===============================================================================
! TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN

!===============================================================================
! 1. DIMENSION DES MACROS TABLEAUX IA ET RA :

!  L'utilisateur peut avoir a modifier la taille des tableaux
!    d'entiers et de reels ici : LONGIA et LONGRA respectivement

!  Le nombre d'entiers LONGIA et le nombre de reels LONGRA
!    dependent des options de calcul, du type d'elements,
!    d es caracteristiques du maillage (2d, 3d, hybride,
!    non conforme ...) et du nombre de variables, par exemple.
!    En k-epsilon, si on note NCEL le nombre de cellules du
!    maillage, on peut la plupart du temps utiliser la majoration
!    grossiere suivante : LONGIA =  45*NCEL et LONGRA = 220*NCEL
!    En Rij-epsilon, une majoration de 20% pourra etre appliquee.
!    Ces valeurs sont relativement elevees pour prendre en compte
!    les maillages 2D qui comprennent de nombreuses faces de bord :
!    une formule plus precise mais plus complexe serait necessaire.
!    Pour les gros cas 3D, une estimation plus precise est donnee
!    par : LONGIA =  25*NCEL et LONGRA = 120*NCEL.

!  Si LONGIA et LONGRA sont laisses a 0, alors ces valeurs sont
!    remplies de maniere automatique

!===============================================================================

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     DIMENSION DU MACRO TABLEAU D'ENTIERS IA

longia = 0

!     DIMENSION DU MACRO TABLEAU DE REELS RA

longra = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN

!===============================================================================
! 2. DIMENSIONS DES TABLEAUX UTILISATEUR ITUSER ET RTUSER
!===============================================================================

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_DEBUT

!     DIMENSION DU TABLEAU D'ENTIERS ITUSER

nituse = 0

!     DIMENSION DU TABLEAU DE REELS RTUSER

nrtuse = 0

! CODE_FOURNI_COMME_EXEMPLE_A_ADAPTER_PAR_L_UTILISATEUR_FIN


!----
! FORMATS
!----



return
end
