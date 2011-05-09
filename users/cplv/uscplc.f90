!-------------------------------------------------------------------------------

!VERS


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

subroutine uscplc &
!================

 ( idbia0 , idbra0 ,                                              &
   nvar   , nscal  , nphas  ,                                     &
   maxelt , lstelt ,                                              &
   icodcl , itrifb , itypfb , izfppp ,                            &
   ia     ,                                                       &
   dt     , rtp    , rtpa   , propce , propfa , propfb ,          &
   coefa  , coefb  , rcodcl ,                                     &
   w1     , w2     , w3     , w4     , w5     , w6     , coefu  , &
   ra     )

!===============================================================================
! FONCTION :
! --------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN COUPLE CHARBON PULVERISE :
!   --------------------------------------------------------------

!    ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE

!      COMBUSTION EULERIENNE DE CHARBON PULVERISE ET
!      TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON

!    REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
!    (ICODCL,RCODCL) POUR LES VARIABLES INCONNUES
!    PENDANT DE USCLIM.F




! INTRODUCTION
! ============

! On donne ici les conditions aux limites par face de bord.

! L'identification des faces de bord concernees se fait grace
! a la commande GETFBR.

!  GETFBR(CHAINE,NLELT,LSTELT) :
!  - CHAINE est une chaine de caractere fournie par l'utilisateur
!    qui donne les criteres de selection
!  - NLTELT est renvoye par la commande. C'est un entier qui
!    correspond au nombre de faces de bord trouveees repondant au
!    critere
!  - LSTELT est renvoye par la commande. C'est un tableau d'entiers
!    de taille NLTELT donnant la liste des faces de bord trouvees
!    repondant au critere.

!  CHAINE peut etre constitue de :
!  - references de couleurs (ex. : 1, 8, 26, ...
!  - references de groupes (ex. : entrees, groupe1, ...)
!  - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
!  Ces criteres peuvent etre combines par des operateurs logiques
!  (AND et OR) et des parentheses
!  ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
!  les faces de bord de couleur 1, appartenant aux groupes 'groupe2'
!  ou 'groupe3' et de coordonnee Y inferieure a 1.


! TYPE DE CONDITIONS AUX LIMITES
! ==============================

! On peut affecter les conditions aux limites de deux manieres.


!    Pour les conditions "standard" :
!    --------------------------------

!     (entree, sortie libre, paroi, symetrie), on donne un code
!     stocke dans le tableau ITYPFB (dimensionne au nombre de
!     faces de bord,nombre de phases). Ce code sera ensuite
!     utilise par un sous-programme non utilisateur pour affecter
!     les conditions correspondantes (les scalaires, en
!     particulier, recevront alors les conditions de la phase a
!     laquelle ils sont associes). Ainsi :

!     Code      |  Type de bord
!     -------------------------
!      IENTRE   |   Entree
!      ISOLIB   |   Sortie libre
!      ISYMET   |   Symetrie
!      IPAROI   |   Paroi (lisse)
!      IPARUG   |   Paroi rugueuse

!     Les entiers IENTRE, ISOLIB, ISYMET, IPAROI, IPARUG
!     sont affectes par ailleurs (include param.h). Leur valeur
!     est superieure ou egale a 1 et
!         inferieure ou egale a NTYPMX (valeur fixee dans paramx.h)



!     En outre, il faut donner certaines valeurs :


!     - Entree (plus precisement entree/sortie a debit impose, car le debit
!               peut etre impose sortant) :

!       -> Conditions de Dirichlet sur les variables
!         autres que la pression obligatoire si le flux est entrant,
!         optionnelle si le flux est sortant (le code affecte flux nul
!         si aucun Dirichlet n'est specifie) ; ainsi
!         en face IFAC, sur la variable IVAR : RCODCL(IFAC,IVAR,1)


!     - Paroi lisse : (= solide impermeable, avec frottement lisse)

!       -> Valeur de la vitesse de paroi defilante s'il y a lieu
!         en face IFAC, RCODCL(IFAC,IU,1)
!                       RCODCL(IFAC,IV,1)
!                       RCODCL(IFAC,IW,1)
!       -> Code specifique et valeur de la temperature imposee
!         en paroi s'il y a lieu :
!         en face IFAC, ICODCL(IFAC,IVAR)   = 5
!                       RCODCL(IFAC,IVAR,1) = Temperature imposee
!       -> Code specifique et valeur du flux imposee
!         en paroi s'il y a lieu :
!         en face IFAC, ICODCL(IFAC,IVAR)   = 3
!                       RCODCL(IFAC,IVAR,3) = Flux impose
!                                        =
!        Noter que la condition par defaut pour les scalaires
!         (hors k et epsilon) est un Neumann homogene


!     - Paroi rugueuse : (= solide impermeable, avec frottement rugueux)

!       -> Valeur de la vitesse de paroi defilante s'il y a lieu
!         en face IFAC, RCODCL(IFAC,IU,1)
!                       RCODCL(IFAC,IV,1)
!                       RCODCL(IFAC,IW,1)
!       -> Valeur de la hauteur de rugosite dynamique a specifier dans
!                       RCODCL(IFAC,IU,3) (valeur pour IV et IW non utilisee)
!       -> Code specifique et valeur de la temperature imposee
!         en paroi rugueuse s'il y a lieu :
!         en face IFAC, ICODCL(IFAC,IVAR)   = 6
!                       RCODCL(IFAC,IVAR,1) = Temperature imposee
!                       RCODCL(IFAC,IVAR,3) = Hauteur de rugosite thermique
!       -> Code specifique et valeur du flux imposee
!         en paroi s'il y a lieu :
!         en face IFAC, ICODCL(IFAC,IVAR)   = 3
!                       RCODCL(IFAC,IVAR,3) = Flux impose
!                                        =
!        Noter que la condition par defaut pour les scalaires
!         (hors k et epsilon) est un Neumann homogene

!     - Symetrie (= paroi impermeable avec glissement) :

!       -> Rien a preciser


!     - Sortie libre (plus precisement entree/sortie libre a pression imposee)

!       -> Rien a preciser pour la pression et la vitesse
!          Pour les scalaires et grandeurs turbulentes, une valeur de Dirichlet
!            peut etre specifiee de maniere optionnelle. Le comportement est le
!            suivant :
!              * la pression est toujours traitee en Dirichlet
!              * si flux de masse entrant :
!                  on retient la vitesse a l'infini
!                  condition de Dirichlet pour les scalaires et grandeurs
!                    turbulentes (ou flux nul si l'utilisateur n'a pas
!                    specifie de Dirichlet)
!                 si flux de masse sortant :
!                    on impose un flux nul sur la vitesse, les scalaires et
!                    les grandeurs turbulentes

!       Noter que la pression sera recalee a P0
!           sur la premiere face de sortie libre trouvee


!    Pour les conditions "non standard" :
!    ------------------------------------

!     Autres que (entree, sortie libre, paroi, symetrie), on donne
!      - d'une part, pour chaque face :
!        -> une valeur de ITYPFB admissible
!           ie superieure ou egale a 1 et inferieure ou egale a
!           NTYPMX (voir sa valeur dans paramx.h).
!           Les valeurs predefinies dans paramx.h
!           IENTRE, ISOLIB, ISYMET, IPAROI, IPARUG sont dans cet
!           intervalle et il est preferable de ne pas affecter
!           inconsidrement et par hasard a ITYPFB la valeur
!           d'un de ces entiers. Pour eviter cela, on peut
!           utiliser IINDEF si l'on souhaite eviter de verifier
!           les valeurs dans paramx.h. IINDEF est une valeur
!           admissible a laquelle n'est attachee aucune condition
!           limite predefinie.
!           Noter que le tableau ITYPFB est
!           reinitialise a chaque pas de temps a la valeur
!           non admissible 0. Si on oublie de modifier ITYPFB pour
!           une face, le code s'arretera.

!      - et d'autre part pour chaque face et chaque variable :
!        -> un code     ICODCL(IFAC,IVAR)
!        -> trois reels RCODCL(IFAC,IVAR,1)
!                       RCODCL(IFAC,IVAR,2)
!                       RCODCL(IFAC,IVAR,3)
!     La valeur de ICODCL est prise parmi les suivantes :
!       1 : Dirichlet      (utilisable pour toute variable)
!       3 : Neumann        (utilisable pour toute variable)
!       4 : Symetrie       (utilisable uniquement pour la vitesse et
!                                   les composantes du tenseur Rij)
!       5 : Paroi lisse    (utilisable pour toute variable sauf la
!                                                         pression)
!       6 : Paroi rugueuse (utilisable pour toute variable sauf la
!                                                         pression)
!       9 : Sortie libre   (utilisable uniquement pour la vitesse)
!     Les valeurs des 3 reels RCODCL sont les suivantes
!      RCODCL(IFAC,IVAR,1) :
!         Dirichlet sur la variable        si ICODCL(IFAC,IVAR)=  1
!         valeur en paroi (defilmnt, temp) si ICODCL(IFAC,IVAR)=  5
!         La dimension de RCODCL(IFAC,IVAR,1) est celle de la
!           variable resolue : ex U (vitesse en m/s),
!                                 T (temperature en degres)
!                                 H (enthalpie en J/kg)
!                                 F (scalaire passif en -)
!      RCODCL(IFAC,IVAR,2) :
!         coefficient d'echange "exterieur" (entre la valeur
!                          imposee et la valeur au bord du domaine)
!                          RINFIN = infini par defaut
!         Pour les vitesses U,             en kg/(m2 s) :
!           RCODCL(IFAC,IVAR,2) =            (VISCL+VISCT) / D
!         Pour la  pression P,             en  s/m          :
!           RCODCL(IFAC,IVAR,2) =                       DT / D
!         Pour les temperatures T,         en Watt/(m2 degres) :
!           RCODCL(IFAC,IVAR,2) = CP*(VISCLS+VISCT/SIGMAS) / D
!         Pour les enthalpies H,           en kg /(m2 s) :
!           RCODCL(IFAC,IVAR,2) =    (VISCLS+VISCT/SIGMAS) / D
!         Pour les autres scalaires F      en :
!           RCODCL(IFAC,IVAR,2) =    (VISCLS+VISCT/SIGMAS) / D
!              (D a la dimension d'une distance en m)

!      RCODCL(IFAC,IVAR,3) si ICODCL(IFAC,IVAR)<>6 :
!        Densite de flux (< 0 si gain, n normale orientee vers l'exterieur)
!                         si ICODCL(IFAC,IVAR)= 3
!         Pour les vitesses U,             en kg/(m s2) = J :
!           RCODCL(IFAC,IVAR,3) =           -(VISCL+VISCT) * (GRAD U).n
!         Pour la  pression P,             en kg/(m2 s)     :
!           RCODCL(IFAC,IVAR,3) =                      -DT * (GRAD P).n
!         Pour les temperatures T,         en Watt/m2       :
!           RCODCL(IFAC,IVAR,3) =-CP*(VISCLS+VISCT/SIGMAS) * (GRAD T).n
!         Pour les enthalpies H,           en Watt/m2       :
!           RCODCL(IFAC,IVAR,3) =   -(VISCLS+VISCT/SIGMAS) * (GRAD H).n
!         Pour les autres scalaires F      en :
!           RCODCL(IFAC,IVAR,3) =   -(VISCLS+VISCT/SIGMAS) * (GRAD F).n

!      RCODCL(IFAC,IVAR,3) SI ICODCL(IFAC,IVAR)=6 :
!        Rugosites de la loi rugueuse
!         Pour les vitesses U, rugosite dynamique
!           RCODCL(IFAC,IVAR,3) = RUGD
!         Pour les autres scalaires, rugosite thermique
!           RCODCL(IFAC,IVAR,3) = RUGT


!      Noter bien que si l'utilisateur affecte une valeur a ITYPFB
!       parmi     IENTRE, ISOLIB, ISYMET, IPAROI, IPARUG
!       et qu'il ne modifie pas ICODCL (valeur nulle par defaut),
!       c'est ITYPFB qui imposera la condition limite.

!      Par contre, si l'utilisateur impose
!        ICODCL(IFAC,IVAR) (non nul),
!        ce sont alors les valeurs de RCODCL qu'il aura fournies
!        qui sont retenues pour la face et la variable consideree
!        (s'il ne precise pas RCODCL, ce sont les valeurs
!        par defaut qui sont retenues pour la face et
!        la variable consideree soit :
!                                 RCODCL(IFAC,IVAR,1) = 0.D0
!                                 RCODCL(IFAC,IVAR,2) = RINFIN
!                                 RCODCL(IFAC,IVAR,3) = 0.D0)
!        En particulier, on peut par exemple
!        -> donner ITYPFB(IFAC,IPHAS) = IPAROI
!        qui impose les conditions de paroi par defaut pour toutes
!        les variables sur la face IFAC,
!        -> et preciser EN OUTRE pour la variable IVAR sur cette
!        face IFAC des conditions paarticulieres en imposant
!        ICODCL(IFAC,IVAR) et les 3 RCODCL.


!      L'utilisateur peut egalement affecter a ITYPFB une valeur
!       non egale a IENTRE, ISOLIB, ISYMET, IPAROI, IPARUG, IINDEF
!       mais superieure ou egale a 1 et inferieure ou egale a
!       NTYPMX (voir les valeurs dans param.h) pour reperer
!       des groupes de couleurs dans d'autres sous programmes
!       qui lui seraient propres et ou ITYPFB serait disponible.
!       Dans ce cas cependant, il faudra
!       imposer les conditions limites en donnant des valeurs a
!       ICODCL et aux trois RCODCL (puisque la valeur de ITYPFB
!       ne sera pas predefinie dans le code).


! REGLES DE COHERENCE
! ===================

!       Quelques regles de coherence entre les codes ICODCL
!         des variables pour les conditions non standard :

!           Les codes des vitesses doivent etre identiques
!           Les codes des Rij      doivent etre identiques
!           Si code (vitesse ou Rij) = 4
!             il faut code (vitesse et Rij) = 4
!           Si code (vitesse ou turbulence) = 5
!             il faut code (vitesse et turbulence) = 5
!           Si code (vitesse ou turbulence) = 6
!             il faut code (vitesse et turbulence) = 6
!           Si code scalaire (hormis pression ou fluctuations) = 5
!             il faut code vitesse = 5
!           Si code scalaire (hormis pression ou fluctuations) = 6
!             il faut code vitesse = 6


! REMARQUES
! ==========

!       Attention : pour imposer un flux (non nul) sur les Rij,
!                   la viscosite a prendre en compte est VISCL
!                   meme si VISCT existe (VISCT=RHO CMU K2/EPSILON)


!       On dispose du tableau de tri des faces de bord au pas
!           de temps precedent (sauf au premier pas de temps, ou
!           ITRIFB n'a pas ete renseigne).
!       Le tableau du type des faces de bord ITYPFB a ete
!           reinitialise avant d'entrer dans le sous programme.



!       Noter comment acceder a certaines variables :

! Valeurs aux cellules
!               Soit        IEL = IFABOR(IFAC)

! * Masse vol                       phase IPHAS, cellule      IEL  :
!                  PROPCE(IEL ,IPPROC(IROM (IPHAS)))
! * Viscosite moleculaire dynamique phase IPHAS, cellule      IEL  :
!                  PROPCE(IEL ,IPPROC(IVISCL(IPHAS)))
! * Viscosite turbulente  dynamique phase IPHAS, cellule      IEL  :
!                  PROPCE(IEL ,IPPROC(IVISCT(IPHAS)))
! * Chaleur specifique              phase IPHAS, cellule      IEL  :
!                  PROPCE(IEL ,IPPROC(ICP   (IPHAS)))
! * Diffusivite lambda           scalaire ISCAL, cellule      IEL  :
!                  PROPCE(IEL ,IPPROC(IVISLS(ISCAL)))

! Valeurs aux faces de bord

! * Masse vol                      phase IPHAS, face de bord IFAC :
!                  PROPFB(IFAC,IPPROB(IROM (IPHAS)))
! * Flux de masse relatif a la variable  IVAR , face de bord IFAC :
!      (i.e. le flux de masse servant a la convection de IVAR)
!                  PROPFB(IFAC,IPPROB(IFLUMA(IVAR )))
! * Pour les autres grandeurs              a la face de bord IFAC :
!      prendre pour approximation la valeur dans la cellule  IEL
!      adjacente i.e. comme plus haut avec IEL = IFABOR(IFAC).

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! nvar             ! e  ! <-- ! nombre total de variables                      !
! nscal            ! e  ! <-- ! nombre total de scalaires                      !
! nphas            ! e  ! <-- ! nombre de phases                               !
! maxelt           !  e ! <-- ! nb max d'elements (cell,fac,fbr)               !
! lstelt(maxelt) te ! --- ! tableau de travail                             !
! icodcl           ! te ! --> ! code de condition limites aux faces            !
!  (nfabor,nvar    !    !     !  de bord                                       !
!                  !    !     ! = 1   -> dirichlet                             !
!                  !    !     ! = 3   -> densite de flux                       !
!                  !    !     ! = 4   -> glissemt et u.n=0 (vitesse)           !
!                  !    !     ! = 5   -> frottemt et u.n=0 (vitesse)           !
!                  !    !     ! = 6   -> rugosite et u.n=0 (vitesse)           !
!                  !    !     ! = 9   -> entree/sortie libre (vitesse          !
!                  !    !     !  entrante eventuelle     bloquee               !
! itrifb(nfabor    ! te ! <-- ! indirection pour tri des faces de brd          !
!  nphas      )    !    !     !                                                !
! itypfb(nfabor    ! te ! --> ! type des faces de bord                         !
!  nphas      )    !    !     !                                                !
! izfppp           ! te ! --> ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! dt(ncelet)       ! tr ! <-- ! pas de temps                                   !
! rtp, rtpa        ! tr ! <-- ! variables de calcul au centre des              !
! (ncelet,*)       !    !     !    cellules (instant courant ou prec)          !
! propce           ! tr ! <-- ! proprietes physiques au centre des             !
! (ncelet,*)       !    !     !    cellules                                    !
! propfa           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfac,*)        !    !     !    faces internes                              !
! propfb           ! tr ! <-- ! proprietes physiques au centre des             !
!  (nfabor,*)      !    !     !    faces de bord                               !
! coefa, coefb     ! tr ! <-- ! conditions aux limites aux                     !
!  (nfabor,*)      !    !     !    faces de bord                               !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
! w1,2,3,4,5,6     ! tr ! --- ! tableaux de travail                            !
!  (ncelet         !    !     !  (calcul du gradient de pression)              !
! coefu            ! tr ! --- ! tab de trav                                    !
!  (nfabor,3)      !    !     !  (calcul du gradient de pression)              !
! ra(*)            ! tr ! --- ! macro tableau reel                             !
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
use pointe
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use period
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0
integer          nvar   , nscal  , nphas

integer          maxelt, lstelt(maxelt)
integer          icodcl(nfabor,nvar)
integer          itrifb(nfabor,nphas), itypfb(nfabor,nphas)
integer          izfppp(nfabor)
integer          ia(*)

double precision dt(ncelet), rtp(ncelet,*), rtpa(ncelet,*)
double precision propce(ncelet,*)
double precision propfa(nfac,*), propfb(nfabor,*)
double precision coefa(nfabor,*), coefb(nfabor,*)
double precision rcodcl(nfabor,nvar,3)
double precision w1(ncelet),w2(ncelet),w3(ncelet)
double precision w4(ncelet),w5(ncelet),w6(ncelet)
double precision coefu(nfabor,ndim)
double precision ra(*)

! VARIABLES LOCALES

integer          idebia, idebra
integer          ifac, iphas, ii
integer          izone
integer          ilelt, nlelt

double precision uref2, d2s3
double precision xkent, xeent

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================
! 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
!       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
!       ET NON CELLE DE LA BIBLIOTHEQUE
!===============================================================================

if(1.eq.1) then
  write(nfecra,9001)
  call csexit (1)
  !==========
endif

 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,&
'@    =========                                               ',/,&
'@     COMBUSTION CHARBON PULVERISE COUPLE AU                 ',/,&
'@     TRANSPORT LAGRANGIEN DES PARTICULES DE CHARBON :       ',/,&
'@     LE SOUS-PROGRAMME UTILISATEUR uscpcl DOIT ETRE COMPLETE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END


!===============================================================================
! 1.  INITIALISATIONS

!===============================================================================

idebia = idbia0
idebra = idbra0

d2s3 = 2.d0/3.d0

!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR LES FACES DE BORD
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LA CONDITION LIMITE

!          IMPOSER ICI LES CONDITIONS LIMITES SUR LES FACES DE BORD

!          INTERVENTION UTLISATEUR

!===============================================================================

  iphas = 1

! ---- Face de type entree correspondant a une entree d'air
!        Par exemple : Air primaire , secondaire ou Air tertiaire

CALL GETFBR('12',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!   Type de condition aux limites pour les variables standard
  itypfb(ifac,iphas) = ientre

!   Numero de zone (on les numerote de 1 a n)
  izone = 1

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

! ------ Pour ces faces d'entree , on est a debit impose

  ientat(izone) = 1
  iqimp(izone)  = 1
!      - Debit en kg/s pour l'air
  qimpat(izone) = 1.46d-03
!      - Temperature en K pour l'air
  timpat(izone) = 400.d0 + tkelvi


! ------ On impose en couleur 12 une entree a debit impose
!        L'utilisateur donne donc ici uniquement
!          la direction du vecteur vitesse

  rcodcl(ifac,iu(iphas),1) = 0.d0
  rcodcl(ifac,iv(iphas),1) = 0.d0
  rcodcl(ifac,iw(iphas),1) = 5.d0

! ------ Traitement de la turbulence

!        La turbulence est calculee par defaut si ICALKE different de 0
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference adaptes a l'entree courante si ICALKE = 1
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference et de l'intensite turvulente
!            adaptes a l'entree courante si ICALKE = 2

!      Choix pour le calcul automatique ICALKE = 1 ou 2
  icalke(izone) = 1
!      Saisie des donnees
  dh(izone)     = 0.1d0
  xintur(izone) = 0.1d0



! Exemple de cas ou ICALKE(IZONE) = 0 : DEBUT
!    Eliminer ces lignes pour la clarte si on a fait le choix ICALKE(IZONE) = 1

  if(icalke(izone).eq.0) then

!         Calcul de k et epsilon en entree (XKENT et XEENT) a partir
!           l'intensite turbulente et de lois standards en conduite
!           circulaire (leur initialisation est inutile mais plus
!           propre)
    uref2 = rcodcl(ifac,iu(iphas),1)**2                           &
           +rcodcl(ifac,iv(iphas),1)**2                           &
           +rcodcl(ifac,iw(iphas),1)**2
    uref2 = max(uref2,1.d-12)
    xkent  = epzero
    xeent  = epzero

    call keenin                                                   &
    !==========
      ( uref2, xintur(izone), dh(izone), cmu, xkappa,             &
        xkent, xeent )

!     (ITYTUR est un indicateur qui vaut ITURB/10)
    if    (itytur(iphas).eq.2) then

      rcodcl(ifac,ik(iphas),1)  = xkent
      rcodcl(ifac,iep(iphas),1) = xeent

    elseif(itytur(iphas).eq.3) then

      rcodcl(ifac,ir11(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir22(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir33(iphas),1) = d2s3*xkent
      rcodcl(ifac,ir12(iphas),1) = 0.d0
      rcodcl(ifac,ir13(iphas),1) = 0.d0
      rcodcl(ifac,ir23(iphas),1) = 0.d0
      rcodcl(ifac,iep(iphas),1)  = xeent

    elseif (iturb(iphas).eq.50) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iep(iphas),1)  = xeent
      rcodcl(ifac,iphi(iphas),1) = d2s3
      rcodcl(ifac,ifb(iphas),1)  = 0.d0

    elseif (iturb(iphas).eq.60) then

      rcodcl(ifac,ik(iphas),1)   = xkent
      rcodcl(ifac,iomg(iphas),1) = xeent/cmu/xkent

    elseif (iturb(iphas).eq.70) then

      rcodcl(ifac,inusa(iphas),1) = cmu*xkent**2/xeent

    endif

  endif
! Exemple de cas ou ICALKE(IZONE) = 0 : FIN

! ------ Traitement des scalaires physiques particulieres
!        Ils sont traites automatiquement


! ------ Traitement des scalaires utilisateurs

  if ( (nscal-nscapp).gt.0 ) then
    do ii = 1, (nscal-nscapp)
      rcodcl(ifac,isca(ii),1) = 1.d0
    enddo
  endif

enddo

! --- On impose en couleur 15 une paroi laterale

CALL GETFBR('15',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          PAROI : DEBIT NUL (FLUX NUL POUR LA PRESSION)
!                  FROTTEMENT POUR LES VITESSES (+GRANDEURS TURB)
!                  FLUX NUL SUR LES SCALAIRES

!   Type de condition aux limites pour les variables standard
  itypfb(ifac,iphas)   = iparoi


!   Numero de zone (on les numerote de 1 a n)
  izone = 2

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

enddo

! --- On impose en couleur 19 une sortie

CALL GETFBR('19',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SORTIE : FLUX NUL VITESSE ET TEMPERATURE, PRESSION IMPOSEE
!            Noter que la pression sera recalee a P0
!                sur la premiere face de sortie libre (ISOLIB)

!   Type de condition aux limites pour les variables standard
  itypfb(ifac,iphas)   = isolib

!   Numero de zone (on les numerote de 1 a n)
  izone = 3

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

enddo


! --- On impose en couleur 14 et 4 une symetrie

CALL GETFBR('14 or 4',NLELT,LSTELT)
!==========

do ilelt = 1, nlelt

  ifac = lstelt(ilelt)

!          SYMETRIES

!   Type de condition aux limites pour les variables standard
  itypfb(ifac,iphas)   = isymet

!   Numero de zone (on les numerote de 1 a n)
  izone = 4

!      - Reperage de la zone a laquelle appartient la face
  izfppp(ifac) = izone

enddo


!----
! FIN
!----

return
end subroutine
