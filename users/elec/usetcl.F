c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS
C
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2008 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
                        SUBROUTINE USETCL
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ITRIFB , ITYPFB , IZFPPP ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     ROUTINE UTILISATEUR POUR PHYSIQUE PARTICULIERE
CFONC                 MODULE ELECTRIQUE
CFONC    (Effet Joule, Arc Electrique, Conduction ionique)
CFONC     REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
CFONC     (ICODCL,RCODCL) POUR LES VARIABLES INCONNUES
CFONC     PENDANT DE USCLIM.F
CFONC
CFONC
CFONC
CFONC     CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE
CFONC     =============================================
CFONC
CFONC
CFONC
CFONC INTRODUCTION
CFONC ============
CFONC
CFONC On donne ici les conditions aux limites par face de bord.
CFONC
CFONC L'identification des faces de bord concernees se fait grace
CFONC a la commande GETFBR.
CFONC
CFONC   GETFBR(CHAINE,NLELT,LSTELT) :
CFONC   - CHAINE est une chaine de caractere fournie par l'utilisateur
CFONC     qui donne les criteres de selection
CFONC   - NLTELT est renvoye par la commande. C'est un entier qui
CFONC     correspond au nombre de faces de bord trouveees repondant au
CFONC     critere
CFONC   - LSTELT est renvoye par la commande. C'est un tableau d'entiers
CFONC     de taille NLTELT donnant la liste des faces de bord trouvees
CFONC     repondant au critere.
CFONC
CFONC   CHAINE peut etre constitue de :
CFONC   - references de couleurs (ex. : 1, 8, 26, ...
CFONC   - references de groupes (ex. : entrees, groupe1, ...)
CFONC   - criteres geometriques (ex. X<0.1, Y>=0.25, ...)
CFONC   Ces criteres peuvent etre combines par des operateurs logiques
CFONC   (AND et OR) et des parentheses
CFONC   ex. : '1 AND (groupe2 OR groupe3) AND Y<1' permettra de recuperer
CFONC   les faces de bord de couleur 1, appartenant aux groupes 'groupe2'
CFONC   ou 'groupe3' et de coordonnee Y inferieure a 1.
CFONC
CFONC
CFONC
CFONC TYPE DE CONDITIONS AUX LIMITES
CFONC ==============================
CFONC
CFONC On peut affecter les conditions aux limites de deux manieres.
CFONC
CFONC
CFONC     Pour les conditions "standard" :
CFONC     --------------------------------
CFONC
CFONC      (entree, sortie libre, paroi, symetrie), on donne un code
CFONC      stocke dans le tableau ITYPFB (dimensionne au nombre de
CFONC      faces de bord,nombre de phases). Ce code sera ensuite
CFONC      utilise par un sous-programme non utilisateur pour affecter
CFONC      les conditions correspondantes (les scalaires, en
CFONC      particulier, recevront alors les conditions de la phase a
CFONC      laquelle ils sont associes). Ainsi :
CFONC
CFONC      Code      |  Type de bord
CFONC      -------------------------
CFONC       IENTRE   |   Entree
CFONC       ISOLIB   |   Sortie libre
CFONC       ISYMET   |   Symetrie
CFONC       IPAROI   |   Paroi (lisse)
CFONC       IPARUG   |   Paroi rugueuse
CFONC
CFONC      Les entiers IENTRE, ISOLIB, ISYMET, IPAROI, IPARUG
CFONC      sont affectes par ailleurs (include param.h). Leur valeur
CFONC      est superieure ou egale a 1 et
CFONC          inferieure ou egale a NTYPMX (valeur fixee dans paramx.h)
CFONC
CFONC
CFONC
CFONC      En outre, il faut donner certaines valeurs :
CFONC
CFONC
CFONC      - Entree (plus precisement entree/sortie a debit impose, car le debit
CFONC                peut etre impose sortant) :
CFONC
CFONC        -> Conditions de Dirichlet sur les variables
CFONC          autres que la pression obligatoire si le flux est entrant,
CFONC          optionnelle si le flux est sortant (le code affecte flux nul
CFONC          si aucun Dirichlet n'est specifie) ; ainsi
CFONC          en face IFAC, sur la variable IVAR : RCODCL(IFAC,IVAR,1)
CFONC
CFONC
CFONC      - Paroi lisse : (= solide impermeable, avec frottement lisse)
CFONC
CFONC        -> Valeur de la vitesse de paroi defilante s'il y a lieu
CFONC          en face IFAC, RCODCL(IFAC,IU,1)
CFONC                        RCODCL(IFAC,IV,1)
CFONC                        RCODCL(IFAC,IW,1)
CFONC        -> Code specifique et valeur de la temperature imposee
CFONC          en paroi s'il y a lieu :
CFONC          en face IFAC, ICODCL(IFAC,IVAR)   = 5
CFONC                        RCODCL(IFAC,IVAR,1) = Temperature imposee
CFONC        -> Code specifique et valeur du flux imposee
CFONC          en paroi s'il y a lieu :
CFONC          en face IFAC, ICODCL(IFAC,IVAR)   = 3
CFONC                        RCODCL(IFAC,IVAR,3) = Flux impose
CFONC                                         =
CFONC         Noter que la condition par defaut pour les scalaires
CFONC          (hors k et epsilon) est un Neumann homogene
CFONC
CFONC
CFONC      - Paroi rugueuse : (= solide impermeable, avec frottement rugueux)
CFONC
CFONC        -> Valeur de la vitesse de paroi defilante s'il y a lieu
CFONC          en face IFAC, RCODCL(IFAC,IU,1)
CFONC                        RCODCL(IFAC,IV,1)
CFONC                        RCODCL(IFAC,IW,1)
CFONC        -> Valeur de la hauteur de rugosite dynamique a specifier dans
CFONC                        RCODCL(IFAC,IU,3) (valeur pour IV et IW non utilisee)
CFONC        -> Code specifique et valeur de la temperature imposee
CFONC          en paroi rugueuse s'il y a lieu :
CFONC          en face IFAC, ICODCL(IFAC,IVAR)   = 6
CFONC                        RCODCL(IFAC,IVAR,1) = Temperature imposee
CFONC                        RCODCL(IFAC,IVAR,3) = Hauteur de rugosite thermique
CFONC        -> Code specifique et valeur du flux imposee
CFONC          en paroi s'il y a lieu :
CFONC          en face IFAC, ICODCL(IFAC,IVAR)   = 3
CFONC                        RCODCL(IFAC,IVAR,3) = Flux impose
CFONC                                         =
CFONC         Noter que la condition par defaut pour les scalaires
CFONC          (hors k et epsilon) est un Neumann homogene
CFONC
CFONC      - Symetrie (= paroi impermeable avec glissement) :
CFONC
CFONC        -> Rien a preciser
CFONC
CFONC
CFONC      - Sortie libre (plus precisement entree/sortie libre a pression imposee)
CFONC
CFONC        -> Rien a preciser pour la pression et la vitesse
CFONC           Pour les scalaires et grandeurs turbulentes, une valeur de Dirichlet
CFONC             peut etre specifiee de maniere optionnelle. Le comportement est le
CFONC             suivant :
CFONC               * la pression est toujours traitee en Dirichlet
CFONC               * si flux de masse entrant :
CFONC                   on retient la vitesse a l'infini
CFONC                   condition de Dirichlet pour les scalaires et grandeurs
CFONC                     turbulentes (ou flux nul si l'utilisateur n'a pas
CFONC                     specifie de Dirichlet)
CFONC                  si flux de masse sortant :
CFONC                     on impose un flux nul sur la vitesse, les scalaires et
CFONC                     les grandeurs turbulentes
CFONC
CFONC        Noter que la pression sera recalee a P0
CFONC            sur la premiere face de sortie libre trouvee
CFONC
CFONC
CFONC     Pour les conditions "non standard" :
CFONC     ------------------------------------
CFONC
CFONC      Autres que (entree, sortie libre, paroi, symetrie), on donne
CFONC       - d'une part, pour chaque face :
CFONC         -> une valeur de ITYPFB admissible
CFONC            ie superieure ou egale a 1 et inferieure ou egale a
CFONC            NTYPMX (voir sa valeur dans paramx.h).
CFONC            Les valeurs predefinies dans paramx.h
CFONC            IENTRE, ISOLIB, ISYMET, IPAROI, IPARUG sont dans cet
CFONC            intervalle et il est preferable de ne pas affecter
CFONC            inconsidrement et par hasard a ITYPFB la valeur
CFONC            d'un de ces entiers. Pour eviter cela, on peut
CFONC            utiliser IINDEF si l'on souhaite eviter de verifier
CFONC            les valeurs dans paramx.h. IINDEF est une valeur
CFONC            admissible a laquelle n'est attachee aucune condition
CFONC            limite predefinie.
CFONC            Noter que le tableau ITYPFB est
CFONC            reinitialise a chaque pas de temps a la valeur
CFONC            non admissible 0. Si on oublie de modifier ITYPFB pour
CFONC            une face, le code s'arretera.
CFONC
CFONC       - et d'autre part pour chaque face et chaque variable :
CFONC         -> un code     ICODCL(IFAC,IVAR)
CFONC         -> trois reels RCODCL(IFAC,IVAR,1)
CFONC                        RCODCL(IFAC,IVAR,2)
CFONC                        RCODCL(IFAC,IVAR,3)
CFONC      La valeur de ICODCL est prise parmi les suivantes :
CFONC        1 : Dirichlet      (utilisable pour toute variable)
CFONC        3 : Neumann        (utilisable pour toute variable)
CFONC        4 : Symetrie       (utilisable uniquement pour la vitesse et
CFONC                                    les composantes du tenseur Rij)
CFONC        5 : Paroi lisse    (utilisable pour toute variable sauf la
CFONC                                                          pression)
CFONC        6 : Paroi rugueuse (utilisable pour toute variable sauf la
CFONC                                                          pression)
CFONC        9 : Sortie libre   (utilisable uniquement pour la vitesse)
CFONC      Les valeurs des 3 reels RCODCL sont les suivantes
CFONC       RCODCL(IFAC,IVAR,1) :
CFONC          Dirichlet sur la variable        si ICODCL(IFAC,IVAR)=  1
CFONC          valeur en paroi (defilmnt, temp) si ICODCL(IFAC,IVAR)=  5
CFONC          La dimension de RCODCL(IFAC,IVAR,1) est celle de la
CFONC            variable resolue : ex U (vitesse en m/s),
CFONC                                  T (temperature en degres)
CFONC                                  H (enthalpie en J/kg)
CFONC                                  F (scalaire passif en -)
CFONC       RCODCL(IFAC,IVAR,2) :
CFONC          coefficient d'echange "exterieur" (entre la valeur
CFONC                           imposee et la valeur au bord du domaine)
CFONC                           RINFIN = infini par defaut
CFONC          Pour les vitesses U,             en kg/(m2 s) :
CFONC            RCODCL(IFAC,IVAR,2) =            (VISCL+VISCT) / D
CFONC          Pour la  pression P,             en  s/m          :
CFONC            RCODCL(IFAC,IVAR,2) =                       DT / D
CFONC          Pour les temperatures T,         en Watt/(m2 degres) :
CFONC            RCODCL(IFAC,IVAR,2) = CP*(VISCLS+VISCT/SIGMAS) / D
CFONC          Pour les enthalpies H,           en kg /(m2 s) :
CFONC            RCODCL(IFAC,IVAR,2) =    (VISCLS+VISCT/SIGMAS) / D
CFONC          Pour les autres scalaires F      en :
CFONC            RCODCL(IFAC,IVAR,2) =    (VISCLS+VISCT/SIGMAS) / D
CFONC               (D a la dimension d'une distance en m)
CFONC
CFONC       RCODCL(IFAC,IVAR,3) si ICODCL(IFAC,IVAR)<>6 :
CFONC         Densite de flux (< 0 si gain, n normale orientee vers l'exterieur)
CFONC                          si ICODCL(IFAC,IVAR)= 3
CFONC          Pour les vitesses U,             en kg/(m s2) = J :
CFONC            RCODCL(IFAC,IVAR,3) =           -(VISCL+VISCT) * (GRAD U).n
CFONC          Pour la  pression P,             en kg/(m2 s)     :
CFONC            RCODCL(IFAC,IVAR,3) =                      -DT * (GRAD P).n
CFONC          Pour les temperatures T,         en Watt/m2       :
CFONC            RCODCL(IFAC,IVAR,3) =-CP*(VISCLS+VISCT/SIGMAS) * (GRAD T).n
CFONC          Pour les enthalpies H,           en Watt/m2       :
CFONC            RCODCL(IFAC,IVAR,3) =   -(VISCLS+VISCT/SIGMAS) * (GRAD H).n
CFONC          Pour les autres scalaires F      en :
CFONC            RCODCL(IFAC,IVAR,3) =   -(VISCLS+VISCT/SIGMAS) * (GRAD F).n
CFONC
CFONC       RCODCL(IFAC,IVAR,3) SI ICODCL(IFAC,IVAR)=6 :
CFONC         Rugosites de la loi rugueuse
CFONC          Pour les vitesses U, rugosite dynamique
CFONC            RCODCL(IFAC,IVAR,3) = RUGD
CFONC          Pour les autres scalaires, rugosite thermique
CFONC            RCODCL(IFAC,IVAR,3) = RUGT
CFONC
CFONC
CFONC       Noter bien que si l'utilisateur affecte une valeur a ITYPFB
CFONC        parmi     IENTRE, ISOLIB, ISYMET, IPAROI, IPARUG
CFONC        et qu'il ne modifie pas ICODCL (valeur nulle par defaut),
CFONC        c'est ITYPFB qui imposera la condition limite.
CFONC
CFONC       Par contre, si l'utilisateur impose
CFONC         ICODCL(IFAC,IVAR) (non nul),
CFONC         ce sont alors les valeurs de RCODCL qu'il aura fournies
CFONC         qui sont retenues pour la face et la variable consideree
CFONC         (s'il ne precise pas RCODCL, ce sont les valeurs
CFONC         par defaut qui sont retenues pour la face et
CFONC         la variable consideree soit :
CFONC                                  RCODCL(IFAC,IVAR,1) = 0.D0
CFONC                                  RCODCL(IFAC,IVAR,2) = RINFIN
CFONC                                  RCODCL(IFAC,IVAR,3) = 0.D0)
CFONC         En particulier, on peut par exemple
CFONC         -> donner ITYPFB(IFAC,IPHAS) = IPAROI
CFONC         qui impose les conditions de paroi par defaut pour toutes
CFONC         les variables sur la face IFAC,
CFONC         -> et preciser EN OUTRE pour la variable IVAR sur cette
CFONC         face IFAC des conditions paarticulieres en imposant
CFONC         ICODCL(IFAC,IVAR) et les 3 RCODCL.
CFONC
CFONC
CFONC       L'utilisateur peut egalement affecter a ITYPFB une valeur
CFONC        non egale a IENTRE, ISOLIB, ISYMET, IPAROI, IPARUG, IINDEF
CFONC        mais superieure ou egale a 1 et inferieure ou egale a
CFONC        NTYPMX (voir les valeurs dans param.h) pour reperer
CFONC        des groupes de couleurs dans d'autres sous programmes
CFONC        qui lui seraient propres et ou ITYPFB serait disponible.
CFONC        Dans ce cas cependant, il faudra
CFONC        imposer les conditions limites en donnant des valeurs a
CFONC        ICODCL et aux trois RCODCL (puisque la valeur de ITYPFB
CFONC        ne sera pas predefinie dans le code).
CFONC
CFONC
CFONC REGLES DE COHERENCE
CFONC ===================
CFONC
CFONC        Quelques regles de coherence entre les codes ICODCL
CFONC          des variables pour les conditions non standard :
CFONC
CFONC            Les codes des vitesses doivent etre identiques
CFONC            Les codes des Rij      doivent etre identiques
CFONC            Si code (vitesse ou Rij) = 4
CFONC              il faut code (vitesse et Rij) = 4
CFONC            Si code (vitesse ou turbulence) = 5
CFONC              il faut code (vitesse et turbulence) = 5
CFONC            Si code (vitesse ou turbulence) = 6
CFONC              il faut code (vitesse et turbulence) = 6
CFONC            Si code scalaire (hormis pression ou fluctuations) = 5
CFONC              il faut code vitesse = 5
CFONC            Si code scalaire (hormis pression ou fluctuations) = 6
CFONC              il faut code vitesse = 6
CFONC
CFONC
CFONC REMARQUES
CFONC ==========
CFONC
CFONC        Attention : pour imposer un flux (non nul) sur les Rij,
CFONC                    la viscosite a prendre en compte est VISCL
CFONC                    meme si VISCT existe (VISCT=RHO CMU K2/EPSILON)
CFONC
CFONC
CFONC        On dispose du tableau de tri des faces de bord au pas
CFONC            de temps precedent (sauf au premier pas de temps, ou
CFONC            ITRIFB n'a pas ete renseigne).
CFONC        Le tableau du type des faces de bord ITYPFB a ete
CFONC            reinitialise avant d'entrer dans le sous programme.
CFONC
CFONC
CFONC
CFONC        Noter comment acceder a certaines variables :
CFONC
CFONC  Valeurs aux cellules
CFONC                Soit        IEL = IFABOR(IFAC)
CFONC
CFONC * Masse vol                       phase IPHAS, cellule      IEL  :
CFONC                   PROPCE(IEL ,IPPROC(IROM (IPHAS)))
CFONC * Viscosite moleculaire dynamique phase IPHAS, cellule      IEL  :
CFONC                   PROPCE(IEL ,IPPROC(IVISCL(IPHAS)))
CFONC * Viscosite turbulente  dynamique phase IPHAS, cellule      IEL  :
CFONC                   PROPCE(IEL ,IPPROC(IVISCT(IPHAS)))
CFONC * Chaleur specifique              phase IPHAS, cellule      IEL  :
CFONC                   PROPCE(IEL ,IPPROC(ICP   (IPHAS)))
CFONC * Diffusivite lambda           scalaire ISCAL, cellule      IEL  :
CFONC                   PROPCE(IEL ,IPPROC(IVISLS(ISCAL)))
CFONC
CFONC  Valeurs aux faces de bord
CFONC
CFONC  * Masse vol                      phase IPHAS, face de bord IFAC :
CFONC                   PROPFB(IFAC,IPPROB(IROM (IPHAS)))
CFONC  * Flux de masse relatif a la variable  IVAR , face de bord IFAC :
CFONC       (i.e. le flux de masse servant a la convection de IVAR)
CFONC                   PROPFB(IFAC,IPPROB(IFLUMA(IVAR )))
CFONC  * Pour les autres grandeurs              a la face de bord IFAC :
CFONC       prendre pour approximation la valeur dans la cellule  IEL
CFONC       adjacente i.e. comme plus haut avec IEL = IFABOR(IFAC).
CFONC
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! NDIM         ! E  !  -> ! DIMENSION DE L'ESPACE                !
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! NFML         ! E  !  -> ! NOMBRE DE FAMILLES D ENTITES         !
CARGU ! NPRFML       ! E  !  -> ! NOMBRE DE PROPRIETESE DES FAMILLES   !
CARGU ! NNOD         ! E  !  -> ! NOMBRE DE SOMMETS                    !
CARGU ! LNDFAC       ! E  !  -> ! LONGUEUR DU TABLEAU NODFAC (OPTIONNEL!
CARGU ! LNDFBR       ! E  !  -> ! LONGUEUR DU TABLEAU NODFBR (OPTIONNEL!
CARGU ! NCELBR       ! E  !  -> ! NOMBRE D'ELEMENTS AYANT AU MOINS UNE !
CARGU !              !    !     ! FACE DE BORD                         !
CARGU ! NVAR         ! E  !  -> ! NOMBRE TOTAL DE VARIABLES            !
CARGU ! NSCAL        ! E  !  -> ! NOMBRE TOTAL DE SCALAIRES            !
CARGU ! NPHAS        ! E  !  -> ! NOMBRE DE PHASES                     !
CARGU ! NIDEVE NRDEVE! E  !  -> ! LONGUEUR DE IDEVEL RDEVEL            !
CARGU ! NITUSE NRTUSE! E  !  -> ! LONGUEUR DE ITUSER RTUSER            !
CARGU ! IFACEL       ! TE !  -> ! ELEMENTS VOISINS D'UNE FACE INTERNE  !
CARGU ! (2, NFAC)    !    !     !                                      !
CARGU ! IFABOR       ! TE !  -> ! ELEMENT  VOISIN  D'UNE FACE DE BORD  !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMFBR       ! TE !  -> ! NUMERO DE FAMILLE D'UNE FACE DE BORD !
CARGU ! (NFABOR)     !    !     !                                      !
CARGU ! IFMCEL       ! TE !  -> ! NUMERO DE FAMILLE D'UNE CELLULE      !
CARGU ! (NCELET)     !    !     !                                      !
CARGU ! IPRFML       ! TE !  -> ! PROPRIETES D'UNE FAMILLE             !
CARGU ! NFML  ,NPRFML!    !     !                                      !
CARGU ! MAXELT       !  E !  -> ! NB MAX D'ELEMENTS (CELL,FAC,FBR)     !
CARGU ! LSTELT(MAXELT) TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! IPNFAC       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFAC)   !    !     !  FACE INTERNE DANS NODFAC (OPTIONNEL)!
CARGU ! NODFAC       ! TE !  -> ! CONNECTIVITE FACES INTERNES/NOEUDS   !
CARGU !   (NFAC+1)   !    !     !  (OPTIONNEL)                         !
CARGU ! IPNFBR       ! TE !  -> ! POSITION DU PREMIER NOEUD DE CHAQUE  !
CARGU !   (LNDFBR)   !    !     !  FACE DE BORD DANS NODFBR (OPTIONNEL)!
CARGU ! NODFBR       ! TE !  -> ! CONNECTIVITE FACES DE BORD/NOEUDS    !
CARGU !   (NFABOR+1) !    !     !  (OPTIONNEL)                         !
CARGU ! ICODCL       ! TE ! <-  ! CODE DE CONDITION LIMITES AUX FACES  !
CARGU !  (NFABOR,NVAR!    !     !  DE BORD                             !
CARGU !              !    !     ! = 1   -> DIRICHLET                   !
CARGU !              !    !     ! = 3   -> DENSITE DE FLUX             !
CARGU !              !    !     ! = 4   -> GLISSEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 5   -> FROTTEMT ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 6   -> RUGOSITE ET U.n=0 (VITESSE) !
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
CARGU ! ITRIFB(NFABOR! TE !  -> ! INDIRECTION POUR TRI DES FACES DE BRD!
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! ITYPFB(NFABOR! TE ! <-  ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! IZFPPP       ! TE ! <-  ! NUMERO DE ZONE DE LA FACE DE BORD    !
CARGU ! (NFABOR)     !    !     !  POUR LE MODULE PHYS. PART.          !
CARGU ! IDEVEL(NIDEVE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE DEVELOPEMT !
CARGU ! ITUSER(NITUSE! TE ! <-> ! TAB ENTIER COMPLEMENTAIRE UTILISATEUR!
CARGU ! IA(*)        ! TR !  -  ! MACRO TABLEAU ENTIER                 !
CARGU ! XYZCEN       ! TR !  -> ! POINT ASSOCIES AUX VOLUMES DE CONTROL!
CARGU ! (NDIM,NCELET !    !     !                                      !
CARGU ! SURFAC       ! TR !  -> ! VECTEUR SURFACE DES FACES INTERNES   !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! SURFBO       ! TR !  -> ! VECTEUR SURFACE DES FACES DE BORD    !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! CDGFAC       ! TR !  -> ! CENTRE DE GRAVITE DES FACES INTERNES !
CARGU ! (NDIM,NFAC)  !    !     !                                      !
CARGU ! CDGFBO       ! TR !  -> ! CENTRE DE GRAVITE DES FACES DE BORD  !
CARGU ! (NDIM,NFABOR)!    !     !                                      !
CARGU ! XYZNOD       ! TR !  -> ! COORDONNES DES NOEUDS (OPTIONNEL)    !
CARGU ! (NDIM,NNOD)  !    !     !                                      !
CARGU ! VOLUME       ! TR !  -> ! VOLUME D'UN DES NCELET ELEMENTS      !
CARGU ! (NCELET      !    !     !                                      !
CARGU ! DT(NCELET)   ! TR !  -> ! PAS DE TEMPS                         !
CARGU ! RTP, RTPA    ! TR !  -> ! VARIABLES DE CALCUL AU CENTRE DES    !
CARGU ! (NCELET,*)   !    !     !    CELLULES (INSTANT COURANT OU PREC)!
CARGU ! PROPCE       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU ! (NCELET,*)   !    !     !    CELLULES                          !
CARGU ! PROPFA       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFAC,*)    !    !     !    FACES INTERNES                    !
CARGU ! PROPFB       ! TR !  -> ! PROPRIETES PHYSIQUES AU CENTRE DES   !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! COEFA, COEFB ! TR ! <-> ! CONDITIONS AUX LIMITES AUX           !
CARGU !  (NFABOR,*)  !    !     !    FACES DE BORD                     !
CARGU ! RCODCL       ! TR ! <-  ! VALEUR DES CONDITIONS AUX LIMITES    !
CARGU !  (NFABOR,NVAR!    !     !  AUX FACES DE BORD                   !
CARGU !              !    !     ! RCODCL(1) = VALEUR DU DIRICHLET      !
CARGU !              !    !     ! RCODCL(2) = VALEUR DU COEF. D'ECHANGE!
CARGU !              !    !     !  EXT. (INFINIE SI PAS D'ECHANGE)     !
CARGU !              !    !     ! RCODCL(3) = VALEUR DE LA DENSITE DE  !
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2 OU      !
CARGU !              !    !     !  HAUTEUR DE RUGOSITE (m) si ICODCL=6 !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! COEFU        ! TR !  -  ! TAB DE TRAV                          !
CARGU !  (NFABOR,3)  !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! RDEVEL(NRDEVE! TR ! <-> ! TAB REEL COMPLEMENTAIRE DEVELOPEMT   !
CARGU ! RTUSER(NRTUSE! TR ! <-> ! TAB REEL COMPLEMENTAIRE UTILISATEUR  !
CARGU ! RA(*)        ! TR !  -  ! MACRO TABLEAU REEL                   !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "paramx.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "ppppar.h"
      INCLUDE "ppthch.h"
      INCLUDE "ppincl.h"
      INCLUDE "elincl.h"
C
C***********************************************************************
C
C ARGUMENTS
C
      INTEGER          IDBIA0 , IDBRA0
      INTEGER          NDIM   , NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          NFML   , NPRFML
      INTEGER          NNOD   , LNDFAC , LNDFBR , NCELBR
      INTEGER          NVAR   , NSCAL  , NPHAS
      INTEGER          NIDEVE , NRDEVE , NITUSE , NRTUSE
C
      INTEGER          IFACEL(2,NFAC) , IFABOR(NFABOR)
      INTEGER          IFMFBR(NFABOR) , IFMCEL(NCELET)
      INTEGER          IPRFML(NFML,NPRFML)
      INTEGER          MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ITRIFB(NFABOR,NPHAS), ITYPFB(NFABOR,NPHAS)
      INTEGER          IZFPPP(NFABOR)
      INTEGER          IDEVEL(NIDEVE), ITUSER(NITUSE), IA(*)
C
      DOUBLE PRECISION XYZCEN(NDIM,NCELET)
      DOUBLE PRECISION SURFAC(NDIM,NFAC), SURFBO(NDIM,NFABOR)
      DOUBLE PRECISION CDGFAC(NDIM,NFAC), CDGFBO(NDIM,NFABOR)
      DOUBLE PRECISION XYZNOD(NDIM,NNOD), VOLUME(NCELET)
      DOUBLE PRECISION DT(NCELET), RTP(NCELET,*), RTPA(NCELET,*)
      DOUBLE PRECISION PROPCE(NCELET,*)
      DOUBLE PRECISION PROPFA(NFAC,*), PROPFB(NFABOR,*)
      DOUBLE PRECISION COEFA(NFABOR,*), COEFB(NFABOR,*)
      DOUBLE PRECISION RCODCL(NFABOR,NVAR,3)
      DOUBLE PRECISION W1(NCELET),W2(NCELET),W3(NCELET)
      DOUBLE PRECISION W4(NCELET),W5(NCELET),W6(NCELET)
      DOUBLE PRECISION COEFU(NFABOR,NDIM)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC  , II     , IPHAS , IEL
      INTEGER          I     , NTF    , NB    , ID , ITROUV
      INTEGER          IZONE
      INTEGER          NBORNE(NBTRMX)
      INTEGER          ILELT , NLELT
C
      DOUBLE PRECISION RNBS2,CAPAEQ
      DOUBLE PRECISION SIR(NELEMX)   ,SII(NELEMX)
      DOUBLE PRECISION SIRB(NBTRMX,6),SIIB(NBTRMX,6)
      DOUBLE PRECISION UR(NBTRMX,6)  ,UI(NBTRMX,6)
      DOUBLE PRECISION SIRT(NBTRMX)  ,SIIT(NBTRMX)
      CHARACTER*200    CHAIN
C
C***********************************************************************
C=======================================================================
C 1.  INITIALISATIONS
C
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
C=======================================================================
C 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
C       ON BOUCLE SUR LES FACES DE BORD
C         ON DETERMINE LA FAMILLE ET SES PROPRIETES
C           ON IMPOSE LA CONDITION LIMITE
C
C          IMPOSER ICI LES CONDITIONS LIMITES SUR LES FACES DE BORD
C
C=======================================================================
C
C 2.1 - Calcul des intensites de courant sur chaque electrode
C       -----------------------------------------------------
C
C   Mise a zero
C
      DO I= 1,NBELEC
        SIR(I) = 0.D0
        SII(I) = 0.D0
      ENDDO
C
      DO NTF= 1,NBTRF
        SIRT(NTF) = 0.D0
        SIIT(NTF) = 0.D0
      ENDDO
C
      IF(NTCABS.LT.(NTPABS+2)) THEN
        DO NTF = 1,NBTRF
          UROFF(NTF) = 0.D0
          UIOFF(NTF) = 0.D0
        ENDDO
      ENDIF
C
C     Boucle sur les faces de bord selectionnees
C
      IPHAS = 1
C
      DO I=1,NBELEC
C
        CHAIN = ' '
        WRITE(CHAIN,3000) IELECC(I)
C
        IF ( IELECT(I).NE. 0 ) THEN
          CALL GETFBR(CHAIN,NLELT,LSTELT)
C         ===========
C
          DO ILELT = 1, NLELT
C
            IFAC = LSTELT(ILELT)
C
            IEL = IFABOR(IFAC)
C
            DO ID=1,NDIMVE
              SIR(I) = SIR(I)
     &             +PROPCE(IEL,IPPROC(IDJR(ID)))*SURFBO(ID,IFAC)
            ENDDO
C
            IF ( IPPMOD(IELJOU) .EQ. 4 ) THEN
              DO ID=1,NDIMVE
                SII(I) = SII(I)
     &               +PROPCE(IEL,IPPROC(IDJI(ID)))*SURFBO(ID,IFAC)
              ENDDO
            ENDIF
C
          ENDDO
C
        ENDIF
C
      ENDDO
C
C
C 2.2 - Definitions des CL sur les tensions de chaque bornes des transfo
C       ----------------------------------------------------------------
C
C  2.2.1 Calcul des intensites de courant sur chaque borne des transfo
C
      DO I=1,NBELEC
        SIRB(IELECT(I),IELECB(I)) = 0.D0
        IF ( IPPMOD(IELJOU) .EQ. 4 ) THEN
          SIIB(IELECT(I),IELECB(I)) = 0.D0
        ENDIF
      ENDDO
C
      DO I=1,NBELEC
        IF ( IELECT(I).NE. 0 ) THEN
          SIRB(IELECT(I),IELECB(I)) = SIRB(IELECT(I),IELECB(I))
     &                               +SIR(I)
          IF ( IPPMOD(IELJOU) .EQ. 4 ) THEN
             SIIB(IELECT(I),IELECB(I)) = SIIB(IELECT(I),IELECB(I))
     &                                  +SII(I)
          ENDIF
        ENDIF
      ENDDO
C
C  2.2.2 Tensions sur chaque bornes des transfo
C
      DO NTF=1,NBTRF
C
C      Primaire et Secondaire en Triangle
C
        IF (IBRPR(NTF) .EQ. 0 .AND. IBRSEC(NTF) .EQ. 0 ) THEN
C
          NBORNE(NTF) = 3
C
          RNBS2 = 3.D0*RNBS(NTF)*RNBS(NTF)
          UR(NTF,1)=  1.154675D0*TENSPR(NTF)/RNBS(NTF)
     &      + (ZR(NTF)*SIRB(NTF,1)-ZI(NTF)*SIIB(NTF,1))/RNBS2
C
          UR(NTF,2)= -0.5773D0*TENSPR(NTF)/RNBS(NTF)
     &      + (ZR(NTF)*SIRB(NTF,2)-ZI(NTF)*SIIB(NTF,2))/RNBS2
          UR(NTF,3)= -0.5773D0*TENSPR(NTF)/RNBS(NTF)
     &      + (ZR(NTF)*SIRB(NTF,3)-ZI(NTF)*SIIB(NTF,3))/RNBS2
C
          UI(NTF,1)=  0.D0
     &      + (ZI(NTF)*SIRB(NTF,1)+ZR(NTF)*SIIB(NTF,1))/RNBS2
          UI(NTF,2)= -1.D0*TENSPR(NTF)/RNBS(NTF)
     &      + (ZI(NTF)*SIRB(NTF,2)+ZR(NTF)*SIIB(NTF,2))/RNBS2
          UI(NTF,3)=  1.D0*TENSPR(NTF)/RNBS(NTF)
     &      + (ZI(NTF)*SIRB(NTF,3)+ZR(NTF)*SIIB(NTF,3))/RNBS2
C
        ELSE
C
          WRITE(NFECRA,*) 'Matrice sur le Transfo a ecrire'
          CALL CSEXIT(1)
C
        ENDIF
      ENDDO
C
C  2.2.3 Intensite de courant totale par transfo
C
      DO NTF=1,NBTRF
        SIRT(NTF) = 0.D0
        IF ( IPPMOD(IELJOU) .EQ. 4 ) THEN
          SIIT(NTF) = 0.D0
        ENDIF
      ENDDO
C
      DO I=1,NBELEC
        IF ( IELECT(I).NE. 0 ) THEN
          SIRT(IELECT(I)) = SIRT(IELECT(I)) + SIR(I)
          IF ( IPPMOD(IELJOU) .EQ. 4 ) THEN
            SIIT(IELECT(I)) = SIIT(IELECT(I)) + SII(I)
          ENDIF
        ENDIF
      ENDDO
C
C  2.2.4 Prise en compte des Offset
C
      CAPAEQ = 3.D0
C
      DO NTF=1,NBTRF
        UROFF(NTF) = UROFF(NTF) + SIRT(NTF)/CAPAEQ
        IF ( IPPMOD(IELJOU) .EQ. 4 ) THEN
          UIOFF(NTF) = UIOFF(NTF) + SIIT(NTF)/CAPAEQ
        ENDIF
      ENDDO
C
C Transfo de reference ===> Offset nul
C
      IF ( NTFREF .GT. 0 ) THEN
        UROFF(NTFREF) = 0.D0
        UIOFF(NTFREF) = 0.D0
      ENDIF
C
      DO NTF=1,NBTRF
        DO NB=1,NBORNE(NTF)
          UR(NTF,NB) = UR(NTF,NB) + UROFF(NTF)
          IF ( IPPMOD(IELJOU) .EQ. 4 ) THEN
            UI(NTF,NB) = UI(NTF,NB) + UIOFF(NTF)
          ENDIF
        ENDDO
      ENDDO
C
C Affichage des UROFF
C
      WRITE(NFECRA,1500)
      DO NTF=1,NBTRF
        WRITE(NFECRA,2000) NTF,UROFF(NTF)
      ENDDO
      WRITE(NFECRA,1501)
C
C  2.2.5 Prise en compte des CL
C
C
C     Boucle sur les faces de bord selectionnees
C
      IPHAS = 1
C
      DO I=1,NBELEC
C
        CHAIN = ' '
        WRITE(CHAIN,3000) IELECC(I)
C
        CALL GETFBR(CHAIN,NLELT,LSTELT)
C       ===========
C
        DO ILELT = 1, NLELT
C
          IFAC = LSTELT(ILELT)
C
          IEL = IFABOR(IFAC)
C
          ITYPFB(IFAC,IPHAS) = IPAROI
C
C     - Numero de zone (choix du numero de l'arc)
          IZONE = I
C
C      - Reperage de la zone a laquelle appartient la face
          IZFPPP(IFAC) = IZONE
C
          IF ( IELECT(I) .NE. 0 ) THEN
            ICODCL(IFAC,ISCA(IPOTR))   = 1
            RCODCL(IFAC,ISCA(IPOTR),1) = UR(IELECT(I),IELECB(I))
C
            IF ( IPPMOD(IELJOU).EQ.4  ) THEN
              ICODCL(IFAC,ISCA(IPOTI))   = 1
              RCODCL(IFAC,ISCA(IPOTI),1) = UI(IELECT(I),IELECB(I))
            ENDIF
C
          ELSE
C
            II = IPOTR
            ICODCL(IFAC,ISCA(II))   = 3
            RCODCL(IFAC,ISCA(II),3) = 0.D0
C
            IF ( IPPMOD(IELJOU).EQ. 4   ) THEN
              II = IPOTI
              ICODCL(IFAC,ISCA(II))   = 3
              RCODCL(IFAC,ISCA(II),3) = 0.D0
            ENDIF
          ENDIF
C
        ENDDO
C
      ENDDO
C
C 3 - Verification dans le cas ou pas de transfo de referenc qu'il existe
C     une paroi avec potentielle nul
C
      IF ( NTFREF .EQ. 0 ) THEN
C
        ITROUV = 0
        DO IFAC = 1, NFABOR
C
          IPHAS = 1
C
          IF ( ITYPFB(IFAC,IPHAS) .EQ. IPAROI ) THEN
C
            IF (ICODCL(IFAC,ISCA(IPOTR)) .EQ. 1 ) THEN
C
              IF ( IPPMOD(IELJOU).EQ.3 ) THEN
C
                IF ( ABS(RCODCL(IFAC,ISCA(IPOTR),1)).LT.1.E-20 ) THEN
                  ITROUV = 1
                ENDIF
C
              ELSE IF ( IPPMOD(IELJOU).EQ.4 ) THEN
C
                IF (ICODCL(IFAC,ISCA(IPOTI)) .EQ. 1 ) THEN
C
                  IF (ABS(RCODCL(IFAC,ISCA(IPOTR),1)).LT.1.E-20
     &           .AND.ABS(RCODCL(IFAC,ISCA(IPOTI),1)).LT.1.E-20 ) THEN
                    ITROUV = 1
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDIF
        ENDDO
C
        IF ( ITROUV .EQ. 0 ) THEN
          WRITE(NFECRA,1000)
          CALL CSEXIT (1)
        ENDIF
C
      ENDIF
C
C----
C FORMATS
C----
C
 1000 FORMAT(1X,' ERREUR VERSION JOULE : ',/,
     &       1X,' ====================   ',/,
     &      10X,' Vous n''avez pas donne de transfo de référence',/,
     &      10X,' et vous n''avez pas de CL a potentiel nul dans le',
     &          ' domaine')
 1500 FORMAT(/,2X,' ** INFORMATIONS SUR LES DIFFERENTS TRANSFO'/,
     &         2X,'    ---------------------------------------'/,/,
     &         1X,'      ---------------------------------',/,
     &         1X,'      Numero du Transfo        UROFF    ',/,
     &         1X,'      ---------------------------------')
 1501 FORMAT(1X,'      ---------------------------------')
 2000 FORMAT(10X,I6,12X,E12.5)
 3000 FORMAT(I7)
C
C----
C FIN
C----
C
      END
c@z
