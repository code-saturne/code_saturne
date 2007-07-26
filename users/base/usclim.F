c@a
c@versb
C-----------------------------------------------------------------------
C
CVERS                  Code_Saturne version 1.3
C                      ------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2007 EDF S.A., France
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
                        SUBROUTINE USCLIM
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ITRIFB , ITYPFB ,
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
CFONC     ROUTINE UTILISATEUR
CFONC     REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
CFONC     (ICODCL,RCODCL) POUR LES VARIABLES INCONNUES
CFONC
CFONC
CFONC
CFONC INTRODUCTION
CFONC ============
CFONC
CFONC On donne ici les conditions aux limites par face de bord.
CFONC
CFONC On balaye les faces de bord et selon un critere on affecte tel
CFONC     ou tel type de conditions aux limites. Dans l'exemple donne
CFONC     ci dessous, c'est la couleur (propriete 1 de la famille)
CFONC     qui permet de distinguer les differents types de bord. On
CFONC     aurait pu aussi travailler avec les coordonnees du centre
CFONC     des faces, mais "c'eut ete moins pratique".
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
CARGU ! MAXELT       !  E !  -> ! NB MAX D'ELEMENTS (CELL,FAC,FBR)     !
CARGU ! LSTELT(MAXELT) TE !  -  ! TABLEAU DE TRAVAIL                   !
CARGU ! NFML  ,NPRFML!    !     !                                      !
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
      INCLUDE "pointe.h"
      INCLUDE "numvar.h"
      INCLUDE "optcal.h"
      INCLUDE "cstphy.h"
      INCLUDE "cstnum.h"
      INCLUDE "entsor.h"
      INCLUDE "parall.h"
      INCLUDE "period.h"
      INCLUDE "ihmpre.h"
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
      INTEGER          IPRFML(NFML,NPRFML), MAXELT, LSTELT(MAXELT)
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ITRIFB(NFABOR,NPHAS), ITYPFB(NFABOR,NPHAS)
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
      INTEGER          IFAC, IEL, II, IVAR, IPHAS
      INTEGER          ILELT, NLELT
      DOUBLE PRECISION UREF2, D2S3
      DOUBLE PRECISION RHOMOY, DH, USTAR2
      DOUBLE PRECISION XINTUR
      DOUBLE PRECISION XKENT, XEENT
C
C***********************************************************************
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_DEBUT
C=======================================================================
C 0.  CE TEST PERMET A L'UTILISATEUR D'ETRE CERTAIN QUE C'EST
C       SA VERSION DU SOUS PROGRAMME QUI EST UTILISEE
C       ET NON CELLE DE LA BIBLIOTHEQUE
C     SI UN FICHIER ISSU DE L'IHM EST UTILISE, LE SOUS-PROGRAMME N'EST
C       PAS NECESSAIREMENT INDISPENSABLE (RETURN DANS LA VERSION DE
C       LA BIBLIOTHEQUE)
C=======================================================================
C
      IF(IIHMPR.EQ.1) THEN
        RETURN
      ELSE
        WRITE(NFECRA,9000)
        CALL CSEXIT (1)
      ENDIF
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,
     &'@    *********                                               ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usclim DOIT ETRE COMPLETE',/,
     &'@                                                            ',/,
     &'@  Le calcul ne sera pas execute.                            ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C
C TEST_A_ENLEVER_POUR_UTILISER_LE_SOUS_PROGRAMME_FIN
C
C=======================================================================
C 1.  INITIALISATIONS
C
C=======================================================================
C
      IDEBIA = IDBIA0
      IDEBRA = IDBRA0
C
      D2S3 = 2.D0/3.D0
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
C
C --- On impose en couleur 2 pour X <= 0.01D0 une entree pour toutes les phases
      CALL GETFBR('2 and X < 0.01',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
        IEL = IFABOR(IFAC)
C
        DO IPHAS = 1, NPHAS
C
          ITYPFB(IFAC,IPHAS) = IENTRE
C
          RCODCL(IFAC,IU(IPHAS),1) = 1.1D0
          RCODCL(IFAC,IV(IPHAS),1) = 1.1D0
          RCODCL(IFAC,IW(IPHAS),1) = 1.1D0
C
          UREF2 = RCODCL(IFAC,IU(IPHAS),1)**2
     &           +RCODCL(IFAC,IV(IPHAS),1)**2
     &           +RCODCL(IFAC,IW(IPHAS),1)**2
          UREF2 = MAX(UREF2,1.D-12)
C
C
C       Exemple de turbulence calculee a partir
C         de formules valables pour une conduite
C
C       On veillera a specifier le diametre hydraulique
C         adapte a l'entree courante.
C
C       On s'attachera egalement a utiliser si besoin une formule
C         plus precise pour la viscosite dynamique utilisee dans le
C         calcul du nombre de Reynolds (en particulier, lorsqu'elle
C         est variable, il peut etre utile de reprendre ici la loi
C         imposee dans USPHYV. On utilise ici par defaut la valeur
C         VISCL0 donnee dans USINI1
C       En ce qui concerne la masse volumique, on dispose directement
C         de sa valeur aux faces de bord (ROMB) et c'est celle que
C         utilise donc ici (elle est en particulier coherente avec
C         le traitement implante dans USPHYV, en cas de masse
C         volumique variable)
C
C         Diametre hydraulique
          DH     = 0.075D0
C
C         Calcul de la vitesse de frottement au carre (USTAR2)
C           et de k et epsilon en entree (XKENT et XEENT) a partir
C           de lois standards en conduite circulaire
C           (leur initialisation est inutile mais plus propre)
          RHOMOY = PROPFB(IFAC,IPPROB(IROM(IPHAS)))
          USTAR2 = 0.D0
          XKENT  = EPZERO
          XEENT  = EPZERO
C
          CALL KEENDB
C         ===========
     &        ( UREF2, DH, RHOMOY, VISCL0(IPHAS), CMU, XKAPPA,
     &          USTAR2, XKENT, XEENT )
C
C     ITYTUR est un indicateur qui vaut ITURB/10
          IF    (ITYTUR(IPHAS).EQ.2) THEN
C
            RCODCL(IFAC,IK(IPHAS),1)  = XKENT
            RCODCL(IFAC,IEP(IPHAS),1) = XEENT
C
          ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
            RCODCL(IFAC,IR11(IPHAS),1) = D2S3*XKENT
            RCODCL(IFAC,IR22(IPHAS),1) = D2S3*XKENT
            RCODCL(IFAC,IR33(IPHAS),1) = D2S3*XKENT
            RCODCL(IFAC,IR12(IPHAS),1) = 0.D0
            RCODCL(IFAC,IR13(IPHAS),1) = 0.D0
            RCODCL(IFAC,IR23(IPHAS),1) = 0.D0
            RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
C
          ELSEIF(ITURB(IPHAS).EQ.50) THEN
C
            RCODCL(IFAC,IK(IPHAS),1)   = XKENT
            RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
            RCODCL(IFAC,IPHI(IPHAS),1) = D2S3
            RCODCL(IFAC,IFB(IPHAS),1)  = 0.D0
C
          ELSEIF(ITURB(IPHAS).EQ.60) THEN
C
            RCODCL(IFAC,IK(IPHAS),1)   = XKENT
            RCODCL(IFAC,IOMG(IPHAS),1) = XEENT/CMU/XKENT
C
          ENDIF
C
C --- On traite les scalaires rattaches a la phase courante
          IF(NSCAL.GT.0) THEN
            DO II = 1, NSCAL
              IF(IPHSCA(II).EQ.IPHAS) THEN
                RCODCL(IFAC,ISCA(II),1) = 1.D0
              ENDIF
            ENDDO
          ENDIF
C
        ENDDO
C
      ENDDO
C
C --- On impose en couleur 3 une entree pour toutes les phases
      CALL GETFBR('3',NLELT,LSTELT)
C     ===========
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
        IEL  = IFABOR(IFAC)
C
        DO IPHAS = 1, NPHAS
C
          ITYPFB(IFAC,IPHAS) = IENTRE
C
          RCODCL(IFAC,IU(IPHAS),1) = 1.1D0
          RCODCL(IFAC,IV(IPHAS),1) = 1.1D0
          RCODCL(IFAC,IW(IPHAS),1) = 1.1D0
C
          UREF2 = RCODCL(IFAC,IU(IPHAS),1)**2
     &           +RCODCL(IFAC,IV(IPHAS),1)**2
     &           +RCODCL(IFAC,IW(IPHAS),1)**2
          UREF2 = MAX(UREF2,1.D-12)
C
C
C       Exemple de turbulence calculee a partir
C         de la donne d'une intensite turbulente
C
C       On veillera a specifier le diametre hydraulique
C         adapte a l'entree courante.
C
C         Diametre hydraulique
          DH     = 0.075D0
C         Intensite turbulente
          XINTUR = 0.02D0
C
C         Calcul de k et epsilon en entree (XKENT et XEENT) a partir
C           l'intensite turbulente et de lois standards en conduite
C           circulaire (leur initialisation est inutile mais plus
C           propre)
          XKENT  = EPZERO
          XEENT  = EPZERO
C
          CALL KEENIN
C         ===========
     &        ( UREF2, XINTUR, DH, CMU, XKAPPA, XKENT, XEENT )
C
C     ITYTUR est un indicateur qui vaut ITURB/10
          IF    (ITYTUR(IPHAS).EQ.2) THEN
C
            RCODCL(IFAC,IK(IPHAS),1)  = XKENT
            RCODCL(IFAC,IEP(IPHAS),1) = XEENT
C
          ELSEIF(ITYTUR(IPHAS).EQ.3) THEN
C
            RCODCL(IFAC,IR11(IPHAS),1) = D2S3*XKENT
            RCODCL(IFAC,IR22(IPHAS),1) = D2S3*XKENT
            RCODCL(IFAC,IR33(IPHAS),1) = D2S3*XKENT
            RCODCL(IFAC,IR12(IPHAS),1) = 0.D0
            RCODCL(IFAC,IR13(IPHAS),1) = 0.D0
            RCODCL(IFAC,IR23(IPHAS),1) = 0.D0
            RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
C
          ELSEIF(ITURB(IPHAS).EQ.50) THEN
C
            RCODCL(IFAC,IK(IPHAS),1)   = XKENT
            RCODCL(IFAC,IEP(IPHAS),1)  = XEENT
            RCODCL(IFAC,IPHI(IPHAS),1) = D2S3
            RCODCL(IFAC,IFB(IPHAS),1)  = 0.D0
C
          ELSEIF(ITURB(IPHAS).EQ.60) THEN
C
            RCODCL(IFAC,IK(IPHAS),1)   = XKENT
            RCODCL(IFAC,IOMG(IPHAS),1) = XEENT/CMU/XKENT
C
          ENDIF
C
C --- On traite les scalaires rattaches a la phase courante
          IF(NSCAL.GT.0) THEN
            DO II = 1, NSCAL
              IF(IPHSCA(II).EQ.IPHAS) THEN
                RCODCL(IFAC,ISCA(II),1) = 1.D0
              ENDIF
            ENDDO
          ENDIF
C
        ENDDO
      ENDDO
C
C --- On impose sur le groupe "sortie" une sortie pour toutes les phases
      CALL GETFBR('sortie',NLELT,LSTELT)
C     ===========
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C          SORTIE : FLUX NUL VITESSE ET TEMPERATURE, PRESSION IMPOSEE
C            Noter que la pression sera recalee a P0
C                sur la premiere face de sortie libre (ISOLIB)
C
        DO IPHAS = 1, NPHAS
          ITYPFB(IFAC,IPHAS)   = ISOLIB
        ENDDO
C
      ENDDO
C
C --- On impose en couleur 5 une paroi pour toutes les phases
      CALL GETFBR('5',NLELT,LSTELT)
C     ===========
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C          PAROI : DEBIT NUL (FLUX NUL POUR LA PRESSION)
C                  FROTTEMENT POUR LES VITESSES (+GRANDEURS TURB)
C                  FLUX NUL SUR LES SCALAIRES
C
        DO IPHAS = 1, NPHAS
          ITYPFB(IFAC,IPHAS)   = IPAROI
        ENDDO
C
C SI DEFILANTE VITESSE U(1) = 1
C         RCODCL(IFAC,IU(1),1) = 1.D0
C SI VITESSE DEFILANTE U = 0, RIEN A FAIRE !
C
        IF(NSCAL.GT.0) THEN
C
C SI TEMPERATURE IMPOSEE A 20 AVEC LOI DE PAROI (SCALAIRE II=1)
C            II = 1
C            ICODCL(IFAC,ISCA(II))   = 5
C            RCODCL(IFAC,ISCA(II),1) = 20.D0
C
C SI TEMPERATURE IMPOSEE A 50 SANS LOI DE PAROI (DIRICHLET SIMPLE)
C               AVEC COEFFICIENT D'ECHANGE DE 8 (SCALAIRE II=2)
C            II = 2
C            ICODCL(IFAC,ISCA(II))   = 1
C            RCODCL(IFAC,ISCA(II),1) = 50.D0
C            RCODCL(IFAC,ISCA(II),2) = 8.D0
C
C SI FLUX IMPOSE A 4.D0 (SCALAIRE II=3)
C            II = 3
C            ICODCL(IFAC,ISCA(II))   = 3
C            RCODCL(IFAC,ISCA(II),3) = 4.D0
C
        ENDIF
      ENDDO
C
C --- On impose en couleur 7 une paroi rugueuse pour toutes les phases
      CALL GETFBR('7',NLELT,LSTELT)
C     ===========
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C          PAROI : DEBIT NUL (FLUX NUL POUR LA PRESSION)
C                  FROTTEMENT RUGUEUX POUR LES VITESSES (+GRANDEURS TURB)
C                  FLUX NUL SUR LES SCALAIRES
C
        DO IPHAS = 1, NPHAS
          ITYPFB(IFAC,IPHAS)   = IPARUG
C
C     Rugosite pour la vitesse : 1cm
          RCODCL(IFAC,IU(IPHAS),3) = 0.01D0
        ENDDO
C
C SI DEFILANTE VITESSE U(1) = 1
C         RCODCL(IFAC,IU(1),1) = 1.D0
C SI VITESSE DEFILANTE U = 0, RIEN A FAIRE !
C
        IF(NSCAL.GT.0) THEN
C
C SI TEMPERATURE IMPOSEE A 20 AVEC LOI DE PAROI RUGUEUSE (SCALAIRE II=1)
C    AVEC RUGOSITE DE 1cm
C            II = 1
C            ICODCL(IFAC,ISCA(II))   = 6
C            RCODCL(IFAC,ISCA(II),1) = 20.D0
C            RCODCL(IFAC,ISCA(II),3) = 0.001D0
C
C SI FLUX IMPOSE A 4.D0 (SCALAIRE II=2)
C            II = 2
C            ICODCL(IFAC,ISCA(II))   = 3
C            RCODCL(IFAC,ISCA(II),3) = 4.D0
C
        ENDIF
C
      ENDDO
C
C --- On impose en couleur 4 une symetrie pour toutes les phases
      CALL GETFBR('4',NLELT,LSTELT)
C     ===========
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C          SYMETRIES
C
        DO IPHAS = 1, NPHAS
          ITYPFB(IFAC,IPHAS)   = ISYMET
        ENDDO
C
      ENDDO
C
C
C    A partir d'ici, les exemples donnes sont destines
C       a une UTILISATION AVANCEE, l'utilisateur etant amene a
C       definir des conbinaisons complexes de conditions aux limites.
C
C
C
C
C
C --- Exemple de conditions aux limites specifiques baties par
C       l'utilisateur, sur la base de conditions de paroi.
C       On impose en couleur 1234 une paroi pour toutes les phases
C       avec, en outre,
C          une condition de Dirichlet sur la vitesse de la phase 1
C          (paroi mobile avec adhérence) et
C          une condition de Dirichlet sur le premier scalaire.
C
      CALL GETFBR('1234',NLELT,LSTELT)
C     ===========
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        DO IPHAS = 1, NPHAS
          ITYPFB(IFAC,IPHAS) = IPAROI
        ENDDO
C
        IPHAS = 1
        ICODCL(IFAC,IU(IPHAS) )  = 1
        RCODCL(IFAC,IU(IPHAS),1) = 1.D0
        RCODCL(IFAC,IU(IPHAS),2) = RINFIN
        RCODCL(IFAC,IU(IPHAS),3) = 0.D0
        ICODCL(IFAC,IV(IPHAS) )  = 1
        RCODCL(IFAC,IV(IPHAS),1) = 0.D0
        RCODCL(IFAC,IV(IPHAS),2) = RINFIN
        RCODCL(IFAC,IV(IPHAS),3) = 0.D0
        ICODCL(IFAC,IW(IPHAS) )  = 1
        RCODCL(IFAC,IW(IPHAS),1) = 0.D0
        RCODCL(IFAC,IW(IPHAS),2) = RINFIN
        RCODCL(IFAC,IW(IPHAS),3) = 0.D0
C
        IVAR = ISCA(1)
        ICODCL(IFAC,IVAR )  = 1
        RCODCL(IFAC,IVAR,1) = 10.D0
        RCODCL(IFAC,IVAR,2) = RINFIN
        RCODCL(IFAC,IVAR,3) = 0.D0
C
      ENDDO
C
C --- Exemple de conditions aux limites specifiques baties de toutes
C       pieces par l'utilisateur, sans definir de type particulier.
C       On impose en couleur 5678 un Neumann homogene pour toutes
C       les variables (quelle que soit la phase).
C
      CALL GETFBR('5678',NLELT,LSTELT)
C     ===========
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C          ATTENTION, la valeur de ITYPFB doit etre affectee a IINDEF
C
        DO IPHAS = 1, NPHAS
          ITYPFB(IFAC,IPHAS) = IINDEF
        ENDDO
C
        DO II = 1, NVAR
          ICODCL(IFAC,II )  = 3
          RCODCL(IFAC,II,1) = 0.D0
          RCODCL(IFAC,II,2) = RINFIN
          RCODCL(IFAC,II,3) = 0.D0
        ENDDO
C
      ENDDO
C
C --- Exemple de conditions aux limites specifiques baties de toutes
C       pieces par l'utilisateur, avec definition d'un type particulier
C       par exemple pour un reperage ulterieur (calcul du flux de masse,
C       impresssions particulieres...)
C       On impose en couleur 6789 un Neumann homogene pour toutes
C       les variables (quelle que soit la phase), sauf le premier
C       scalaire, pour lequel on selectionne un Dirichlet homogene.
C
      CALL GETFBR('6789',NLELT,LSTELT)
C     ===========
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C          ATTENTION, la valeur de ITYPFB doit etre differente de
C              IPAROI, IENTRE, ISYMET, ISOLIB, IINDEF
C              superieure ou egale a 1
C              inferieure ou egale a NTYPMX
C          ces entiers sont definis en include
C
        DO IPHAS = 1, NPHAS
          ITYPFB(IFAC,IPHAS) = 89
        ENDDO
C
        DO II = 1, NVAR
          ICODCL(IFAC,II )  = 3
          RCODCL(IFAC,II,1) = 0.D0
          RCODCL(IFAC,II,2) = RINFIN
          RCODCL(IFAC,II,3) = 0.D0
        ENDDO
C
        ICODCL(IFAC,ISCA(1) )  = 1
        RCODCL(IFAC,ISCA(1),1) = 0.D0
        RCODCL(IFAC,ISCA(1),2) = RINFIN
        RCODCL(IFAC,ISCA(1),3) = 0.D0
C
      ENDDO
C
C----
C FORMATS
C----
C
C----
C FIN
C----
C
      RETURN
      END
c@z
