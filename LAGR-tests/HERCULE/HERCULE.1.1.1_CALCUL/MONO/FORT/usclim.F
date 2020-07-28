c@a
                        SUBROUTINE USCLIM
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , 
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR , 
     &   NVAR   , NSCAL  , NPHAS  , 
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE , 
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ITRIFB , ITYPFB ,
     &   ICOUPL , IDEVEL , ITUSER , IA     , 
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL ,
     &   W1     , W2     , W3     , W4     , W5     , W6     , COEFU  , 
     &   RCOUPL , RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
c@versb
CVERS    CODE_SATURNE VERSION 1.1             COPYRIGHT EDF 1999, 2000 *
c@verse
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
CFONC   Type de bord  Entree  Sortie "9"  Sortie "10"  Symetrie   Paroi
CFONC   Code          IENTRE  ISOR09      ISOR10       ISYMET     IPAROI
CFONC
CFONC      Les entiers IENTRE, ISOR09, ISOR10, ISYMET, IPAROI,
CFONC      sont affectes par ailleurs (include param.h). Leur valeur
CFONC      est superieure ou egale a 1 et
CFONC          inferieure ou egale a NTYPMX (valeur fixee dans paramx.h)
CFONC
CFONC
CFONC
CFONC      En outre, il faut donner certaines valeurs : 
CFONC
CFONC
CFONC      - Entree : 
CFONC
CFONC        -> Conditions de Dirichlet sur les variables 
CFONC          autres que la pression ; ainsi
CFONC          en face IFAC, sur la variable IVAR : RCODCL(IFAC,IVAR,1)
CFONC
CFONC
CFONC      - Paroi  : (= solide impermeable, avec frottement)
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
CFONC
CFONC      - Symetrie (= paroi impermeable avec glissement) :
CFONC
CFONC        -> Rien a preciser
CFONC
CFONC
CFONC      - Sortie "9" 
CFONC               (= si flux de masse entrant : 
CFONC                     on impose un flux nul sur les scalaires
CFONC                                                  (hormis pression)
CFONC                     on retient la vitesse a l'infini         
CFONC                  si flux de masse sortant : 
CFONC                     on impose un flux nul sur les scalaires
CFONC                                                  (hormis pression)
CFONC                     on impose un flux nul sur la vitesse
CFONC                  dans les deux cas : Dirichlet sur la pression
CFONC
CFONC        -> Rien a preciser
CFONC
CFONC
CFONC      - Sortie "10" 
CFONC               (= si flux de masse entrant : 
CFONC                     on impose un Dirichlet sur les scalaires 
CFONC                                                  (hormis pression)
CFONC                     on retient la vitesse a l'infini         
CFONC                  si flux de masse sortant : 
CFONC                     on impose un flux nul sur les scalaires
CFONC                                                  (hormis pression)
CFONC                     on impose un flux nul sur la vitesse
CFONC                  dans les deux cas : Dirichlet sur la pression
CFONC
CFONC
CFONC        -> Conditions de Dirichlet sur les scalaires
CFONC                                                  (hormis pression)
CFONC          en face IFAC, sur la variable IVAR : RCODCL(IFAC,IVAR,1)
CFONC
CFONC          Dans ce cas pour le Rij-epsilon, donner, a
CFONC              defaut de mieux une valeur de Dirichlet 
CFONC              NULLE POUR LES TENSIONS EXTRADIAGONALES.
CFONC
CFONC
CFONC        Noter que la pression sera recalee a P0 
CFONC            sur la premiere face de sortie de type 9 (ISOR09) ou, 
CFONC            s'il n'y en a pas, 
CFONC            sur la premiere face de sortie de type 10 (ISOR10)
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
CFONC            IENTRE, ISOR09, ISOR10, ISYMET, IPAROI, sont dans cet 
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
CFONC        1 : Dirichlet
CFONC        3 : Neumann
CFONC        4 : Symetrie
CFONC        5 : Paroi
CFONC        9 : Sortie 9
CFONC       10 : Sortie 10 (interdit pour la vitesse)
CFONC      Les valeurs des 3 reels RCODCL sont les suivantes
CFONC       RCODCL(IFAC,IVAR,1) : 
CFONC          Dirichlet sur la variable        si ICODCL(IFAC,IVAR)=  1
CFONC          valeur en paroi (defilmnt, temp) si ICODCL(IFAC,IVAR)=  5
CFONC          valeur en sortie "entrante"      si ICODCL(IFAC,IVAR)= 10
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
CFONC       RCODCL(IFAC,IVAR,3) : 
CFONC         Densite de flux (< 0 si gain) 
CFONC                          si ICODCL(IFAC,IVAR)= 3, 9 ou 10
CFONC          Pour les vitesses U,             en kg/(m s2) = J :
CFONC            RCODCL(IFAC,IVAR,3) =            (VISCL+VISCT) * GRAD U
CFONC          Pour la  pression P,             en kg/(m2 s)     :
CFONC            RCODCL(IFAC,IVAR,3) =                       DT * GRAD P
CFONC          Pour les temperatures T,         en Watt/m2       :
CFONC            RCODCL(IFAC,IVAR,3) = CP*(VISCLS+VISCT/SIGMAS) * GRAD T
CFONC          Pour les enthalpies H,           en Watt/m2       :
CFONC            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD H
CFONC          Pour les autres scalaires F      en :
CFONC            RCODCL(IFAC,IVAR,3) =    (VISCLS+VISCT/SIGMAS) * GRAD F
CFONC
CFONC
CFONC
CFONC       Noter bien que si l'utilisateur affecte une valeur a ITYPFB
CFONC        parmi     IENTRE, ISOR09, ISOR10, ISYMET, IPAROI, 
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
CFONC        non egale a IENTRE, ISOR09, ISOR10, ISYMET, IPAROI, IINDEF
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
CFONC            Si code vitesse = 4          
CFONC              il faut code (k et epsilon)   = 3
CFONC            Si code (vitesse ou Rij ou k ou epsilon) = 5
CFONC              il faut code (Vitesse et Rij et k et epsilon) = 5
CFONC            Si code scalaire (hormis pression ou fluctuations) = 5 
CFONC              il faut code Vitesse = 5
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
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
CARGU !              !    !     ! = 10  -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE NON BLOQUEE :   !
CARGU !              !    !     !  PRESCRIRE UNE VALEUR DE DIRICHLET EN!
CARGU !              !    !     !  PREVISION POUR LES SCALAIRES K, EPS,!
CARGU !              !    !     !  SCAL EN PLUS DU NEUMANN USUEL       !
CARGU ! ITRIFB(NFABOR! TE !  -> ! INDIRECTION POUR TRI DES FACES DE BRD!
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! ITYPFB(NFABOR! TE ! <-  ! TYPE DES FACES DE BORD               !
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! ICOUPL(*)    ! TE !  -> ! TAB ENTIER COMPLEMENTAIRE COUPLAGES  !
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
CARGU !              !    !     !  FLUX (NEGATIF SI GAIN) W/m2         !
CARGU !              !    !     ! POUR LES VITESSES (VISTL+VISCT)*GRADU!
CARGU !              !    !     ! POUR LA PRESSION             DT*GRADP!
CARGU !              !    !     ! POUR LES SCALAIRES                   !
CARGU !              !    !     !        CP*(VISCLS+VISCT/SIGMAS)*GRADT!
CARGU ! W1,2,3,4,5,6 ! TR !  -  ! TABLEAUX DE TRAVAIL                  !
CARGU !  (NCELET     !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! COEFU        ! TR !  -  ! TAB DE TRAV                          !
CARGU !  (NFABOR,3)  !    !     !  (CALCUL DU GRADIENT DE PRESSION)    !
CARGU ! RCOUPL(*)    ! TR !  -> ! TABLEAU REEL COMPOSITE COUPLAGES     !
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
      INTEGER          IPNFAC(NFAC+1), NODFAC(LNDFAC)
      INTEGER          IPNFBR(NFABOR+1), NODFBR(LNDFBR)
      INTEGER          ICODCL(NFABOR,NVAR)
      INTEGER          ITRIFB(NFABOR,NPHAS), ITYPFB(NFABOR,NPHAS)
      INTEGER          ICOUPL(*)
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
      DOUBLE PRECISION RCOUPL(*)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES 
C
      INTEGER          IDEBIA, IDEBRA 
      INTEGER          IFAC, IEL, II, ICOUL, IVAR, IPHAS
      INTEGER          IFML
      DOUBLE PRECISION UREF2
      DOUBLE PRECISION RHOMOY, DH, REYNOL, XLAMBD, USTAR2
      DOUBLE PRECISION XINTUR
      DOUBLE PRECISION XKENT, XEENT

      DOUBLE PRECISION UR11,UA,UC,AK0,EP0,RAYON,TETA
      DOUBLE PRECISION R11,R12,R13,R22,R23,R33
C
C***********************************************************************
C
C
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
C     BOUCLE SUR LES FACES DE BORD
C
      DO IFAC = 1, NFABOR
C
C      FAMILLE ET PREMIERE PROPRIETE
C
C        En general, la couleur des faces est la premiere propriete :
C          ICOUL = IPRFML(IFML,1).
C
C        Pour les fichiers slc et tlc, dans la rubrique propriete, 
C          le premier entier est le numero de famille, les suivants 
C          sont les couleurs. L'enveloppe cherche a lire la seconde 
C          couleur si elle existe et la premiere sinon. S'il y a plus 
C          de 2 couleurs, les autres sont ignorees. 
C          Dans le cas de maillages repris de N3S, la deuxieme colonne 
C          est celle qui contient la couleur a utiliser (issue des 
C          codes de conditions aux limites portant sur la composante 
C          x de la vitesse). 
C          Dans le cas des maillages issus de ESTET, les colonnes 
C          sont en general identiques. 
C          Pour les fichiers slc et tlc, on determine donc la couleur 
C          par ICOUL = IPRFML(IFML,1) si on utilise l'enveloppe
C          et par ICOUL = IPRFML(IFML,2) si on lit directement un
C          fichier tlc provenant de N3S
C          et par ICOUL = IPRFML(IFML,1) si on lit directement un 
C          fichier slc provenant de ESTET
C
C
C        Dans le cas plus specifique de maillages contenant des groupes 
C          (ensembles dont l'intersection peut ne pas etre vide), 
C          il est possible que l'on obtienne plusieurs proprietes par 
C          famille. Ainsi, une face appartenant aux groupes 3 et 5 et 
C          portant la couleur 8 sera dans une famille dont les 
C          proprietes seront -3, -5 et 8 (ordre arbitraire). 
C          Il faut donc, si l'on veut savoir si la face IFAC porte 
C          la couleur 8, balayer toutes les proprietes de la 
C          famille : 
C            IFML  = IFMFBR(IFAC   )
C            DO IPROP = 1, NPRFML
C              IF(IPRFML(IFML,IPROP).EQ.8) THEN
C                La face porte la couleur 8 ...
C              ENDIF
C            ENDDO
C
C
C          
        IFML  = IFMFBR(IFAC   )
        ICOUL = IPRFML(IFML,1)
C
        IEL = IFABOR(IFAC)
C
C      ELEMENT ADJACENT A LA FACE DE BORD
C
C
C --- On impose en couleur 2 pour X <= 0.01D0 une entree pour toutes les phases

       IF( ICOUL.EQ.11 .OR. ICOUL.EQ.21 
     & .OR. ICOUL.EQ.7  .OR. ICOUL.EQ.9 ) THEN
C
          DO IPHAS = 1, NPHAS
C
            ITYPFB(IFAC,IPHAS) = IENTRE
C
            RAYON = SQRT( CDGFBO(1,IFAC)*CDGFBO(1,IFAC)
     &                   +CDGFBO(2,IFAC)*CDGFBO(2,IFAC) )
C
            CALL INLET(UR11,UA,UC,AK0,EP0,RAYON,R11,R12,R13,R22,R23,R33)
C
            IF ( CDGFBO(1,IFAC) .NE. 0.D0 .AND.
     &           CDGFBO(2,IFAC) .NE. 0.D0       ) THEN
              TETA = atan(CDGFBO(2,IFAC)/CDGFBO(1,IFAC))
            ELSE IF ( CDGFBO(1,IFAC) .EQ. 0.D0 ) THEN
              TETA = 0.
            ELSE
              TETA = 3.14159 /2.
            ENDIF
C
            RCODCL(IFAC,IU(IPHAS),1) = UR11*COS(TETA)-UC*SIN(TETA)
            RCODCL(IFAC,IV(IPHAS),1) = UR11*SIN(TETA)+UC*COS(TETA)
            RCODCL(IFAC,IW(IPHAS),1) = UA
C
            RCODCL(IFAC,IR11(IPHAS),1)  = 2.D0/3.D0 * AK0
            RCODCL(IFAC,IR12(IPHAS),1)  = 0.D0
            RCODCL(IFAC,IR13(IPHAS),1)  = 0.D0
            RCODCL(IFAC,IR22(IPHAS),1)  = 2.D0/3.D0 * AK0
            RCODCL(IFAC,IR23(IPHAS),1)  = 0.D0
            RCODCL(IFAC,IR33(IPHAS),1)  = 2.D0/3.D0 * AK0
            RCODCL(IFAC,IEP(IPHAS),1)   = EP0
 
C
          ENDDO
C
C --- On impose en couleur 19 une sortie
        ELSE IF( ICOUL.EQ.19 .OR. ICOUL.EQ.10 ) THEN
C
C          SORTIE : FLUX NUL VITESSE ET TEMPERATURE, PRESSION IMPOSEE
C            Noter que la pression sera recalee a P0 
C                sur la premiere face de sortie de type 9 (ISOR09) ou, 
C                s'il n'y en a pas, 
C                sur la premiere face de sortie de type 10 (ISOR10)
C
          DO IPHAS = 1, NPHAS
            ITYPFB(IFAC,IPHAS)   = ISOR09
          ENDDO
C
C --- On impose en couleur 15 une paroi pour toutes les phases
        ELSE IF( ICOUL.EQ.5 .OR. ICOUL.EQ.8 .OR. ICOUL.EQ.15 ) THEN
C
C          PAROI : DEBIT NUL (FLUX NUL POUR LA PRESSION)
C                  FROTTEMENT POUR LES VITESSES (+GRANDEURS TURB)
C                  FLUX NUL SUR LES SCALAIRES 
C
          DO IPHAS = 1, NPHAS
            ITYPFB(IFAC,IPHAS)   = IPAROI
          ENDDO
C
C --- On impose en couleur 4 une symetrie pour toutes les phases
        ELSEIF( ICOUL.EQ.14 .OR. ICOUL.EQ.4 
     &     .OR. ICOUL.EQ.2  .OR. ICOUL.EQ.3 ) THEN
C
C          SYMETRIES
C
          DO IPHAS = 1, NPHAS
            ITYPFB(IFAC,IPHAS)   = ISYMET
          ENDDO
C
        ENDIF
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
