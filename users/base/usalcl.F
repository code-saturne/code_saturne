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
                        SUBROUTINE USALCL
C                       *****************
C     --------------------------------------------------------------
     & ( IDBIA0 , IDBRA0 , ITRALE ,
     &   NDIM   , NCELET , NCEL   , NFAC   , NFABOR , NFML   , NPRFML ,
     &   NNOD   , LNDFAC , LNDFBR , NCELBR ,
     &   NVAR   , NSCAL  , NPHAS  ,
     &   NIDEVE , NRDEVE , NITUSE , NRTUSE ,
     &   IFACEL , IFABOR , IFMFBR , IFMCEL , IPRFML , MAXELT , LSTELT ,
     &   IPNFAC , NODFAC , IPNFBR , NODFBR ,
     &   ICODCL , ITYPFB , IALTYB , IMPALE ,
     &   IDEVEL , ITUSER , IA     ,
     &   XYZCEN , SURFAC , SURFBO , CDGFAC , CDGFBO , XYZNOD , VOLUME ,
     &   DT     , RTP    , RTPA   , PROPCE , PROPFA , PROPFB ,
     &   COEFA  , COEFB  , RCODCL , XYZNO0 , DEPALE ,
     &   RDEVEL , RTUSER , RA     )
C     --------------------------------------------------------------
C***********************************************************************
C FONCTION :
C --------
c@foncb
CFONC
CFONC     ROUTINE UTILISATEUR
CFONC     REMPLISSAGE DU TABLEAU DE CONDITIONS AUX LIMITES
CFONC       (IALTYB,ICODCL,RCODCL) POUR LA VITESSE DE MAILLAGE
CFONC
CFONC     POSSIBILITE DE FIXER DIRECTEMENT LE DEPLACEMENT DES
CFONC       NOEUDS
CFONC
CFONC
CFONC
CFONC     CE SOUS PROGRAMME UTILISATEUR EST OBLIGATOIRE EN ALE
CFONC     ====================================================
CFONC
CFONC
CFONC
CFONC INTRODUCTION
CFONC ============
CFONC
CFONC On donne ici les conditions aux limites par face de bord, pour
CFONC     pour la vitesse de maillage uniquement (pour les autres
CFONC     variables, se reporter a l'Interface Graphique ou a usclim)
CFONC
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
CFONC     Pour les conditions "standards" :
CFONC     --------------------------------
CFONC
CFONC Trois types de conditions "standards" sont disponibles. On
CFONC      specifie pour chaque face le type de condition choisi.
CFONC      Le choix se fait en remplissant le tableau IALTYB.
CFONC
CFONC
CFONC * IALTYB(IFAC) = IBFIXE : la face IFAC correspond a un bord
CFONC      de maillage fixe. Une condition de Dirichlet nulle sera
CFONC      affectee automatiquement a la vitesse de maillage. En
CFONC      outre, le deplacement des noeuds sera automatiquement
CFONC      impose a 0 (cf. plus base sur IMPALE), sauf si
CFONC      l'utilisateur a modifie la condition pour au moins une
CFONC      des conposantes de la vitesse de maillage (ICODCL modifie,
CFONC      cf. conditions non standards)
CFONC
CFONC
CFONC * IALTYB(IFAC) = IGLISS : la face IFAC est une face a maillage
CFONC      glissant. La composante normale a la face de la vitesse
CFONC      de maillage sera forcee a 0, les autres composantes seront
CFONC      traitees en Neumann homogene (pour les faces non alignees
CFONC      avec les axes, la condition de Neumann homogene est
CFONC      seulement partiellement implicitee), comme pour les
CFONC      conditions de symetrie de la vitesse fluide.
CFONC
CFONC
CFONC * IALTYB(IFAC) = IVIMPO : la face IFAC a une vitesse de maillage
CFONC      imposee. Cette vitesse est specifiee en remplissant le
CFONC      tableau RCODCL correspondant :
CFONC      RCODCL(IFAC,IUMA,1) = vitesse de maillage selon X
CFONC      RCODCL(IFAC,IVMA,1) = vitesse de maillage selon Y
CFONC      RCODCL(IFAC,IWMA,1) = vitesse de maillage selon Z
CFONC      Les composantes non specifiees de RCODCL(.,I.MA,1) seront mises
CFONC      a 0 (seules les composantes non nulles sont donc a specifier)
CFONC
CFONC
CFONC
CFONC     Pour les conditions "non standards" :
CFONC     ------------------------------------
CFONC
CFONC      Autres que fixe, glissement ou vitesse imposee, on donne
CFONC
CFONC       - pour chaque face et chaque composante IVAR=IUMA, IVMA, IWMA :
CFONC         -> un code     ICODCL(IFAC,IVAR)
CFONC         -> trois reels RCODCL(IFAC,IVAR,1)
CFONC                        RCODCL(IFAC,IVAR,2)
CFONC                        RCODCL(IFAC,IVAR,3)
CFONC      La valeur de ICODCL est prise parmi les suivantes :
CFONC        1 : Dirichlet
CFONC        3 : Neumann
CFONC        4 : Symetrie
CFONC      Les valeurs des 3 reels RCODCL sont les suivantes
CFONC
CFONC       RCODCL(IFAC,IVAR,1) :
CFONC          Dirichlet sur la variable si ICODCL(IFAC,IVAR)=  1
CFONC          La dimension de RCODCL(IFAC,IVAR,1) est en m/s
CFONC
CFONC       RCODCL(IFAC,IVAR,2) :
CFONC          coefficient d'echange "exterieur" en kg/(m2 s) (entre la
CFONC                     valeur imposee et la valeur au bord du domaine)
CFONC                     RINFIN = infini par defaut
CFONC            RCODCL(IFAC,IVAR,2) =            (VISCMA) / D
CFONC               (D a la dimension d'une distance en m, VISCMA est
CFONC                la viscosite de maillage)
CFONC          NB : la definition de RCODCL(.,.,2) est calquee sur le cas
CFONC               des autres variables standards de l'ecoulement. Il n'a pas
CFONC               grand sens physique dans le cas de la vitesse de maillage.
CFONC
CFONC       RCODCL(IFAC,IVAR,3) :
CFONC         Densite de flux en kg/(m s2) = J si ICODCL(IFAC,IVAR)= 3
CFONC                        (< 0 si gain, n normale orientee vers l'exterieur)
CFONC            RCODCL(IFAC,IVAR,3) =           -(VISCMA) * (GRAD Um).n
CFONC             (ou Um represente la vitesse de maillage)
CFONC          NB : la definition de RCODCL(.,.,3) est calquee sur le cas
CFONC               des autres variables standards de l'ecoulement.
CFONC               RCODCL(.,.,3) = 0.D0 permet de specifier un Neumann
CFONC               homogene sur la vitesse de maillage. Toute autre valeur
CFONC               aura une signification physique plus discutable.
CFONC
CFONC
CFONC
CFONC       Noter bien que si l'utilisateur affecte une valeur a IALTYB
CFONC        parmi IBFIXE, IGLISS, IVIMPO,
CFONC        et qu'il ne modifie pas ICODCL (valeur nulle par defaut),
CFONC        c'est ITYPFB qui imposera la condition limite.
CFONC
CFONC       Par contre, si l'utilisateur impose
CFONC         ICODCL(IFAC,IVAR) (non nul),
CFONC         ce sont alors les valeurs de ICODCL et RCODCL qu'il aura
CFONC         fournies qui sont retenues pour la face et la composante
CFONC         consideree (s'il ne precise pas RCODCL, ce sont les valeurs
CFONC         par defaut qui sont retenues pour la face et
CFONC         la variable consideree soit :
CFONC                                  RCODCL(IFAC,IVAR,1) = 0.D0
CFONC                                  RCODCL(IFAC,IVAR,2) = RINFIN
CFONC                                  RCODCL(IFAC,IVAR,3) = 0.D0)
CFONC
CFONC
CFONC       Si l'utilisateur specifie de lui-meme la valeur de ICODCL
CFONC          et de RCODCL pour TOUTES les composantes de la vitesse
CFONC          de maillage, il n'est pas necessaire de specifier IALTYB
CFONC          pour la face concernee (sa valeur ne sera pas utilisee)
CFONC
CFONC
CFONC REGLE DE COHERENCE
CFONC ==================
CFONC
CFONC        Si une condition de symetrie (ICODCL=4) a ete imposee sur
CFONC           une des composantes, elle doit l'etre sur toutes les
CFONC           composantes.
CFONC
CFONC
CFONC DEPLACEMENT FORCE DES NOEUDS
CFONC ============================
CFONC
CFONC Pour plus de precision dans le deplacement du maillage, on peut
CFONC    aussi forcer directement le deplacement de certains noeuds,
CFONC    internes ou de bord. Pour cela on remplit les tableaux DEPALE
CFONC    et IMPALE :
CFONC    DEPALE(INOD,1) = deplacement du noeud INOD dans la direction X
CFONC    DEPALE(INOD,2) = deplacement du noeud INOD dans la direction Y
CFONC    DEPALE(INOD,3) = deplacement du noeud INOD dans la direction Z
CFONC    Ce deplacement s'entend comme le deplacement absolu du noeud
CFONC    a partir de sa position dans le maillage d'origine.
CFONC    IMPALE(INOD) = 1 indique que le noeud INOD est a deplacement
CFONC    impose (IMPALE est initialise a 0 ; si sa valeur est laissee
CFONC    nulle, DEPALE ne sera pas pris en compte).
CFONC
CFONC Lors de la mise a jour du maillage, les noeuds tels que IMPALE=1
CFONC    (internes ou de bord) ne seront pas deplaces a partir de la
CFONC    vitesse de maillage mais directement a partir de DEPALE.
CFONC Dans le cas ou tous les noeuds d'une face sont a deplacement
CFONC    impose, il n'est pas necessaire de remplir les tableaux de
CFONC    conditions aux limites de vitesse de maillage pour cette faces,
CFONC    ils seront ecrases :
CFONC     -> ICODCL sera mis a 1 (Dirichlet)
CFONC     -> RCODCL sera mis a la valeur moyenne des vitesses des noeuds
CFONC        de la face (vitesses calculees a partir de DEPALE)
CFONC
CFONC Dans le cas de faces specifies a maillage fixes (IALTYB(IFAC)=IBFIXE),
CFONC    tous les noeuds de la face sont automatiquement mis en
CFONC    deplacement impose, avec une valeur nulle de DEPALE.
CFONC
CFONC
CFONC  CARACTERISTIQUES DES NOEUDS
CFONC  ===========================
CFONC Le nombre total de noeuds est stocke dans la variable NNOD.
CFONC Les coordonnees des noeuds sont accessibles par le tableau
CFONC    XYZNOD(3,NNOD).
CFONC Les coordonnees des noeuds dans le maillage initiale sont
CFONC    accessibles par le tableau XYZNO0(3,NNOD).
CFONC
CFONC Le reperage des noeuds est possible a partir des faces internes et
CFONC    de bord, grace aux tableaux IPNFAC, NODFAC, IPNFBR, NODFBR.
CFONC    NODFAC (NODFBR) contient sequentiellement la liste des noeuds de
CFONC      toutes les faces internes (de bord).
CFONC    IPNFAC (IPNFBR) contient pour chaque face interne (de bord) le
CFONC      numero de la premiere case de NODFAC (NODFBR) lui correspondant.
CFONC
CFONC Par exemple, pour recuperer sequentiellement tous les noeuds de la
CFONC    face interne IFAC, il suffit d'ecrire la boucle suivante :
CFONC    DO II = IPNFAC(IFAC), IPNFAC(IFAC+1)-1 <- indices des elements
CFONC                                              de NODFAC correspondant
CFONC                                              a IFAC
CFONC      INOD = NODFAC(II)                    <- recuperation du numero
CFONC                                              du IIeme noeud de la face
CFONC                                              IFAC
CFONC      ...
CFONC    ENDDO
CFONC
CFONC
CFONC  INFLUENCE SUR LES CONDITIONS AUX LIMITES DE VITESSE FLUIDE EN PAROI
CFONC  ===================================================================
CFONC Dans le cas de faces de paroi pour le fluide (ITYPFB=IPAROI ou IPARUG),
CFONC   l'influence de la vitesse de maillage depend de sa signification
CFONC   physique.
CFONC   En effet, dans le cas d'une structure immergee par exemple, le
CFONC   mouvement des faces de paroi de la structure correspond a un mouvement
CFONC   physique et doit donc entrainer le fluide.
CFONC   Dans le cas d'un piston au contraire, les bords lateraux sont des
CFONC   parois ou le mouvement des noeuds n'a aucune signification physique et
CFONC   ne doit pas entrainer le fluide.
CFONC   Dans tous les cas, la vitesse de maillage normale a la face est prise
CFONC   en compte (u.n=w.n a la face), c'est le traitement de la vitesse
CFONC   tangentielle qui differe suivant les cas.
CFONC
CFONC Par defaut, Code_Saturne gere la relation entre vitesse fluide et vitesse
CFONC   de maillage pour les parois de la maniere suivante :
CFONC   - Si IALTYB(IFAC) = IBFIXE, la vitesse de maillage est nulle, il n'y a
CFONC     pas de probleme (et si la paroi est defilante, l'utilisateur le
CFONC     specifiera dans l'interface ou dans usclim).
CFONC   - Si IALTYB(IFAC) = IGLISS, la vitesse tangentielle de maillage n'est
CFONC     pas prise en compte dans les conditions aux limites de paroi pour le
CFONC     fluide (et si la paroi est defilante, l'utilisateur le specifiera
CFONC     dans l'interface ou dans usclim).
CFONC   - Si IALTYB(IFAC) = IVIMPO, la vitesse tangentielle de maillage est
CFONC     prise en compte comme une vitesse de defilement dans les conditions
CFONC     aux limites de paroi pour le fluide, sauf si une vitesse de
CFONC     defilement de paroi a ete specifiee par l'utilisateur dans l'interface
CFONC     ou dans usclim (auquel cas c'est cette vitesse qui est consideree).
CFONC   - Si IMPALE(INOD) = 1 pour tous les noeuds d'une face, la vitesse
CFONC     tangentielle de maillage deduite de ce deplacement sera prise en
CFONC     compte comme une vitesse de defilement dans les conditions
CFONC     aux limites de paroi pour le fluide, sauf si une vitesse de
CFONC     defilement de paroi a ete specifiee par l'utilisateur dans l'interface
CFONC     ou dans usclim (auquel cas c'est cette vitesse qui est consideree).
CFONC
CFONC Pour les autres types de conditions aux limites pour le fluide (IENTRE et
CFONC   ISOLIB), la vitesse de maillage n'a pas d'influence.
CFONC
CFONC Dans le cas de conditions aux limites non standards, c'est a l'utilisateur
CFONC   de gerer directement la relation entre les conditions sur la vitesse de
CFONC   maillage et celles sur la vitesse fluide (les conditions aux limites du
CFONC   fluide pouvant etre modifiees dans cette routine).
CFONC
CFONC
CFONC
CFONC IDENTIFICATION DES FACES DE BORD
CFONC ================================
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
c@fonce
C-----------------------------------------------------------------------
C                             ARGUMENTS
c@argub
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! IDBIA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS IA !
CARGU ! IDBRA0       ! E  !  -> ! NUMERO DE LA 1ERE CASE LIBRE DANS RA !
CARGU ! ITRALE       ! E  !  -> ! NUMERO D'ITERATION POUR L'ALE        !
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
CARGU !              !    !     ! = 9   -> ENTREE/SORTIE LIBRE (VITESSE!
CARGU !              !    !     !  ENTRANTE EVENTUELLE     BLOQUEE     !
CARGU ! ITYPFB(NFABOR! TE !  -> ! TYPE DES FACES DE BORD POUR LE FLUIDE!
CARGU !  NPHAS      )!    !     !                                      !
CARGU ! IALTYB(NFABOR! TE ! <-  ! TYPE DES FACES DE BORD POUR LA       !
CARGU !              ! TE !     !                  VITESSE DE MAILLAGE !
CARGU ! IMPALE(NNOD) ! TE !  -> ! INDICATEUR DE DELACEMENT IMPOSE      !
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
CARGU ! DEPALE(NNOD,3! TR !  -> ! DEPLACEMENT AUX NOEUDS               !
CARGU ! XYZNO0(3,NNOD! TR !  -> ! COORDONNEES NOEUDS MAILLAGE INITIAL  !
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
      INTEGER          IDBIA0 , IDBRA0 , ITRALE
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
      INTEGER          ITYPFB(NFABOR,NPHAS), IALTYB(NFABOR)
      INTEGER          IMPALE(NNOD)
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
      DOUBLE PRECISION DEPALE(NNOD,3), XYZNO0(3,NNOD)
      DOUBLE PRECISION RDEVEL(NRDEVE), RTUSER(NRTUSE), RA(*)
C
C VARIABLES LOCALES
C
      INTEGER          IDEBIA, IDEBRA
      INTEGER          IFAC, IEL, II
      INTEGER          INOD
      INTEGER          ILELT, NLELT
C
      DOUBLE PRECISION DELTA, DELTAA
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
      WRITE(NFECRA,9000)
      CALL CSEXIT (1)
C
 9000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET LORS DE L''ENTREE DES COND. LIM.      ',/,
     &'@    *********                                               ',/,
     &'@     LA METHODE ALE A ETE ENCLENCHEE                        ',/,
     &'@     LE SOUS-PROGRAMME UTILISATEUR usalcl DOIT ETRE COMPLETE',/,
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
C     Calcul du deplacement au temps courant
      DELTAA = SIN(3.141596D0*(NTCABS-1)/50.D0)
      DELTA  = SIN(3.141596D0*NTCABS/50.D0)
C
C --- On impose en couleur 4 une vitesse de deplacement des faces
C
      CALL GETFBR('4',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
C     ELEMENT ADJACENT A LA FACE DE BORD
        IEL = IFABOR(IFAC)
C
        IALTYB(IFAC) = IVIMPO
        RCODCL(IFAC,IUMA,1) = 0.D0
        RCODCL(IFAC,IVMA,1) = 0.D0
        RCODCL(IFAC,IWMA,1) = (DELTA-DELTAA)/DT(IEL)
C
      ENDDO
C
C --- On impose en couleur 5 un deplacement des noeuds
C
      CALL GETFBR('5',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        DO II = IPNFBR(IFAC), IPNFBR(IFAC+1)-1
          INOD = NODFBR(II)
          IF (IMPALE(INOD).EQ.0) THEN
            DEPALE(INOD,1) = 0.D0
            DEPALE(INOD,2) = 0.D0
            DEPALE(INOD,3) = DELTA
            IMPALE(INOD) = 1
          ENDIF
        ENDDO
C
      ENDDO
C
C --- On impose en couleur 6 un glissement
C
      CALL GETFBR('6',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IALTYB(IFAC) = IGLISS
C
      ENDDO
C
C --- On impose ailleurs une condition de bord fixe
C
      CALL GETFBR( 'not (4 or 5 or 6)',NLELT,LSTELT)
C     ===========
C
      DO ILELT = 1, NLELT
C
        IFAC = LSTELT(ILELT)
C
        IALTYB(IFAC) = IBFIXE
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
