c@a
c@versb
C-----------------------------------------------------------------------
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
C                              pointe.h
C***********************************************************************
C
C... GEOMETRIE
C      ACCESSIBLES DIRECTEMENT DANS IA, RA
C
C Pointeur Dimension       Description
C IICELB ! NCELBR        ! REPERAGE DES CELLULES DE BORD
C ISRFAN ! NFAC          ! SURFACE DES FACES INTERNES
C ISRFBN ! NFABOR        ! SURFACE DES FACES DE BORD
C IDIST  ! NFAC          ! IJ.Nij
C IDISTB ! NFABOR        ! EQUIVALENT DE IDIST AU BORD
C IPOND  ! NFAC          ! PONDERATION (Aij=POND Ai+(1-POND)Aj)
C IDIJPF ! NFAC*NDIM     ! VECTEUR I'J'
C IDIIPB ! NFAC*NDIM     ! EQUIVALENT DE IDIIPF AU BORD
C IDOFIJ ! NFAC*NDIM     ! VECTEUR OF A LA FACE IJ
C
      INTEGER           IICELB ,
     &                  ISRFAN , ISRFBN , IDIST  , IDISTB , IPOND ,
     &                  IDIJPF , IDIIPB , IDOFIJ
      COMMON / IPGEOM / IICELB ,
     &                  ISRFAN , ISRFBN , IDIST  , IDISTB , IPOND ,
     &                  IDIJPF , IDIIPB , IDOFIJ
C
C
C... AUXILIAIRES INDEPENDANTS DU NOMBRE DE PHASES
C      ACCESSIBLES DIRECTEMENT DANS IA, RA
C
C Pointeur Dimension       Description
C ICOCG  ! NCELET*9      ! STOCKAGE POUR GRADIENT
C ICOCGB ! NCELBR*9      ! STOCKAGE POUR GRADIENT BORD
C ICOCI  ! NCELET*9      ! STOCKAGE POUR GRADIENT SI INIT. PAR MC
C ICOCIB ! NCELBR*9      ! STOCKAGE POUR GRADIENT BORD SI INIT. PAR MC
C ITPUCO !    -          ! VARIABLES DU COUPLAGE U-P
C IDIPAR ! NCELET        ! DISTANCE A LA FACE DE TYPE 5 (PHASE 1) LA
C                            PLUS PROCHE
C IYPPAR ! NCELET        ! YPLUS ASSOCIE (LES only)
C IFORBR ! NFABOR*3      ! EFFORTS AUX BORDS (SI POSTTRAITE)
C IIDFST ! NFABOR        ! TABLEAU D'INDIRECTION POUR LES STRUCTURES
C                          MOBILES EN ALE
C
C... PARAMETRES DU MODULE THERMIQUE 1D
C NFPT1D !               ! NB DE FACES DE BORD AVEC MODULE THERMIQUE 1D
C INPPT1 ! NFPT1D        ! NOMBRE DE MAILLES DANS LA PAROI
C IEPPT1 ! NFPT1D        ! EPAISSEUR DE LA PAROI
C IRGPT1 ! NFPT1D        ! RAISON DU MAILLAGE
C IIFPT1 ! NFPT1D        ! NUMERO DE LA FACE
C ITPPT1 ! NFPT1D        ! TEMPERATURE DE PAROI
C IICLT1 ! NFPT1D        ! TYPE DE CONDITION LIMITE
C ITEPT1 ! NFPT1D        ! TEMPERATURE EXTERIEURE
C IHEPT1 ! NFPT1D        ! COEFFICIENT D ECHANGE EXTERIEUR
C IFEPT1 ! NFPT1D        ! FLUX THERMIQUE EXTERIEUR
C IXLMT1 ! NFPT1D        ! DIFFUSIVITE THERMIQUE
C IRCPT1 ! NFPT1D        ! RHO*CP
C IDTPT1 ! NFPT1D        ! PAS DE TEMPS
C
      INTEGER           ICOCG  , ICOCGB , ICOCI  , ICOCIB ,
     &                  ITPUCO , IDIPAR , IYPPAR , IFORBR , IIDFST ,
     &                  NFPT1D , NMXT1D , INPPT1 , IIFPT1 , IICLT1 ,
     &                  IEPPT1 , IRGPT1 , ITPPT1 ,
     &                  ITEPT1 , IHEPT1 , IFEPT1 ,
     &                  IXLMT1 , IRCPT1 , IDTPT1

      COMMON / IPAUX0 / ICOCG  , ICOCGB , ICOCI  , ICOCIB ,
     &                  ITPUCO , IDIPAR , IYPPAR , IFORBR , IIDFST ,
     &                  NFPT1D , NMXT1D , INPPT1 , IIFPT1 , IICLT1 ,
     &                  IEPPT1 , IRGPT1 , ITPPT1 ,
     &                  ITEPT1 , IHEPT1 , IFEPT1 ,
     &                  IXLMT1 , IRCPT1 , IDTPT1
C
C... AUXILIAIRES ACCESSIBLES DIRECTEMENT DANS IA, RA
C     TOUS CES TABLEAUX SONT (Dimension,NPHAS)
C
C Pointeur Dimension       Description
C IITYPF ! NFABOR        ! TYPE DES FACES DE BORD
C IITRIF ! NFABOR        ! INDIRECTION POUR TRI FACES DE BORD
C IYPLBR ! NFABOR        ! YPLUS BORD (SI POST-TRAITE)
C IISYMP ! NFABOR        ! ZERO POUR ANNULER LE FLUX DE MASSE
C        !               !   (SYMETRIES ET PAROIS AVEC CL COUPLEES)
C        !               ! UN SINON
C IIFAPA ! NCELET        ! NUMERO DE FACE DE BORD 5 LA PLUS PROCHE
C NCEPDC !               ! NOMBRE DE CELLULES AVEC PDC
C NCKPDC !               ! DIMENSION TENSEUR PDC
C IICEPD ! NCEPDC        ! NUMERO DES CELLULES AVEC PerteDeCharge
C ICKUPD ! NCEPDC,NCKPDC ! VALEUR DES COEFF DE PDC
C NCETSM !               ! NOMBRE DE CELLULES AVEC TSM
C IICESM ! NCETSM        ! NUMERO DES CELLULES AVEC TSMasse
C IITPSM ! NCETSM        ! TYPE DE TSM
C ISMACE ! NCETSM        ! VALEUR DE TSM
C IS2KW  ! NCELET        ! STOCKAGE DE 2 Sij.Sij EN K-OMEGA
C IDVUKW ! NCELET        ! STOCKAGE DE DIVU EN K-OMEGA (EN MEME TEMPS QUE S2KW)
C
      INTEGER           IITYPF        , IITRIF        ,
     &                  IISYMP        , IYPLBR        ,
     &                  IIFAPA(NPHSMX), NCEPDC(NPHSMX), NCKPDC(NPHSMX),
     &                  IICEPD(NPHSMX), ICKUPD(NPHSMX), NCETSM(NPHSMX),
     &                  IICESM(NPHSMX), IITPSM(NPHSMX), ISMACE(NPHSMX),
     &                  IS2KW (NPHSMX), IDVUKW(NPHSMX)
      COMMON / IPOSUP / IITYPF        , IITRIF        ,
     &                  IISYMP        , IYPLBR        ,
     &                  IIFAPA        , NCEPDC        , NCKPDC        ,
     &                  IICEPD        , ICKUPD        , NCETSM        ,
     &                  IICESM        , IITPSM        , ISMACE        ,
     &                  IS2KW         , IDVUKW
C
C
C... AUXILIAIRES ACCESSIBLES DIRECTEMENT DANS IA, RA
C    POUR LA PERODICITE
C
C Pointeur Dimension                 |   Description
C IDUDXY ! (NCELET-NCEL,3,3,NPHAS)   ! SAUVEGARDE DU GRADIENT DE LA
C        !                           ! VITESSE EN CAS DE ROTATION
C IWDUDX ! (NCELET-NCEL,3,3,NPHAS)   ! TABLEAU DE TRAVAIL LIE A DUDXYZ
C IDRDXY ! (NCELET-NCEL,6,3,NPHAS)   ! SAUVEGARDE DU GRADIENT DE RIJ
C        !                           ! EN CAS DE ROTATION
C IWDRDX ! (NCELET-NCEL,6,3,NPHAS)   ! TABLEAU DE TRAVAIL LIE A DRDXYZ
C
      INTEGER            IDUDXY,IDRDXY,IWDUDX,IWDRDX
      COMMON / IPAUX1 /  IDUDXY,IDRDXY,IWDUDX,IWDRDX
C
C
C FIN
c@z
