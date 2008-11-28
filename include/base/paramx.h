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
C                              paramx.h
C***********************************************************************
C
C
C
C
C
C                         =================
C                         =================
C
C                             ATTENTION
C
C                         =================
C                         =================
C
C
C
C
C
C
C              LA MODIFICATION DES PARAMETRES CI DESSOUS
C
C
C
C                           EST INTERDITE
C
C                         =================
C                         =================
C
C
C
C
C
C
C
C
C       Elle demande la recompilation de la totalite de la bibliotheque
C         operation qui ne peut etre effectuee que si l'on dispose de
C         la totalite des sources.
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C
C PARAMETRES DIVERS
C =================
C
C NPHSMX : NOMBRE MAX DE PHASES
C          (keep it coherent with CS_NPHSMX in cs_perio.h)
C NSCAMX : NOMBRE MAX DE SCALAIRES
C NVARMX : NOMBRE MAX DE VARIABLES =
C          NOMBRE MAX DE SCALAIRES + 12 (U,V,W,P,Rij,E,ALP)*NPHSMX
C NPRCMX : NOMBRE MAX DE PROPRIETES PHYSIQUES AUX CELLULES (TOTAL) =
C          NSCAMX (Lambda) + 7 (RHO,CP,VISCL,VISCT,COU,FOU,IPRTOT) NPHSMX
C                          + 4 (ESTIM) NPHSMX
C NPRFMX : NOMBRE MAX DE PROPRIETES PHYSIQUES AUX FACES INTERNES =
C          NSCAMX (Flumas) + 2*NPHSMX(Flumas,ALP)
C NPRBMX : NOMBRE MAX DE PROPRIETES PHYSIQUES AUX FACES DE BORD =
C          NSCAMX (Flumab) + 3*NPHSMX(Flumab,ALP, ROMB)
C NPROMX : NOMBRE MAX DE PROPRIETES PHYSIQUES TOUT CONFONDU
C          Majore par NPRCMX+NPRFMX+NPRBMX
C NGRDMX : NOMBRE MAX DE GRANDEURS =
C          NVARMX + NPROMX
C NSMAMX : NOMBRE MAX DE CASES POUR LES TABLEAUX TERMES SOURCE DE MASSE
C          NVARMX + NPHSMX pour SMACEL
C NVPPMX : NOMBRE DE VARIABLES POUR AFFICHAGES
C          NGRDMX + 20 (20 couvre DT, TPUCOU, et une marge de 16 ...)
C
      INTEGER   NPHSMX, NSCAMX, NVARMX, NPRCMX, NPRFMX, NPRBMX, NPROMX
      INTEGER   NGRDMX, NSMAMX, NVPPMX
      PARAMETER(NPHSMX=1)
      PARAMETER(NSCAMX=200)
      PARAMETER(NVARMX=NSCAMX+12*NPHSMX)
      PARAMETER(NPRCMX=NSCAMX+11*NPHSMX)
      PARAMETER(NPRFMX=NSCAMX+ 2*NPHSMX)
      PARAMETER(NPRBMX=NSCAMX+ 3*NPHSMX)
      PARAMETER(NPROMX=NPRCMX+ NPRFMX+NPRBMX)
      PARAMETER(NGRDMX=NVARMX+ NPROMX)
      PARAMETER(NSMAMX=NVARMX+ NPHSMX)
      PARAMETER(NVPPMX=NGRDMX+20)
C
C NUSHMX : NOMBRE MAX DE FICHIERS UTILISATEUR POUR HISTORIQUES
      INTEGER    NUSHMX
      PARAMETER(NUSHMX=16)
C
C NUSRMX : NOMBRE MAX DE FICHIERS UTILISATEUR
      INTEGER    NUSRMX
      PARAMETER(NUSRMX=10)
C
C NCAPTM : NOMBRE MAX DE SONDES (POUR HISTORIQUES)
C          Voir le format associe dans ecrhis
      INTEGER    NCAPTM
      PARAMETER(NCAPTM=100)
C
C NTYPMX NOMBRE DE TYPES DE CONDITIONS AUX LIMITES POSSIBLES
C
      INTEGER    NTYPMX
      PARAMETER(NTYPMX=200)
C
      INTEGER    IINDEF, IENTRE, ISOLIB, ISYMET, IPAROI,
     &   IPARUG, IESICF, ISSPCF, ISOPCF, IERUCF, IEQHCF
C
      PARAMETER(IINDEF=1, IENTRE=2, ISOLIB=3, ISYMET=4, IPAROI=5,
     & IPARUG=6, IESICF=7, ISSPCF=8, ISOPCF=9, IERUCF=10, IEQHCF=11)
C
C NESTMX : NOMBRE MAX D'ESTIMATEURS
C  IESPRE, IESDER, IESCOR, IESTOT : Numeros
      INTEGER    NESTMX
      PARAMETER (NESTMX=4)
      INTEGER    IESPRE  , IESDER  , IESCOR  , IESTOT
      PARAMETER (IESPRE=1, IESDER=2, IESCOR=3, IESTOT=4)
C
C NBMOMX : NOMBRE MAX DE MOYENNES (MOMENTS) CALCULE
C NDGMOX : DEGRE MAX DES MOMENTS
      INTEGER    NBMOMX, NDGMOX
      PARAMETER (NBMOMX = 50, NDGMOX = 5)
C
C IPST* : SELECTION POST TRAITEMENT AUTOMATIQUE BORD : VOIR IPSTDV
C
      INTEGER    IPSTYP  , IPSTCL  , IPSTFT, IPSTFO
      PARAMETER (IPSTYP=2, IPSTCL=3, IPSTFT=5, IPSTFO=7)
C
C CONDITIONS AUX LIMITES POSSIBLES POUR LA VITESSE DE MAILLAGE EN ALE
C
      INTEGER    IBFIXE, IGLISS, IVIMPO
      PARAMETER(IBFIXE=1, IGLISS=2, IVIMPO=3 )
C
C NOMBRE DE STRUCTURES MAX EN ALE
C
      INTEGER NSTRMX
      PARAMETER (NSTRMX=20)
C
C NOMBRE DE STRUCTURES MAX EN ALE ET COUPLAGE CODE_ASTER
C
      INTEGER NASTMX
      PARAMETER (NASTMX=20)
C
c@z
