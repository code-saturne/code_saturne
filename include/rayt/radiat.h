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
C                              radiat.h
C***********************************************************************
C
C
C-->  IIRAYO = 1 rayonnement, =0 pas de rayonnement
C-->  IRAYON = 1 DOM de la phase 2 P-1 de la phase, =0 pas de rayonnement de la phase
C-->  NPHAST = nombre des phases qui rayonnent
C-->  NPHASC = nombre des phases qui rayonnent augmente eventuellement
C              du nombre de classe (Charbon)
C-->  IRAPHA = numero des phases qui rayonnent
C-->  IIMPAR = 0,1,2 niveau d'impression du calcul des temperatures de paroi
C-->  IIMLUM = 0,1,2 niveau d'impression de la resolution luminance
C-->  IMODAK = 1 calcul du coefficient d'absorption a l'aide de Modak
C            = 0 on n'utilise pas Modak
C
      INTEGER           IIRAYO ,
     &                  IRAYON(NPHSMX)  ,
     &                  NPHAST , NPHASC ,
     &                  IRAPHA(NPHSMX)  ,
     &                  IIMPAR ,
     &                  IIMLUM ,
     &                  IMODAK
      COMMON / IIIRAY / IIRAYO ,
     &                  IRAYON          ,
     &                  NPHAST , NPHASC ,
     &                  IRAPHA          ,
     &                  IIMPAR ,
     &                  IIMLUM ,
     &                  IMODAK
C
C
C
C--> pointeur dans le macrotableau RA :
C
C                       ITSRE --> Terme source explicite
C                       ITSRI --> Terme source implicite
C                       IQX,IQY,IQZ --> Composantes du vecteur densite de flux radiatif
C                       IABS --> part d'absorption dans le terme source explicite
C                       IEMI --> part d'emission dans le terme source explicite
C                       ICAK --> coefficient d'absorption
C
      INTEGER           ITSRE ,
     &                  ITSRI ,
     &                  IQX   ,
     &                  IQY   ,
     &                  IQZ   ,
     &                  IABS  ,
     &                  IEMI  ,
     &                  ICAK
C
      COMMON / IPRAYO / ITSRE ,
     &                  ITSRI ,
     &                  IQX   ,
     &                  IQY   ,
     &                  IQZ   ,
     &                  IABS  ,
     &                  IEMI  ,
     &                  ICAK
C
C--> pointeur dans le macrotableau RA :
C                       ITPARO --> temperature de paroi
C                       IQINCI --> densite de flux incident radiatif
C                       IXLAM  --> conductivite thermique de la paroi
C                       IEPA   --> epaisseur de la paroi
C                       IEPS   --> emissivite de la paroi
C                       IFNET  --> Flux Net radiatif
C                       IFCONV --> Flux Convectif
C                       IHCONV --> Coef d'echange fluide
C
      INTEGER           ITPARO ,
     &                  IQINCI ,
     &                  IXLAM  ,
     &                  IEPA   ,
     &                  IEPS   ,
     &                  IFNET  ,
     &                  IFCONV ,
     &                  IHCONV
C
      COMMON / IMRAYO / ITPARO ,
     &                  IQINCI ,
     &                  IXLAM  ,
     &                  IEPA   ,
     &                  IEPS   ,
     &                  IFNET  ,
     &                  IFCONV ,
     &                  IHCONV
C
C
C--> XNP1MX : pour le modele P-1,
C     pourcentage de cellules pour lesquelles on admet que l'epaisseur
C     optique depasse l'unite bien que ce ne soit pas souhaitable
C
      DOUBLE PRECISION  XNP1MX
      COMMON / RRAYP1 / XNP1MX
C
C--> ISTPP1 : pour le modele P-1,
C     indicateur d'arret mis a 1 dans ppcabs si le pourcentage de cellules
C     pour lesquelles l'epaisseur optique depasse l'unite est superieur a
C     XNP1MX  (on s'arrete a la fin du pas de temps)
C
      INTEGER           ISTPP1
      COMMON / IRAYP1 / ISTPP1
C
C--> IDIVER =0 1 ou 2 suivant le calcul du terme source explicite
C
      INTEGER           IDIVER
      COMMON / IKRAYO / IDIVER
C
C--> parametre sur le nombre de directions de discretisation angulaire
C
      INTEGER     NDIRS8
      PARAMETER ( NDIRS8 = 16 )
C
C--> suite de calcul (0 : non, 1 : oui)
C
      INTEGER           ISUIRD
      COMMON / ISRAYO / ISUIRD
C
C--> frequence de passage dans le module (=1 si tous les pas de temps)
C
      INTEGER           NFREQR
      COMMON / IFRAYO / NFREQR
C
C--> nombre de bandes spectrales
C
      INTEGER NBANDE
      COMMON / IBANDE / NBANDE
C
C--> nombre de directions de discretisation angulaire
C
      INTEGER NDIREC
      COMMON / IDIREC / NDIREC
C
C
C--> Informations sur les zones frontieres
C
C NBZRDM Nombre max. de  zones frontieres
C NOZRDM Numero max. des zones frontieres
C
      INTEGER    NBZRDM
      PARAMETER (NBZRDM=2000)
      INTEGER    NOZRDM
      PARAMETER (NOZRDM=2000)
C
C IIZFRD Pointeur dans IA sur IZFRAD pour reperage des zones
C          frontieres associees aux faces de bord
C
      INTEGER           IIZFRD
      COMMON / IFRORD / IIZFRD
C
C NZFRAD Nombre de zones de bord (sur le proc courant)
C ILZRAY Liste des numeros de zone de bord (du proc courant)
C NOZARM Numero de zone de bord atteint max
C   exemple zones 1 4 2 : NZFRAD=3,NOZARM=4
C
      INTEGER           NOZARM(NPHSMX), NZFRAD(NPHSMX),
     &                  ILZRAD(NBZRDM,NPHSMX)
      COMMON / IZONRD / NOZARM        , NZFRAD        ,
     &                  ILZRAD
C
C
C--> Types de condition pour les temperatures de paroi :
C       ITPIMP Profil de temperature imposee
C       IPGRNO Parois grises ou noires
C       IPREFL Parois reflechissante
C       IFGRNO Flux de conduction impose dans la paroi
C                   ET paroi non reflechissante (EPS non nul)
C       IFREFL Flux de conduction impose dans la paroi
C                   ET paroi reflechissante     (EPS = 0)
C
      INTEGER   ITPIMP   , IPGRNO   , IPREFL   , IFGRNO   , IFREFL
      PARAMETER(ITPIMP=1 , IPGRNO=21, IPREFL=22, IFGRNO=31, IFREFL=32)
C
C
C--> sortie postprocessing P0
C    NBRAYP : nombre max de sorties cellules
C    NBRAYF : nombre max de sorties facettes de bord
C
      INTEGER     NBRAYP,NBRAYF
      PARAMETER ( NBRAYP = 5 , NBRAYF = 8 )
C
      CHARACTER*80      NBRVAP(NBRAYP,NPHSMX) , NBRVAF(NBRAYF,NPHSMX)
      COMMON / AENRAY / NBRVAP , NBRVAF
C
      INTEGER           IRAYVP(NBRAYP,NPHSMX) , IRAYVF(NBRAYF,NPHSMX)
      COMMON / IENRAY / IRAYVP , IRAYVF
C
      INTEGER           ITSRAY , IQRAYP , IABSP , IEMIP , ICAKP
      PARAMETER (ITSRAY=1 , IQRAYP=2 , IABSP=3 , IEMIP=4 , ICAKP=5)
C
      INTEGER           ITPARP , IQINCP , IXLAMP , IEPAP  ,
     &                  IEPSP  , IFNETP , IFCONP , IHCONP
      PARAMETER (ITPARP=1 , IQINCP=2 , IXLAMP=3 , IEPAP=4  ,
     &           IEPSP=5  , IFNETP=6 , IFCONP=7 , IHCONP=8)
C
c@z
