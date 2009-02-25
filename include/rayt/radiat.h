!-------------------------------------------------------------------------------

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

!                              radiat.h
!===============================================================================


!-->  IIRAYO = 1 rayonnement, =0 pas de rayonnement
!-->  IRAYON = 1 DOM de la phase 2 P-1 de la phase, =0 pas de rayonnement de la phase
!-->  NPHAST = nombre des phases qui rayonnent
!-->  NPHASC = nombre des phases qui rayonnent augmente eventuellement
!              du nombre de classe (Charbon)
!-->  IRAPHA = numero des phases qui rayonnent
!-->  IIMPAR = 0,1,2 niveau d'impression du calcul des temperatures de paroi
!-->  IIMLUM = 0,1,2 niveau d'impression de la resolution luminance
!-->  IMODAK = 1 calcul du coefficient d'absorption a l'aide de Modak
!            = 0 on n'utilise pas Modak

integer           iirayo ,                                        &
                  irayon(nphsmx)  ,                               &
                  nphast , nphasc ,                               &
                  irapha(nphsmx)  ,                               &
                  iimpar ,                                        &
                  iimlum ,                                        &
                  imodak
common / iiiray / iirayo ,                                        &
                  irayon          ,                               &
                  nphast , nphasc ,                               &
                  irapha          ,                               &
                  iimpar ,                                        &
                  iimlum ,                                        &
                  imodak



!--> pointeur dans le macrotableau RA :

!                       ITSRE --> Terme source explicite
!                       ITSRI --> Terme source implicite
!                       IQX,IQY,IQZ --> Composantes du vecteur densite de flux radiatif
!                       IABS --> part d'absorption dans le terme source explicite
!                       IEMI --> part d'emission dans le terme source explicite
!                       ICAK --> coefficient d'absorption

integer           itsre ,                                         &
                  itsri ,                                         &
                  iqx   ,                                         &
                  iqy   ,                                         &
                  iqz   ,                                         &
                  iabs  ,                                         &
                  iemi  ,                                         &
                  icak

common / iprayo / itsre ,                                         &
                  itsri ,                                         &
                  iqx   ,                                         &
                  iqy   ,                                         &
                  iqz   ,                                         &
                  iabs  ,                                         &
                  iemi  ,                                         &
                  icak

!--> pointeur dans le macrotableau RA :
!                       ITPARO --> temperature de paroi
!                       IQINCI --> densite de flux incident radiatif
!                       IXLAM  --> conductivite thermique de la paroi
!                       IEPA   --> epaisseur de la paroi
!                       IEPS   --> emissivite de la paroi
!                       IFNET  --> Flux Net radiatif
!                       IFCONV --> Flux Convectif
!                       IHCONV --> Coef d'echange fluide

integer           itparo ,                                        &
                  iqinci ,                                        &
                  ixlam  ,                                        &
                  iepa   ,                                        &
                  ieps   ,                                        &
                  ifnet  ,                                        &
                  ifconv ,                                        &
                  ihconv

common / imrayo / itparo ,                                        &
                  iqinci ,                                        &
                  ixlam  ,                                        &
                  iepa   ,                                        &
                  ieps   ,                                        &
                  ifnet  ,                                        &
                  ifconv ,                                        &
                  ihconv


!--> XNP1MX : pour le modele P-1,
!     pourcentage de cellules pour lesquelles on admet que l'epaisseur
!     optique depasse l'unite bien que ce ne soit pas souhaitable

double precision  xnp1mx
common / rrayp1 / xnp1mx

!--> ISTPP1 : pour le modele P-1,
!     indicateur d'arret mis a 1 dans ppcabs si le pourcentage de cellules
!     pour lesquelles l'epaisseur optique depasse l'unite est superieur a
!     XNP1MX  (on s'arrete a la fin du pas de temps)

integer           istpp1
common / irayp1 / istpp1

!--> IDIVER =0 1 ou 2 suivant le calcul du terme source explicite

integer           idiver
common / ikrayo / idiver

!--> parametre sur le nombre de directions de discretisation angulaire

integer     ndirs8
parameter ( ndirs8 = 16 )

!--> suite de calcul (0 : non, 1 : oui)

integer           isuird
common / israyo / isuird

!--> frequence de passage dans le module (=1 si tous les pas de temps)

integer           nfreqr
common / ifrayo / nfreqr

!--> nombre de bandes spectrales

integer nbande
common / ibande / nbande

!--> nombre de directions de discretisation angulaire

integer ndirec
common / idirec / ndirec


!--> Informations sur les zones frontieres

! NBZRDM Nombre max. de  zones frontieres
! NOZRDM Numero max. des zones frontieres

integer    nbzrdm
parameter (nbzrdm=2000)
integer    nozrdm
parameter (nozrdm=2000)

! IIZFRD Pointeur dans IA sur IZFRAD pour reperage des zones
!          frontieres associees aux faces de bord

integer           iizfrd
common / ifrord / iizfrd

! NZFRAD Nombre de zones de bord (sur le proc courant)
! ILZRAY Liste des numeros de zone de bord (du proc courant)
! NOZARM Numero de zone de bord atteint max
!   exemple zones 1 4 2 : NZFRAD=3,NOZARM=4

integer           nozarm(nphsmx), nzfrad(nphsmx),                 &
                  ilzrad(nbzrdm,nphsmx)
common / izonrd / nozarm        , nzfrad        ,                 &
                  ilzrad


!--> Types de condition pour les temperatures de paroi :
!       ITPIMP Profil de temperature imposee
!       IPGRNO Parois grises ou noires
!       IPREFL Parois reflechissante
!       IFGRNO Flux de conduction impose dans la paroi
!                   ET paroi non reflechissante (EPS non nul)
!       IFREFL Flux de conduction impose dans la paroi
!                   ET paroi reflechissante     (EPS = 0)

integer   itpimp   , ipgrno   , iprefl   , ifgrno   , ifrefl
parameter(itpimp=1 , ipgrno=21, iprefl=22, ifgrno=31, ifrefl=32)


!--> sortie postprocessing P0
!    NBRAYP : nombre max de sorties cellules
!    NBRAYF : nombre max de sorties facettes de bord

integer     nbrayp,nbrayf
parameter ( nbrayp = 5 , nbrayf = 8 )

character*80      nbrvap(nbrayp,nphsmx) , nbrvaf(nbrayf,nphsmx)
common / aenray / nbrvap , nbrvaf

integer           irayvp(nbrayp,nphsmx) , irayvf(nbrayf,nphsmx)
common / ienray / irayvp , irayvf

integer           itsray , iqrayp , iabsp , iemip , icakp
parameter (itsray=1 , iqrayp=2 , iabsp=3 , iemip=4 , icakp=5)

integer           itparp , iqincp , ixlamp , iepap  ,             &
                  iepsp  , ifnetp , ifconp , ihconp
parameter (itparp=1 , iqincp=2 , ixlamp=3 , iepap=4  ,            &
           iepsp=5  , ifnetp=6 , ifconp=7 , ihconp=8)

