!-------------------------------------------------------------------------------

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

!                              paramx.h
!===============================================================================





!                         =================
!                         =================

!                             ATTENTION

!                         =================
!                         =================






!              LA MODIFICATION DES PARAMETRES CI DESSOUS



!                           EST INTERDITE

!                         =================
!                         =================








!       Elle demande la recompilation de la totalite de la bibliotheque
!         operation qui ne peut etre effectuee que si l'on dispose de
!         la totalite des sources.

















! PARAMETRES DIVERS
! =================

! NPHSMX : NOMBRE MAX DE PHASES
!          (keep it coherent with CS_NPHSMX in cs_perio.h)
! NSCAMX : NOMBRE MAX DE SCALAIRES
! NVARMX : NOMBRE MAX DE VARIABLES =
!          NOMBRE MAX DE SCALAIRES + 12 (U,V,W,P,Rij,E,ALP)*NPHSMX
! NPRCMX : NOMBRE MAX DE PROPRIETES PHYSIQUES AUX CELLULES (TOTAL) =
!          NSCAMX (Lambda) + 7 (RHO,CP,VISCL,VISCT,COU,FOU,IPRTOT) NPHSMX
!                          + 4 (ESTIM) NPHSMX
! NPRFMX : NOMBRE MAX DE PROPRIETES PHYSIQUES AUX FACES INTERNES =
!          NSCAMX (Flumas) + 2*NPHSMX(Flumas,ALP)
! NPRBMX : NOMBRE MAX DE PROPRIETES PHYSIQUES AUX FACES DE BORD =
!          NSCAMX (Flumab) + 3*NPHSMX(Flumab,ALP, ROMB)
! NPROMX : NOMBRE MAX DE PROPRIETES PHYSIQUES TOUT CONFONDU
!          Majore par NPRCMX+NPRFMX+NPRBMX
! NGRDMX : NOMBRE MAX DE GRANDEURS =
!          NVARMX + NPROMX
! NSMAMX : NOMBRE MAX DE CASES POUR LES TABLEAUX TERMES SOURCE DE MASSE
!          NVARMX + NPHSMX pour SMACEL
! NVPPMX : NOMBRE DE VARIABLES POUR AFFICHAGES
!          NGRDMX + 20 (20 couvre DT, TPUCOU, et une marge de 16 ...)

integer   nphsmx, nscamx, nvarmx, nprcmx, nprfmx, nprbmx, npromx
integer   ngrdmx, nsmamx, nvppmx
parameter(nphsmx=1)
parameter(nscamx=200)
parameter(nvarmx=nscamx+12*nphsmx)
parameter(nprcmx=nscamx+11*nphsmx)
parameter(nprfmx=nscamx+ 2*nphsmx)
parameter(nprbmx=nscamx+ 3*nphsmx)
parameter(npromx=nprcmx+ nprfmx+nprbmx)
parameter(ngrdmx=nvarmx+ npromx)
parameter(nsmamx=nvarmx+ nphsmx)
parameter(nvppmx=ngrdmx+20)

! NUSHMX : NOMBRE MAX DE FICHIERS UTILISATEUR POUR HISTORIQUES
integer    nushmx
parameter(nushmx=16)

! NUSRMX : NOMBRE MAX DE FICHIERS UTILISATEUR
integer    nusrmx
parameter(nusrmx=10)

! NCAPTM : NOMBRE MAX DE SONDES (POUR HISTORIQUES)
!          Voir le format associe dans ecrhis
integer    ncaptm
parameter(ncaptm=100)

! NTYPMX NOMBRE DE TYPES DE CONDITIONS AUX LIMITES POSSIBLES

integer    ntypmx
parameter(ntypmx=200)

integer    iindef, ientre, isolib, isymet, iparoi,                &
   iparug, iesicf, isspcf, isopcf, ierucf, ieqhcf

parameter(iindef=1, ientre=2, isolib=3, isymet=4, iparoi=5,       &
 iparug=6, iesicf=7, isspcf=8, isopcf=9, ierucf=10, ieqhcf=11)

! NESTMX : NOMBRE MAX D'ESTIMATEURS
!  IESPRE, IESDER, IESCOR, IESTOT : Numeros
integer    nestmx
parameter (nestmx=4)
integer    iespre  , iesder  , iescor  , iestot
parameter (iespre=1, iesder=2, iescor=3, iestot=4)

! NBMOMX : NOMBRE MAX DE MOYENNES (MOMENTS) CALCULE
! NDGMOX : DEGRE MAX DES MOMENTS
integer    nbmomx, ndgmox
parameter (nbmomx = 50, ndgmox = 5)

! IPST* : SELECTION POST TRAITEMENT AUTOMATIQUE BORD : VOIR IPSTDV

integer    ipstyp  , ipstcl  , ipstft, ipstfo
parameter (ipstyp=2, ipstcl=3, ipstft=5, ipstfo=7)

! CONDITIONS AUX LIMITES POSSIBLES POUR LA VITESSE DE MAILLAGE EN ALE

integer    ibfixe, igliss, ivimpo
parameter(ibfixe=1, igliss=2, ivimpo=3 )

! NOMBRE DE STRUCTURES MAX EN ALE

integer nstrmx
parameter (nstrmx=20)

! NOMBRE DE STRUCTURES MAX EN ALE ET COUPLAGE CODE_ASTER

integer nastmx
parameter (nastmx=20)

! DIMENSIONS MAXIMALES DES TABLEAUX CONTENANT LES BORNES PAR THREAD

integer nthrd1, nthrd2
parameter (nthrd1=4, nthrd2=16)

