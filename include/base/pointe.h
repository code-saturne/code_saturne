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

!                              pointe.h
!===============================================================================

!... GEOMETRIE
!      ACCESSIBLES DIRECTEMENT DANS IA, RA

! Pointeur Dimension       Description
! IICELB ! NCELBR                  ! REPERAGE DES CELLULES DE BORD
! ISRFAN ! NFAC                    ! SURFACE DES FACES INTERNES
! ISRFBN ! NFABOR                  ! SURFACE DES FACES DE BORD
! IDIST  ! NFAC                    ! IJ.Nij
! IDISTB ! NFABOR                  ! EQUIVALENT DE IDIST AU BORD
! IPOND  ! NFAC                    ! PONDERATION (Aij=POND Ai+(1-POND)Aj)
! IDIJPF ! NFAC*NDIM               ! VECTEUR I'J'
! IDIIPB ! NFAC*NDIM               ! EQUIVALENT DE IDIIPF AU BORD
! IDOFIJ ! NFAC*NDIM               ! VECTEUR OF A LA FACE IJ

integer           iicelb ,                                        &
                  isrfan , isrfbn , idist  , idistb , ipond ,     &
                  idijpf , idiipb , idofij
common / ipgeom / iicelb ,                                        &
                  isrfan , isrfbn , idist  , idistb , ipond ,     &
                  idijpf , idiipb , idofij


!... AUXILIAIRES INDEPENDANTS DU NOMBRE DE PHASES
!      ACCESSIBLES DIRECTEMENT DANS IA, RA

! Pointeur Dimension       Description
! ICOCG  ! NCELET*9                ! STOCKAGE POUR GRADIENT
! ICOCGB ! NCELBR*9                ! STOCKAGE POUR GRADIENT BORD
! ICOCI  ! NCELET*9                ! STOCKAGE POUR GRADIENT SI INIT. PAR MC
! ICOCIB ! NCELBR*9                ! STOCKAGE POUR GRADIENT BORD SI INIT. PAR MC
! ITPUCO !    -                    ! VARIABLES DU COUPLAGE U-P
! IDIPAR ! NCELET                  ! DISTANCE A LA FACE DE TYPE 5 (PHASE 1) LA
!                            PLUS PROCHE
! IYPPAR ! NCELET                  ! YPLUS ASSOCIE (LES only)
! IFORBR ! NFABOR*3                ! EFFORTS AUX BORDS (SI POSTTRAITE)
! IIDFST ! NFABOR                  ! TABLEAU D'INDIRECTION POUR LES STRUCTURES
!                          MOBILES EN ALE

!... PARAMETRES DU MODULE THERMIQUE 1D
! NFPT1D !                         ! NB DE FACES DE BORD AVEC MODULE THERMIQUE 1D
! INPPT1 ! NFPT1D                  ! NOMBRE DE MAILLES DANS LA PAROI
! IEPPT1 ! NFPT1D                  ! EPAISSEUR DE LA PAROI
! IRGPT1 ! NFPT1D                  ! RAISON DU MAILLAGE
! IIFPT1 ! NFPT1D                  ! NUMERO DE LA FACE
! ITPPT1 ! NFPT1D                  ! TEMPERATURE DE PAROI
! IICLT1 ! NFPT1D                  ! TYPE DE CONDITION LIMITE
! ITEPT1 ! NFPT1D                  ! TEMPERATURE EXTERIEURE
! IHEPT1 ! NFPT1D                  ! COEFFICIENT D ECHANGE EXTERIEUR
! IFEPT1 ! NFPT1D                  ! FLUX THERMIQUE EXTERIEUR
! IXLMT1 ! NFPT1D                  ! DIFFUSIVITE THERMIQUE
! IRCPT1 ! NFPT1D                  ! RHO*CP
! IDTPT1 ! NFPT1D                  ! PAS DE TEMPS

integer           icocg  , icocgb , icoci  , icocib ,             &
                  itpuco , idipar , iyppar , iforbr , iidfst ,    &
                  nfpt1d , nmxt1d , inppt1 , iifpt1 , iiclt1 ,    &
                  ieppt1 , irgpt1 , itppt1 ,                      &
                  itept1 , ihept1 , ifept1 ,                      &
                  ixlmt1 , ircpt1 , idtpt1

common / ipaux0 / icocg  , icocgb , icoci  , icocib ,             &
                  itpuco , idipar , iyppar , iforbr , iidfst ,    &
                  nfpt1d , nmxt1d , inppt1 , iifpt1 , iiclt1 ,    &
                  ieppt1 , irgpt1 , itppt1 ,                      &
                  itept1 , ihept1 , ifept1 ,                      &
                  ixlmt1 , ircpt1 , idtpt1

!... AUXILIAIRES ACCESSIBLES DIRECTEMENT DANS IA, RA
!     TOUS CES TABLEAUX SONT (Dimension,NPHAS)

! Pointeur Dimension       Description
! IITYPF ! NFABOR                  ! TYPE DES FACES DE BORD
! IITRIF ! NFABOR                  ! INDIRECTION POUR TRI FACES DE BORD
! IYPLBR ! NFABOR                  ! YPLUS BORD (SI POST-TRAITE)
! IISYMP ! NFABOR                  ! ZERO POUR ANNULER LE FLUX DE MASSE
!        !                         !   (SYMETRIES ET PAROIS AVEC CL COUPLEES)
!        !                         ! UN SINON
! IIFAPA ! NCELET                  ! NUMERO DE FACE DE BORD 5 LA PLUS PROCHE
! NCEPDC !                         ! NOMBRE DE CELLULES AVEC PDC
! IICEPD ! NCEPDC                  ! NUMERO DES CELLULES AVEC PerteDeCharge
! ICKUPD ! (NCEPDC,6)              ! VALEUR DES COEFF DE PDC
! NCETSM !                         ! NOMBRE DE CELLULES AVEC TSM
! IICESM ! NCETSM                  ! NUMERO DES CELLULES AVEC TSMasse
! IITPSM ! NCETSM                  ! TYPE DE TSM
! ISMACE ! NCETSM                  ! VALEUR DE TSM
! IS2KW  ! NCELET                  ! STOCKAGE DE 2 Sij.Sij EN K-OMEGA
! IDVUKW ! NCELET                  ! STOCKAGE DE DIVU EN K-OMEGA (EN MEME TEMPS QUE S2KW)

integer           iitypf        , iitrif        ,                 &
                  iisymp        , iyplbr        ,                 &
                  iifapa(nphsmx), ncepdc(nphsmx),                 &
                  iicepd(nphsmx), ickupd(nphsmx), ncetsm(nphsmx), &
                  iicesm(nphsmx), iitpsm(nphsmx), ismace(nphsmx), &
                  is2kw (nphsmx), idvukw(nphsmx)
common / iposup / iitypf        , iitrif        ,                 &
                  iisymp        , iyplbr        ,                 &
                  iifapa        , ncepdc        ,                 &
                  iicepd        , ickupd        , ncetsm        , &
                  iicesm        , iitpsm        , ismace        , &
                  is2kw         , idvukw


!... AUXILIAIRES ACCESSIBLES DIRECTEMENT DANS IA, RA
!    POUR LA PERODICITE

! Pointeur Dimension                 |   Description
! IDUDXY ! (NCELET-NCEL,3,3,NPHAS)             ! SAUVEGARDE DU GRADIENT DE LA
!        !                                     ! VITESSE EN CAS DE ROTATION
! IWDUDX ! (NCELET-NCEL,3,3,NPHAS)             ! TABLEAU DE TRAVAIL LIE A DUDXYZ
! IDRDXY ! (NCELET-NCEL,6,3,NPHAS)             ! SAUVEGARDE DU GRADIENT DE RIJ
!        !                                     ! EN CAS DE ROTATION
! IWDRDX ! (NCELET-NCEL,6,3,NPHAS)             ! TABLEAU DE TRAVAIL LIE A DRDXYZ

integer            idudxy,idrdxy,iwdudx,iwdrdx
common / ipaux1 /  idudxy,idrdxy,iwdudx,iwdrdx


! FIN
