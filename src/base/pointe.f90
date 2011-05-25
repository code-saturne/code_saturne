!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2011 EDF S.A., France

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

! Module for pointer variables

module pointe

  !=============================================================================

  use paramx

  !=============================================================================

  !... Auxiliaires independants du nombre de phases
  !      accessibles directement dans ia, ra

  ! Pointeur Dimension       Description
  ! icocg  ! ncelet*9                ! stockage pour gradient
  ! icocgb ! ncelbr*9                ! stockage pour gradient bord
  ! icoci  ! ncelet*9                ! stockage pour gradient si init. par mc
  ! icocib ! ncelbr*9                ! stockage pour gradient bord si init. par
  ! itpuco !    -                    ! mc variables du couplage U-P
  ! idipar ! ncelet                  ! distance a la face de type 5 (phase 1) la
  !                                    plus proche
  ! iyppar ! ncelet                  ! yplus associe (LES only)
  ! iforbr ! nfabor*3                ! efforts aux bords (si posttraite)
  ! iidfst ! nfabor                  ! tableau d'indirection pour les structures
  !                                    mobiles EN ALE

  !... Parametres du module thermique 1D
  ! nfpt1d !                         ! nb faces de bord avec module thermique 1D
  ! inppt1 ! nfpt1d                  ! nombre de mailles dans la paroi
  ! ieppt1 ! nfpt1d                  ! epaisseur de la paroi
  ! irgpt1 ! nfpt1d                  ! raison du maillage
  ! iifpt1 ! nfpt1d                  ! numero de la face
  ! itppt1 ! nfpt1d                  ! temperature de paroi
  ! iiclt1 ! nfpt1d                  ! type de condition limite
  ! itept1 ! nfpt1d                  ! temperature exterieure
  ! ihept1 ! nfpt1d                  ! coefficient d'echange exterieur
  ! ifept1 ! nfpt1d                  ! flux thermique exterieur
  ! ixlmt1 ! nfpt1d                  ! diffusivite thermique
  ! ircpt1 ! nfpt1d                  ! rho*Cp
  ! idtpt1 ! nfpt1d                  ! pas de temps

  integer, save :: icocg  , icocgb , icoci  , icocib ,             &
                   itpuco , idipar , iyppar , iforbr , iidfst ,    &
                   nfpt1d , nmxt1d , inppt1 , iifpt1 , iiclt1 ,    &
                   ieppt1 , irgpt1 , itppt1 ,                      &
                   itept1 , ihept1 , ifept1 ,                      &
                   ixlmt1 , ircpt1 , idtpt1

  !... Auxiliaires accessibles directement dans ia, ra
  !     tous ces tableaux sont (dimension)

  ! Pointeur Dimension       Description
  ! iitypf ! nfabor                  ! type des faces de bord
  ! iitrif ! nfabor                  ! indirection pour tri faces de bord
  ! iyplbr ! nfabor                  ! yplus bord (si post-traite)
  ! iuetbo ! nfabor                  ! uetbor  bord (si LES+VanDriest)
  !        !                         ! ou si le modele de depot est actif
  ! iisymp ! nfabor                  ! zero pour annuler le flux de masse
  !        !                         !   (symetries et parois avec cl couplees)
  !        !                         ! un sinon
  ! iifapa ! ncelet                  ! numero de face de bord 5 la plus proche
  ! ncepdc !                         ! nombre de cellules avec pdc
  ! iicepd ! ncepdc                  ! numero des cellules avec pertedecharge
  ! ickupd ! (ncepdc,6)              ! valeur des coeff de pdc
  ! ncetsm !                         ! nombre de cellules avec tsm
  ! iicesm ! ncetsm                  ! numero des cellules avec tsmasse
  ! iitpsm ! ncetsm                  ! type de tsm
  ! ismace ! ncetsm                  ! valeur de tsm
  ! is2kw  ! ncelet                  ! stockage de 2 Sij.Sij en k-omega
  ! idvukw ! ncelet                  ! stockage de divu en k-omega (en meme
  !                                    temps que s2kw)

  integer, save :: iitypf        , iitrif        ,                 &
                   iisymp        , iyplbr        , iuetbo        , &
                   iifapa, ncepdc,                 &
                   iicepd, ickupd, ncetsm, &
                   iicesm, iitpsm, ismace, &
                   is2kw , idvukw

  !... Auxiliaires accessibles directement dans ia, ra
  !    pour la perodicite

  ! Pointeur Dimension           |   Description
  ! idudxy ! (ncelet-ncel,3,3)   ! sauvegarde du gradient de la
  !        !                     ! vitesse en cas de rotation
  ! iwdudx ! (ncelet-ncel,3,3)   ! tableau de travail lie a dudxyz
  ! idrdxy ! (ncelet-ncel,6,3)   ! sauvegarde du gradient de rij
  !        !                     ! en cas de rotation
  ! iwdrdx ! (ncelet-ncel,6,3)   ! tableau de travail lie a drdxyz

  integer, save :: idudxy, idrdxy, iwdudx, iwdrdx

  !=============================================================================

end module pointe
