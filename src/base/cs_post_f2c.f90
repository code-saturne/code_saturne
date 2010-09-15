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

subroutine pstcwr &
!=================

 ( numgep , nomcas , nomrep , nomfmt , optfmt , indmod , ntchr )

!===============================================================================
! FONCTION :
! --------

! CREATION D'UN "WRITER" A PARTIR DES DONNEES FOURNIES PAR LA
! COUCHE FORTRAN : ENCAPSULATION COUCHE C POUR LA TRANSMISSION
! DES LONGUEURS DES CHAINES DE CARACTERES

! UN WRITER CORRESPOND AU CHOIX D'UN NOM DE CAS, DE REPERTOIRE,
! ET DE FORMAT, AINSI QU'UN INDICATEUR PRECISANT SI LES MAILLAGES
! ASSOCIES DOIVENT DEPENDRE OU NON DU TEMPS, ET LA FREQUENCE
! DE SORTIE PAR DEFAUT POUR LES VARIABLES ASSOCIEES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! numgep           ! e  ! <-- ! identificateur du gestionnaire                 !
!                  !    !     ! (< 0 pour gestionnaire reserve,                !
!                  !    !     !  > 0 pour gestionnaire utilisateur)            !
! nomcas           ! a  ! <-- ! nom du cas associe                             !
! nomrep           ! a  ! <-- ! nom du repertoire associe                      !
! nomfmt           ! a  ! <-- ! nom de format associe                          !
! optfmt           ! e  ! <-- ! options associees au format                    !
! indmod           ! e  ! <-- ! 0 : maillages figes                            !
!                  !    !     ! 1 : maillages deformables                      !
!                  !    !     ! 2 : maillages modifiables                      !
! ntchr            ! e  ! <-- ! frequence de sortie par defaut                 !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

character*32     nomcas , nomfmt
character*96     nomrep , optfmt
integer          numgep , indmod , ntchr

! Local variables

integer          lnmcas , lnmrep , lnmfmt , lopfmt

!===============================================================================

lnmcas = len(nomcas)
lnmrep = len(nomrep)
lnmfmt = len(nomfmt)
lopfmt = len(optfmt)

call pstcw1 (numgep, nomcas, nomrep, nomfmt, optfmt,              &
!==========
             lnmcas, lnmrep, lnmfmt, lopfmt,                      &
             indmod, ntchr)

return

end subroutine
subroutine pstcma &
!=================

 ( nummai , nommai , indgrp ,                                     &
   nbrcel , nbrfac , nbrfbr , lstcel , lstfac , lstfbr )

!===============================================================================
! FONCTION :
! --------

! CREATION D'UN MAILLAGE DE POST TRAITEMENT A PARTIR DES DONNEES
! FOURNIES PAR LA COUCHE FORTRAN : ENCAPSULATION COUCHE C
! POUR LA TRANSMISSION DES LONGUEURS DES CHAINES DE CARACTERES

! LES LISTES DE CELLULES OU FACES A EXTRAIRE SONT TRIEES EN SORTIE,
! QU'ELLES LE SOIENT DEJA EN ENTREE OU NON.

! LA LISTE DES CELLULES ASSOCIEES N'EST NECESSAIRE QUE SI LE NOMBRE
! DE CELLULES A EXTRAIRE EST STRICTEMENT SUPERIEUR A 0 ET INFERIEUR
! AU NOMBRE DE CELLULES DU MAILLAGE.

! LES LISTES DE FACES NE SONT PRISES EN COMPTE QUE SI LE NOMBRE DE
! CELLULES A EXTRAIRE EST NUL ; SI LE NOMBRE DE FACES DE BORD A
! EXTRAIRE EST EGAL AU NOMBRE DE FACES DE BORD DU MAILLAGE GLOBAL,
! ET LE NOMBRE DE FACES INTERNES A EXTRAIRE EST NUL, ALORS ON
! EXTRAIT PAR DEFAUT LE MAILLAGE DE BORD, ET LA LISTE DES FACES DE
! BORD ASSOCIEES N'EST DONC PAS NECESSAIRE.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nummai           ! e  ! <-- ! identificateur du maillage                     !
!                  !    !     ! (< 0 pour maillage reserve,   ,                !
!                  !    !     !  > 0 pour maillage utilisateur)                !
! nommai           ! a  ! <-- ! nom du maillage associe                        !
! indgrp           ! e  ! <-- ! 1 to add group information, or O               !
! nbrcel           ! e  ! <-- ! nombre de cellules associees                   !
! nbrfac           ! e  ! <-- ! nombre de faces internes associees             !
! nbrfbr           ! e  ! <-- ! nombre de faces de bord associees              !
! lstcel           ! e  ! <-- ! liste des cellules associees                   !
!                  ! e  !     ! (inutile si nbrcel >= ncel)                    !
! lstfac           ! e  ! <-- ! liste des faces internes associees             !
! lstfbr           ! e  ! <-- ! liste des faces de bord associees              !
!                  ! e  !     ! (inutile si    nbrfbr = nfabor                 !
!                  ! e  !     !             et nbrfac = 0     )                !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

character*32     nommai
integer          nummai, indgrp, nbrcel, nbrfac, nbrfbr

integer          lstcel(nbrcel), lstfac(nbrfac), lstfbr(nbrfbr)

! Local variables

integer          lnmmai

!===============================================================================

lnmmai = len(nommai)

call pstcm1 (nummai, nommai, lnmmai, indgrp,                      &
!==========
             nbrcel, nbrfac, nbrfbr, lstcel, lstfac, lstfbr)

return

end subroutine
subroutine psteva &
!================

 ( nummai , nomvar , dimvar , ientla , ivarpr , ntcabs , ttcabs , &
   varcel , varfac , varfbo )

!===============================================================================
! FONCTION :
! --------

! ECRITURE D'UN CHAMP DE POST TRAITEMENT ASSOCIE AUX CELLULES
! OU FACES D'UN MAILLAGE A PARTIR DES DONNEES FOURNIES PAR LA
! COUCHE FORTRAN :
! ENCAPSULATION COUCHE C POUR LA TRANSMISSION DES LONGUEURS DES
! CHAINES DE CARACTERES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nummai           ! a  ! <-- ! numero du maillage associe                     !
! nomvar           ! e  ! <-- ! nom de la variable associee                    !
! dimvar           ! e  ! <-- ! 1 pour scalaire, 3 pour vecteur                !
! ientla           ! e  ! <-- ! si vecteur, 1 si valeurs entrelacees           !
!                  !    !     ! (x1, y1, z1, x2, y2, ..., yn, zn),             !
!                  !    !     ! 0 sinon (x1, x2, ...xn, y1, y2, ...)           !
! ivarpr           ! e  ! <-- ! 1 si variable definie sur maillage             !
!                  !    !     ! "parent", 0 si variable restreinte             !
!                  !    !     ! au maillage nummai                             !
! ntcabs           ! e  ! <-- ! numero de pas de temps (-1 pour une            !
!                  !    !     ! variable independante du temps)                !
! ttcabs           ! r  ! <-- ! temps physique associe                         !
! varcel(*)        ! r  ! <-- ! valeurs aux cellules associees                 !
! varfac(*)        ! r  ! <-- ! valeurs aux faces internes associees           !
! varfbo(*)        ! r  ! <-- ! valeurs aux faces de bord associees            !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

character*32     nomvar
integer          nummai, dimvar, ientla, ivarpr, ntcabs
double precision ttcabs, varcel(*), varfac(*), varfbo(*)

! Local variables

integer          lnmvar

!===============================================================================

lnmvar = len(nomvar)

call pstev1 (nummai, nomvar, lnmvar, dimvar, ientla, ivarpr,      &
!==========
             ntcabs, ttcabs, varcel, varfac, varfbo)

return

end subroutine

subroutine pstsnv &
!================

 ( nomvar , nomva2 , nomva3 )

!===============================================================================
! FONCTION :
! --------

! SUPPRESSION DU CARACTERE X, x, OU 1 D'UNE CHAINE DE CARACTERES
! FORTRAN SI LES CHAINES COMPAREES SONT IDENTIQUES AU DERNIER
! CARACTERE PRES, RESPECTIVEMENT Y, y, OU 2 ET Z, z, OU 3

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nomvar           ! e  ! <-- ! nom de la variable associee                    !
! nomva2           ! e  ! <-- ! nom de la variable 2 associee                  !
! nomva3           ! e  ! <-- ! nom de la variable 3 associee                  !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

!===============================================================================

implicit none

! Arguments

character*32     nomvar, nomva2, nomva3

! Local variables

integer          ii, jj
integer          lnmvar, lnmva2, lnmva3

!===============================================================================

lnmvar = len(nomvar)
lnmva2 = len(nomva2)
lnmva3 = len(nomva3)

if ((lnmvar .eq. lnmva2) .and. (lnmvar .eq. lnmva3)) then

  do 10 ii = lnmvar, 1, -1
    if (     nomvar(ii:ii) .ne. ' '                               &
        .or. nomva2(ii:ii) .ne. ' '                               &
        .or. nomva3(ii:ii) .ne. ' ') then
      goto 20
    endif
 10     continue

 20     continue

  if (ii .gt. 1) then

    jj = ii

    ! Handle the case where the next-to-last character changes, such
    ! as with VelocityX1, VelocityX2, ... in case of a calculation
    ! with multiple phases.

    if (      (ii .gt. 2)                                         &
        .and. (nomvar(ii:ii) .eq. nomva2(ii:ii))                  &
        .and. (nomvar(ii:ii) .eq. nomva3(ii:ii))) then
      ii = jj-1
    endif

    ! Remove the character related to the spatial axis

    if (      nomvar(ii:ii) .eq. 'X'                              &
        .and. nomva2(ii:ii) .eq. 'Y'                              &
        .and. nomva3(ii:ii) .eq. 'Z') then
      nomvar(ii:ii) = ' '
    else if (      nomvar(ii:ii) .eq. 'x'                         &
             .and. nomva2(ii:ii) .eq. 'y'                         &
             .and. nomva3(ii:ii) .eq. 'z') then
      nomvar(ii:ii) = ' '
    else if (      nomvar(ii:ii) .eq. 'U'                         &
             .and. nomva2(ii:ii) .eq. 'V'                         &
             .and. nomva3(ii:ii) .eq. 'W') then
      nomvar(ii:ii) = ' '
    else if (      nomvar(ii:ii) .eq. 'u'                         &
             .and. nomva2(ii:ii) .eq. 'v'                         &
             .and. nomva3(ii:ii) .eq. 'w') then
      nomvar(ii:ii) = ' '
    else if (      nomvar(ii:ii) .eq. '1'                         &
             .and. nomva2(ii:ii) .eq. '2'                         &
             .and. nomva3(ii:ii) .eq. '3') then
      nomvar(ii:ii) = ' '
    endif

    ! If the next-to last character was removed, the last one must be shifted.

    if (ii .eq. jj+1) then
      nomvar(ii:ii) = nomvar(jj:jj)
      nomvar(jj:jj) = ' '
    endif

  endif

endif

return

end subroutine

