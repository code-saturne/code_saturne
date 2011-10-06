!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2011 EDF S.A.
!
! This program is free software; you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation; either version 2 of the License, or (at your option) any later
! version.
!
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
! details.
!
! You should have received a copy of the GNU General Public License along with
! this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
! Street, Fifth Floor, Boston, MA 02110-1301, USA.

!-------------------------------------------------------------------------------

subroutine ouestu &
!================

 ( nfecra , ndim   , npoint ,                                     &
   ierror ,                                                       &
   px     , py     , pz     ,                                     &
   qx     , qy     , qz     ,                                     &
   cdgx   , cdgy   , cdgz   ,                                     &
   celx   , cely   , celz   ,                                     &
   itypf7 , iconf7 , xyznod ,                                     &
   indian )


!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!  reperage d'une particule par rapport a une face :
!  soient P le point de depart de la particule
!  et Q le point d'arrivee ; on obtient :

!  INDIAN =  0 le rayon PQ ne sort pas de la cellule par cette face
!  INDIAN = -1 meme cellule
!  INDIAN =  1 sortie de la cellule par cette face

!  ATTENTION : POUR UTILISER CES SUBROUTINE DE REPERAGE, IL EST
!  NECESSAIRE QUE LE TYPE DOUBLE PRECISION SOIT CODE SUR 8 OCTETS,
!  i.e. DOUBLE PRECISION = REAL*8

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nfecra           ! e  ! <-- ! unite du fichier de sortie listing             !
! ndim             ! e  ! <-- ! dimension de l'espace (=3)                     !
! npoint           ! e  ! <-- ! nombre de noeuds                               !
! ierror           ! e  ! --> ! indicateur d'erreur                            !
! px,py,pz         ! r  ! <-- ! point de depart de la particule                !
! qx,qy,qz         ! r  ! <-- ! point de d'arrive de la particule              !
! cdgx,..,cdgz     ! r  ! <-- ! centre de gravite de la face                   !
! celx,..,celz     ! r  ! <-- ! centre de gravite de la cellule                !
!                  !    !     !   contenant le point de depart                 !
! itypf7           ! e  ! <-- ! nombre de points support de la face            !
! iconf7(itypf7    ! e  ! <-- ! numero des points support                      !
! xyznod           ! tr ! <-- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
! indian           ! e  ! --> ! indicateur d'orientation                       !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use cstnum
use period

!===============================================================================

implicit none

! Arguments

integer          nfecra , ndim , npoint
integer          indian , ierror
integer          itypf7 , iconf7(itypf7)

double precision px , py , pz
double precision qx , qy , qz
double precision cdgx , cdgy , cdgz
double precision celx , cely , celz
double precision xyznod(ndim,npoint)

! Local variables

integer          in , is , isign , il , it , ii
integer          ijklug , ittour , ipturb , isensf , ipos

double precision valmax
double precision cr1x , cr1y , cr1z
double precision cr2x , cr2y , cr2z

!===============================================================================


!===============================================================================
! -1.  MACRO DE DEBUGGAGE DEVELOPPEUR
!===============================================================================

!            ATTENTION INTERVENTION DEVELOPPEUR UNIQUEMENT.

!     Si cette macro est vrai elle permet les impressions
!     dans le listing du detail des erreurs de reperage.

!     PAR DEFAUT : DEBUG_OUESTU = 0


#define DEBUG_OUESTU 0

!===============================================================================
! 0. Initialisations
!===============================================================================

ierror = 0

ittour = 0

indian = 0

ijklug = 0

ipturb = 0

! Calcul de ValMax en local :  VALEUR MAX POUR L'ARRONDI

valmax = epzero

! ---> Face
do is = 1,itypf7

  ii = iconf7(is)
  do il = 1,3
    valmax = max( valmax,abs( xyznod(il,ii) ) )
  enddo

enddo

! ---> Centre de gravite de la face
valmax = max(valmax,abs(cdgx))
valmax = max(valmax,abs(cdgy))
valmax = max(valmax,abs(cdgz))

! ---> Centre de gravite de la cellule de depart
valmax = max(valmax,abs(celx))
valmax = max(valmax,abs(cely))
valmax = max(valmax,abs(celz))

! ---> Point de depart de la particule
valmax = max(valmax,abs(px))
valmax = max(valmax,abs(py))
valmax = max(valmax,abs(pz))

! ---> Point d'arrivee de la particule
valmax = max(valmax,abs(qx))
valmax = max(valmax,abs(qy))
valmax = max(valmax,abs(qz))

!===============================================================================
! 1. Position relative de P et Q : les 2 points sont-ils confondus ?
!===============================================================================

call coloca                                                       &
!==========
  ( valmax ,                                                      &
    px     , py     , pz     ,                                    &
    qx     , qy     , qz     ,                                    &
    ipos   )

!--> Si P et Q sont confondus, la particule est dans le meme cellule

if (ipos.eq.1) then
  indian = -1
  return
endif

!===============================================================================
! 2. Verification de l'orientation de la face
!===============================================================================

!   Ici on verifie dans quel sens sont lus les points
!   de la face par rapport a la cellule dans laquelle
!   se trouve la particule (au point P)
!   "bien" orientee : ISENSF = 1
!   "mal" orientee : ISENSF = -1

!--> Coordonnees du premier point support S(1) de la face

ii = iconf7(1)
cr1x = xyznod(1,ii)
cr1y = xyznod(2,ii)
cr1z = xyznod(3,ii)

!--> Coordonnees du deuxieme point support S(2) de la face

ii = iconf7(2)
cr2x = xyznod(1,ii)
cr2y = xyznod(2,ii)
cr2z = xyznod(3,ii)

!--> Orientation de PgS(1)S(2)

call coturn                                                       &
!==========
  ( valmax ,                                                      &
    px     , py    , pz    ,                                      &
    cdgx   , cdgy  , cdgz  ,                                      &
    cr1x   , cr1y  , cr1z  ,                                      &
    cr2x   , cr2y  , cr2z  ,                                      &
    isensf , ipturb        )

!--> Si l'orientation precedente a echouee, on en essaie une autre,
!    mais elle est dangereuse dans le cas des cellules convaves.

!    En periodicite, on force ISENSF a zero a cause d'un probleme
!    potentiel de reperage des particules dans le cas contraire
!    (FIXME: determiner la source du probleme)
if (iperio.eq.1) isensf = 0

if (isensf.eq.0) then

!--> Orientation de GgS(1)S(2)

  call coturn                                                     &
  !==========
  ( valmax ,                                                      &
    celx   , cely  , celz  ,                                      &
    cdgx   , cdgy  , cdgz  ,                                      &
    cr1x   , cr1y  , cr1z  ,                                      &
    cr2x   , cr2y  , cr2z  ,                                      &
    isensf , ipturb        )

endif

if (isensf.eq.0) then
#if DEBUG_OUESTU
  write(nfecra,1010)
#endif
  ierror = 1
  return
endif

!--------
! FORMATS
!--------

#if DEBUG_OUESTU
 1010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,&
'@     LA RECHERCHE DE L''ORIENTATION DE LA FACE A ECHOUEE    ',/,&
'@                                                            ',/,&
'@  La particule est perdue.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#endif

!===============================================================================
! 2. Test sur le 1er point support de la face
!===============================================================================

!-->orientation de PQgS(1)

call coturn                                                       &
!==========
  ( valmax ,                                                      &
    px     , py    , pz    ,                                      &
    qx     , qy    , qz    ,                                      &
    cdgx   , cdgy  , cdgz  ,                                      &
    cr1x   , cr1y  , cr1z  ,                                      &
    isign  , ipturb        )

isign = isign * isensf

if (isign.eq.0) then
#if DEBUG_OUESTU
  write(nfecra,1020)
#endif
  ierror = 1
  return
endif

!--------
! FORMATS
!--------

#if DEBUG_OUESTU
 1020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,&
'@     L''ORIENTATION DE PQgS(1) A ECHOUEE                    ',/,&
'@                                                            ',/,&
'@     Cause probable : les points PQgS(1) sont coplanaires.  ',/,&
'@                                                            ',/,&
'@  La particule est perdue.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#endif

!===============================================================================
! 3. Boucle sur tous les points supports S(i) de la face
!    On cherche dans quel est le triangle construit sur les sommets et
!    le CDG de la face, qui est perce par la droite PQ
!===============================================================================

do is = 2, itypf7

!--> Coordonnees des points supports

  ii = iconf7(is)
  cr1x = xyznod(1,ii)
  cr1y = xyznod(2,ii)
  cr1z = xyznod(3,ii)

!--> Orientation de PQgS(i) avec 2 =< i =< ITYPF7

  call coturn                                                     &
  !==========
    ( valmax ,                                                    &
      px     , py     , pz     ,                                  &
      qx     , qy     , qz     ,                                  &
      cdgx   , cdgy   , cdgz   ,                                  &
      cr1x   , cr1y   , cr1z   ,                                  &
      in     , ipturb        )

  in = in * isensf

  if (in.eq.0) then
#if DEBUG_OUESTU
    write(nfecra,1030) is
#endif
    ierror = 1
    return
  endif

!--> Si inversion de l'orientation de PQgS(i) alors orientation
!    de PQS(i-1)S(i)

  if (isign.eq.-in) then

!--> Le triangle gS(i-1)S(i) est t-il traverse par le rayon PQ ?
!    si oui, on repere le triangle en question

    if (isign.eq.1) ijklug = is

    isign = in

!--> Orientation de PQS(i-1)S(i)

    ii = iconf7(is-1)
    cr2x = xyznod(1,ii)
    cr2y = xyznod(2,ii)
    cr2z = xyznod(3,ii)

    call coturn                                                   &
    !==========
      ( valmax ,                                                  &
        px     , py     , pz     ,                                &
        qx     , qy     , qz     ,                                &
        cr2x   , cr2y   , cr2z   ,                                &
        cr1x   , cr1y   , cr1z   ,                                &
        it     , ipturb          )

    it = it * isensf

    if (it.eq.0) then
#if DEBUG_OUESTU
      write(nfecra,1040)
#endif
      ierror = 1
      return
    endif

    ittour = ittour + it

  endif

enddo

if (ittour.ne.-2 .and. ittour.ne.2 .and. ittour.ne.0) then
#if DEBUG_OUESTU
  write(nfecra,2010) ittour
#endif
  ierror = 1
  return
else if ((ittour.eq.-2 .or. ittour.eq.2) .and. ijklug.eq.0) then
#if DEBUG_OUESTU
  write(nfecra,2020)
#endif
  ierror = 1
  return
endif

!--> Si la droite PQ ne traverse pas la face, ou s'elle rentre dans
!    la cellule par cette face, on retourne dans LAGCEL

if (ittour.eq.0 .or. ittour.eq.-2) return

!--------
! FORMATS
!--------

#if DEBUG_OUESTU
 1030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,&
'@     L''ORIENTATION DE PQgS(i) A ECHOUEE, i = ',I2           ,/,&
'@                                                            ',/,&
'@     Cause probable : les points PQgS(i) sont coplanaires.  ',/,&
'@                                                            ',/,&
'@  La particule est perdue.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 1040 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,&
'@     L''ORIENTATION DE PQS(i-1)S(i) A ECHOUEE               ',/,&
'@                                                            ',/,&
'@     Cause Probable : points PQS(i-1)S(i) coplanaires.      ',/,&
'@                                                            ',/,&
'@  La particule est perdue.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,&
'@     L''INDICE DE PASSAGE DOIT ESTRE EGALE A -2, 0 ou 2     ',/,&
'@     IL EST CALCULE A : ',I2                                 ,/,&
'@                                                            ',/,&
'@  La particule est perdue.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,&
'@     LA DETERMINATION DU TRIANGLE DE PASSAGE                ',/,&
'@     DE LA PARTICULE A ECHOUEE                              ',/,&
'@                                                            ',/,&
'@  La particule est perdue.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#endif

!===============================================================================
! 4. Position relative des points d'arrive et de depart et de la face
!===============================================================================

!--> Dans la cas ou la droite PQ construite sur les points d'arrive Q
!    et de depart P de la particule, sort de la cellule par la face
!    courante (ie ITTOUR = 2), on veut savoir si P et Q se trouve
!    de part et d'autre de la face, ou du meme cote :

!    ATTENTION : test perturbe (IPTURB=1) : si Q est sur la face,
!                la detection fonctionnera (a moins que QgS(i-1)s(i)
!                ne soient coplanaires).

ipturb = 1

ii = iconf7(ijklug-1)
cr1x = xyznod(1,ii)
cr1y = xyznod(2,ii)
cr1z = xyznod(3,ii)

ii = iconf7(ijklug)
cr2x = xyznod(1,ii)
cr2y = xyznod(2,ii)
cr2z = xyznod(3,ii)

!--> Orientation de QgS(i-1)S(i)

    call coturn                                                   &
    !==========
      ( valmax ,                                                  &
        qx     , qy     , qz     ,                                &
        cdgx   , cdgy   , cdgz   ,                                &
        cr1x   , cr1y   , cr1z   ,                                &
        cr2x   , cr2y   , cr2z   ,                                &
        isign  , ipturb          )

if (isign.eq.0) then
#if DEBUG_OUESTU
  write(nfecra,3010)
#endif
  ierror = 1
  return
endif

!--> Position relative entre P Q et la face

indian = -isign * isensf

!--------
! FORMATS
!--------

#if DEBUG_OUESTU
 3010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : A L''EXECUTION DU MODULE LAGRANGIEN         ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@     ERREUR DANS LE REPERAGE D''UNE PARTICULE :             ',/,&
'@     LA POSITION RELATIVE DE Q PAR RAPPORT A P              ',/,&
'@     A LA FACE A ECHOUEE                                    ',/,&
'@                                                            ',/,&
'@     Cause Probable : points QgS(i-1)S(i) coplanaires.      ',/,&
'@                                                            ',/,&
'@  La particule est perdue.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#endif

!----
! FIN
!----

end subroutine
