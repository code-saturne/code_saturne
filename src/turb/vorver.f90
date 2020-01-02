!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2020 EDF S.A.
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

subroutine vorver &
 ( nfabor , iappel )

!===============================================================================
!  FONCTION  :
!  ---------

! VERIFICATION DES PARAMETRES DE CALCUL DE LA METHODE DES VORTEX
!   APRES INTERVENTION UTILISATEUR
!   (COMMONS)
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! iappel           ! e  ! <-- ! indique les donnes a verifier                  !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use cstnum
use entsor
use optcal
use vorinc
use parall

!===============================================================================

implicit none

! Arguments

integer          nfabor , iappel

! Local variables

integer          iok, nmax , ii, jj
double precision norm1, norm2, crosp

!===============================================================================

iok = 0

!===============================================================================
! 1. IAPPEL = 1
!===============================================================================

if(iappel.eq.1) then

! --- Nombre d'entree

  if (nnent.gt.nentmx.or.nnent.le.0) then
    write(nfecra,1000) nentmx, nnent
    iok = iok + 1
  endif

! --- Nombre de vortex

  do ii = 1, nnent
    if(nvort(ii).le.0) then
      write(nfecra,2000) ii, nvort(ii)
      iok = iok + 1
    endif
  enddo

! --- Reperage des entrees

  nmax = 0
  do ii = 1, nfabor
    if(irepvo(ii).lt.0.or.irepvo(ii).gt.nnent) then
      write(nfecra,2200) nnent,irepvo(ii)
      iok = iok + 1
    else
      nmax = max(nmax,irepvo(ii))
    endif
  enddo
  if(irangp.ge.0) then
    call parcmx(nmax)
  endif
  if(nmax.eq.0) then
    write(nfecra,2200) nnent,nmax
    iok = iok + 1
  endif

!===============================================================================
! 2. IAPPEL = 2
!===============================================================================

elseif(iappel.eq.2) then

! --- Dimensions des Tableaux

  if (nnent.gt.nentmx.or.nnent.le.0) then
    write(nfecra,1000) nentmx, nnent
    iok = iok + 1
  endif

  nmax = 0

  do ii = 1, nnent
    if(icas(ii).eq.4) then
      ndat(ii) = 1
    endif
  enddo

  do ii = 1, nnent
    if(ndat(ii).le.0) then
      write(nfecra,1100) ndatmx, ndat(ii), ii
      iok = iok + 1
    else
      nmax = max(nmax,ndat(ii))
    endif
    if(nmax.gt.ndatmx) then
      write(nfecra,1100) ndatmx, nmax, ii
      iok = iok + 1
    endif
  enddo

! --- Nombre de vortex

  do ii = 1, nnent
    if(nvort(ii).le.0) then
      write(nfecra,2000) ii, nvort(ii)
      iok = iok + 1
    endif
  enddo

! --- Suite de calcul

  if(isuivo.ne.0.and.isuivo.ne.1) then
    write(nfecra,2100) isuivo
    iok = iok + 1
  endif

! --- Reperage des entrees

  nmax = 0
  do ii = 1, nfabor
    if(irepvo(ii).lt.0.or.irepvo(ii).gt.nnent) then
      write(nfecra,2200) nnent,irepvo(ii)
      iok = iok + 1
    else
      nmax = max(nmax,irepvo(ii))
    endif
  enddo
  if(irangp.ge.0) then
    call parcmx(nmax)
  endif
  if(nmax.eq.0) then
    write(nfecra,2200) nnent,nmax
    iok = iok + 1
  endif

! --- Definition du cas

  do ii = 1, nnent
    if(icas(ii).ne.1.and.icas(ii).ne.2.and.                       &
         icas(ii).ne.3.and.icas(ii).ne.4) then
      write(nfecra,3000) ii, icas(ii)
      iok = iok + 1
    endif
  enddo

! --- Repere locale

  do ii = 1, nnent
    if (icas(ii).ne.4) then
      norm1 = dir1(1,ii)**2.d0+dir1(2,ii)**2.d0+dir1(3,ii)**2.d0
      norm2 = dir1(1,ii)**2.d0+dir1(2,ii)**2.d0+dir1(3,ii)**2.d0
      crosp = (dir1(2,ii)*dir2(3,ii)-dir1(3,ii)*dir2(2,ii))**2.d0 &
            + (dir1(1,ii)*dir2(3,ii)-dir1(3,ii)*dir2(1,ii))**2.d0 &
            + (dir1(1,ii)*dir2(2,ii)-dir1(2,ii)*dir2(1,ii))**2.d0
      if(abs(norm1*norm2*crosp).lt.epzero) then
        iok = iok + 1
        write(nfecra,3100)ii,sqrt(norm1),sqrt(norm2),sqrt(crosp)
      endif
      if(abs(norm1-1.d0).gt.epzero.or.                            &
         abs(norm2-1.d0).gt.epzero.and.                           &
         abs(norm1*norm2*crosp).gt.epzero) then
        write(nfecra,3150) ii, sqrt(norm1),sqrt(norm2)
      endif
    endif
  enddo

! --- Conditions aux limites

  do ii = 11, nnent
    if(icas(ii).eq.1) then
      do jj = 1, 4
        if(iclvor(jj,ii).ne.1                                     &
             .and.iclvor(jj,ii).ne.2                              &
             .and.iclvor(jj,ii).ne.3) then
          write(nfecra,3200) iclvor(jj,ii), jj, ii
          iok = iok + 1
        endif
      enddo
      if((iclvor(1,ii).eq.3.and.iclvor(3,ii).ne.3).or.            &
           (iclvor(1,ii).ne.3.and.iclvor(3,ii).eq.3).or.          &
           (iclvor(2,ii).eq.3.and.iclvor(4,ii).ne.3).or.          &
           (iclvor(2,ii).ne.3.and.iclvor(4,ii).eq.3)) then
        write(nfecra,3250) iclvor(jj,ii), jj, ii
        iok = iok + 1
      endif
    endif
  enddo

! --- Dimension carateristiques

  do ii = 1, nnent
    if(icas(ii).eq.1.and.(llz(ii).le.0.d0.or.                     &
         lly(ii).le.0.d0)) then
      write(nfecra,3400) llz(ii), lly(ii), ii
      iok = iok + 1
    endif
    if(icas(ii).eq.2.and.lld(ii).le.0.d0) then
      write(nfecra,3500) lld(ii), ii
      iok = iok + 1
    endif
  enddo

! --- Parametres pour la duree de vie, la taille, et le deplacement des vortex

  do ii = 1, nnent
    if(itlivo(ii).ne.1.and.itlivo(ii).ne.2) then
      write(nfecra,4000) itlivo(ii), ii
      iok = iok + 1
    endif
  enddo

  do ii = 1, nnent
    if(itlivo(ii).eq.1.and.tlimvo(ii).le.0.d0) then
      write(nfecra,4100) tlimvo(ii), ii
      iok = iok + 1
    endif
  enddo

  do ii = 1, nnent
    if(isgmvo(ii).ne.1.and.isgmvo(ii).ne.2.and.                   &
         isgmvo(ii).ne.3) then
      write(nfecra,4200) isgmvo(ii), ii
      iok = iok + 1
    endif
  enddo

  do ii = 1, nnent
    if(isgmvo(ii).eq.1.and.xsgmvo(ii).le.0.d0) then
      write(nfecra,4300) xsgmvo(ii), ii
      iok = iok + 1
    endif
  enddo

  do ii = 1, nnent
    if(idepvo(ii).ne.1.and.idepvo(ii).ne.2.and.                   &
         idepvo(ii).ne.3) then
      write(nfecra,4400) idepvo(ii), ii
      iok = iok + 1
    endif
  enddo

  do ii = 1, nnent
    if(idepvo(ii).eq.1.and.ud(ii).le.0.d0) then
      write(nfecra,4500) ud(ii), ii
!     IOK = IOK + 1
    endif
  enddo

! --- Donnees utilisateur

  do ii = 1, nnent
    if(icas(ii).eq.4.and.udebit(ii).le.0.d0) then
      write(nfecra,5000) ii, udebit(ii)
      iok = iok + 1
    endif
  enddo

  do ii = 1, nnent
    if(icas(ii).eq.4.and.kdebit(ii).le.0.d0) then
      write(nfecra,5100) ii, kdebit(ii)
      iok = iok + 1
    endif
  enddo

  do ii = 1, nnent
    if(icas(ii).eq.4.and.edebit(ii).le.0.d0) then
      write(nfecra,5200) ii, edebit(ii)
      iok = iok + 1
    endif
  enddo

endif
!===============================================================================
! 3. FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    NNENT DOIT ETRE UN ENTIER POSITIF INFERIEUR A ',I10      ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort ou augmenter la taille de NENTMX dans     ',/,&
'@  vortex.h                                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    NDAT DOIT ETRE UN ENTIER POSITIF INFERIEUR A ',I10       ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort ou agmenter la taille de NDAT dans        ',/,&
'@  vortex.h                                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE NOMBRE DE VORTEX NVORT A L''ENTREE',I10               ,/,&
'@    DOIT ETRE UN ENTIER POSITIF                             ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ISUIVO DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    IREPVO INDIQUE LE NUMERO DE L''ENTREE                   ',/,&
'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A ',I10  ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ICAS INDIQUE LE CAS TRAITE                              ',/,&
'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 4     ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES VECTEURS DIR1 ET DIR2 NE PEUVENT ETRE NI NULS       ',/,&
'@    NI COLINEAIRES                                          ',/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@    ||DIR1|| VAUT ',E14.5                                    ,/,&
'@    ||DIR2|| VAUT ',E14.5                                    ,/,&
'@    ||DIR1 ^ DIR2 || VAUT ',E14.5                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3150 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES VECTEURS DIR1 ET DIR2 NE SONT PAS NORMES            ',/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@    ||DIR1|| VAUT ',E14.5                                    ,/,&
'@    ||DIR2|| VAUT ',E14.5                                    ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ICLVOR INDIQUE LE TYPE DE CONDITION AUX LIMITES         ',/,&
'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 4     ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@    SUR LE COTE ',I10                                        ,/,&
'@    DE L''ENTREE ',I10                                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3250 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VERFIFIER LA COHERENCE DES CONDITION AUX LIMITES        ',/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@    NOMBRE IMPAIR DE CONDITIONS PERIODIQUES                 ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3400 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LLY ET LLZ SONT LES DIMENSIONS DE L''ENTREE             ',/,&
'@    ELLES VALENT ICI ',E14.5, ' ET ',E14.5                   ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3500 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LLD EST LE DIAMETRE DE LA CONDUITE                      ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ITLIVO INDIQUE LE TYPE DE MODELE UTILISE POUR CALCULER  ',/,&
'@    LA DUREE DE VIE DES VORTEX                              ',/,&
'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 2     ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    TLIMVO EST LA DUREE DE VIE MAXIMALE DES VORTEX          ',/,&
'@    ELLE VAUT ICI ',E14.5                                    ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ISGMVO INDIQUE LE TYPE DE MODELE UTILISE POUR CALCULER  ',/,&
'@    LA TAILLE DES VORTEX                                    ',/,&
'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 3     ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    XSGMVO EST LA TAILLE MAXIMALE DES VORTEX                ',/,&
'@    ELLE VAUT ICI ',E14.5                                    ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4400 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    IDEPVO INDIQUE LE TYPE DE MODELE UTILISE POUR CALCULER  ',/,&
'@    LE DEPLACEMENT DES VORTEX                               ',/,&
'@    IL S AGIT D UN ENTIER POSITIF INFERIEUR OU EGAL A 3     ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4500 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@    UD EST LA VITESSE MAXIMALE DE DEPLCAMENT DES VORTEX     ',/,&
'@    ELLE VAUT ICI ',E14.5                                    ,/,&
'@    A L''ENTREE ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LA VITESSE DEBITANTE UDEBIT A L''ENTREE ',I10            ,/,&
'@    VAUT ',E14.5                                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    L ENERGIE CINETIQUE KDEBIT A L''ENTREE ',I10             ,/,&
'@    VAUT ',E14.5                                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LA DISSIPATION EDEBIT A L''ENTREE ',I10                  ,/,&
'@    VAUT ',E14.5                                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    NNENT MUST BE A POSITIVE INTEGER LOWER THAN',I10         ,/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort or increase the size of NENTMX in vortex.h  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    NDAT MUST BE A POSITIVE INTEGER LOWER THAN ',I10         ,/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort or increase the size of NDAT in vortex.h    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE NUMBER OF VORTICES NVORT AT THE INLET ',I10          ,/,&
'@    MUST BE A POSITIVE INTEGER                              ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    ISUIVO MUST BE AN INTEGER EQUAL TO 0 OR 1               ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    IREPVO GIVES THE INLET NUMBER                           ',/,&
'@    IT IS A POSITIVE INTEGER LOWER THAN ',I10                ,/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    ICAS GIVES THE CHOSEN CASE                              ',/,&
'@    IT IS A POSITIVE INTEGER LOWER THAN 4                   ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE VECTORS DIR1 AND DIR2 CANNOT BE NULL NOR COLINEAR   ',/,&
'@    AT THE INLET ',I10                                       ,/,&
'@    ||DIR1|| IS ',E14.5                                      ,/,&
'@    ||DIR2|| IS ',E14.5                                      ,/,&
'@    ||DIR1 ^ DIR2 || IS ',E14.5                              ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3150 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE VECTORS DIR1 AND DIR2 ARE NOT NORMED                ',/,&
'@    AT THE INLET ',I10                                       ,/,&
'@    ||DIR1|| IS ',E14.5                                      ,/,&
'@    ||DIR2|| IS ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  The calculation will be run.                              ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    ICLVOR GIVES THE TYPE OF THE BOUNDARY CONDITIONS        ',/,&
'@    IT IS A POSITIVE INTEGER LOWER THAN 4                   ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@    ON THE SIDE ',I10                                        ,/,&
'@    OF THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3250 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    VERIFY THE COHERENCY OF THE BOUNDARY CONDITIONS         ',/,&
'@    AT THE INLET ',I10                                       ,/,&
'@    ODD NUMBER OF PERIODIC CONDITIONS                       ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3400 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    LLY AND LLZ ARE THE INLET DIMENSIONS                    ',/,&
'@    THEIR VALUES ARE ',E14.5, ' AND ',E14.5                  ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3500 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    LLD IS THE PIPE DIAMETER                                ',/,&
'@    ITS VALUE IS ',E14.5                                     ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    ITLIVO GIVES THE MODEL TYPE USED TO COMPUTE THE         ',/,&
'@    VORTICES LIFE TIME                                      ',/,&
'@    IT IS A POSITIVE INTEGER LOWER THAN OR EQUAL TO 2       ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    TLIMVO IS THE MAXIMUM LIFE TIME OF THE VORTICES         ',/,&
'@    ITS VALUE IS ',E14.5                                     ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    ISGMVO GIVES THE MODEL TYPE USED TO COMPUTE THE         ',/,&
'@    VORTICES SIZE                                           ',/,&
'@    IT IS A POSITIVE INTEGER LOWER THAN OR EQUAL TO 3       ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    XSGMVO IS THE MAXIMUM SIZE OF THE VORTICES              ',/,&
'@    ITS VALUE IS ',E14.5                                     ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4400 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    IDEPVO GIVES THE MODEL TYPE USED TO COMPUTE THE         ',/,&
'@    VORTICES DISPLACEMENT                                   ',/,&
'@    IT IS A POSITIVE INTEGER LOWER THAN OR EQUAL TO 3       ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4500 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING:       IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    UD IS THE MAXIMUM DISPLACEMENT VELOCITY OF THE VORTICES ',/,&
'@    ITS VALUE IS ',E14.5                                     ,/,&
'@    AT THE INLET ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will be run.                              ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE INLET VELOCITY UDEBIT AT THE INLET ',I10             ,/,&
'@    IS ',E14.5                                               ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE INLET KINETIC ENERGY KDEBIT AT THE INLET ',I10       ,/,&
'@    IS ',E14.5                                               ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE INLET DISSIPATION EDEBIT AT THE INLET ',I10          ,/,&
'@    IS ',E14.5                                               ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#endif

!===============================================================================
! 3. SORTIE ET IMPRESSIONS FINALES
!===============================================================================

if(iok.gt.0) then
  write(nfecra,9999)iok
  call csexit (1)
else
  write(nfecra,9998)
endif

#if defined(_CS_LANG_FR)

 9998 format(                                                           &
'                                                             ',/,&
' Pas d''erreur detectee lors de la verification des donnees  ',/,&
'    pour la methode des vortex (usvort).                     ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LES PARAMETRES DE CALCUL SONT INCOHERENTS OU INCOMPLETS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@  Verifier usvort.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 9998 format(                                                           &
'                                                             ',/,&
' No error detected during the verification of the parameters ',/,&
'    for the vortex method (usvort).                          ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE CALCULATION PARAMETERS ARE INCOHERENT OR INCOMPLETE ',/,&
'@                                                            ',/,&
'@  The calculation will not be run (',I10,' errors).         ',/,&
'@                                                            ',/,&
'@  Refer to the previous warnings for further information.   ',/,&
'@  Verify usvort.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

return
end subroutine
