!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2014 EDF S.A.
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

subroutine numvec &
!================

 ( ncelet , ncel   , nfac   , nfabor ,                            &
   lregis , irveci , irvecb ,                                     &
   ifacel , ifabor ,                                              &
   inumfi , inumfb , iworkf , ismbs  )

!===============================================================================
! FONCTION :
! ---------

! CALCUL D UNE TABLE DE RENUMEROTATION DES FACES INTERNES ET DE BORD

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac  /nfabor    ! e  ! <-- ! nombre total de faces internes/de brd          !
! lregis           ! e  ! --> ! longueur de registre vectoriel                 !
! irveci           ! e  ! --> ! indicateur vectorisation face intern           !
! irvecb           ! e  ! --> ! indicateur vectorisation face bord             !
! ifacel           ! te ! <-- ! no des elts voisins d'une face intern          !
! ifabor           ! te ! <-- ! no de l'elt voisin d'une face de bord          !
! nfabor  )        !    !     !                                                !
! inumfi(nfac)     ! te ! --- ! table de renum des faces internes              !
! inumfb(nfabor    ! te ! --- ! table de renum des faces de bord               !
! iworkf(*         ! te ! --- ! tab de trav de dim max(nfac,nfabor)            !
! ismbs (ncelet    ! te ! --- ! tab de trav pour assemblage scalaire           !
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
use entsor
use parall

!===============================================================================

implicit none

! Arguments

integer          ncelet, ncel, nfac, nfabor
integer          lregis, irveci, irvecb
integer          ifacel(2,nfac),ifabor(nfabor)
integer          inumfi(nfac), inumfb(nfabor)
integer          iworkf(*), ismbs(ncelet)

! Local variables

integer          irelii, nregii, irelib, nregib
integer          iloop, imodav, iregip, iregic, jregic
integer          ilast, inext
integer          iiech
integer          ii, jj, iok, ifac, ifac1, nfanp1
integer          ibloc, iel, ireg, ilig, nfamax, itmp


!===============================================================================

! Initialize variables to avoid compiler warnings

iiech = 0
irelii = 0
irelib = 0

!===============================================================================
! 1. SORTIE IMMEDIATE SI L'ON N'A PAS UN CALCULATEUR VECTORIEL
!===============================================================================

if (lregis .eq. 1) then
  ivecti = 0
  ivectb = 0
  irveci = 0
  irvecb = 0
  return
endif

!===============================================================================
! 2. INITIALISATIONS COMMUNES
!===============================================================================

! --- Numerotation
do ifac = 1, nfac
  inumfi(ifac) = ifac
enddo
do ifac = 1, nfabor
  inumfb(ifac) = ifac
enddo


!===============================================================================
! 3. RANGEMENT DES FACES INTERNES
!     (pour le raisonnement, on place le reliquat a la fin)
!===============================================================================

! --- Si l'utilisateur a indique IVECTI = 0, il ne souhaite pas
!       vectoriser, sinon, IVECTI a ete initialise a -1

if(ivecti.eq.0) goto 400

! --- Indicateur de vectorisation possible = 1
ivecti = 0

! --- Determination du reliquat et du nbre de registres complets
irelii = mod(nfac,lregis)
nregii = nfac/lregis

! --- Compteur de boucles
iloop = 0

! --- Debut de la boucle externe
  100 continue
iloop = iloop + 1

! --- IMODAV = 1 si on a echange un element et un element precedent
!                 dans la table INUMFI
imodav = 0

! --- Registre precedent

iregic = 0

! --- On parcourt les elements de INUMFI

do jj = 1, nfac

! --- Registre courant et position dans le registre
  iregip = iregic
  iregic = (jj-1)/lregis+1
  jregic = mod(jj-1,lregis)+1

! --- On teste entre ILAST, debut du registre,
!       et la position courante

!     En prenant le cas le plus penalisant entre
!       reliquat au debut et reliquat a la fin : reliquat au debut :

  if(iregic.eq.1) then
    ilast = 1
  elseif (jregic.le.irelii) then
    ilast = (iregic-2)*lregis+irelii+1
  else
    ilast = (iregic-1)*lregis+1
  endif

!     Avec reliquat a la fin :

!        ILAST = (IREGIC-1)*LREGIS+1

! --- On echange a partir de INEXT, debut du registre suivant

!     En prenant le cas le plus penalisant entre reliquat au debut et
!       reliquat a la fin : reliquat au debut :

  if ((iregic.eq.nregii.and.jregic.gt.irelii).or.                 &
      (iregic.eq.nregii+1)                     ) then
    inext = 1
  elseif (jregic.gt.irelii) then
    inext = iregic*lregis+irelii+1
  else
    inext = iregic*lregis+1
  endif

!     Sinon, reliquat a la fin :

!        IF ((IREGIC.EQ.NREGII.AND.IRELII.EQ.0) .OR.
!     &       IREGIC.EQ.NREGII+1                   ) THEN
!          INEXT = 1
!        ELSE
!          INEXT = IREGIC*LREGIS+1
!        ENDIF

  if(iregic.ne.iregip) iiech = inext-1

! --- Compteur pour ne pas echanger avec tous les elements de INUMFI
!       plus de n fois
  ibloc = 0


! --- Test avec tous les elements precedents depuis ILAST
!      IIECH indique avec quel element de INUMFI on echange
!      IMODAV indique qu'on modifie un element deja vu
!      IBLOC indique qu'on a vu tous les elements et qu'il faut
!              melanger (il n'y a pas de solution)

 200    continue

  ifac = inumfi(jj)
  do ii = ilast, jj-1

    if ( (ifacel(1,inumfi(ii)).eq.ifacel(1,ifac)).or.             &
         (ifacel(2,inumfi(ii)).eq.ifacel(2,ifac))   ) then

      iiech           = iiech+1

      if(iiech.gt.nfac) then
        iiech = 1
        ibloc = ibloc + 1
      endif
      if (iiech.lt.jj) imodav = 1
      if (ibloc.ge.2) then
        ibloc = -1
        goto 450
      endif

      itmp             = inumfi(iiech )
      inumfi(iiech )   = inumfi(jj    )
      inumfi(jj    )   = itmp

      goto 200
    endif
  enddo
enddo

! --- Si on n'a pas touche aux elements precedant le courant,
!       ca a marche
if(imodav.eq.0) then
  ivecti = 1
  goto 400
endif

! --- On melange s'il n'y a pas de solution ou si on a boucle 10 fois
 450  continue
if (iloop.le.100.and.(mod(iloop,10).eq.0.or.ibloc.eq.-1)) then
  do ii = 1, (nfac-4)/2, 2
    jj = nfac-ii+1
    itmp             = inumfi(ii    )
    inumfi(ii    )   = inumfi(jj    )
    inumfi(jj    )   = itmp
  enddo
endif


! --- Et on recommence
if(iloop.le.100) goto 100

 400  continue

!===============================================================================
! 4. RANGEMENT DES FACES DE BORD
!===============================================================================

! --- Si l'utilisateur a indique IVECTB = 0, il ne souhaite pas
!       vectoriser, sinon, IVECTB a ete initialise a -1

if(ivectb.eq.0) goto 900

! --- Indicateur de vectorisation possible = 1
ivectb = 0

! --- Determination du reliquat et du nbre de registres complets
irelib = mod(nfabor,lregis)
nregib = nfabor/lregis

! --- Nombre max de faces de bord
!      Si > NREGIB : il n'y a pas de solution
do iel = 1, ncel
  ismbs(iel) = 0
enddo
do ifac = 1, nfabor
  ii = ifabor(ifac)
  ismbs(ii) = ismbs(ii) + 1
enddo
nfamax = 0
nfanp1 = 0
do iel = 1, ncel
  nfamax = max(nfamax,ismbs(iel))
  if(ismbs(iel).eq.nregib+1) nfanp1 = nfanp1 + 1
enddo
if ( nfamax.gt.nregib+1.or.                                       &
    (nfamax.eq.nregib+1.and.nfanp1.gt.irelib)) then
  goto 900
endif

! --- On classe par nombre de faces de bord du voisin decroissant
!          et numero de voisin decroissant

do ifac = 1, nfabor
  iel = ifabor(ifac)
  ifabor(ifac) = iel + ncel*ismbs(iel)
enddo
call ordita(nfabor,ifabor,iworkf)
!==========
do ifac = 1, nfabor
  iel = mod(ifabor(ifac)-1,ncel)+1
  ifabor(ifac) = ifabor(ifac) - ncel*ismbs(iel)
enddo

! --- On distribue les faces dans les registres
do ifac = 1, nfabor
  if(ifac.le.irelib*(nregib+1)) then
    ireg=mod(ifac-1,nregib+1)+1
    ilig=(ifac-1)/(nregib+1)+1
    ii = (ireg-1)*lregis+ilig
  else
    ifac1 = ifac-irelib*(nregib+1)
    ireg=mod(ifac1-1,nregib)+1
    ilig=(ifac1-1)/nregib+1+irelib
    ii = (ireg-1)*lregis+ilig
  endif
  inumfb(ii)=iworkf(ifac)
enddo
ivectb=1

  900 continue

!===============================================================================
! 5. VERIFICATIONS
!===============================================================================

! -----> Verif que toutes les faces se retrouvent une et une seule fois
!       dans INUMFB et INUMFI

if(ivecti.eq.1) then

  call ordita(nfac  ,inumfi,iworkf)
  !==========

  iok = 0
  do ii = 1, nfac
    if(inumfi(iworkf(ii)).ne.nfac-ii+1) iok = iok - 1
  enddo
  if (iok.ne.0) then
    write(nfecra,1100)iok,1101
    ivecti = 0
  endif

endif

if(ivectb.eq.1) then

  call ordita(nfabor,inumfb,iworkf)
  !==========

  iok = 0
  do ii = 1, nfabor
    if(inumfb(iworkf(ii)).ne.nfabor-ii+1) iok = iok - 1
  enddo
  if (iok.ne.0) then
    write(nfecra,1200)iok,1201
    ivectb = 0
  endif

endif


! -----> Test classique en balayant les faces precedentes

! --- Faces internes

if(ivecti.eq.1) then

  iok = 0
  do jj = 1, nfac

! --- Registre courant et position dans le registre
    iregic = (jj-1)/lregis+1
    jregic = mod(jj-1,lregis)+1

! --- On teste entre ILAST, debut du registre,
!       et la position courante

!     En prenant le cas le plus penalisant entre
!       reliquat au debut et reliquat a la fin : reliquat au debut :

    if(iregic.eq.1) then
      ilast = 1
    elseif (jregic.le.irelii) then
      ilast = (iregic-2)*lregis+irelii+1
    else
      ilast = (iregic-1)*lregis+1
    endif

!     Avec reliquat a la fin :

!          ILAST = (IREGIC-1)*LREGIS+1


! --- Test avec tous les elements precedents depuis ILAST

    do ii = ilast, jj-1
      ifac = inumfi(jj)
      if((ifacel(1,inumfi(ii)).eq.ifacel(1,ifac)).or.             &
         (ifacel(2,inumfi(ii)).eq.ifacel(2,ifac))    )then
        iok = iok - 1
      endif
    enddo
  enddo

  if(iok.ne.0) then
    write(nfecra,1100)iok,1102
    ivecti = 0
  endif

endif


! --- Faces de bord

if(ivectb.eq.1) then

  iok = 0
  do jj = 1, nfabor

! --- Registre courant et position dans le registre
    iregic = (jj-1)/lregis+1
    jregic = mod(jj-1,lregis)+1

! --- On teste entre ILAST, debut du registre,
!       et la position courante

!     En prenant le cas le plus penalisant entre
!       reliquat au debut et reliquat a la fin : reliquat au debut :

    if(iregic.eq.1) then
      ilast = 1
    elseif (jregic.le.irelib) then
      ilast = (iregic-2)*lregis+irelib+1
    else
      ilast = (iregic-1)*lregis+1
    endif

!     Avec reliquat a la fin :

!          ILAST = (IREGIC-1)*LREGIS+1


! --- Test avec tous les elements precedents depuis ILAST

    do ii = ilast, jj-1
      ifac = inumfb(jj)
      if (ifabor(inumfb(ii)).eq.ifabor(ifac))then
        iok = iok - 1
      endif
    enddo
  enddo

  if(iok.ne.0) then
    write(nfecra,1200)iok,1202
    ivectb = 0
  endif

endif

!===============================================================================
! 6. INDICATEURS
!===============================================================================

irveci = ivecti
irvecb = ivectb

!===============================================================================
! 7. FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)

 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RENUMEROTATION DES FACES                    ',/,&
'@    =========                                               ',/,&
'@     PROBLEME DE RENUMEROTATION DES FACES INTERNES          ',/,&
'@     ',I10   ,' OCCURRENCES DE L''EVENEMENT ',I5             ,/,&
'@                                                            ',/,&
'@  Le calcul peut etre execute.                              ',/,&
'@  La vectorisation restera non forcee                       ',/,&
'@                                                            ',/,&
'@  Signaler ce message a l''assistance.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : RENUMEROTATION DES FACES                    ',/,&
'@    =========                                               ',/,&
'@     PROBLEME DE RENUMEROTATION DES FACES DE BORD           ',/,&
'@     ',I10   ,' OCCURRENCES DE L''EVENEMENT ',I5             ,/,&
'@                                                            ',/,&
'@  Le calcul peut etre execute.                              ',/,&
'@  La vectorisation restera non forcee                       ',/,&
'@                                                            ',/,&
'@  Signaler ce message a l''assistance.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: FACES RENUMBERING                              ',/,&
'@    ========                                                ',/,&
'@     PROBLEM IN THE INTERIOR FACES RENUMBERING              ',/,&
'@     ',I10   ,' OCCURRENCES OF THE EVENT ',I5                ,/,&
'@                                                            ',/,&
'@  The calculation will be run.                              ',/,&
'@  The vectorization will not be forced.                     ',/,&
'@                                                            ',/,&
'@  Please contact the support about this message.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: FACES RENUMBERING                              ',/,&
'@    ========                                                ',/,&
'@     PROBLEM IN THE BOUNDARY FACES RENUMBERING              ',/,&
'@     ',I10   ,' OCCURRENCES OF THE EVENT ',I5                ,/,&
'@                                                            ',/,&
'@  The calculation will be run.                              ',/,&
'@  The vectorization will not be forced.                     ',/,&
'@                                                            ',/,&
'@  Please contact the support about this message.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

return
end subroutine
