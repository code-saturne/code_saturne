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

subroutine strini &
!================

 ( idbia0 , idbra0 , ifinia , ifinra ,                            &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   dt     , rdevel , rtuser ,                                     &
   ra     )

!===============================================================================
! FONCTION :
! ----------

! INITILISATION DES DONNEES DES STRUCTURES MOBILES EN ALE

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ifinia           ! i  ! --> ! number of first free position in ia (at exit)  !
! ifinra           ! i  ! --> ! number of first free position in ra (at exit)  !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ncelbr           ! i  ! <-- ! number of cells with faces on boundary         !
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! nphas            ! i  ! <-- ! number of phases                               !
! nideve, nrdeve   ! i  ! <-- ! sizes of idevel and rdevel arrays              !
! nituse, nrtuse   ! i  ! <-- ! sizes of ituser and rtuser arrays              !
! ifacel(2, nfac)  ! ia ! <-- ! interior faces -> cells connectivity           !
! ifabor(nfabor)   ! ia ! <-- ! boundary faces -> cells connectivity           !
! ifmfbr(nfabor)   ! ia ! <-- ! boundary face family numbers                   !
! ifmcel(ncelet)   ! ia ! <-- ! cell family numbers                            !
! iprfml           ! ia ! <-- ! property numbers per family                    !
!  (nfml, nprfml)  !    !     !                                                !
! ipnfac           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfac+1)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! <-- ! connectivite faces internes/noeuds             !
!   (lndfac)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! <-- ! position du premier noeud de chaque            !
!   (nfabor+1)     !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! <-- ! connectivite faces de bord/noeuds              !
!   (lndfbr)       !    !     !  (optionnel)                                   !
! idevel(nideve)   ! ia ! <-> ! integer work array for temporary development   !
! ituser(nituse)   ! ia ! <-> ! user-reserved integer work array               !
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! ra ! <-- ! cell centers                                   !
!  (ndim, ncelet)  !    !     !                                                !
! surfac           ! ra ! <-- ! interior faces surface vectors                 !
!  (ndim, nfac)    !    !     !                                                !
! surfbo           ! ra ! <-- ! boundary faces surface vectors                 !
!  (ndim, nfabor)  !    !     !                                                !
! cdgfac           ! ra ! <-- ! interior faces centers of gravity              !
!  (ndim, nfac)    !    !     !                                                !
! cdgfbo           ! ra ! <-- ! boundary faces centers of gravity              !
!  (ndim, nfabor)  !    !     !                                                !
! xyznod           ! ra ! <-- ! vertex coordinates (optional)                  !
!  (ndim, nnod)    !    !     !                                                !
! volume(ncelet)   ! ra ! <-- ! cell volumes                                   !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rdevel(nrdeve)   ! ra ! <-> ! real work array for temporary development      !
! rtuser(nrtuse)   ! ra ! <-> ! user-reserved real work array                  !
! ra(*)            ! ra ! --- ! main real work array                           !
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
use optcal
use cstphy
use entsor
use pointe
use albase
use alstru
use alaste
use ihmpre
use parall
use period

!===============================================================================

implicit none

! Arguments

integer          idbia0 , idbra0 , ifinia , ifinra
integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr , ncelbr
integer          nideve , nrdeve , nituse , nrtuse

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)
integer          idevel(nideve), ituser(nituse)
integer          ia(*)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac), surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac), cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod), volume(ncelet)
double precision dt(ncelet)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! Local variables

integer          idebia, idebra, itrav
integer          ifac  , istr, icompt, ii
integer          mbstru, mbaste
integer          maxelt, ifnia2, ils

integer          jj, inod
integer          ilstfa, indast, iidflo, iidnlo
integer          nbpdt, nbssit, ihi, chro

double precision delta, pdt, tt

!===============================================================================


!===============================================================================
! 1. INITIALISATION
!===============================================================================

idebia = idbia0
idebra = idbra0

do istr = 1, nstrmx
  dtstr(istr) = dt(1)
enddo

do istr = 1, nastmx
  do ii = 1, 3
    asddlf(ii,istr) = 1
    asddlc(ii,istr) = 1
  enddo
enddo

!     NBSTRU et NBASTE valent -999 si le calcul n'est pas un calcul
!       suite ou s'il est une suite d'un calcul sans ALE
mbstru = nbstru
mbaste = nbaste

!     En ALE on met IHISTR et ISYNCP a 1 par defaut
!       (remis a zero sans structure)
ihistr = 1
isyncp = 1

!===============================================================================
! 2.  RESERVATION DU TABLEAU IDFSTR
!===============================================================================

iidfst = idebia
ifinia = iidfst + nfabor

maxelt = max(ncelet,nfac,nfabor)
ils    = ifinia
ifnia2 = ils + maxelt

CALL IASIZE('STRINI',IFNIA2)
!==========

do ifac = 1, nfabor
  ia(iidfst+ifac-1) = 0
enddo


!===============================================================================
! 2.  REMPLISSAGE DE IDFSTR PAR L'UTILISATEUR
!===============================================================================

! 2.1 STRUCTURES INTERNES :
! -----------------------

if (iihmpr.eq.1) then

  call uistr1 &
  !==========
( nfabor,                  &
  ia(iidfst),              &
  aexxst, bexxst, cfopre,  &
  ihistr,                  &
  xstp, xstreq, xpstr )

endif

call usstr1                                                       &
!==========
 ( ifnia2 , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(iidfst),                &
   idevel , ituser , ia     ,                                     &
   aexxst , bexxst , cfopre ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   xstp   , xpstr  , xstreq ,                                     &
   rdevel , rtuser , ra     )

! 2.2 STRUCTURES EXTERNES : COUPLAGE CODE_SATURNE / CODE_ASTER
! -----------------------

call usaste                                                       &
!==========
 ( ifnia2 , idebra ,                                              &
   ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr , ncelbr ,                            &
   nideve , nrdeve , nituse , nrtuse ,                            &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml , maxelt , ia(ils), &
   ipnfac , nodfac , ipnfbr , nodfbr , ia(iidfst),                &
   idevel , ituser , ia     ,                                     &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod , volume , &
   rdevel , rtuser , ra     )


!===============================================================================
! 3.  CALCUL DE NBSTRU ET NBASTE
!===============================================================================

! 3.1 STRUCTURES INTERNES :
! -----------------------

nbstru = 0
do ifac = 1, nfabor
  if (ia(iidfst+ifac-1).gt.nbstru) nbstru = ia(iidfst+ifac-1)
enddo

if (irangp.ge.0) call parcmx(nbstru)
                 !==========

if (nbstru.gt.nstrmx) then
  write(nfecra,4000)
  call csexit(1)
endif

!     On compare NBSTRU a la valeur eventuelle anterieure

if (mbstru.gt.-999) then
  if (nbstru.ne.mbstru) then
    write(nfecra,4001)mbstru,nbstru
    call csexit(1)
  endif
endif

! 3.2 STRUCTURES EXTERNES : COUPLAGE CODE_SATURNE / CODE_ASTER
! -----------------------

nbaste = 0
do ifac = 1, nfabor
  if (-ia(iidfst+ifac-1).gt.nbaste) nbaste = -ia(iidfst+ifac-1)
enddo

if (irangp.ge.0) call parcmx(nbaste)
                 !==========

if (nbaste.gt.nastmx) then
  write(nfecra,4002)
  call csexit(1)
endif

!     On compare NBASTE a la valeur eventuelle anterieure

if (mbaste.gt.-999) then
  if (nbaste.ne.mbaste) then
    write(nfecra,4003)mbaste,nbaste
    call csexit(1)
  endif
endif


!===============================================================================
! 5.  CALCUL ET ENVOI A CODE_ASTER DES PARAMETRES GEOMETRIQUES
!===============================================================================

if (nbaste.gt.0) then

  itrav  = ifinia
  ifnia2 = itrav + nnod
  CALL IASIZE('STRINI',IFNIA2)

  do inod = 1, nnod
     ia(itrav + inod-1) = 0
  enddo

  nbfast = 0
  nbnast = 0

!       Calcul du nombre de faces et noeuds couples avec Code_Aster
  do ifac = 1, nfabor
    istr = ia(iidfst+ifac-1)
    if (istr.lt.0) then
      nbfast = nbfast + 1
      do ii = ipnfbr(ifac), ipnfbr(ifac+1)-1
        inod = nodfbr(ii)
        ia(itrav + inod-1) = istr
      enddo
    endif
  enddo
  do inod = 1, nnod
    if (ia(itrav+inod-1).gt.0) nbnast = nbnast + 1
  enddo

  ilstfa = ifinia
  iidflo = ilstfa + nbfast
  iidnlo = iidflo + nbfast
  ifnia2 = iidnlo + nbnast
  CALL IASIZE('STRINI',IFNIA2)

  indast = 0
  do ifac = 1, nfabor
    istr = ia(iidfst+ifac-1)
    if (istr.lt.0) then
      indast = indast + 1
      ia(ilstfa + indast-1) = ifac
      ia(iidflo + indast-1) = -istr
    endif
  enddo
  nbfast = indast

  indast = 0
  do inod = 1, nnod
    istr = ia(itrav+inod-1)
    if (istr.lt.0) then
      indast = indast + 1
      ia(iidnlo + indast-1) = -istr
    endif
  enddo
  nbnast = indast

!       Recuperation des parametres commun du couplage
  call astpar                                                     &
  !==========
 ( ntmabs, nalimx, epalim, isyncp, ntchr, ttpabs, dtref )

!       Envoi des donnees geometriques a Code_Aster
  call astgeo                                                     &
  !==========
 ( nbfast, nbnast, ia(ilstfa), ia(iidflo), ia(iidnlo), almax )

endif


!===============================================================================
! 6.  MESSAGES D'INFORMATION SUR LE COUPLAGE
!===============================================================================

!     Valeur par defaut et verifiction de IHISTR
if (nbstru.eq.0) ihistr = 0
if (nbaste.eq.0) isyncp = 0

icompt = 0
do ii = 2, nvppmx
  if(ihisvr(ii,1).ne.0) icompt = icompt+1
enddo
if( (icompt.eq.0.or.ncapt.eq.0) .and. ihistr.eq.0 ) then
  nthist = -1
endif

if (ihistr.ne.0 .and. ihistr.ne.1) then
  write(nfecra,1000)ihistr
  call csexit(1)
endif
if (isyncp.ne.0 .and. isyncp.ne.1) then
  write(nfecra,1002)isyncp
  call csexit(1)
endif

!     Si NBSTRU=0 et NBASTE=0, on desalloue IDFSTR et on passe NALIMX a 1
!       si necessaire
if (nbstru.gt.0) then
  write(nfecra,2010) nbstru,alpnmk,betnmk,gamnmk,ihistr
else
  write(nfecra,2000) nbstru
endif
if (nbaste.gt.0) then
  write(nfecra,2012) nbaste,isyncp
else
  write(nfecra,2002) nbaste
endif
if (nbstru.eq.0.and.nbaste.eq.0) then
  ifinia = idebia
  ifinra = idebra
  if (nalimx.gt.1) then
    write(nfecra,2001)
    nalimx = 1
  endif
else if (nbstru.gt.0) then
  if (nalimx.eq.1) then
    write(nfecra,2020) aexxst, bexxst, cfopre
  else
    cfopre = 1.d0
    write(nfecra,2030) nalimx, epalim
  endif
endif
write(nfecra,3000)

!----
! FORMATS
!----

#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR D''ECRITURE DES FICHIERS HISTORIQUES DES  ',/,&
'@      STRUCTURES MOBILES NE PEUT VALOIR QUE 0 OU 1.         ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes dans usstru.               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR D''ECRITURE AU MEME INSTANT DES SORTIES   ',/,&
'@      CODE_SATURNE ET CODE_ASTER NE PEUT VALOIR QUE 0 OU 1. ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@    LA PERIODE D''IMPRESSION EST DEFINIE PAR LA VARIABLE :  ',/,&
'@      NTCHR                                                 ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes dans usstru.               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
    /,'TTES PHASES  : MODE COUPLAGE DE STRUCTURES NON ACTIVE  ',/,&
      '                 NBSTRU = ',I10                         ,/)
 2001 format(                                                           &
      '               NALIMX INUTILE ET POSITIONNE A 1        ',/)
 2002 format(                                                           &
    /,'TTES PHASES  : MODE COUPLAGE CODE_ASTER NON ACTIVE     ',/,&
      '                 NBASTE = ',I10                         ,/)
 2010 format(                                                           &
    /,'TTES PHASES  : MODE COUPLAGE DE STRUCTURES ACTIVE      ',/,&
      '                 AVEC NBSTRU = ',I10   ,' STRUCTURE(S) ',/,&
      '                                                       ',/,&
      '               COEFFICIENTS DE NEWMARK :               ',/,&
      '                 ALPNMK = ',E12.4                       ,/,&
      '                 BETNMK = ',E12.4                       ,/,&
      '                 GAMNMK = ',E12.4                       ,/,&
      '                                                       ',/,&
      '               FICHIERS HISTORIQUES DES STRUCTURES :   ',/,&
      '                 IHISTR = ',I4,' ( 1 : active)         ',/)
 2012 format(                                                           &
    /,'TTES PHASES  : MODE COUPLAGE CODE_ASTER ACTIVE         ',/,&
      '                 AVEC NBASTE = ',I10   ,' STRUCTURE(S) ',/,&
      '                                                       ',/,&
      '               IMPRESSIONS  DES SORTIES AUX MEMES      ',/,&
      '               INSTANTS DES DEUX CODES :               ',/,&
      '                 ISYNCP = ',I4,' ( 1 : active)         ',/)
 2020 format(                                                           &
    /,'TTES PHASES  : SCHEMA DE COUPLAGE EXPLICITE ACTIVE     ',/,&
      '                                                       ',/,&
      '               COEFFICIENTS DU SCHEMA :                ',/,&
      '                 AEXXST = ',E12.4                       ,/,&
      '                 BEXXST = ',E12.4                       ,/,&
      '                 CFOPRE = ',E12.4                       ,/)
 2030 format(                                                           &
    /,'TTES PHASES  : SCHEMA DE COUPLAGE IMPLICITE ACTIVE     ',/,&
      '                                                       ',/,&
      '               NB DE SOUS-ITERATIONS MAX. : ',I10       ,/,&
      '               SEUIL DE CONVERGENCE       : ',E12.4     ,/)

 3000 format(                                                           &
'-------------------------------------------------------------',/)
 4000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES        ',/,&
'@                INTERNES                                    ',/,&
'@                                                            ',/,&
'@    Le nombre de structures definies est superieur au nombre',/,&
'@    maximum autorise NSTRMX :                               ',/,&
'@      Nombre de structures definies      : ',I10             ,/,&
'@      Nombre de structures autorisees    : ',I10             ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Reduire le nombre de structure                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES MOBILES',/,&
'@                INTERNES                                    ',/,&
'@                                                            ',/,&
'@    Le nombre de structures internes definies est           ',/,&
'@      different de celui du calcul precedent :              ',/,&
'@      Nombre de structures calcul precedent : ',I10          ,/,&
'@      Nombre de structures calcul actuel    : ',I10          ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire ou la specification',/,&
'@      des structures dans usstru.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES        ',/,&
'@                              EXTERNES (COUPLAGE CODE_ASTER)',/,&
'@                                                            ',/,&
'@    Le nombre de structures definies est superieur au nombre',/,&
'@    maximum autorise NASTMX :                               ',/,&
'@      Nombre de structures definies      : ',I10             ,/,&
'@      Nombre de structures autorisees    : ',I10             ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Reduire le nombre de structure                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA DEFINITION DES STRUCTURES MOBILES',/,&
'@                             EXTERNES  (COUPLAGE CODE_ASTER)',/,&
'@                                                            ',/,&
'@    Le nombre de structures externes definies est           ',/,&
'@      different de celui du calcul precedent :              ',/,&
'@      Nombre de structures calcul precedent : ',I10          ,/,&
'@      Nombre de structures calcul actuel    : ',I10          ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire ou la specification',/,&
'@      des structures dans usaste.                           ',/,&
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
'@    THE TIME MONITORING FILES INDICATOR FOR THE MOBILE      ',/,&
'@      STRUCTURES CAN ONLY TAKE THE VALUES 0 OR 1.           ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify the parameters given in usstru.                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE DATA SPECIFICATION                ',/,&
'@    ========                                                ',/,&
'@    THE INDICATOR OF SYNCHRONIZED OUTPUTS FOR CODE_SATURNE  ',/,&
'@      AND CODE_ASTER CAN ONLY TAKE THE VALUES 0 OR 1.       ',/,&
'@    ITS VALUE IS ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify the parameters given in usstru.                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2000 format(                                                           &
    /,'ALL PHASES: COUPLING MODE FOR STRUCTURES NOT ACTIVATED ',/,&
      '              NBSTRU = ',I10                            ,/)
 2001 format(                                                           &
      '            NALIMX USELESS AND SET TO 1                ',/)
 2002 format(                                                           &
    /,'ALL PHASES: CODE_ASTER COUPLING MODE NOT ACTIVATED     ',/,&
      '              NBASTE = ',I10                            ,/)
 2010 format(                                                           &
    /,'ALL PHASES: COUPLING MODE FOR STRUCTURES ACTIVATED     ',/,&
      '              WITH NBSTRU = ',I10   ,' STRUCTURE(S)    ',/,&
      '                                                       ',/,&
      '            NEWMARK COEFFICIENTS:                      ',/,&
      '              ALPNMK = ',E12.4                          ,/,&
      '              BETNMK = ',E12.4                          ,/,&
      '              GAMNMK = ',E12.4                          ,/,&
      '                                                       ',/,&
      '            MONITORING FILES FOR STRUCTURES:           ',/,&
      '                 IHISTR = ',I4,' ( 1 : activated)      ',/)
 2012 format(                                                           &
    /,'ALL PHASES: CPDE_ASTER COUPLING MODE ACTIVATED         ',/,&
      '              WITH NBASTE = ',I10   ,' STRUCTURE(S)    ',/,&
      '                                                       ',/,&
      '            SYNCHRONIZED OUTPUTS FOR BOTH CODES:       ',/,&
      '                 ISYNCP = ',I4,' ( 1 : activated)      ',/)
 2020 format(                                                           &
    /,'ALL PHASES: EXPLICIT SCHEME FOR COUPLING ACTIVATED     ',/,&
      '                                                       ',/,&
      '            SCHEME COEFFICIENTS:                       ',/,&
      '              AEXXST = ',E12.4                          ,/,&
      '              BEXXST = ',E12.4                          ,/,&
      '              CFOPRE = ',E12.4                          ,/)
 2030 format(                                                           &
    /,'ALL PHASES: IMPLICIT SCHEME FOR COUPING ACTIVATED      ',/,&
      '                                                       ',/,&
      '            NB OF MAX INNER ITERATIONS : ',I10          ,/,&
      '            CONVERGENCE THRESHOLD      : ',E12.4        ,/)

 3000 format(                                                           &
'-------------------------------------------------------------',/)
 4000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE INTERNAL STRUCTURES SPECIFICATION ',/,&
'@                                                            ',/,&
'@    The number of defined structures is greater than the    ',/,&
'@      allowed maximum NSTRMX:                               ',/,&
'@      Number of defined structures: ',I10                    ,/,&
'@      Number of allowed structures: ',I10                    ,/,&
'@                                                            ',/,&
'@    The calculation will not be run.                        ',/,&
'@                                                            ',/,&
'@    Decrease the number of structures                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE INTERNAL MOBILE STRUCTURES        ',/,&
'@             SPECIFICATION                                  ',/,&
'@                                                            ',/,&
'@    The number of defined structures is different from the  ',/,&
'@      previous calculation:                                 ',/,&
'@      Number of structures previous calculation: ',I10       ,/,&
'@      Number of structures current  calculation: ',I10       ,/,&
'@                                                            ',/,&
'@    The calculation will not be run.                        ',/,&
'@                                                            ',/,&
'@    Verify the auxilliary restart file or the structures    ',/,&
'@      specifications in usstru.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE EXTERNAL MOBILE STRUCTURES        ',/,&
'@             SPECIFICATION (CODE_ASTER COUPLING)            ',/,&
'@                                                            ',/,&
'@    The number of defined structures is greater than the    ',/,&
'@      allowed maximum NASTMX:                               ',/,&
'@      Number of defined structures: ',I10                    ,/,&
'@      Number of allowed structures: ',I10                    ,/,&
'@                                                            ',/,&
'@    The calculation will not be run.                        ',/,&
'@                                                            ',/,&
'@    Decrease the number of structures                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE EXTERNAL MOBILE STRUCTURES        ',/,&
'@             SPECIFICATION (CODE_ASTER COUPLING)            ',/,&
'@                                                            ',/,&
'@    The number of defined structures is different from the  ',/,&
'@      previous calculation:                                 ',/,&
'@      Number of structures previous calculation: ',I10       ,/,&
'@      Number of structures current  calculation: ',I10       ,/,&
'@                                                            ',/,&
'@    The calculation will not be run.                        ',/,&
'@                                                            ',/,&
'@    Verify the auxilliary restart file or the structures    ',/,&
'@      specifications in usstru.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----

end subroutine
