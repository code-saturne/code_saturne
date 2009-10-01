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

subroutine strhis &
!================

 ( idbia0 , idbra0 , ncelet , ncel ,                              &
   nideve , nrdeve , nituse , nrtuse , modhis ,                   &
   idevel , ituser , ia     ,                                     &
   rdevel , rtuser , ra )

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE D'ECRITURE DES HISTORIQUES

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
! idbia0           ! e  ! <-- ! numero de la 1ere case libre dans ia           !
! idbra0           ! e  ! <-- ! numero de la 1ere case libre dans ra           !
! ncelet           ! e  ! <-- ! nombre d'elements halo compris                 !
! ncel             ! e  ! <-- ! nombre d'elements actifs                       !
! nideve nrdeve    ! e  ! <-- ! longueur de idevel rdevel                      !
! nituse nrtuse    ! e  ! <-- ! longueur de ituser rtuser                      !
! modhis           ! e  ! <-- ! indicateur valant 0,1 ou 2                     !
!                  !    !               ! 1,2 = ecriture intermediaire, finale |
! idevel(nideve    ! te ! <-- ! tab entier complementaire developemt           !
! ituser(nituse    ! te ! <-- ! tab entier complementaire utilisateur          !
! ia(*)            ! tr ! --- ! macro tableau entier                           !
! rdevel(nrdeve    ! tr ! <-- ! tab reel complementaire developemt             !
! rtuser(nrtuse    ! tr ! <-- ! tab reel complementaire utilisateur            !
! ra               ! tr !  -- ! tableau des reels                              !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "entsor.h"
include "optcal.h"
include "parall.h"
include "alstru.h"

!===============================================================================

! Arguments

integer          idbia0, idbra0
integer          ncelet, ncel
integer          nideve , nrdeve , nituse , nrtuse
integer          modhis
integer          idevel(nideve), ituser(nituse), ia(*)
double precision rdevel(nrdeve), rtuser(nrtuse), ra(*)

! VARIABLES LOCALES

integer          nbname
parameter        (nbname=12)
character        nomfic*300, nenvar*300
character*80     namevr(nbname)
integer          ii, jj, istr, ii1, ii2, lpos, inam1, inam2
integer          nbpdte, jtcabs
integer          idebia, idebra
double precision xtcabs
double precision vartmp(nstrmx)

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
! 0. INITIALISATIONS LOCALES
!===============================================================================

!     Seul le processeur 0 ecrit, on n'ecrit rien si pas de structure
if (irangp.gt.0 .or. nbstru.le.0) return

ipass = ipass + 1

idebia = idbia0
idebra = idbra0

namevr(1 ) = "deplacement x"
namevr(2 ) = "deplacement y"
namevr(3 ) = "deplacement z"
namevr(4 ) = "vitesse x"
namevr(5 ) = "vitesse y"
namevr(6 ) = "vitesse z"
namevr(7 ) = "acceleration x"
namevr(8 ) = "acceleration y"
namevr(9 ) = "acceleration z"
namevr(10) = "force x"
namevr(11) = "force y"
namevr(12) = "force z"

!--> Il n'y a pas eu d'historiques
if(ipass.eq.1.and.modhis.eq.2) return

!===============================================================================
! 2. OUVERTURE DU FICHIER DE STOCKAGE hist.tmp
!===============================================================================

if(ipass.eq.1) then
  NOMFIC = ' '
  nomfic = emphis
  call verlon ( nomfic,ii1,ii2,lpos)
  !==========

  NOMFIC(II2+1:II2+11) = 'histstr.tmp'
  ii2 = ii2+11
  open ( unit=impsth(1), file=nomfic (ii1:ii2),                   &
         STATUS='UNKNOWN', FORM='UNFORMATTED',                    &
         ACCESS='SEQUENTIAL')
endif

!===============================================================================
! 3. ECRITURE DES RESULTATS dans le FICHIER DE STOCKAGE
!===============================================================================

if(modhis.eq.0.or.modhis.eq.1) then

  write(impsth(1)) ntcabs, ttcabs, (xstr  (1,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (xstr  (2,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (xstr  (3,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (xpstr (1,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (xpstr (2,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (xpstr (3,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (xppstr(1,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (xppstr(2,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (xppstr(3,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (forstr(1,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (forstr(2,ii),ii=1,nbstru)
  write(impsth(1)) ntcabs, ttcabs, (forstr(3,ii),ii=1,nbstru)

endif


!===============================================================================
! 4. EN CAS DE SAUVEGARDE INTERMEDIAIRE OU FINALE,
!    TRANSMISSION DES INFORMATIONS DANS LES DIFFERENTS FICHIERS
!===============================================================================

! On sauve aussi au premier passage pour permettre une
!     verification des le debut du calcul

if(modhis.eq.1.or.modhis.eq.2.or.ipass.eq.1) then

!       --> nombre de pas de temps enregistres

  if(modhis.eq.2) then
    nbpdte = ipass - 1
  else
    nbpdte = ipass
  endif

!       --> ecriture un fichier par variable
  do ii = 1, nbname

!         --> nom du fichier
    NOMFIC = ' '
    nomfic = emphis
    call verlon ( nomfic,ii1,ii2,lpos)
    !==========
    NENVAR = 'str_'//NAMEVR(II)
    call verlon(nenvar,inam1,inam2,lpos)
    !==========
    call undscr(inam1,inam2,nenvar)
    !==========
    nomfic(ii2+1:ii2+lpos) = nenvar(inam1:inam2)
    ii2 = ii2+lpos
    NOMFIC(II2+1:II2+1) = '.'
    ii2 = ii2+1
    nenvar = exthis
    call verlon(nenvar,inam1,inam2,lpos)
    !==========
    call undscr(inam1,inam2,nenvar)
    !==========
    nomfic(ii2+1:ii2+lpos) = nenvar(inam1:inam2)
    ii2 = ii2+lpos
!         --> ouverture
    open ( unit=impsth(2), file=nomfic (ii1:ii2),                 &
         STATUS='UNKNOWN', FORM='FORMATTED',                      &
         ACCESS='SEQUENTIAL')
!         --> entete
    write(impsth(2),1000)namevr(ii),nbstru,nbpdte,nbstru+2
    write(impsth(2),2000) (istr,istr=1,nbstru)
    write(impsth(2),2005)
!         --> impression des matrices de masse
    write(impsth(2),2001) ((xmstru(1,jj,istr),jj=1,3),            &
         istr=1,nbstru)
    write(impsth(2),2004) ((xmstru(2,jj,istr),jj=1,3),            &
         istr=1,nbstru)
    write(impsth(2),2004) ((xmstru(3,jj,istr),jj=1,3),            &
         istr=1,nbstru)
    write(impsth(2),2005)
!         --> impression des matrices d'amortissement
    write(impsth(2),2002) ((xcstru(1,jj,istr),jj=1,3),            &
         istr=1,nbstru)
    write(impsth(2),2004) ((xcstru(2,jj,istr),jj=1,3),            &
         istr=1,nbstru)
    write(impsth(2),2004) ((xcstru(3,jj,istr),jj=1,3),            &
         istr=1,nbstru)
    write(impsth(2),2005)
!         --> impression des matrices de raideur
    write(impsth(2),2003) ((xkstru(1,jj,istr),jj=1,3),            &
         istr=1,nbstru)
    write(impsth(2),2004) ((xkstru(2,jj,istr),jj=1,3),            &
         istr=1,nbstru)
    write(impsth(2),2004) ((xkstru(3,jj,istr),jj=1,3),            &
         istr=1,nbstru)
    write(impsth(2),2005)
!         --> impression de fin de section
    write(impsth(2),2006)


!         --> boucle sur les differents enregistrements
!             et les variables
    rewind(impsth(1))
    do jj = 1, nbpdte
      do ii1 = 1, nbname
        read(impsth(1))                                           &
             jtcabs, xtcabs, (vartmp(istr),istr=1,nbstru)
        if(ii1.eq.ii)                                             &
             write(impsth(2),3000)                                &
             jtcabs, xtcabs, (vartmp(istr),istr=1,nbstru)
      enddo
    enddo

!         --> fermeture fichier
    close(impsth(2))

  enddo

endif


!----
! FORMATS
!----
#if defined(_CS_LANG_FR)

 1000 format(                                                           &
'# ---------------------------------------------------',/,  &
'#      FICHIER HISTORIQUE EN TEMPS                   ',/,  &
'#      VARIABLE ',A30                                 ,/,  &
'# ---------------------------------------------------',/,  &
'# ',/,                                                     &
'#      NOMBRE DE STRUCTURES      : ',I8,/,                 &
'# ',/,                                                     &
'#      NOMBRE D''ENREGISTREMENTS  : ',I8,/,                &
'# ',/,                                                     &
'# ',/,                                                     &
'# COLONNE 1        : NUMERO DU PAS DE TEMPS',/,            &
'#         2        : TEMPS PHYSIQUE (ou No pas de temps*DTREF',/,&
'#                               en pas de temps non uniforme)',/,&
'#         3 A ',I4,' : VALEUR POUR CHAQUE STRUCTURE          ',/,&
'# ',/,                                                     &
'# ---------------------------------------------------',/,  &
'# ')
 2000 format('# STRUCTURE    |',20(21X,I3,22X,'|'))
 2001 format('# MASSE        |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2002 format('# AMORTISSEMENT|',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2003 format('# RAIDEUR      |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2004 format('#              |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2005 format('#')
 2006 format(                                                           &
'# (dans le cas ou les caracteristiques des structures sont   ',/,&
'#  variables, les valeurs ci-dessus correspondent au dernier ',/,&
'#  pas de temps d''ecriture dans le fichier',/,            &
'#',/,                                                      &
'# ---------------------------------------------------',/,  &
'# ')
 3000 format ( 1x,i7,1x,21(1x,e14.7))

#else

 1000 format (                                                          &
'# ---------------------------------------------------',/,  &
'#      TIME MONITORING FILE                          ',/,  &
'#      VARIABLE ',A30                                 ,/,  &
'# ---------------------------------------------------',/,  &
'# ',/,                                                     &
'#      NUMBER OF STRUCTURES  : ',I8,/,                     &
'# ',/,                                                     &
'#      NUMBER OF RECORDS     : ',I8,/,                     &
'# ',/,                                                     &
'# ',/,                                                     &
'# COLUMN 1        : TIME STEP NUMBER ',/,                  &
'#        2        : PHYSICAL TIME (or Nb of time steps*DTREF ',/,&
'#                                with non uniform time step)',/, &
'#        3 TO ',I4,' : VALUE FOR EACH STRUCTURE              ',/,&
'# ',/,                                                     &
'# ---------------------------------------------------',/,  &
'# ')
 2000 format('# STRUCTURE    |',20(21X,I3,22X,'|'))
 2001 format('# MASS         |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2002 format('# DAMPING      |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2003 format('# STIFNESS     |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2004 format('#              |',20(1X,G14.7,1X,G14.7,1X,G14.7,' |'))
 2005 format('#')
 2006 format(                                                           &
'# (in the case where the structures characteristics are      ',/,&
'#  variables, these values correspond to the latest time step',/,&
'#  written in the file',/,                                 &
'#',/,                                                      &
'# ---------------------------------------------------',/,  &
'# ')
 3000 format ( 1x,i7,1x,21(1x,e14.7))

#endif

return
end subroutine


