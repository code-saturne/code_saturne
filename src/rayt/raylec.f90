!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2013 EDF S.A.
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

subroutine raylec &
!================

 ( ndim   , ncelet , ncel   , nfac   , nfabor ,                   &
   propce , propfb )

!===============================================================================
! Purpose:
! --------

!   SOUS-PROGRAMME DU MODULE RAYONNEMENT :
!   --------------------------------------

!         Lecture du fichier suite au 1er passage


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
! propfb(nfabor, *)! ra ! <-- ! physical properties at boundary face centers   !
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use pointe
use entsor
use cstphy
use ppppar
use ppthch
use cpincl
use ppincl
use radiat

!===============================================================================

implicit none

! Arguments

integer          ndim   , ncelet , ncel   , nfac   , nfabor

double precision propce(ncelet,*)
double precision propfb(nfabor,*)
!
! Local variables

character        rubriq*64
character        cphase*2
character        ficsui*32
integer          iok

integer          jphast
integer          ncelok , nfaiok , nfabok , nsomok
integer          ierror , irtyp  , itysup , nbval
integer          ilecec , nberro , ivers
integer          impamr

!===============================================================================
! 0 - GESTION MEMOIRE
!===============================================================================

!===============================================================================
! 1. LECTURE DU FICHIER SUITE
!===============================================================================

if (isuird.eq.1) then

!  ---> Ouverture

    write(nfecra,6000)

!     (ILECEC=1:lecture)
    ilecec = 1
    ficsui = 'radiative_transfer'
    call opnsui(ficsui,len(ficsui),ilecec,impamr,ierror)
    !==========
    if (ierror.ne.0) then
      write(nfecra,9011) ficsui
      call csexit (1)
    endif

    write(nfecra,6010)


!  ---> Type de fichier suite
!        Pourrait porter le numero de version si besoin.
!        On ne se sert pas de IVERS pour le moment

    itysup = 0
    nbval  = 1
    irtyp  = 1
    RUBRIQ = 'version_fichier_suite_rayonnement'
    call lecsui(impamr,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ivers,ierror)

    if (ierror.ne.0) then
      write(nfecra,9200)ficsui
      call csexit (1)
    endif


!  ---> Tests

    iok = 0

!     Dimensions des supports

    call tstsui(impamr,ncelok,nfaiok,nfabok,nsomok)
    !==========
    if (ncelok.eq.0) then
      write(nfecra,9210)
      iok = iok + 1
    endif
    if (nfabok.eq.0) then
      write(nfecra,9211)
      iok = iok + 1
    endif


!  ---> Pour test ulterieur si pb : arret

    nberro = 0

!  ---> Lecture des donnees

!     Aux faces de bord

      itysup = 3
      nbval  = 1
      irtyp  = 2

      RUBRIQ = 'tparoi_fb'
      call lecsui(impamr,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           propfb(1,ipprob(itparo)),ierror)
      nberro=nberro+ierror

      RUBRIQ = 'qincid_fb'
      call lecsui(impamr,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           propfb(1,ipprob(iqinci)),ierror)
      nberro=nberro+ierror

      RUBRIQ = 'hfconv_fb'
      call lecsui(impamr,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           propfb(1,ipprob(ihconv)),ierror)
      nberro=nberro+ierror

     RUBRIQ = 'flconv_fb'
      call lecsui(impamr,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           propfb(1,ipprob(ifconv)),ierror)
      nberro=nberro+ierror


!     Aux cellules

     itysup = 1
      nbval  = 1
      irtyp  = 2

      RUBRIQ = 'rayimp_ce'
      call lecsui(impamr,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           propce(1,ipproc(itsri(1))),ierror)
      nberro=nberro+ierror

      RUBRIQ = 'rayexp_ce'
      call lecsui(impamr,rubriq,len(rubriq),itysup,nbval,irtyp,   &
           propce(1,ipproc(itsre(1))),ierror)
      nberro=nberro+ierror

      RUBRIQ = 'luminance'
      call lecsui(impamr,rubriq,len(rubriq),itysup,nbval,irtyp,   &
                  propce(1,ipproc(ilumin)),ierror)
      nberro=nberro+ierror


!  ---> Si pb : arret

    if(nberro.ne.0) then
      write(nfecra,9100)
      call csexit (1)
    endif

    write(nfecra,6011)

!  ---> Fermeture du fichier suite

    call clssui(impamr,ierror)

    if (ierror.ne.0) then
      write(nfecra,8011) ficsui
    endif

    write(nfecra,6099)

! Fin detection suite rayonnement
endif

!--------
! FORMATS
!--------

 6000 FORMAT (   3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT   ',/,&
           3X,'   ------------------------------------------  ',/,&
           3X,' Lecture d''un fichier suite                   '  )
 6010 FORMAT (   3X,'   Debut de la lecture                         '  )
 6011 FORMAT (   3X,'   Fin   de la lecture                         '  )
 6099 FORMAT (   3X,' Fin de la lecture du fichier suite            ',/,&
'                                                             ',/,&
'-------------------------------------------------------------',/)

 8011 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,&
'@    =========   RAYONNEMENT                                 ',/,&
'@                                                            ',/,&
'@    Probleme sur le fichier de nom (',A13,')                ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9011 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========   RAYONNEMENT                                 ',/,&
'@      ERREUR A L''OUVERTURE DU FICHIER SUITE                ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier l''existence et le nom (',A13,') du            ',/,&
'@        fichier suite dans le repertoire de travail.        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/)
 9200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                    RAYONNEMENT',/,&
'@      TYPE DE FICHIER INCORRECT                             ',/,&
'@                                                            ',/,&
'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,&
'@      suite rayonnement.                                    ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        a un fichier suite rayonnement.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========   RAYONNEMENT                                 ',/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    Le nombre de cellules a ete modifie                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9211 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========   RAYONNEMENT                                 ',/,&
'@      DONNEES AMONT ET ACTUELLES INCOHERENTES               ',/,&
'@                                                            ',/,&
'@    Le nombre de faces de bord a ete modifie                ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9100 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION: ARRET A LA LECTURE DU FICHIER SUITE          ',/,&
'@    =========   RAYONNEMENT                                 ',/,&
'@      ERREUR LORS DE LA LECTURE DES DONNEES                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! End
!----

return

end subroutine
