!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2015 EDF S.A.
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

subroutine raylec
!================

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
!__________________!____!_____!________________________________________________!

!     Type: i (integer), r (real), s (string), a (array), l (logical),
!           and composite types (ex: ra real array)
!     mode: <-- input, --> output, <-> modifies data, --- work array
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use, intrinsic :: iso_c_binding

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
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

character        rubriq*64
character        ficsui*32
integer          iok

integer          ierror , itysup , nbval
integer          nberro , ivers
integer          ival(1)

logical(kind=c_bool) :: ncelok, nfaiok, nfabok, nsomok

type(c_ptr) :: rp

!===============================================================================
! 1. LECTURE DU FICHIER SUITE
!===============================================================================

if (isuird.eq.1) then

  ! Ouverture

  write(nfecra,6000)

  ficsui = 'radiative_transfer'
  call restart_create(ficsui, '', 0, rp)

  write(nfecra,6010)

  ! Type de fichier suite
  !    Pourrait porter le numero de version si besoin.
  !    On ne se sert pas de IVERS pour le moment

  itysup = 0
  nbval  = 1
  rubriq = 'version_fichier_suite_rayonnement'
  call restart_read_section_int_t(rp,rubriq,itysup,nbval,ival,ierror)
  ivers = ival(1)

  if (ierror.ne.0) then
    write(nfecra,9200)ficsui
    call csexit (1)
  endif

  !  Tests

  iok = 0

  ! Dimensions des supports

  call restart_check_base_location(rp,ncelok,nfaiok,nfabok,nsomok)
  if (ncelok.eqv..false.) then
    write(nfecra,9210)
    iok = iok + 1
  endif
  if (nfabok.eqv..false.) then
    write(nfecra,9211)
    iok = iok + 1
  endif

  ! ---> Pour test ulterieur si pb : arret

  nberro = 0

  ! Aux faces de bord

  call restart_read_field_vals(rp, itparo, 0, ierror)
  nberro=nberro+ierror

  call restart_read_field_vals(rp, iqinci, 0, ierror)
  nberro=nberro+ierror

  call restart_read_field_vals(rp, ihconv, 0, ierror)
  nberro=nberro+ierror

  call restart_read_field_vals(rp, ifconv, 0, ierror)
  nberro=nberro+ierror

  ! Aux cellules

  call restart_read_field_vals(rp, iprpfl(itsri(1)), 0, ierror)
  nberro=nberro+ierror

  call restart_read_field_vals(rp, iprpfl(itsre(1)), 0, ierror)
  nberro=nberro+ierror

  call restart_read_field_vals(rp, iprpfl(ilumin), 0, ierror)
  nberro=nberro+ierror

  !---> Si pb : arret

  if (nberro.ne.0) then
    write(nfecra,9100)
    call csexit (1)
  endif

  write(nfecra,6011)

  ! ---> Fermeture du fichier suite

  call restart_destroy(rp)

  write(nfecra,6099)

  ! Fin detection suite rayonnement
endif

!--------
! Formats
!--------

 6000 format (   3X,'** INFORMATIONS SUR LE MODULE DE RAYONNEMENT   ',/,&
           3X,'   ------------------------------------------  ',/,&
           3X,' Lecture d''un fichier suite                   '  )
 6010 format (   3X,'   Debut de la lecture                         '  )
 6011 format (   3X,'   Fin   de la lecture                         '  )
 6099 format (   3X,' Fin de la lecture du fichier suite            ',/,&
'                                                             ',/,&
'-------------------------------------------------------------',/)

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
