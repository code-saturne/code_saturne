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

subroutine cs_coal_prop &
!======================
 ( ipropp , ipppst )

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES POSITIONS DES VARIABLES D'ETAT SELON
!         COMBUSTION CHARBON PULVERISE
!   (DANS VECTEURS PROPCE)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ipropp           ! e  ! <-- ! numero de la derniere propriete                !
!                  !    !     !  (les proprietes sont dans propce)             !
! ipppst           ! e  ! <-- ! pointeur indiquant le rang de la               !
!                  !    !     !  derniere grandeur definie aux                 !
!                  !    !     !  cellules (rtp,propce...) pour le              !
!                  !    !     !  post traitement                               !
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
use dimens
use numvar
use optcal
use cstphy
use entsor
use cstnum
use ppppar
use ppthch
use ppcpfu
use coincl
use cpincl
use ppincl
use ihmpre
use cs_coal_incl
use field

!===============================================================================

implicit none

! Arguments

integer          ipropp, ipppst

! Local variables

integer          iprop, ige , icla, iprop2
integer          f_id, itycat, ityloc, idim1
integer          keyccl
integer          iopchr

logical          ilved, iprev, inoprv

character*80     f_name


!===============================================================================

itycat = FIELD_INTENSIVE + FIELD_PROPERTY
ityloc = 1 ! variables defined on cells
idim1  = 1
ilved  = .false.   ! not interleaved by default
iprev  = .true.    ! variables have previous value
inoprv = .false.   ! variables have no previous value
iopchr = 1         ! Postprocessing level for variables

! Key id of the coal scalar class
call field_get_key_id("scalar_class", keyccl)

! ---> Definition des pointeurs relatifs aux variables d'etat

iprop = ipropp

!    Phase continue (melange gazeux)
iprop   = iprop + 1
itemp1  = iprop
iprop   = iprop + 1
irom1  = iprop
do ige = 1,ngazg
! ---- Cf. definition de NGAZE dans cs_coal_readata
  iprop     = iprop + 1
  iym1(ige) = iprop
enddo
iprop = iprop + 1
immel = iprop

if ( ieqnox .eq. 1 ) then
  iprop  = iprop + 1
  ighcn1 = iprop

  iprop  = iprop + 1
  ighcn2 = iprop

  iprop  = iprop + 1
  ignoth = iprop

  iprop  = iprop + 1
  ignh31 = iprop

  iprop  = iprop + 1
  ignh32 = iprop

  iprop  = iprop + 1
  ifhcnd = iprop

  iprop  = iprop + 1
  ifhcnc = iprop

  iprop  = iprop + 1
  ifnh3d = iprop

  iprop  = iprop + 1
  ifnh3c = iprop

  iprop  = iprop + 1
  ifnohc = iprop

  iprop  = iprop + 1
  ifnonh = iprop

  iprop  = iprop + 1
  ifnoch = iprop

  iprop  = iprop + 1
  ifnoth = iprop

  iprop  = iprop + 1
  icnohc = iprop

  iprop  = iprop + 1
  icnonh = iprop

  iprop  = iprop + 1
  ifhcnr = iprop

  iprop  = iprop + 1
  icnorb = iprop

  iprop  = iprop + 1
  igrb   = iprop

endif

iprop2 = iprop

!   Phase dispersee (classes de particules)
do icla = 1, nclacp
  iprop        = iprop2 + icla
  itemp2(icla) = iprop
  iprop        = iprop2 + 1*nclacp + icla
  ix2(icla)    = iprop
  iprop        = iprop2 + 2*nclacp + icla
  irom2(icla)  = iprop
  iprop        = iprop2 + 3*nclacp + icla
  idiam2(icla) = iprop
  iprop        = iprop2 + 4*nclacp + icla
  igmdch(icla) = iprop
  iprop        = iprop2 + 5*nclacp + icla
  igmdv1(icla) = iprop
  iprop        = iprop2 + 6*nclacp + icla
  igmdv2(icla) = iprop
  iprop        = iprop2 + 7*nclacp + icla
  igmhet(icla) = iprop

  if (i_coal_drift.eq.1) then
    write(f_name,'(a6,i2.2)')'Age_CP' ,icla
    call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
    call field_set_key_str(f_id, keylbl, f_name)
    ! Set the index of the scalar class in the field structure
    call field_set_key_int(f_id, keyccl, icla)
    ! For post-processing
    call field_set_key_int(f_id, keyvis, iopchr)
    ! For log in the listing
    call field_set_key_int(f_id, keylog, 1)
  endif

  if (ihtco2 .eq. 1) then
    iprop        = iprop2 + 8*nclacp + icla
    ighco2(icla) = iprop
    if ( ihth2o .eq. 1 ) then
      iprop        = iprop2 + 9*nclacp + icla
      ighh2o(icla) = iprop
      if ( ippmod(iccoal) .ge. 1 ) then
        iprop        = iprop2 + 10*nclacp + icla
        igmsec(icla) = iprop
      endif
    else
      if ( ippmod(iccoal) .ge. 1 ) then
        iprop        = iprop2 + 9*nclacp + icla
        igmsec(icla) = iprop
      endif
    endif
  else
    if ( ihth2o .eq. 1 ) then
      iprop        = iprop2 + 8*nclacp + icla
      ighh2o(icla) = iprop
      if ( ippmod(iccoal) .ge. 1 ) then
        iprop        = iprop2 + 9*nclacp + icla
        igmsec(icla) = iprop
      endif
    else
      if ( ippmod(iccoal) .ge. 1 ) then
        iprop        = iprop2 + 8*nclacp + icla
        igmsec(icla) = iprop
      endif
    endif
  endif
enddo

if (i_coal_drift.eq.1) then
  icla = -1
  f_name = 'Age_Gas'
  call field_create(f_name, itycat, ityloc, idim1, ilved, inoprv, f_id)
  call field_set_key_str(f_id, keylbl, f_name)
  ! Set the index of the scalar class in the field structure
  call field_set_key_int(f_id, keyccl, icla)

  ! For post-processing
  call field_set_key_int(f_id, keyvis, iopchr)
  ! For log in the listing
  call field_set_key_int(f_id, keylog, 1)
endif

!
! Bilan : C , O , H
!
iprop     = iprop + 1
ibcarbone = iprop
iprop     = iprop + 1
iboxygen  = iprop
iprop     = iprop + 1
ibhydrogen= iprop

! ---- Nb de variables algebriques (ou d'etat)
!         propre a la physique particuliere NSALPP
!         total NSALTO

nsalpp = iprop - ipropp
nsalto = iprop

! ----  On renvoie IPROPP au cas ou d'autres proprietes devraient
!         etre numerotees ensuite

ipropp = iprop

! ---> Positionnement dans le tableau PROPCE
!      et reperage du rang pour le post-traitement

iprop         = nproce

!    Phase continue (melange gazeux)
iprop           = iprop + 1
ipproc(itemp1)  = iprop
ipppst          = ipppst + 1
ipppro(iprop)   = ipppst

iprop           = iprop + 1
ipproc(irom1)   = iprop
ipppst          = ipppst + 1
ipppro(iprop)   = ipppst

do ige = 1, (ngaze-2*ncharb)
! ---- Cf. definition de NGAZE dans cs_coal_readata
  iprop                 = iprop + 1
  ipproc(iym1(ige))     = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
enddo

iprop                 = iprop + 1
ipproc(immel)         = iprop
ipppst                = ipppst + 1
ipppro(iprop)         = ipppst

!
if ( ieqnox .eq. 1 ) then
  iprop                 = iprop + 1
  ipproc(ighcn1)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ighcn2)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ignoth)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
! Positionnement des variables EXP4 et EXP5 dans le tableau de propriete
  iprop                 = iprop + 1
  ipproc(ignh31)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ignh32)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
! Affichage des termes source
  iprop                 = iprop + 1
  ipproc(ifhcnd)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ifhcnc)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ifnh3d)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ifnh3c)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ifnohc)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ifnonh)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ifnoch)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ifnoth)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(icnohc)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(icnonh)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(ifhcnr)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(icnorb)        = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
!
  iprop                 = iprop + 1
  ipproc(igrb)          = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst
endif

iprop2 = iprop

!   Phase dispersee (classes de particules)
do icla = 1, nclacp

  iprop                 = iprop2 + icla
  ipproc(itemp2(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 1*nclacp + icla
  ipproc(ix2(icla))     = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 2*nclacp + icla
  ipproc(irom2(icla))   = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 3*nclacp + icla
  ipproc(idiam2(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 4*nclacp + icla
  ipproc(igmdch(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 5*nclacp + icla
  ipproc(igmdv1(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 6*nclacp + icla
  ipproc(igmdv2(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iprop2 + 7*nclacp + icla
  ipproc(igmhet(icla))  = iprop
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  if ( ihtco2 .eq. 1 ) then
    iprop                 = iprop2 + 8*nclacp + icla
    ipproc(ighco2(icla))  = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    if ( ihth2o .eq. 1 ) then
      iprop                 = iprop2 + 9*nclacp + icla
      ipproc(ighh2o(icla))  = iprop
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst

      if ( ippmod(iccoal) .eq. 1 ) then
        iprop                 = iprop2 + 10*nclacp + icla
        ipproc(igmsec(icla))  = iprop
        ipppst                = ipppst + 1
        ipppro(iprop)         = ipppst
      endif

    else
      if ( ippmod(iccoal) .eq. 1 ) then
        iprop                 = iprop2 + 9*nclacp + icla
        ipproc(igmsec(icla))  = iprop
        ipppst                = ipppst + 1
        ipppro(iprop)         = ipppst
      endif
    endif

  else
    if ( ihth2o .eq. 1 ) then
      iprop                 = iprop2 + 8*nclacp + icla
      ipproc(ighh2o(icla))  = iprop
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst

      if ( ippmod(iccoal) .eq. 1 ) then
        iprop                 = iprop2 + 9*nclacp + icla
        ipproc(igmsec(icla))  = iprop
        ipppst                = ipppst + 1
        ipppro(iprop)         = ipppst
      endif

    else
      if ( ippmod(iccoal) .eq. 1 ) then
        iprop                 = iprop2 + 8*nclacp + icla
        ipproc(igmsec(icla))  = iprop
        ipppst                = ipppst + 1
        ipppro(iprop)         = ipppst
      endif
    endif
  endif

enddo
!
! Bilan C , O , H
!
iprop              = iprop  + 1
ipproc(ibcarbone)  = iprop
ipppst             = ipppst + 1
ipppro(iprop)      = ipppst
!
iprop              = iprop  + 1
ipproc(iboxygen)   = iprop
ipppst             = ipppst + 1
ipppro(iprop)      = ipppst
!
iprop              = iprop  + 1
ipproc(ibhydrogen) = iprop
ipppst             = ipppst + 1
ipppro(iprop)      = ipppst
!
nproce = iprop

!   - Interface Code_Saturne
!     ======================
!     Construction de l'indirection entre la numerotation du noyau et XML
if (iihmpr.eq.1) then
  call uicppr &
  !==========
 ( nclacp, nsalpp, ippmod, iccoal, ipppro,            &
   ipproc, ieqnox, ihtco2, ihth2o, itemp1,            &
   irom1, iym1, ighcn1, ighcn2, ignoth,               &
   ignh31, ignh32, ifhcnd, ifhcnc, ifnh3d, ifnh3c,    &
   ifnohc, ifnonh, ifnoch, ifnoth, icnohc, icnonh,    &
   ifhcnr, icnorb, igrb,   immel,                     &
   itemp2, ix2, irom2, idiam2, igmdch, igmdv1,        &
   igmdv2, igmhet, ighco2, ighh2o, igmsec,            &
   ibcarbone, iboxygen, ibhydrogen)
endif

!--------
! Formats
!--------

!----
! End
!----

return
end subroutine
