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

subroutine enslag &
!================

 ( idbia0 , idbra0 ,                                              &
   nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   nfin   , iforce ,                                              &
   itepa  ,                                                       &
   ettp   , tepa   , ra)

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!    SAUVEGARDE DES TRAJECTOIRES AU FORMAT ENSIGHT GOLD

!    1. AU PREMIER PASSAGE : OUVERTURE  D'UN FICHIER
!       SCRATCH ET DEBUT D'ECRITURE SEQUENTIELLE

!    2. A CHAQUE PASSAGE : ECRITURE DANS LE TEMPORAIRE

!    3. AU DERNIER PASSAGE : ECRITURE SUR LE FICHIER FINAL
!       AU BON FORMAT POUR ENSIGHT

!    REMARQUE : on dispose de 15 fichiers de sortie,
!               les unites logiques sont defini dans INIINI.F
!               et on y accede par IMPLA5(1) a IMPLA5(15)

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! nfin             ! e  ! <-- ! nfin = 1 si dernier pas de temps               !
!                  !    !     ! nfin = 0 sinon                                 !
! iforce           ! e  ! <-- ! force l'ecriture si = numero de la             !
!                  !    !     !   particule courante                           !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax, nivep   !    !     !   (cellule de la particule, ...)               !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax, nvp)   !    !     !   aux particules                               !
!                  !    !     !   etape courante ou precedente                 !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax, nvep)   !    !     !   (poids statistiques, ...)                    !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail

!===============================================================================

implicit none

!==============================================================================
! Common blocks
!==============================================================================

include "paramx.h"
include "entsor.h"
include "lagpar.h"
include "lagran.h"

!==============================================================================

! Arguments

integer          idbia0 , idbra0
integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          nfin   , iforce
integer          itepa(nbpmax, nivep)
double precision ettp(nbpmax, nvp) , tepa(nbpmax, nvep)
double precision ra(*)

! Local variables

integer          idebia , idebra
integer          nl , np

integer          numl
integer          nume(nliste)
integer          ii1 , ii2 , lpos , nlmax , ios
integer          npt , ipt , lmax
integer          ix, iy, iz, ii
integer          iu1l , iv1l  , iw1l
integer          iu2l , iv2l  , iw2l
integer          itpl , idml  , itel  , impl
integer          ihpl , idckl , imchl , imckl
integer          ifinia, ifinra

double precision xpl , ypl , zpl
double precision u1l , v1l , w1l
double precision u2l , v2l , w2l
double precision tpl , dml , tel , mpl
double precision hpl , dckl , mchl , mckl
character        fich*80
character        name*80

integer ipass
data    ipass /0/
save    ipass

!==============================================================================

!===============================================================================
! -1.  GESTION MEMOIRE
!===============================================================================

idebia = idbia0
idebra = idbra0

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

!--> Compteur de passages

ipass = ipass+1

nlmax = 0

!===============================================================================
! 1. OUVERTURE DU FICHIER TAMPON AU PREMIER PASSAGE
!===============================================================================

if (ipass.eq.1) then

  open(impla3, file='scratch3.lag',                               &
       status='unknown', form='unformatted', access='sequential')

  do nl = 1, nbvis
    nplist(nl) = 0
    list0(nl) = liste(nl)
  enddo

endif

!===============================================================================
! 2. ECRITURE SEQUENTIELLE DES INFOS TRAJECTOIRES
!    POUR CHAQUE PARTICULE A VISUALISER SUIVANT LA FREQUENCE
!===============================================================================

if ((mod(ipass-1, nvisla).eq.0 .or. iforce.gt.0) .and. nfin.eq.0) then

  do nl = 1, nbvis

    np = liste(nl)

    if ((np.ge.1 .and. iforce.eq.0) .or. (iforce.eq.np)) then

      ! sortie du domaine ?
      if (itepa(np, jisor).gt.0) then

        !--->incrementation du nombre d'enregistrement pour la particule NP :
        nplist(nl) = nplist(nl)+1

        if (nplist(nl).gt.nlmax) nlmax = nplist(nl)

        !--->numero de liste :
        write(impla3) nl

        !--->coordonnees de la particule NP :
        write(impla3) ettp(np, jxp), ettp(np, jyp), ettp(np, jzp)

        !--->vitesse du fluide vu :
        if (ivisv1.eq.1) then
          write(impla3) ettp(np, juf), ettp(np, jvf), ettp(np, jwf)
        endif

        !--->vitesse de la particule :
        if (ivisv2.eq.1) then
          write(impla3) ettp(np, jup), ettp(np, jvp), ettp(np, jwp)
        endif

        !--->temps de sejour :
        if (ivistp.eq.1) then
          write(impla3) tepa(np, jrtsp)
        endif

        !--->diametre :
        if (ivisdm.eq.1) then
            write(impla3) ettp(np, jdp)
        endif

        !--->masse :
        if (ivismp.eq.1) then
          write(impla3) ettp(np, jmp)
        endif

        !--->temperature :
        if (iviste.eq.1) then
          write(impla3) ettp(np, jtp)
        endif

        !--->Specifique charbon :
        ! Temperature
        if (ivishp.eq.1) then
          write(impla3) ettp(np, jhp)
        endif
        ! Diametre du coeur retrecisant
        if (ivisdk.eq.1) then
          write(impla3) tepa(np, jrdck)
        endif
        ! Masse charbon reactif
        if (ivisch.eq.1) then
          write(impla3) ettp(np, jmch)
        endif
        ! Masse de coke
        if (ivisck.eq.1) then
          write(impla3) ettp(np, jmck)
        endif

      endif

    endif

  enddo

endif

!===============================================================================
! 3.1 FIN DU CALCUL (NFIN=1) : OUVERTURE DES FICHIERS RESU
!===============================================================================

if (nfin.eq.1) then

  NAME = ' '
  NAME = 'trajectoire'

  !  0) ouverture du fichier .ensight.CASE :

  fich = name
  call verlon(fich, ii1, ii2, lpos)
  fich(ii2+1:ii2+14) = '.ensight.CASE'
  ii2 = ii2 + 14
  open(unit=impla2, file=fich (ii1:ii2),                        &
       status='unknown', form='formatted', access='sequential', &
       iostat=ios, err=99)
  rewind(unit=impla2, err=99)

  rewind(impla2)
  write(impla2, 5010)
  write(impla2, 5011)
  write(impla2, 5012)

  !  1) ouverture du fichier .ensight.geom + entete fichier case(suite)

  fich = name
  call verlon(fich, ii1, ii2, lpos)
  fich(ii2+1:ii2+13) = '.ensight.geom'
  ii2 = ii2 + 13

  write(impla2, 5013) fich (ii1:ii2)
  write(impla2, 5014)

  open(unit=impla1, file=fich (ii1:ii2),                        &
       status='unknown', form='formatted', access='sequential', &
       iostat=ios, err=99)

  rewind(unit=impla1, err=99)

  !  2) ouverture du fichier .vitflu + entete fichier case(suite)

  if (ivisv1.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.vitflu'
    ii2 = ii2 + 7
    open(unit=impla5(1), file=fich (ii1:ii2),                    &
        status='unknown', form='formatted', access='sequential', &
        iostat=ios, err=99)
    rewind(unit=impla5(1), err=99)

    write(impla2, 5015) fich (ii1:ii2)

  endif

  !  3) ouverture du fichier .vitpar + entete fichier case(suite)

  if (ivisv2.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.vitpar'
    ii2 = ii2 + 7
    open(unit=impla5(2), file=fich (ii1:ii2),                    &
         status='unknown', form='formatted', access='sequential', &
         iostat=ios, err=99)
    rewind(unit=impla5(2), err=99)

    write(impla2, 5016) fich (ii1:ii2)

  endif

!  4) ouverture du fichier .tpssej + entete fichier case(suite)

  if (ivistp.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.tpssej'
    ii2 = ii2 + 7
    open(unit=impla5(3), file=fich (ii1:ii2),                     &
         status='unknown', form='formatted', access='sequential', &
         iostat=ios, err=99)
    rewind(unit=impla5(3), err=99)

    write(impla2, 5017) fich (ii1:ii2)

  endif

!  5) ouverture du fichier .diamet + entete fichier case(suite)

  if (ivisdm.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.diamet'
    ii2 = ii2 + 7
    open(unit=impla5(4), file=fich (ii1:ii2),                     &
         status='unknown', form='formatted', access='sequential', &
         iostat=ios, err=99)
    rewind(unit=impla5(4), err=99)

    write(impla2, 5018) fich (ii1:ii2)

  endif

!  6) ouverture du fichier .masse + entete fichier case(suite)

  if (ivismp.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.masse'
    ii2 = ii2 + 7
    open(unit=impla5(5), file=fich (ii1:ii2),                     &
         status='unknown', form='formatted', access='sequential', &
         iostat=ios, err=99)
    rewind(unit=impla5(5), err=99)

    write(impla2, 5019) fich (ii1:ii2)

  endif

!  7) ouverture du fichier .temper + entete fichier case(suite)

  if (iviste.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.temper'
    ii2 = ii2 + 7
    open(unit=impla5(6), file=fich (ii1:ii2),                     &
         status='unknown', form='formatted', access='sequential', &
         iostat=ios, err=99)
    rewind(unit=impla5(6), err=99)

    write(impla2, 5020) fich (ii1:ii2)

  endif

!  8) ouverture du fichier .tempch + entete fichier case(suite)

  if (ivishp.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.tempch'
    ii2 = ii2 + 7
    open(unit=impla5(7), file=fich (ii1:ii2),                     &
         status='unknown', form='formatted', access='sequential', &
         iostat=ios, err=99)
    rewind(unit=impla5(7), err=99)

    write(impla2, 5021) fich (ii1:ii2)

  endif

!  9) ouverture du fichier .dck + entete fichier case(suite)

  if (ivisdk.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.dck'
    ii2 = ii2 + 7
    open(unit=impla5(8), file=fich (ii1:ii2),                     &
         status='unknown', form='formatted', access='sequential', &
         iostat=ios, err=99)
    rewind(unit=impla5(8), err=99)

    write(impla2, 5022) fich (ii1:ii2)

  endif

!  10) ouverture du fichier .mch + entete fichier case(suite)

  if (ivisch.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.mch'
    ii2 = ii2 + 7
    open(unit=impla5(9), file=fich (ii1:ii2),                     &
         status='unknown', form='formatted', access='sequential', &
         iostat=ios, err=99)
    rewind(unit=impla5(9), err=99)

    write(impla2, 5023) fich (ii1:ii2)

  endif

!  11) ouverture du fichier .mch + entete fichier case(suite)

  if (ivisck.eq.1) then
    fich = name
    call verlon(fich, ii1, ii2, lpos)
    fich(ii2+1:ii2+7) = '.mck'
    ii2 = ii2 + 7
    open(unit=impla5(10), file=fich (ii1:ii2),                    &
         status='unknown', form='formatted', access='sequential', &
         iostat=ios, err=99)
    rewind(unit=impla5(10), err=99)

    write(impla2, 5024) fich (ii1:ii2)

  endif

!===============================================================================
! 3.2 FIN DU CALCUL (NFIN=1): ECRITURE SUR FICHIERS FINAUX
!                             ET FERMETURE DU TEMPORAIRE
!                             OU SAUVEGARDE INTERMEDIAIRE
!===============================================================================

! 1) ON COMPTE LE NOMBRE TOTAL D'ENREGISTREMENTS

  npt = 0
  nlmax = 0
  lmax = 0
  do nl = 1, nbvis
    npt = npt + nplist(nl)
    lmax = max(lmax, nplist(nl))
    if (nplist(nl).gt.nlmax) nlmax = nplist(nl)
    nume(nl) = 0
  enddo

! 2) Allocation Memoire

  ifinia = idebia
  ix     = idebra
  iy     = ix+lmax
  iz     = iy+lmax
  ifinra = iz +lmax

  if (ivisv1.eq.1) then
    iu1l   = ifinra
    iv1l   = iu1l  +lmax
    iw1l   = iv1l  +lmax
    ifinra = iw1l  +lmax
  endif
  if (ivisv2.eq.1) then
    iu2l   = ifinra
    iv2l   = iu2l  +lmax
    iw2l   = iv2l  +lmax
    ifinra = iw2l  +lmax
  endif
  if (ivistp.eq.1) then
    itpl = ifinra
    ifinra = itpl +lmax
  endif
  if (ivisdm.eq.1) then
    idml = ifinra
    ifinra = idml +lmax
  endif
  if (ivismp.eq.1) then
    impl = ifinra
    ifinra = impl +lmax
  endif
  if (iviste.eq.1) then
    itel = ifinra
    ifinra = itel +lmax
  endif
  if (ivishp.eq.1) then
    ihpl = ifinra
    ifinra = ihpl +lmax
  endif
  if (ivisdk.eq.1) then
    idckl= ifinra
    ifinra = idckl +lmax
  endif
  if (ivisch.eq.1) then
    imchl = ifinra
    ifinra = imchl +lmax
  endif
  if (ivisck.eq.1) then
    imckl = ifinra
    ifinra = imckl +lmax
  endif

  call rasize('enslag', ifinra)
  !==========

! 3) ON REMPLIT LES ENTETES DES FICHIERS : geo + variable

  rewind(impla1)

  write(impla1, 3000)
  write(impla1, 3001)
  write(impla1, 3002)
  write(impla1, 3003)

  if (ivisv1.eq.1) then
    rewind(impla5(1))
    write(impla5(1), 4000)
  endif

  if (ivisv2.eq.1) then
    rewind(impla5(2))
    write(impla5(2), 4001)
  endif

  if (ivistp.eq.1) then
    rewind(impla5(3))
    write(impla5(3), 4002)
  endif

  if (ivisdm.eq.1) then
    rewind(impla5(4))
    write(impla5(4), 4003)
  endif

  if (ivismp.eq.1) then
    rewind(impla5(5))
    write(impla5(5), 4004)
  endif

  if (iviste.eq.1) then
    rewind(impla5(6))
    write(impla5(6), 4005)
  endif

  if (ivishp.eq.1) then
    rewind(impla5(7))
    write(impla5(7), 4006)
  endif

  if (ivisdk.eq.1) then
    rewind(impla5(8))
    write(impla5(8), 4007)
  endif

  if (ivisch.eq.1) then
    rewind(impla5(9))
    write(impla5(9), 4008)
 endif

  if (ivisck.eq.1) then
    rewind(impla5(10))
    write(impla5(10), 4009)
  endif

  do nl = 1, nbvis

    np = liste(nl)

    if (itepa(np, jisor).gt.0) then

      rewind(impla3)

      ipt = 0
      do ii=1, npt

        read(impla3) numl
        read(impla3) xpl, ypl, zpl
        if (ivisv1.eq.1) read(impla3) u1l, v1l, w1l
        if (ivisv2.eq.1) read(impla3) u2l, v2l, w2l
        if (ivistp.eq.1) read(impla3) tpl
        if (ivisdm.eq.1) read(impla3) dml
        if (ivismp.eq.1) read(impla3) mpl
        if (iviste.eq.1) read(impla3) tel
        if (ivishp.eq.1) read(impla3) hpl
        if (ivisdk.eq.1) read(impla3) dckl
        if (ivisch.eq.1) read(impla3) mchl
        if (ivisck.eq.1) read(impla3) mckl

        if (numl .eq. nl) then

          ipt = ipt+1

          ra(ix+ipt-1) = xpl
          ra(iy+ipt-1) = ypl
          ra(iz+ipt-1) = zpl
          if (ivisv1.eq.1) then
            ra(iu1l+ipt-1) = u1l
            ra(iv1l+ipt-1) = v1l
            ra(iw1l+ipt-1) = w1l
          endif
          if (ivisv2.eq.1) then
            ra(iu2l+ipt-1) = u2l
            ra(iv2l+ipt-1) = v2l
            ra(iw2l+ipt-1) = w2l
          endif
          if (ivistp.eq.1) then
            ra(itpl+ipt-1) = tpl
          endif
          if (ivisdm.eq.1) then
            ra(idml+ipt-1) = dml
          endif
          if (ivismp.eq.1) then
            ra(impl+ipt-1) = mpl
          endif
          if (iviste.eq.1) then
            ra(itel+ipt-1) = tel
          endif
          if (ivishp.eq.1) then
            ra(ihpl+ipt-1) = hpl
          endif
          if (ivisdk.eq.1) then
            ra(idckl+ipt-1) = dckl
          endif
          if (ivisch.eq.1) then
            ra(imchl+ipt-1) = mchl
          endif
          if (ivisck.eq.1) then
            ra(imckl+ipt-1) = mckl
          endif

        endif
      enddo

!  Ecriture Fichier Geometrie

      write(impla1, 3010)
      write(impla1, 1010) nl
      write(impla1, 3004) list0(nl)
      write(impla1, 3005)
      write(impla1, 1010) ipt
      do ii=1, ipt
        write(impla1, 1030) ra(ix+ii-1)
      enddo
      do ii=1, ipt
        write(impla1, 1030) ra(iy+ii-1)
      enddo
      do ii=1, ipt
        write(impla1, 1030) ra(iz+ii-1)
      enddo
      write(impla1, 3006)
      if (ipt.eq.0) then
        write(impla1, 1010) 0
      else
        write(impla1, 1010) ipt-1
      endif
      do ii=1, ipt-1
        write(impla1, 1020) ii, ii+1
      enddo

!  Ecriture Fichiers Variables

      if (ivisv1.eq.1) then

        write(impla5(1), 3010)
        write(impla5(1), 1010) nl
        write(impla5(1), 3005)
        do ii=1, ipt
          write(impla5(1), 1030) ra(iu1l+ii-1)
        enddo
        do ii=1, ipt
          write(impla5(1), 1030) ra(iv1l+ii-1)
        enddo
        do ii=1, ipt
          write(impla5(1), 1030) ra(iw1l+ii-1)
        enddo
      endif

      if (ivisv2.eq.1) then
        write(impla5(2), 3010)
        write(impla5(2), 1010) nl
        write(impla5(2), 3005)
        do ii=1, ipt
          write(impla5(2), 1030) ra(iu2l+ii-1)
        enddo
        do ii=1, ipt
          write(impla5(2), 1030) ra(iv2l+ii-1)
        enddo
        do ii=1, ipt
          write(impla5(2), 1030) ra(iw2l+ii-1)
        enddo
      endif

      if (ivistp.eq.1) then
        write(impla5(3), 3010)
        write(impla5(3), 1010) nl
        write(impla5(3), 3005)
        do ii=1, ipt
          write(impla5(3), 1030) ra(itpl+ii-1)
        enddo
      endif
      if (ivisdm.eq.1) then
        write(impla5(4), 3010)
        write(impla5(4), 1010) nl
        write(impla5(4), 3005)
        do ii=1, ipt
          write(impla5(4), 1030) ra(idml+ii-1)
        enddo
      endif
      if (ivismp.eq.1) then
        write(impla5(5), 3010)
        write(impla5(5), 1010) nl
        write(impla5(5), 3005)
        do ii=1, ipt
          write(impla5(5), 1030) ra(impl+ii-1)
        enddo
      endif
      if (iviste.eq.1) then
        write(impla5(6), 3010)
        write(impla5(6), 1010) nl
        write(impla5(6), 3005)
        do ii=1, ipt
          write(impla5(6), 1030) ra(itel+ii-1)
        enddo
      endif
      if (ivishp.eq.1) then
        write(impla5(7), 3010)
        write(impla5(7), 1010) nl
        write(impla5(7), 3005)
        do ii=1, ipt
          write(impla5(7), 1030) ra(ihpl+ii-1)
        enddo
      endif
      if (ivisdk.eq.1) then
        write(impla5(8), 3010)
        write(impla5(8), 1010) nl
        write(impla5(8), 3005)
        do ii=1, ipt
          write(impla5(8), 1030) ra(idckl+ii-1)
        enddo
      endif
      if (ivisch.eq.1) then
        write(impla5(9), 3010)
        write(impla5(9), 1010) nl
        write(impla5(9), 3005)
        do ii=1, ipt
          write(impla5(9), 1030) ra(imchl+ii-1)
        enddo
      endif
      if (ivisck.eq.1) then
        write(impla5(10), 3010)
        write(impla5(10), 1010) nl
        write(impla5(10), 3005)
        do ii=1, ipt
          write(impla5(10), 1030) ra(imckl+ii-1)
        enddo
      endif

    endif

  enddo

  close(impla1)
  close(impla2)
  close(impla3)
  if (ivisv1.eq.1) close(impla5(1))
  if (ivisv2.eq.1) close(impla5(2))
  if (ivistp.eq.1) close(impla5(3))
  if (ivisdm.eq.1) close(impla5(4))
  if (ivismp.eq.1) close(impla5(5))
  if (iviste.eq.1) close(impla5(6))
  if (ivishp.eq.1) close(impla5(7))
  if (ivisdk.eq.1) close(impla5(8))
  if (ivisch.eq.1) close(impla5(9))
  if (ivisck.eq.1) close(impla5(10))

endif


return

   99 continue
write (nfecra, 9999) fich (ii1:ii2), ios
call csexit(1)

!--------
! FORMATS
!--------

 1010 format (i10)
 1020 format (i10, i10)
 1030 format (e12.5)

 3000 format('geometrie trajectoire')
 3001 format('au format ensight6 : .case')
 3002 format('node id assign')
 3003 format('element id assign')
 3004 format('trajectoire', I10)
 3005 format('coordinates')
 3006 format('bar2')

 3010 format('part')

 4000 format('vitesse fluide vu')
 4001 format('vitesse particules')
 4002 format('temps de sejour')
 4003 format('diametre')
 4004 format('masse')
 4005 format('temerature')
 4006 format('temperature')
 4007 format('dck')
 4008 format('mch')
 4009 format('mck')

 5010 format('FORMAT')
 5011 format('type: ensight gold')
 5012 format('GEOMETRY')
 5013 format('model: ', A)
 5014 format('VARIABLE')
 5015 format('vector per node: vitesse_fluide_vu       ', A)
 5016 format('vector per node: vitesse_particules      ', A)
 5017 format('scalar per node: temps_de_sejour         ', A)
 5018 format('scalar per node: diametre                ', A)
 5019 format('scalar per node: masse                   ', A)
 5020 format('scalar per node: temperature             ', A)
 5021 format('scalar per node: temperature             ', A)
 5022 format('scalar per node: dck                     ', A)
 5023 format('scalar per node: mch                     ', A)
 5024 format('scalar per node: mck                     ', A)

 9999 format(                                                       &
'@                                                            ', /, &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@', /, &
'@                                                            ', /, &
'@ @@ ATTENTION : ARRET A L''EXECUTION DU MODULE LAGRANGIEN   ', /, &
'@    =========                                               ', /, &
'@    ERREUR D''OUVERTURE SUR LE FICHIER : ', A                , /, &
'@    AVEC UN IOSTAT EGAL A : ', I6                            , /, &
'@    (ENSLAG)                                                ', /, &
'@                                                            ', /, &
'@  Verifier les numero de fichiers utilises par le Lagrangien', /, &
'@                                                            ', /, &
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@', /, &
'@                                                            ', /)

!----
! FIN
!----

end subroutine
