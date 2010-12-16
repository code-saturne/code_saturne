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

subroutine laghis &
!================

 ( idbia0 , idbra0 , ndim   , ncelet , ncel ,                     &
   modhis , nvlsta ,                                              &
   ia     ,                                                       &
   xyzcen , volume , statis , stativ ,                            &
   ra     )

!===============================================================================
!  FONCTION  :
!  ---------

! ROUTINE D'ECRITURE DES HISTORIQUES POUR LE LAGRANGIEN

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! idbia0           ! i  ! <-- ! number of first free position in ia            !
! idbra0           ! i  ! <-- ! number of first free position in ra            !
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! modhis           ! e  ! <-- ! indicateur valant 0,1 ou 2                     !
!                  !    !               ! 1,2 = ecriture intermediaire, finale |
! ia(*)            ! ia ! --- ! main integer work array                        !
! xyzcen           ! tr ! <-- ! point associes aux volumes de control          !
! (ndim,ncelet)    !    !     !                                                !
! ra               ! tr !  -- ! tableau des reels                              !
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
use numvar
use entsor
use parall
use cstnum
use optcal
use lagpar
use lagran

!===============================================================================

implicit none

! Arguments

integer          idbia0, idbra0
integer          ndim, ncelet, ncel
integer          modhis, nvlsta
integer          ia(*)
double precision xyzcen(ndim,ncelet) , volume(ncelet)
double precision statis(ncelet,nvlsta)
double precision stativ(ncelet,nvlsta-1)
double precision ra(*)

! Local variables

character        nomfic*300, nenvar*300
integer          ii, ii1, ii2, lpos, inam1, inam2, lng
integer          icap,ncap,ipp,ipp2,nbpdte, jtcabs
integer          idmoyd, idebia, idebra, ifinia ,ifinra
integer          iel   , ivarl
integer          nbcap(nvppmx)
integer          ipas  , ilpd1 , il , ilfv1 , icla , ilts1
integer          iokhis, iok2

double precision xtcabs,xyztmp(3)
double precision varcap(ncaptm)
double precision dmoy

! NOMBRE DE PASSAGES DANS LA ROUTINE

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 0. INITIALISATIONS LOCALES
!===============================================================================

ipass = ipass + 1

idebia = idbia0
idebra = idbra0

! Test : Si il n'y a pas de capteur ====> On ne fait rien

if(ncapt.le.0) return

!===============================================================================
! 2. OUVERTURE DU FICHIER DE STOCKAGE histla.tmp
!===============================================================================

if(ipass.eq.1 .and. irangp.le.0) then
  NOMFIC = ' '
  nomfic = emphis
  call verlon ( nomfic,ii1,ii2,lpos)
  !==========

  NOMFIC(II2+1:II2+10) = 'histla.tmp'
  ii2 = ii2+10
  open ( unit=impli1, file=nomfic (ii1:ii2),                      &
         STATUS='UNKNOWN', FORM='UNFORMATTED',                    &
         ACCESS='SEQUENTIAL')
endif

!===============================================================================
! 3. ECRITURE DES RESULTATS dans le FICHIER DE STOCKAGE
!===============================================================================

if(modhis.eq.0.or.modhis.eq.1) then

  ifinia = idebia
  idmoyd = idebra
  ifinra = idmoyd + ncelet
  call rasize('laghis',ifinra)
  !==========

  do ipas  = 1,1+nbclst

!   Moyenne

  do il = 1,nvlsta

    ivarl = (ipas-1)*nvlsta+il
    ilpd1 = (ipas-1)*nvlsta+ilpd
    ilfv1 = (ipas-1)*nvlsta+ilfv
    icla  = ipas -1

! Pour l'instant on fait des chrono sur toutes les variables Stat. Lag.
! et sur tout les capteurs

   if ( ihslag(ivarl).ge. 1 ) then

     if ( ivarl .ne. ilpd1 .and. ivarl .ne. ilfv1 ) then
       do iel=1,ncel
         if ( statis(iel,ilpd1) .gt. seuil ) then
           ra(idmoyd+iel-1) = statis(iel,ivarl)                   &
                             /statis(iel,ilpd1)
         else
           ra(idmoyd+iel-1) = 0.d0
         endif
       enddo

     else if (ivarl.eq.ilpd1) then
       do iel=1,ncel
         ra(idmoyd+iel-1) = statis(iel,ivarl)
       enddo
     else
       do iel=1,ncel
         if (npst.gt.0) then
           ra(idmoyd+iel-1) = statis(iel,ivarl)                   &
                             /(dble(npst)*volume(iel))
         else
           ra(idmoyd+iel-1) = 0.d0
         endif
       enddo
     endif

     do icap = 1, ncapt
       if (irangp.lt.0) then
         varcap(icap) = ra(idmoyd+nodcap(icap)-1)
       else
         call parhis(nodcap(icap), ndrcap(icap),                  &
         !==========
                     ra(idmoyd), varcap(icap))
       endif
     enddo
     ncap = ncapt

     if (irangp.le.0) then
       write(impli1) ntcabs, ttcabs, (varcap(icap),               &
                                           icap=1,ncap)


     endif

   endif
  enddo

!   Variance

  do il = 1,nvlsta-1

   ivarl = (ipas-1)*nvlsta+il
   ilpd1 = (ipas-1)*nvlsta+ilpd
   ilfv1 = (ipas-1)*nvlsta+ilfv
   ilts1 = (ipas-1)*nvlsta+ilts
   icla  = ipas -1

! Pour l'instant on fait des chrono sur toutes les variables Stat. Lag.
! et sur tout les capteurs

   if ( ihslag(ivarl).eq. 2 ) then
     do iel = 1, ncel

      if ( ivarl.ne.ilfv ) then
        if ( statis(iel,ilpd1).gt.seuil ) then
          ra(idmoyd+iel-1) = stativ(iel,ivarl)/statis(iel,ilpd1)  &
                       -( statis(iel,ivarl)/statis(iel,ilpd1)     &
                         *statis(iel,ivarl)/statis(iel,ilpd1) )
        else
          ra(idmoyd+iel-1) = zero
        endif
      else
        if ( statis(iel,ilpd1).gt.seuil .and. npst.gt.0 ) then
          dmoy = statis(iel,ivarl)                                &
                /(dble(npst)*volume(iel))
          ra(idmoyd+iel-1) = stativ(iel,ivarl)                    &
                             /( dble(npst) * volume(iel) )        &
                             -dmoy*dmoy

        else if ( statis(iel,ilpd1).gt.seuil .and.                &
                  iplas.ge.idstnt                  ) then
          dmoy =  statis(iel,ivarl) / volume(iel)
          ra(idmoyd+iel-1) = stativ(iel,ilfv) / volume(iel)       &
                            -dmoy*dmoy
        else
          ra(idmoyd+iel-1) = zero
        endif
      endif
      ra(idmoyd+iel-1) = sqrt( max(zero,ra(idmoyd+iel-1)))
    enddo

    do icap = 1, ncapt
      if (irangp.lt.0) then
        varcap(icap) = ra(idmoyd+nodcap(icap)-1)
      else
        call parhis(nodcap(icap), ndrcap(icap),                   &
        !==========
                    ra(idmoyd), varcap(icap))
      endif
    enddo
    ncap = ncapt

    if (irangp.le.0) then
      write(impli1) ntcabs, ttcabs, (varcap(icap),                &
                                         icap=1,ncap)
    endif

   endif
  enddo

  enddo

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

!       --> nombre de capteur par variable
  do ipp = 1, 2*nvlsta
    nbcap(ipp) = ncapt
  enddo

!       --> ecriture un fichier par variable


  do ipas  = 1,1+nbclst

   do ipp = 1, 2*nvlsta

    if ( ipp .le. nvlsta) then
      ivarl = (ipas-1)*nvlsta+il
      ilpd1 = (ipas-1)*nvlsta+ilpd
      ilfv1 = (ipas-1)*nvlsta+ilfv
      ilts1 = (ipas-1)*nvlsta+ilts
      icla  = ipas -1
    else
      ivarl = (ipas-1)*nvlsta+(il-nvlsta)
      ilpd1 = (ipas-1)*nvlsta+ilpd
      ilfv1 = (ipas-1)*nvlsta+ilfv
      ilts1 = (ipas-1)*nvlsta+ilts
      icla  = ipas -1
    endif

    iokhis = 0
    if (ipp.le.nvlsta) then
      if (ihslag(ipp).ge.1) iokhis = 1
    else
      if ((ipp-nvlsta).ne.ilpd .and. ihslag(ipp-nvlsta).eq.2) iokhis = 1
    endif

    if (iokhis.eq.1) then

      if(irangp.le.0) then
!           --> nom du fichier
        NOMFIC = ' '
        nomfic = emphis
        call verlon (nomfic,ii1,ii2,lpos)
        !==========

        if ( ipas.eq.1 ) then
          if ( ipp.le.nvlsta ) then
            nenvar = nomlag(ipp)
          else
            nenvar = nomlav(ipp-nvlsta)
          endif
        else
          if ( ipp.le.nvlsta ) then
            WRITE(NENVAR,'(A8,A4,I3)')                            &
                     NOMLAG(IPP),'_grp',ICLA
          else
            WRITE(NENVAR,'(A8,A4,I3)')                            &
                     NOMLAV(IPP-NVLSTA),'_grp',ICLA
          endif
        endif
        call verlon(nenvar,inam1,inam2,lpos)
        !==========
        call undscr(inam1,inam2,nenvar)
        !==========
        NOMFIC(II2+1:II2+4)='Lag_'
        nomfic(ii2+4+1:ii2+4+lpos) = nenvar(inam1:inam2)
        ii2 = ii2+4+lpos
        NOMFIC(II2+1:II2+1) = '.'
        ii2 = ii2+1
        nenvar = exthis
        call verlon(nenvar,inam1,inam2,lpos)
        !==========
        call undscr(inam1,inam2,nenvar)
        !==========
        nomfic(ii2+1:ii2+lpos) = nenvar(inam1:inam2)
        ii2 = ii2+lpos
!             --> ouverture
        open ( unit=impli2, file=nomfic (ii1:ii2),                &
               STATUS='UNKNOWN', FORM='FORMATTED',                &
               ACCESS='SEQUENTIAL')
!             --> entete
        write(impli2,100)
        write(impli2,101)
        write(impli2,102) nomlag(ipp)
        write(impli2,100)
        write(impli2,103)
        write(impli2,104)
        write(impli2,103)
      endif

      do ii=1,ncapt
        if (irangp.lt.0 .or.                                      &
          irangp.eq.ndrcap(ii)) then
          xyztmp(1) = xyzcen(1,nodcap(ii))
          xyztmp(2) = xyzcen(2,nodcap(ii))
          xyztmp(3) = xyzcen(3,nodcap(ii))
        endif
        if (irangp.ge.0) then
          lng = 3
          call parbcr(ndrcap(ii), lng , xyztmp)
            !==========
        endif
        if(irangp.le.0) then
          write(impli2,105) ii,                                   &
                            xyztmp(1), xyztmp(2), xyztmp(3)
        endif
      enddo

      if(irangp.le.0) then

        write(impli2,103)
        write(impli2,106) nbpdte
        write(impli2,103)

        write(impli2,103)
        write(impli2,107)
        write(impli2,103)

        write(impli2,100)
        write(impli2,103)

!             --> boucle sur les differents enregistrements
!                et les variables
        rewind(impli1)
        do ii = 1, nbpdte
          do ipp2 = 1, 2*nvlsta

            iok2 = 0
            if (ipp2.le.nvlsta) then
              if (ihslag(ipp2).ge.1) iok2 = 1
            else
              if (ihslag(ipp2-nvlsta).eq.2) iok2 = 1
            endif

            if (iok2.eq.1) then

              read(impli1) jtcabs, xtcabs, (varcap(icap),icap=1,nbcap(ipp2))

              if (ipp2.eq.ipp) then
                write(impli2,1000) &
                     jtcabs, xtcabs, (varcap(icap),icap=1,nbcap(ipp))
              endif

            endif

          enddo
        enddo

!         --> fermeture fichier
        close(impli2)

      endif

    endif

  enddo

  enddo

endif

!===============================================================================
! 5. EN CAS DE SAUVEGARDE FINALE, DESTRUCTION DU TMP
!===============================================================================

!MO      IF(MODHIS.EQ.2) THEN
!MO        NOMFIC = ' '
!MO        NOMFIC = EMPHIS
!MO        CALL VERLON ( NOMFIC,II1,II2,LPOS)
!MOC !==========
!MO        NOMFIC(II2+1:II2+8) = 'hist.tmp'
!MO        II2 = II2+8
!MO        LPOS = LPOS+8
!MOC
!MO        NENVAR = ' '
!MO        NENVAR = 'rm '
!MO        NENVAR(4:4+LPOS) = NOMFIC ( II1:II2 )
!MO        CALL SYSTEM(NENVAR)
!MOC !==========
!MO      ENDIF

!===============================================================================
! 6. AFFICHAGES
!===============================================================================

 100  FORMAT ('# ---------------------------------------------------')
 101  FORMAT ('#      FICHIER HISTORIQUE EN TEMPS')
 102  FORMAT ('#      VARIABLE    ',A16)
 103  FORMAT ('# ')
 104  FORMAT ('#      POSITION DES CAPTEURS (colonne)')
 105  FORMAT ('# ',I6,')',3(1X,E14.7))
 106  FORMAT ('#      NOMBRE D''ENREGISTREMENTS :',I7)
 107  format (                                                          &
'# COLONNE 1       : NUMERO DU PAS DE TEMPS ',/,            &
'#         2       : TEMPS PHYSIQUE (ou No pas de temps*DTREF ',/,&
'#                               en pas de temps non uniforme)',/,&
'#         3 A 100 : VALEUR AUX CAPTEURS')
 1000 format ( 1(1x,i7,1x),101(1x,e14.7))

return

!----
! FIN
!----

end subroutine
