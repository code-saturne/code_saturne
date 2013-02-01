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

subroutine enswaf &
!================

 ( nbpmax , nvp    , nvp1   , nvep   , nivep  ,                   &
   nfin   ,                                                       &
   itepa  ,                                                       &
   ettp   , tepa   )

!===============================================================================
! FONCTION :
! ----------

!   SOUS-PROGRAMME DU MODULE LAGRANGIEN :
!   -------------------------------------

!   Ecriture des fichiers pour Ensight7 au format CASE pour la
!   visualisation des deplacements des particules et de variables
!   associees.

!   La visualisation des deplacement et le choix des variables
!   associees est realise dans le sous-programme USLAG1.

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbpmax           ! e  ! <-- ! nombre max de particulies autorise             !
! nvp              ! e  ! <-- ! nombre de variables particulaires              !
! nvp1             ! e  ! <-- ! nvp sans position, vfluide, vpart              !
! nvep             ! e  ! <-- ! nombre info particulaires (reels)              !
! nivep            ! e  ! <-- ! nombre info particulaires (entiers)            !
! nfin             ! e  ! <-- ! nfin = 1 si dernier pas de temps               !
!                  !    !     ! nfin = 0 sinon                                 !
! itepa            ! te ! <-- ! info particulaires (entiers)                   !
! (nbpmax,nivep    !    !     !   (cellule de la particule,...)                !
! ettp             ! tr ! <-- ! tableaux des variables liees                   !
!  (nbpmax,nvp)    !    !     !   aux particules                               !
!                  !    !     !   etape courante ou precedente                 !
! tepa             ! tr ! <-- ! info particulaires (reels)                     !
! (nbpmax,nvep)    !    !     !   (poids statistiques,...)                     !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!==============================================================================
! Module files
!==============================================================================

use paramx
use entsor
use lagpar
use lagran

!==============================================================================

implicit none

! Arguments

integer          nbpmax , nvp    , nvp1   , nvep  , nivep
integer          nfin
integer          itepa(nbpmax,nivep)

double precision ettp(nbpmax,nvp)
double precision tepa(nbpmax,nvep)

! Local variables

integer          npt , ipt
integer          np , nl
integer          ii1 , ii2 , lpos , n1 , n2

character        fich*80 , name*80 , entet*80

double precision, allocatable, dimension(:,:) :: trav

integer          ipwaf
data             ipwaf /0/
save             ipwaf

!===============================================================================

!===============================================================================
! 0. GESTION MEMOIRE
!===============================================================================

! Allocate a work array
allocate(trav(nbpmax,3))

!===============================================================================
! 1. Initialisations
!===============================================================================

if (nfin.eq.0) ipwaf = ipwaf + 1

if (ipwaf.eq.1) itlag = 0

FICH = ' '
FICH = 'Displacement'
call verlon (fich,ii1,ii2,lpos)
entet = fich(ii1:ii2)

!===============================================================================
! 2. ENREGISTREMENTS des deplacement.geom====
!===============================================================================

!-->Faut-il enregistrer ?

if ( (mod(ipwaf-1,nvisla).eq.0 .and. nfin.eq.0)   .or.            &
     (nfin.eq.1 .and. mod(ipwaf-1,nvisla).ne.0) ) then

!-->Nombre de particules a visualisees encore presentent dans le domaine

  npt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) npt = npt + 1
    endif
  enddo

!-->Y a t-il encore des particules a visualiser ?

  if (npt.eq.0) goto 100

!-->Nombre d'enregistrements et incrementation du temps physique

  if (itlag.le.9999) then
    itlag = itlag + 1
    timlag(itlag) = ttclag
  else
    write(nfecra,9000) itlag
    goto 100
  endif

!-->Ouverture des fichiers type deplacement.geo0001

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+4) = '.geo'
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  ii2 = ii2 + 5
  fich(ii2:ii2+lpos) = name(n1:n2)

  ii2 = ii2 + lpos
  open ( impla1, file=fich(ii1:ii2),                              &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture de l'entete

  WRITE(IMPLA1,'(A)') 'geometrie deplacement'
  WRITE(IMPLA1,'(A)') 'au format ensight6'
  WRITE(IMPLA1,'(A)') 'node id given'
  WRITE(IMPLA1,'(A)') 'element id given'
  WRITE(IMPLA1,'(A)') 'coordinates'
  WRITE(IMPLA1,'(I8)') NPT

!-->Ecriture des points

  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        write(impla1,'(I8,3E12.5)') &
                  np,               &
        ettp(np,jxp),               &
        ettp(np,jyp),               &
        ettp(np,jzp)
      endif
    endif
  enddo

!-->Ecriture de la geometrie Ensight

  WRITE(IMPLA1,'(A)') 'part   1'
  WRITE(IMPLA1,'(A)') 'Displacements'
  WRITE(IMPLA1,'(A)') 'point'
  WRITE(IMPLA1,'(I8)') NPT

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        write(impla1,'(2I8)') ipt , np
      endif
    endif
  enddo
  close(impla1)

else

  if (nfin.eq.0) return
  goto 100

endif

!===============================================================================
! 3. Ecriture de deplacement.tpssej0001
!===============================================================================

if (ivistp.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = tepa(np,jrtsp)
      endif
    endif
 enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+7) = '.tpssej'
  ii2 = ii2 + 7
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
  close(impla1)

endif

!===============================================================================
! 4. Ecriture de deplacement.temper0001
!===============================================================================

if (iviste.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = ettp(np,jtp)
      endif
    endif
  enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+7) = '.temper'
  ii2 = ii2 + 7
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
  close(impla1)

endif

!===============================================================================
! 5. Ecriture de deplacement.diamet0001
!===============================================================================

if (ivisdm.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = ettp(np,jdp)
      endif
    endif
  enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+7) = '.diamet'
  ii2 = ii2 + 7
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
  close(impla1)

endif

!===============================================================================
! 6. Ecriture de deplacement.massep0001
!===============================================================================

if (ivismp.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = ettp(np,jmp)
      endif
    endif
  enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+7) = '.massep'
  ii2 = ii2 + 7
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
  close(impla1)

endif

!===============================================================================
! 7. Charbon : Ecriture de deplacement.temp_ch0001
!===============================================================================

if (ivishp.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = ettp(np,jhp)
      endif
    endif
  enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+7) = '.tempch'
  ii2 = ii2 + 7
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
  close(impla1)

endif

!===============================================================================
! 8. Charbon : Ecriture de deplacement.dck0001
!===============================================================================

if (ivisdk.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = tepa(np,jrdck)
      endif
    endif
  enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+4) = '.dck'
  ii2 = ii2 + 4
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
  close(impla1)

endif

!===============================================================================
! 9. Charbon : Ecriture de deplacement.mch0001
!===============================================================================

if (ivisch.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = ettp(np,jmch)
      endif
    endif
  enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+4) = '.mch'
  ii2 = ii2 + 4
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
  close(impla1)

endif

!===============================================================================
! 10. Charbon : Ecriture de deplacement.mck0001
!===============================================================================

if (ivisck.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = ettp(np,jmck)
      endif
    endif
  enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+4) = '.mck'
  ii2 = ii2 + 4
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( REAL(TRAV(NP,1)), NP=1,NPT )
  close(impla1)

endif

!===============================================================================
! 11. Ecriture de deplacement.vitflu0001
!===============================================================================

if (ivisv1.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = ettp(np,juf)
        trav(ipt,2) = ettp(np,jvf)
        trav(ipt,3) = ettp(np,jwf)
      endif
    endif
  enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+7) = '.vitflu'
  ii2 = ii2 + 7
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( (REAL(TRAV(NP,NL)),NL=1,3),NP=1,NPT )
  close(impla1)

endif

!===============================================================================
! 12. Ecriture de deplacement.vitpar0001
!===============================================================================

if (ivisv2.eq.1) then

  ipt = 0
  do nl = 1,nbvis
    np = liste(nl)
    if (np.ge.1) then
      if (itepa(np,jisor).ne.0) then
        ipt = ipt + 1
        trav(ipt,1) = ettp(np,jup)
        trav(ipt,2) = ettp(np,jvp)
        trav(ipt,3) = ettp(np,jwp)
      endif
    endif
  enddo

  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+7) = '.vitpar'
  ii2 = ii2 + 7
  WRITE (NAME,'(I4.4)') ITLAG
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)

  open ( impla1, file=fich(ii1:ii2+lpos),                         &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )

!-->Ecriture

  WRITE(IMPLA1,'(A)') FICH(II1:II2+LPOS)
  WRITE(IMPLA1,'(6E12.5)') ( (REAL(TRAV(NP,NL)),NL=1,3),NP=1,NPT )
  close(impla1)

endif

! Free memory
deallocate(trav)

!===============================================================================
! 13. Ecriture du deplacement.case au dernier passage
!===============================================================================

 100  continue


  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  NAME = ' '
  NAME = '.CASE'
  call verlon (name,n1,n2,lpos)
  fich(ii2+1:ii2+lpos) = name(n1:n2)
  ii2 = ii2 + lpos
  open ( unit=impla1, file=fich (ii1:ii2),                        &
         STATUS='UNKNOWN', FORM='FORMATTED', ACCESS='SEQUENTIAL' )
  rewind ( unit=impla1 )

  WRITE(IMPLA1,'(A)') 'FORMAT'
  WRITE(IMPLA1,'(A)') 'type:     ensight'
  WRITE(IMPLA1,'(A)') 'GEOMETRY'
  FICH = ' '
  fich = entet
  call verlon (fich,ii1,ii2,lpos)
  FICH(II2+1:II2+8) = '.geo****'
  call verlon (fich,ii1,ii2,lpos)
  NAME = ' '
  NAME = 'model:                 1    '
  name(29:29+ii2-ii1+1) = fich(ii1:ii2)
  call verlon (name,ii1,ii2,lpos)
  WRITE(IMPLA1,'(A)') NAME(II1:II2)

  WRITE(IMPLA1,'(A)') 'VARIABLE'
  FICH = ' '
  fich = entet
  call verlon (fich,n1,n2,lpos)

! Rem : les trois lignes suivantes sont pour eviter une erreur
!       de lecture fichier .CASE lors de sa lecture par ensight
!       s'il n'y a aucune VARIABLE a voir.

!       NAME = 'constant per case :    1    constant     1.0'
!       CALL VERLON (NAME,II1,II2,LPOS)
!       WRITE(IMPLA1,'(A)') NAME(II1:II2)

  if (ivistp.eq.1) then
    NAME = 'scalar per node:      1    Dis_Resid_Time  '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 2
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.tpssej****'
    ii2 = ii2 + 11
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  if (iviste.eq.1) then
    NAME = 'scalar per node:      1    Dis_Temp      '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 6
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.temper****'
    ii2 = ii2 + 11
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  if (ivisdm.eq.1) then
    NAME = 'scalar per node:      1    Dis_Diameter         '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 9
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.diamet****'
    ii2 = ii2 + 11
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  if (ivismp.eq.1) then
    NAME = 'scalar per node:      1    Dis_Mass            '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 12
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.massep****'
    ii2 = ii2 + 11
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  if (ivishp.eq.1) then
    NAME = 'scalar per node:      1    tempch           '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 11
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.tempch****'
    ii2 = ii2 + 11
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  if (ivisdk.eq.1) then
    NAME = 'scalar per node:      1    dck              '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 14
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.dck****'
    ii2 = ii2 + 8
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  if (ivisch.eq.1) then
    NAME = 'scalar per node:      1    mch              '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 14
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.mch****'
    ii2 = ii2 + 8
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  if (ivisck.eq.1) then
    NAME = 'scalar per node:      1    mck              '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 14
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.mck****'
    ii2 = ii2 + 8
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  if (ivisv1.eq.1) then
    NAME = 'vector per node:      1    Dis_Fluid_Velo   '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 3
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.vitflu****'
    ii2 = ii2 + 11
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  if (ivisv2.eq.1) then
    NAME = 'vector per node:      1    Dis_Part_Velo   '
    call verlon (name,ii1,ii2,lpos)
    ii2 = ii2 + 3
    name(ii2+1:ii2+n2)=fich(n1:n2)
    ii2 = ii2 + n2
    call verlon (name,ii1,ii2,lpos)
    NAME(II2+1:II2+11) = '.vitpar****'
    ii2 = ii2 + 11
    WRITE(IMPLA1,'(A)') NAME(II1:II2)
  endif

  WRITE(IMPLA1,'(A)') 'TIME'
  WRITE(IMPLA1,'(A)') 'time set:               1'
  FICH = ' '
  FICH = 'number of steps:'
  WRITE(NAME,'(I4)') ITLAG
  call verlon(name,n1,n2,lpos)
  !==========
  fich(25+1:25+lpos) = name(n1:n2)
  WRITE(IMPLA1,'(A)') FICH(1:25+LPOS)
  WRITE(IMPLA1,'(A)') 'filename start number:  1'
  WRITE(IMPLA1,'(A)') 'filename increment:     1'
  WRITE(IMPLA1,'(A)') 'time values:'
  do nl = 1, itlag
    WRITE(IMPLA1,'(E12.5)') TIMLAG(NL)
  enddo

  close(impla1)


return

!-------
! FORMAT
!-------

 9000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A L''EXECUTION DU MODULE LAGRANGIEN  ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    LE NOMBRE D''ENREGISTREMENTS TEMPORELS DEMANDES POUR    ',/,&
'@      LE POST-PROCESSING EN MODE DEPLACEMENT DEPASSE        ',/,&
'@      LE MAXIMUM ADMISSIBLE.                                ',/,&
'@                                                            ',/,&
'@    LE NOMBRE DE PAS DE TEMPS DEMANDE EST DE : ',I10         ,/,&
'@      LE MAXIMUM ADMISSIBLE EST 9999                        ',/,&
'@                                                            ',/,&
'@  Le calcul continue, mais les enregistrements sont arretes.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

end subroutine
