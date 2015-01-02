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

subroutine elreca &
!================

 ( dt     , propce )

!===============================================================================
! FONCTION :
! ----------

! ROUTINE PHYSIQUE PARTICULIERE POUR LE MODULE ELECTRIQUE

!             CALCULS DU COEFFICIENT DE RECALAGE
!               POUR LES VARIABLES ELECTIQUES
!             RECALAGE DES VARIABLES ELECTRIQUES
!               EN FONCTION DE CE COEFFICIENT

!-------------------------------------------------------------------------------
!ARGU                             ARGUMENTS
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use entsor
use optcal
use cstphy
use cstnum
use parall
use period
use ppppar
use ppthch
use ppincl
use elincl
use mesh
use field

!===============================================================================

implicit none

! Arguments

double precision dt(ncelet)
double precision propce(ncelet,*)

! Local variables

integer          iel    , ifac
integer          ipcefj , ipcdc3
integer          ipdcrp , idimve, idir

double precision somje , coepoa , coefav , coepot
double precision dtj   , dtjm   , delhsh , cdtj , cpmx

logical          ok

double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cvar_scalt
double precision, dimension(:), pointer :: cvar_potr, cvar_poti

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================


!===============================================================================
! 2.  ARC ELECTRIQUE
!===============================================================================


if ( ippmod(ielarc).ge.1 ) then

! 2.1 :  cas general
! ===============================

  if ( modrec .eq. 1) then

!       CALCUL DU COEFFICIENT DE RECALAGE
!       -------------------------------

!  Calcul de l'integrale sur le Volume de J.E
!     (c'est forcement positif ou nul)

    ipcefj = ipproc(iefjou)
    somje = 0.d0
    do iel = 1, ncel
      somje = somje+propce(iel,ipcefj)*volume(iel)
    enddo

    if(irangp.ge.0) then
      call parsom (somje)
    endif

    coepot = couimp*dpot/max(somje,epzero)

    coepoa = coepot

!  On impose COEPOT >= 0.75 et COEPOT <= 1.5

    if ( coepot .gt. 1.50d0 ) coepot = 1.50d0
    if ( coepot .lt. 0.75d0 ) coepot = 0.75d0

    write(nfecra,1000)coepoa,coepot

 1000     format(/,                                               &
 ' Courant impose/Courant= ',E14.5,', Coeff. recalage= ',E14.5)

  else if ( modrec .eq. 2) then

! 2.2 : 2eme exemple : Autre methode de recalage
! ==============================================
!    Ceci est un cas particulier et doit etre adapte en fonction
!    du cas et du maillage (intervenir aussi dans uselcl)

    ipcefj = ipproc(iefjou)
    call uielrc(izreca, crit_reca)

!   Calcul de l'intensite du courant d'arc
!   --------------------------------------
!      Calcul de l'integrale de J sur une surface plane
!      perpendiculaire a l'axe de l'arc

!   ATTENTION : changer la valeur des tests sur CDGFAC(3,IFAC)
!                en fonction du maillage

    ipcdc3 = ipproc(idjr(idreca))
    elcou  = 0.d0
    do ifac = 1, nfac
      if (izreca(ifac) .gt. 0) then
        ok = .true.
        do idir = 1, 3
          if (abs(surfac(idir, ifac)) .gt. 0.d0 .and. idir.ne.idreca) then
            ok = .false.
          endif
        enddo
        if (ok .eqv. .true.) then
          iel = ifacel(1,ifac)
          elcou = elcou + propce(iel,ipcdc3) * surfac(idreca,ifac)
        endif
      endif
    enddo

    if(irangp.ge.0) then
      call parsom (elcou)
    endif

    if ( abs(elcou).ge.1.d-06 ) then
      elcou=abs(elcou)
    else
      elcou=0.d0
    endif

    coepoa = 1.d0
    if(elcou.ne.0.d0) coepoa = couimp/elcou
    coepot = coepoa

    write(nfecra,*) ' ELCOU = ',ELCOU
  endif

  ! Limitation de l'evolution de l'enthalpie et recalage EM
  if ( modrec .eq. 1 .or. modrec .eq. 2 ) then

    dtj = 1.d15
    dtjm = dtj
    delhsh = 0.d0
    cdtj= 20.d0

    call field_get_val_s(icrom, crom)
    call field_get_val_s(ivarfl(isca(iscalt)), cvar_scalt)

    do iel = 1, ncel
      if (crom(iel).ne.0.d0)                               &
        delhsh =  propce(iel,ipcefj) * dt(iel) /crom(iel)

      if(delhsh.ne.0.d0) then
        dtjm = cvar_scalt(iel)/delhsh
      else
        dtjm = dtj
      endif
      dtjm=abs(dtjm)
      dtj =min(dtj,dtjm)
    enddo
    if(irangp.ge.0) then
      call parmin (dtj)
    endif
    WRITE(NFECRA,*) ' DTJ = ',DTJ

    cpmx= sqrt(cdtj*dtj)
    coepot=cpmx
    if(ntcabs.gt.2) then
      if(coepoa.ge.1.05d0) then
        coepot=cpmx
      else
        coepot=coepoa
      endif
    endif

    write(nfecra,1008) cpmx,coepoa,coepot
    write(nfecra,1009) dpot*coepot

!   RECALAGE DES VARIABLES ELECTRIQUES
!   ----------------------------------

!   Valeur de DPOT
!   --------------

    dpot = dpot*coepot

!   Potentiel Electrique (on pourrait eviter ; c'est pour le post)
!   --------------------

    call field_get_val_s(ivarfl(isca(ipotr)), cvar_potr)
    do iel = 1, ncel
      cvar_potr(iel) = cvar_potr(iel)*coepot
    enddo


!   Densite de courant (sert pour A et pour jXB)
!   ------------------

    if(ippmod(ielarc).ge.1 ) then
      do idimve = 1, ndimve
        do iel = 1, ncel
          ipdcrp = ipproc(idjr(idimve))
          propce(iel,ipdcrp) = propce(iel,ipdcrp) * coepot
        enddo
      enddo
    endif

!   Effet Joule (sert pour H au pas de temps suivant)
!   -----------

    ipcefj = ipproc(iefjou)
    do iel = 1, ncel
      propce(iel,ipcefj) = propce(iel,ipcefj)*coepot**2
    enddo
  endif
endif

!===============================================================================
! 3.  EFFET JOULE
!===============================================================================

if ( ippmod(ieljou).ge.1 ) then

! 3.1  CALCUL DU COEFFICIENT DE RECALAGE
! --------------------------------------

!  Calcul de l'integrale sur le Volume de J.E
!     (c'est forcement positif ou nul)

  ipcefj = ipproc(iefjou)
  somje = 0.d0
  do iel = 1, ncel
    somje = somje+propce(iel,ipcefj)*volume(iel)
  enddo

  if(irangp.ge.0) then
    call parsom (somje)
  endif

  coepot = sqrt(puisim/max(somje,epzero))

  coefav = coepot

!  On impose COEF >= 0.75 et COEF <= 1.5

  if ( coepot .gt. 1.50d0 ) coepot = 1.50d0
  if ( coepot .lt. 0.75d0 ) coepot = 0.75d0

  write(nfecra,2000)coefav,coejou
 2000   format(/,                                                 &
 ' Puissance impose/Somme jE= ',E14.5,', Coeff. recalage= ',E14.5)


! 3.2  RECALAGE DES VARIABLES JOULE
! ---------------------------------

!       Valeur de DPOT (au cas ou utile)
!       --------------

  dpot = dpot*coepot

!       Coefficient correcteur COEJOU cumule
!       ------------------------------------

  coejou = coejou*coepot

!       Potentiel Electrique (on pourrait eviter ; c'est pour le post)
!       --------------------

  if ( ippmod(ieljou).ne.3 .and. ippmod(ieljou).ne.4 ) then
    call field_get_val_s(ivarfl(isca(ipotr)), cvar_potr)
    do iel = 1, ncel
      cvar_potr(iel) = cvar_potr(iel)*coepot
    enddo
  endif

!      Potentiel complexe (on pourrait eviter ; c'est pour le post)
!      -----------------

  if ( ippmod(ieljou).eq.2 ) then
    call field_get_val_s(ivarfl(isca(ipoti)), cvar_poti)
    do iel = 1, ncel
      cvar_poti(iel) = cvar_poti(iel)*coepot
    enddo
  endif

!      Effet Joule (sert pour H au pas de temps suivant)
!      -----------

  ipcefj = ipproc(iefjou)
  do iel = 1, ncel
    propce(iel,ipcefj) = propce(iel,ipcefj)*coepot**2
  enddo

endif

!--------
! FORMATS
!--------

 1008  format(/,' Cpmx   = ',E14.5,/,                             &
          ' COEPOA = ',E14.5,/,                             &
          ' COEPOT = ',E14.5)

 1009  format(/,' Dpot recale     = ',E14.5)

!----
! FIN
!----

return
end subroutine
