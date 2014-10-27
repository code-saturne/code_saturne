!-------------------------------------------------------------------------------

!VERS

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

subroutine uselrc &
!================

 ( nvar   , nscal  ,                                              &
   dt     )

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
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! itypsm           ! te ! <-- ! type de source de masse pour les               !
! (ncesmp,nvar)    !    !     !  variables (cf. cs_user_mass_source_terms)     !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! smacel           ! tr ! <-- ! valeur des variables associee a la             !
! (ncesmp,*   )    !    !     !  source de masse                               !
!                  !    !     !  pour ivar=ipr, smacel=flux de masse           !
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

integer          nvar   , nscal

double precision dt(ncelet)

! Local variables

integer          iel    , ifac   , iutile
integer          idimve , jaiex  , ifcsig

double precision somje , coepoa , coefav , coepot
double precision emax  , aiex   , amex
double precision rayo  , econs  , z1     , z2   , posi
double precision dtj   , dtjm   , delhsh , cdtj , cpmx
double precision xelec , yelec  , zelec, diff

double precision, allocatable, dimension(:) :: w1
double precision, dimension(:), pointer :: crom
double precision, dimension(:), pointer :: cscalt, cpotr
double precision, dimension(:), pointer :: cpro_sig, efjou, djr3
double precision, dimension(:,:), pointer :: djr

!===============================================================================

! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_START
!===============================================================================

if (1.eq.1) return

!===============================================================================
! TEST_TO_REMOVE_FOR_USE_OF_SUBROUTINE_END

!===============================================================================

!===============================================================================
! 1. INITIALISATION
!===============================================================================

call field_get_val_s(icrom, crom)
call field_get_val_s(ivarfl(isca(iscalt)), cscalt)
call field_get_val_s(ivarfl(isca(ipotr)), cpotr)

!===============================================================================
! 2.  ARC ELECTRIQUE
!===============================================================================

if ( ippmod(ielarc).ge.1 ) then

! 2.1 :  exemple : cas avec claquage
! =======================================
!    Ceci est un cas particulier et doit etre adapte en fonction
!    du cas et du maillage (intervenir aussi dans uselcl)


!        Utilisation d'une rampe d'intensite
!        -----------------------------------

    if ( ntcabs.le.200 ) then
      couimp = 200.d0
    endif

    if ( ntcabs.gt.200.and.ntcabs.le.400 ) then
      couimp = 200.d0 + 2 * (ntcabs-200)
    endif

    if ( ntcabs.gt.400 ) then
      couimp = 600.d0
    endif

!        UTILISANT D'UN CLAQUAGE AUTO
!        ----------------------------

    if(ntcabs.le.400.or.ntcabs.eq.ntpabs+1) iclaq = 0

    econs = 1.5d5
    jaiex = 0

!        ON REPERE SI IL Y A CLAQUAGE ET SI OUI OU
!        -----------------------------------------

    if(ntcabs.ge.400 .and. iclaq .eq. 0 ) then

      ! Allocate a work array
      allocate(w1(ncelet))

      amex = 1.d30
      aiex = -1.d30
      emax = 0.d0

!     les composantes du champ electrique : J/SIGMA

      call field_get_val_v(iprpfl(idjr(1)), djr)

      call field_get_key_int (ivarfl(isca(ipotr)), kivisl, ifcsig)
      call field_get_val_s(ifcsig, cpro_sig)

      do iel = 1, ncel

        xelec = djr(1,iel)/cpro_sig(iel)
        yelec = djr(2,iel)/cpro_sig(iel)
        zelec = djr(3,iel)/cpro_sig(iel)

!       Calcul du champ E
        w1(iel) = sqrt ( xelec**2 + yelec**2 + zelec**2 )
        amex =  min(amex,w1(iel))
        aiex =  max(aiex,w1(iel))
      enddo

      if(irangp.ge.0) then
        call parmin (amex)
        call parmax (aiex)
      endif
!
      write(nfecra,*) 'Min et Max de E : amex, aiex = ',amex,aiex

! Si le champ E max depasse la valeur seuil imposé, on claque à l'endroit du max de E
      if(aiex .ge. econs) then
        iclaq = 1
        ntdcla = ntcabs

!Initialisation des variables
        xclaq = 1.d-8
        yclaq = 1.d-8
        zclaq = 1.d-8
        diff  = 0.d0
        write(nfecra,*) '0000 xclaq, yclaq, zclaq = ',xclaq,yclaq,zclaq
!
        do iel = 1, ncel
          diff = aiex - w1(iel)

! Pour le multiprocessing, il ne doit y avoir le claquage que sur un seul processeur
! Pour le verifier, taper ARG_CS_OUTPUT = "--logp 1" dans le runcase
          if(diff .le. 1.d-6) then
            emax  =  w1(iel)
            xclaq =  xyzcen(1,iel)
            yclaq =  xyzcen(2,iel)
            zclaq =  xyzcen(3,iel)
            write(nfecra,*) '0011 xclaq, yclaq, zclaq = ',xclaq,yclaq,zclaq
          endif
        enddo

        call parmax (emax)
        call parmax (zclaq)

! Transfert des bonnes valeurs de x,y,zclaq entre les processeurs.
! On compare abs(xclaq) entre tous les processeurs
! Si valeur négative, parmin se charge de transférer le signe
! Attention : tous les processeurs doivent recevoir le signal PARMAX/PARMIN
! pour se synchroniser sinon calcul sans fin
!

        if(irangp .ge. 0) then
          write(nfecra,*) '1111 xclaq, yclaq, zclaq =',xclaq,yclaq,zclaq
          write(nfecra,*) 'diff =', diff
!
          if (xclaq .gt. 1.d-7) then
            call parmax (xclaq)
            call parmin (xclaq)
          elseif (xclaq .lt. -1.d-7) then
            xclaq=abs(xclaq)
            call parmax (xclaq)
            xclaq=-xclaq
            call parmin (xclaq)
          else
            call parmax (xclaq)
            call parmin (xclaq)
          endif
!
          if(yclaq .gt. 1.d-7) then
            call parmax (yclaq)
            call parmin (yclaq)
          elseif (yclaq .lt. -1.d-7) then
            yclaq=abs(yclaq)
            call parmax (yclaq)
            yclaq=-yclaq
            call parmin (yclaq)
          else
            call parmax (yclaq)
            call parmin (yclaq)
          endif
!
        endif
!
        write(nfecra,*) 'claquage : ntdcla, emax ', ntcabs, emax
        write(nfecra,*) 'xclaq, yclaq, zclaq = ',xclaq,yclaq,zclaq
      endif
      ! Free memory
      deallocate(w1)
    endif

!        SI IL Y A CLAQUAGE : ON IMPOSE COLONNE CHAUDE DU CENTRE VERS
!        LE POINT DE CLAQUAGE
!        =============================================================

    if(iclaq .eq. 1) then
      if(ntcabs.le.ntdcla+30) then
        z1 = zclaq - 3.d-4
        if(z1.le.0.d0) z1 = 0.d0
        z2 = zclaq + 3.d-4
        if(z2.ge.2.d-2) z2 = 2.d-2

        do iel = 1, ncel

          if( xyzcen(3,iel).ge.z1 .and. xyzcen(3,iel).le.z2) then
            rayo = sqrt((xclaq*xyzcen(1,iel)-yclaq*xyzcen(2,iel)  &
                 /sqrt(xclaq**2+yclaq**2))**2+(xyzcen(3,iel)      &
                 -zclaq)**2)
            posi=xclaq*xyzcen(1,iel)
            if( rayo.le.5d-4 .and. posi.ge.0d0 ) then
              cscalt(iel) = 8.d7
            endif
          endif
        enddo
      else
        iclaq = 0
      endif
    endif

!        Calcul de l'integrale sur le Volume de J.E
!        -----------------------------------
!        (c'est forcement positif ou nul)

    call field_get_val_s(iprpfl(iefjou), efjou)
    somje = 0.d0
    do iel = 1, ncel
      somje = somje+efjou(iel)*volume(iel)
    enddo

    if(irangp.ge.0) then
      call parsom (somje)
    endif

    if (somje .ne. 0) then
      coepot = couimp*dpot/max(somje,epzero)
    endif
    write(nfecra,1001) couimp,dpot,somje

!        Calcul de l'intensite du courant d'arc
!        --------------------------------------
!          Calcul de l'integrale de J sur une surface plane
!          perpendiculaire a l'axe de l'arc

!       ATTENTION : changer la valeur des tests sur CDGFAC(3,IFAC)
!                   en fonction du maillage

    call field_get_val_s(iprpfl(idjr(3)), djr3)
    elcou = 0.d0
    do ifac = 1, nfac
      if( abs(surfac(1,ifac)).le.1.d-8 .and. abs(surfac(2,ifac)).le.1.d-8 &
           .and. cdgfac(3,ifac) .gt. 0.05d-2                              &
           .and. cdgfac(3,ifac) .lt. 0.08d-2 ) then
        iel = ifacel(1,ifac)
        elcou = elcou + djr3(iel) * surfac(3,ifac)
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
    if(elcou.ne.0.d0) coepoa = couimp/elcou
    coepot = coepoa

    WRITE(NFECRA,*) ' ELCOU = ',ELCOU

    dtj = 1.d15
    dtjm =dtj
    delhsh = 0.d0
    cdtj= 2.0d2

    do iel = 1, ncel
      if (crom(iel).ne.0.d0)                                   &
           delhsh =  efjou(iel) * dt(iel) / crom(iel)

      if(delhsh.ne.0.d0) then
        dtjm= cscalt(iel)/delhsh
      else
        dtjm= dtj
      endif
      dtjm=abs(dtjm)
      dtj =min(dtj,dtjm)
    enddo

    if(irangp.ge.0) then
      call parmin (dtj)
    endif

    cpmx= sqrt(cdtj*dtj)
    coepot=cpmx
    if(ntcabs.gt.3) then
      if(coepoa.ge.1.05d0) then
        coepot=cpmx
      else
        coepot=coepoa
      endif
    endif

    write(nfecra,1008)cpmx,coepoa,coepot
    write(nfecra,1009)elcou,dpot*coepot

!        RECALAGE DES VARIABLES ELECTRIQUES
!        ----------------------------------

!         Valeur de DPOT
!         --------------

    dpot = dpot*coepot

!         Potentiel Electrique (on pourrait eviter ; c'est pour le post)
!         --------------------

    do iel = 1, ncel
      cpotr(iel) = cpotr(iel)*coepot
    enddo


!      Densite de courant (sert pour A et pour jXB)
!      ------------------

    if(ippmod(ielarc).ge.1 ) then
      call field_get_val_v(iprpfl(idjr(1)), djr)
      do idimve = 1, ndimve
        do iel = 1, ncel
          djr(idimve,iel) =  djr(idimve,iel) * coepot
        enddo
      enddo
    endif

!      Effet Joule (sert pour H au pas de temps suivant)
!      -----------

    call field_get_val_s(iprpfl(iefjou), efjou)
    do iel = 1, ncel
      efjou(iel) = efjou(iel)*coepot**2
    enddo

endif

!--------
! FORMATS
!--------

 1001  format(/, ' Courant impose= ',E14.5, /,                    &
              ' Dpot= ',E14.5,/,                            &
              ' Somje= ',E14.5)

 1008  format(/,' Cpmx   = ',E14.5,/,                             &
          ' COEPOA = ',E14.5,/,                             &
          ' COEPOT = ',E14.5)

 1009  format(/,' Courant calcule = ',E14.5,/,                    &
          ' Dpot recale     = ',E14.5)

!----
! FIN
!----

return
end subroutine uselrc
