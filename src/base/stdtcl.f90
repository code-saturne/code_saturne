!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2019 EDF S.A.
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

subroutine stdtcl &
!================

 ( nbzfmx , nozfmx ,                                              &
   iqimp  , icalke , qimp   , dh     , xintur ,                   &
   itypfb , iznfbr ,                                              &
   rcodcl )

!===============================================================================
! FONCTION :
! --------

!    CONDITIONS AUX LIMITES AUTOMATIQUES

!          EN STANDARD


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nbzfmx           ! e  ! <-- ! nb max de zones de faces de bord               !
! nozfmx           ! e  ! <-- ! numero max de zones de faces de bord           !
! itypfb           ! ia ! <-- ! boundary face types                            !
! iznfbr           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !                                                !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!                  !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/turb_schmidt)*gradt    !
!__________________!____!_____!________________________________________________!

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only: nvar
use entsor
use parall
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          nozfmx
integer          nbzfmx

integer          iqimp(nozfmx), icalke(nozfmx)
integer          itypfb(nfabor)
integer          iznfbr(nfabor)

double precision qimp(nozfmx), dh(nozfmx), xintur(nozfmx)
double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          ifac, izone, ifvu, izonem
integer          nozapm, nzfppp
integer          icke, ii, iel, iok
double precision qisqc, viscla, uref2, rhomoy, dhy, xiturb
double precision, dimension(:), pointer :: brom
double precision, dimension(:), pointer :: viscl
integer, allocatable, dimension(:) :: ilzfbr
double precision, allocatable, dimension(:) :: qcalc

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================

allocate(ilzfbr(nbzfmx))
allocate(qcalc(nozfmx))


call field_get_val_s(ibrom, brom)
call field_get_val_s(iviscl, viscl)

! --> On construit une liste des numeros des zones frontieres.
!           (liste locale au processeur, en parallele)
nzfppp = 0
do ifac = 1, nfabor
  ifvu = 0
  do ii = 1, nzfppp
    if (ilzfbr(ii).eq.iznfbr(ifac)) then
      ifvu = 1
    endif
  enddo
  if(ifvu.eq.0) then
    nzfppp = nzfppp + 1
    if(nzfppp.le.nbzfmx) then
      ilzfbr(nzfppp) = iznfbr(ifac)
    else
      write(nfecra,1001) nbzfmx
      write(nfecra,1002)(ilzfbr(ii),ii=1,nbzfmx)
      call csexit (1)
      !==========
    endif
  endif
enddo

! ---> Plus grand numero de zone

izonem = 0
do ii = 1, nzfppp
  izone = ilzfbr(ii)
  izonem = max(izonem,izone)
enddo
if(irangp.ge.0) then
  call parcmx(izonem)
  !==========
endif
nozapm = izonem

#if defined(_CS_LANG_FR)

 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS LES CONDITIONS AUX LIMITES    ',/,&
'@    =========                                               ',/,&
'@                ARRET DANS LE SOUS-PROGRAMME STDTCL         ',/,&
'@                                                            ',/,&
'@  Le nombre maximal de zones frontieres qui peuvent etre    ',/,&
'@    definies par l''utilisateur est NBZPPM = ',I10           ,/,&
'@    Il a ete depasse.                                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les conditions aux limites.                      ',/,&
'@                                                            ',/,&
'@  Les NBZPPM premieres zones frontieres                     ',/,&
'@    portent ici les numeros suivants :                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(i10)

#else

 1001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS LES CONDITIONS AUX LIMITES    ',/,&
'@    =========                                               ',/,&
'@                ARRET DANS LE SOUS-PROGRAMME STDTCL         ',/,&
'@                                                            ',/,&
'@  The maximum number of boundary zones which can be defined ',/,&
'@    by the user is NBZPPM = ',I10                            ,/,&
'@    It has been exceeded.                                   ',/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify the boundary conditions.                           ',/,&
'@                                                            ',/,&
'@  The first NBZPPM boundary zones have the following number:',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1002 format(i10)

#endif

!===============================================================================
! 2.  ECHANGES EN PARALLELE POUR LES DONNEES UTILISATEUR
!===============================================================================

!  En realite on pourrait eviter cet echange en modifiant usebuc et en
!    demandant a l'utilisateur de donner les grandeurs dependant de la
!    zone hors de la boucle sur les faces de bord : les grandeurs
!    seraient ainsi disponibles sur tous les processeurs. Cependant,
!    ca rend le sous programme utilisateur un peu plus complique et
!    surtout, si l'utilisateur le modifie de travers, ca ne marche pas.
!  On suppose que toutes les gandeurs fournies sont positives, ce qui
!    permet d'utiliser un max pour que tous les procs les connaissent.
!    Si ce n'est pas le cas, c'est plus complique mais on peut s'en tirer
!    avec un max quand meme.

if (irangp.ge.0) then
  call parrmx(nozapm, qimp)
  !==========
  call parimx(nozapm, iqimp)
  !==========
endif

!===============================================================================
! 3.  SI IQIMP = 1 : CORRECTION DES VITESSES (EN NORME) POUR CONTROLER
!                    LES DEBITS IMPOSES
!     SI IQIMP = 0 : CALCUL DE QIMP

!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!===============================================================================


! --- Debit calcule

do izone = 1,nozfmx
  qcalc(izone) = 0.d0
enddo
do ifac = 1,nfabor
  izone = iznfbr(ifac)
  if (izone .gt. 0) then
    if (iqimp(izone).eq.2) then
      qcalc(izone) = qcalc(izone) -                               &
         ( rcodcl(ifac,iu,1)*surfbo(1,ifac) +              &
           rcodcl(ifac,iv,1)*surfbo(2,ifac) +              &
           rcodcl(ifac,iw,1)*surfbo(3,ifac) )
    else
      qcalc(izone) = qcalc(izone) - brom(ifac) *           &
         ( rcodcl(ifac,iu,1)*surfbo(1,ifac) +              &
           rcodcl(ifac,iv,1)*surfbo(2,ifac) +              &
           rcodcl(ifac,iw,1)*surfbo(3,ifac) )
    endif
  endif
enddo
if(irangp.ge.0) then
  call parrsm(nozapm,qcalc)
endif
do izone = 1, nozapm
  if ( iqimp(izone).eq.0 ) then
    qimp(izone) = qcalc(izone)
  endif
enddo

! --- Correction des vitesses en norme

iok = 0
do ii = 1, nzfppp
  izone = ilzfbr(ii)
  if (izone .gt. 0) then
    if ( iqimp(izone).eq.1 .or.  iqimp(izone).eq.2) then
      if(qcalc(izone).lt.epzero) then
        write(nfecra,2001)izone,iqimp(izone),qcalc(izone)
        iok = iok + 1
      endif
    endif
  endif
enddo
if(iok.ne.0) then
  call csexit (1)
  !==========
endif
do ifac = 1, nfabor
  izone = iznfbr(ifac)
  if (izone .gt. 0) then
    if ( iqimp(izone).eq.1 .or.  iqimp(izone).eq.2) then
      qisqc = qimp(izone)/qcalc(izone)
      rcodcl(ifac,iu,1) = rcodcl(ifac,iu,1)*qisqc
      rcodcl(ifac,iv,1) = rcodcl(ifac,iv,1)*qisqc
      rcodcl(ifac,iw,1) = rcodcl(ifac,iw,1)*qisqc
    endif
  endif
enddo

#if defined(_CS_LANG_FR)

 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : PROBLEME DANS LES CONDITIONS AUX LIMITES    ',/,&
'@    =========                                               ',/,&
'@                ARRET DANS LE SOUS-PROGRAMME STDTCL         ',/,&
'@                                                            ',/,&
'@  Le debit est impose sur la zone IZONE = ', I10             ,/,&
'@    puisque                IQIMP(IZONE) = ', I10             ,/,&
'@  Or, sur cette zone, le produit RHO D S integre est nul :  ',/,&
'@    il vaut                             = ',E14.5            ,/,&
'@    (D est la direction selon laquelle est impose le debit).',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les donnees dans l''interface et en particulier  ',/,&
'@    - que le vecteur  RCODCL(IFAC,IU,1),             ',/,&
'@                      RCODCL(IFAC,IV,1),             ',/,&
'@                      RCODCL(IFAC,IW,1) qui determine',/,&
'@      la direction de la vitesse est non nul et n''est pas  ',/,&
'@      uniformement perpendiculaire aux face d''entree       ',/,&
'@    - que la surface de l''entree n''est pas nulle (ou que  ',/,&
'@      le nombre de faces de bord dans la zone est non nul)  ',/,&
'@    - que la masse volumique n''est pas nulle               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: PROBLEM IN THE BOUNDARY CONDITIONS             ',/,&
'@    ========                                                ',/,&
'@                ABORT IN THE SUBROUTINE STDTCL              ',/,&
'@                                                            ',/,&
'@  The flow is imposed on the zone IZONE = ', I10             ,/,&
'@    since                  IQIMP(IZONE) = ', I10             ,/,&
'@  But, on this zone, the integrated product RHO D S is zero:',/,&
'@    its value is                        = ',E14.5            ,/,&
'@    (D is the direction along which is imposed the flow).   ',/,&
'@                                                            ',/,&
'@  The calculation will not run.                             ',/,&
'@                                                            ',/,&
'@  Verify the data in the interface and particularly         ',/,&
'@    - that the vector RCODCL(IFAC,IU,1),             ',/,&
'@                      RCODCL(IFAC,IV,1),             ',/,&
'@                      RCODCL(IFAC,IW,1) which gives  ',/,&
'@      the velocity direction is non null and not uniformly  ',/,&
'@      perpendicular to the inlet faces                      ',/,&
'@    - that the inlet surface is not zero (or that the number',/,&
'@      of boundary faces within the zone is not zero)        ',/,&
'@    - that the density is not zero                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!===============================================================================
! 4.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LA TURBULENCE
!===============================================================================


do ifac = 1, nfabor

  izone = iznfbr(ifac)

  if (izone .gt. 0) then

    if (     itypfb(ifac).eq.ientre  &
        .or. itypfb(ifac).eq.i_convective_inlet) then

! ----  Traitement automatique de la turbulence

      if ( icalke(izone).ne.0 ) then

!       La turbulence est calculee par defaut si ICALKE different de 0
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference adaptes a l'entree courante si ICALKE = 1
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference et de l'intensite turvulente
!            adaptes a l'entree courante si ICALKE = 2

        uref2 = rcodcl(ifac,iu,1)**2                       &
              + rcodcl(ifac,iv,1)**2                       &
              + rcodcl(ifac,iw,1)**2
        uref2 = max(uref2,epzero)
        rhomoy = brom(ifac)
        iel    = ifabor(ifac)
        viscla = viscl(iel)
        icke   = icalke(izone)
        dhy    = dh(izone)
        xiturb = xintur(izone)
        if (icke.eq.1) then
          call turbulence_bc_inlet_hyd_diam(ifac,                            &
                                            uref2, dhy, rhomoy, viscla,      &
                                            rcodcl)
        else if (icke.eq.2) then
          call turbulence_bc_inlet_turb_intensity(ifac,                      &
                                                  uref2, xiturb, dhy,        &
                                                  rcodcl)
        endif

      endif

    endif

    if (itypfb(ifac).eq.iparoi) then
        ! condition automatique en paroi pour alpha
        if (iturb.eq.32) then
          rcodcl(ifac,ial,1)  = 0.d0
        endif
    endif

  endif

enddo

! Free memory
deallocate(ilzfbr)
deallocate(qcalc)

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
