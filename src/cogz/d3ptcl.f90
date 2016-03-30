!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2016 EDF S.A.
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

subroutine d3ptcl &
!================

 ( itypfb , izfppp ,                                              &
   rcodcl )

!===============================================================================
! FONCTION :
! --------

!    CONDITIONS AUX LIMITES AUTOMATIQUES

!           COMBUSTION GAZ CHIMIE 3 POINTS


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvarcl) !    !     !  aux faces de bord                             !
!                  !    !     ! rcodcl(1) = valeur du dirichlet                !
!                  !    !     ! rcodcl(2) = valeur du coef. d'echange          !
!                  !    !     !  ext. (infinie si pas d'echange)               !
!                  !    !     ! rcodcl(3) = valeur de la densite de            !
!                  !    !     !  flux (negatif si gain) w/m2 ou                !
!                  !    !     !  hauteur de rugosite (m) si icodcl=6           !
!                  !    !     ! pour les vitesses (vistl+visct)*gradu          !
!                  !    !     ! pour la pression             dt*gradp          !
!                  !    !     ! pour les scalaires                             !
!                  !    !     !        cp*(viscls+visct/sigmas)*gradt          !
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

!===============================================================================
! Module files
!===============================================================================

! Arguments

use paramx
use numvar
use optcal
use cstphy
use cstnum
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use field
!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          izfppp(nfabor)

double precision rcodcl(nfabor,nvarcl,3)

! Local variables

integer          igg, ifac, izone, mode
integer          ii, iel, ifue, ioxy, iok
integer          icke
double precision qisqc, viscla, d2s3, uref2, rhomoy, dhy, xiturb
double precision ustar2, xkent, xeent
double precision qcalc(nozppm)
double precision coefg(ngazgm)
double precision, dimension(:), pointer ::  brom
double precision, dimension(:), pointer :: viscl

!===============================================================================
! 0.  INITIALISATIONS
!===============================================================================


call field_get_val_s(ibrom, brom)
call field_get_val_s(iprpfl(iviscl), viscl)

d2s3 = 2.d0/3.d0

do igg = 1, ngazgm
  coefg(igg) = zero
enddo


!===============================================================================
! 1.  ECHANGES EN PARALLELE POUR LES DONNEES UTILISATEUR
!===============================================================================

!  En realite on pourrait eviter cet echange en modifiant usd3pc et en
!    demandant a l'utilisateur de donner les grandeurs dependant de la
!    zone hors de la boucle sur les faces de bord : les grandeurs
!    seraient ainsi disponibles sur tous les processeurs. Cependant,
!    ca rend le sous programme utilisateur un peu plus complique et
!    surtout, si l'utilisateur le modifie de travers, ca ne marche pas.
!  On suppose que toutes les gandeurs fournies sont positives, ce qui
!    permet d'utiliser un max pour que tous les procs les connaissent.
!    Si ce n'est pas le cas, c'est plus complique mais on pourrait
!    s'en tirer avec un max quand meme.

if(irangp.ge.0) then
  call parmax(tinfue)
  !==========
  call parmax(tinoxy)
  !==========
  call parrmx(nozapm,qimp  )
  !==========
  call parimx(nozapm,iqimp )
  !==========
  call parimx(nozapm,ientox)
  !==========
  call parimx(nozapm,ientfu)
  !==========
endif

!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES POUR LES SCALAIRES
!
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================

!       ON DETERMINE LA FAMILLE ET SES PROPRIETES
!       ON IMPOSE LES CONDITIONS AUX LIMITES
!===============================================================================

! Fuel and Oxidant inlet are checked
ifue = 0
ioxy = 0
do ii = 1, nzfppp
  izone = ilzppp(ii)
  if (ientfu(izone).eq.1) then
    ifue = 1
  elseif(ientox(izone).eq.1) then
    ioxy = 1
  endif
enddo
if(irangp.ge.0) then
  call parcmx(ifue)
  call parcmx(ioxy)
endif

! Fuel inlet at TINFUE: HINFUE calculation
if (ifue.eq.1) then
  coefg(1) = 1.d0
  coefg(2) = zero
  coefg(3) = zero
  mode    = -1
  call cothht                                                   &
  !==========
      ( mode   , ngazg , ngazgm  , coefg  ,                     &
        npo    , npot   , th     , ehgazg ,                     &
        hinfue , tinfue )
endif

! Oxidant inlet at TINOXY: HINOXY calculation
if (ioxy.eq.1) then
  coefg(1) = zero
  coefg(2) = 1.d0
  coefg(3) = zero
  mode    = -1
  call cothht                                                   &
  !==========
      ( mode   , ngazg , ngazgm  , coefg  ,                     &
        npo    , npot   , th     , ehgazg ,                     &
        hinoxy , tinoxy  )
endif


do ifac = 1, nfabor

  izone = izfppp(ifac)

  if (itypfb(ifac).eq.ientre.or.itypfb(ifac).eq.i_convective_inlet) then

    ! Fuel inlet at TINFUE
    if (ientfu(izone).eq.1) then

      ! Mean mixture fraction
      rcodcl(ifac,isca(ifm),1)   = 1.d0

      ! Mixture fraction variance
      rcodcl(ifac,isca(ifp2m),1) = 0.d0

      ! Mixture enthalpy
      if (ippmod(icod3p).eq.1) then
        rcodcl(ifac,isca(iscalt),1) = hinfue
      endif

      ! Soot model
      if (isoot.ge.1) then
        rcodcl(ifac,isca(ifsm),1) = 0.d0
        rcodcl(ifac,isca(inpm),1) = 0.d0
      endif

      ! Density
      brom(ifac) = pther/(cs_physical_constants_r*tinfue/wmolg(1))

    ! Oxydant inlet at TINOXY
    elseif (ientox(izone).eq.1) then

      ! Mean mixture fraction
      rcodcl(ifac,isca(ifm),1)   = 0.d0

      ! Mixture fraction variance
      rcodcl(ifac,isca(ifp2m),1) = 0.d0

      ! Mixture enthalpy
      if ( ippmod(icod3p).eq.1 ) then
        rcodcl(ifac,isca(iscalt),1) = hinoxy
      endif

      ! Soot model
      if (isoot.ge.1) then
        rcodcl(ifac,isca(ifsm),1) = 0.d0
        rcodcl(ifac,isca(inpm),1) = 0.d0
      endif

      ! Density
      brom(ifac) = pther/(cs_physical_constants_r*tinoxy/wmolg(2))

    endif

  endif

enddo

!===============================================================================
! 3.  SI IQIMP = 1 : CORRECTION DES VITESSES (EN NORME) POUR CONTROLER
!                    LES DEBITS IMPOSES
!     SI IQIMP = 0 : CALCUL DE QIMP

!     ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                   =========================
!===============================================================================

! --- Debit calcule

do izone = 1, nozppm
  qcalc(izone) = 0.d0
enddo
do ifac = 1, nfabor
  izone = izfppp(ifac)
  qcalc(izone) = qcalc(izone) - brom(ifac) *             &
     ( rcodcl(ifac,iu,1)*surfbo(1,ifac) +                  &
       rcodcl(ifac,iv,1)*surfbo(2,ifac) +                  &
       rcodcl(ifac,iw,1)*surfbo(3,ifac) )
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
  izone = ilzppp(ii)
  if ( iqimp(izone).eq.1 ) then
    if(qcalc(izone).lt.epzero) then
      write(nfecra,2001)izone,iqimp(izone),qcalc(izone)
      iok = iok + 1
    endif
  endif
enddo
if(iok.ne.0) then
  call csexit (1)
  !==========
endif
do ifac = 1, nfabor
  izone = izfppp(ifac)
  if ( iqimp(izone).eq.1 ) then
    qisqc = qimp(izone)/qcalc(izone)
    rcodcl(ifac,iu,1) = rcodcl(ifac,iu,1)*qisqc
    rcodcl(ifac,iv,1) = rcodcl(ifac,iv,1)*qisqc
    rcodcl(ifac,iw,1) = rcodcl(ifac,iw,1)*qisqc
  endif
enddo

 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MODULE PHYSIQUES PARTICULIERES              ',/,&
'@    =========                                               ',/,&
'@    PROBLEME DANS LES CONDITIONS AUX LIMITES                ',/,&
'@                                                            ',/,&
'@  Le debit est impose sur la zone IZONE = ', I10             ,/,&
'@    puisque                IQIMP(IZONE) = ', I10             ,/,&
'@  Or, sur cette zone, le produit RHO D S integre est nul :  ',/,&
'@    il vaut                             = ',E14.5            ,/,&
'@    (D est la direction selon laquelle est impose le debit).',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usd3pc, et en particulier                        ',/,&
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

!===============================================================================
! 4.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES POUR LA TURBULENCE
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!         ON IMPOSE LES CONDITIONS AUX LIMITES POUR LA TURBULENCE
!         (pour n'importe quel modele)
!===============================================================================

do ifac = 1, nfabor

  izone = izfppp(ifac)

  if (itypfb(ifac).eq.ientre.or.itypfb(ifac).eq.i_convective_inlet) then

    ! La turbulence est calculee par defaut si ICALKE different de 0
    !    - soit a partir du diametre hydraulique, d'une vitesse
    !      de reference adaptes a l'entree courante si ICALKE = 1
    !    - soit a partir du diametre hydraulique, d'une vitesse
    !      de reference et de l'intensite turvulente
    !      adaptes a l'entree courante si ICALKE = 2

    if ( icalke(izone).ne.0 ) then

      uref2 = rcodcl(ifac,iu,1)**2                         &
            + rcodcl(ifac,iv,1)**2                         &
            + rcodcl(ifac,iw,1)**2
      uref2 = max(uref2,epzero)
      rhomoy = brom(ifac)
      iel    = ifabor(ifac)
      viscla = viscl(iel)
      icke   = icalke(izone)
      dhy    = dh(izone)
      xiturb = xintur(izone)
      ustar2 = 0.d0
      xkent = epzero
      xeent = epzero
      if (icke.eq.1) then
        call keendb ( uref2, dhy, rhomoy, viscla, cmu, xkappa, &
          ustar2, xkent, xeent )
      else if (icke.eq.2) then
        call keenin ( uref2, xiturb, dhy, cmu, xkappa, xkent, xeent )
      endif

      if (itytur.eq.2) then

        rcodcl(ifac,ik,1)  = xkent
        rcodcl(ifac,iep,1) = xeent

      elseif (itytur.eq.3) then

        rcodcl(ifac,ir11,1) = d2s3*xkent
        rcodcl(ifac,ir22,1) = d2s3*xkent
        rcodcl(ifac,ir33,1) = d2s3*xkent
        rcodcl(ifac,ir12,1) = 0.d0
        rcodcl(ifac,ir13,1) = 0.d0
        rcodcl(ifac,ir23,1) = 0.d0
        rcodcl(ifac,iep,1)  = xeent

      elseif (iturb.eq.50) then

        rcodcl(ifac,ik,1)   = xkent
        rcodcl(ifac,iep,1)  = xeent
        rcodcl(ifac,iphi,1) = d2s3
        rcodcl(ifac,ifb,1)  = 0.d0

      elseif (iturb.eq.60) then

        rcodcl(ifac,ik,1)   = xkent
        rcodcl(ifac,iomg,1) = xeent/cmu/xkent

      elseif(iturb.eq.70) then

        rcodcl(ifac,inusa,1) = cmu*xkent**2/xeent

      endif

    endif

  endif

enddo

!----
! FORMATS
!----


!----
! FIN
!----

return
end subroutine
