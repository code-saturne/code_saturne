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

subroutine ebutcl &
!================

 ( itypfb , izfppp ,                                              &
   rcodcl )

!===============================================================================
! FONCTION :
! --------

!    CONDITIONS AUX LIMITES AUTOMATIQUES

!           COMBUSTION GAZ MODELE EBU


!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! itypfb           ! ia ! <-- ! boundary face types                            !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! rcodcl           ! tr ! --> ! valeur des conditions aux limites              !
!  (nfabor,nvar    !    !     !  aux faces de bord                             !
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
integer          icke, ii, iel, iok
double precision qisqc, viscla, d2s3, uref2, rhomoy, dhy, xiturb
double precision ustar2, xkent, xeent, hgazf , tgazf, hgazb, tgazb
double precision qcalc(nozppm), hgent(nozppm)
double precision coefg(ngazgm)
double precision, dimension(:), pointer ::  brom
double precision, dimension(:), pointer :: viscl

!===============================================================================
!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================


call field_get_val_s(ibrom, brom)
call field_get_val_s(iprpfl(iviscl), viscl)

d2s3 = 2.d0/3.d0

do igg = 1, ngazgm
  coefg(igg) = 0
enddo

!===============================================================================
! 1.  ECHANGES EN PARALLELE POUR LES DONNEES UTILISATEUR
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

if(irangp.ge.0) then
  call parrmx(nozapm,qimp  )
  !==========
  call parrmx(nozapm,fment )
  !==========
  call parrmx(nozapm,tkent )
  !==========
  call parimx(nozapm,iqimp )
  !==========
  call parimx(nozapm,ientgf)
  !==========
  call parimx(nozapm,ientgb)
  !==========
endif

!===============================================================================
! 2.  SI IQIMP = 1 : CORRECTION DES VITESSES (EN NORME) POUR CONTROLER
!                    LES DEBITS IMPOSES
!     SI IQIMP = 0 : CALCUL DE QIMP

!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
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
'@  Verifier usebuc, et en particulier                        ',/,&
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
! 4.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LA TURBULENCE
!===============================================================================


do ifac = 1, nfabor

  izone = izfppp(ifac)

!      ELEMENT ADJACENT A LA FACE DE BORD

  if ( itypfb(ifac).eq.ientre ) then

! ----  Traitement automatique de la turbulence

    if ( icalke(izone).ne.0 ) then

!       La turbulence est calculee par defaut si ICALKE different de 0
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference adaptes a l'entree courante si ICALKE = 1
!          - soit a partir du diametre hydraulique, d'une vitesse
!            de reference et de l'intensite turvulente
!            adaptes a l'entree courante si ICALKE = 2

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
        call keendb                                               &
        !==========
        ( uref2, dhy, rhomoy, viscla, cmu, xkappa,                &
          ustar2, xkent, xeent )
      else if (icke.eq.2) then
        call keenin                                               &
        !==========
        ( uref2, xiturb, dhy, cmu, xkappa, xkent, xeent )
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

!===============================================================================
! 3.  VERIFICATION DES DONNEES POUR LA FRACTION DE MELANGE
!                              ET LA TEMPERATURE DES GAZ FRAIS
!    (modele ebu)
!===============================================================================

! --- FRMEL et TGF (on n'en veut qu'un : on prend le max)
!     EBU nominal est a f homogene
!     On se limite pour l'instant a une temperature
!       des gaz frais identiques

frmel = 0.d0
tgf   = 0.d0
do ifac = 1, nfabor
  if ( itypfb(ifac).eq.ientre ) then
    izone = izfppp(ifac)
    if ( ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.1 ) then
      frmel = max(fment(izone),frmel)
    endif
    if (ientgf(izone).eq.1) then
      tgf = max(tkent(izone),tgf)
    endif
  endif
enddo

if(irangp.ge.0) then
  call parmax(frmel)
  call parmax(tgf  )
endif

! Attention, ici on modifie FMENT et TKENT sur les zones
!  presentes sur le proc local. Ca suffit pour le traitement qui suit.
do ifac = 1, nfabor
  izone = izfppp(ifac)
  if ( itypfb(ifac).eq.ientre ) then
    if ( ippmod(icoebu).eq.0 .or. ippmod(icoebu).eq.1 ) then
      fment(izone) = frmel
    endif
    if (ientgf(izone).eq.1) then
      tkent(izone) = tgf
    endif
  endif
enddo

!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LES SCALAIRES
!    (modele ebu)
!===============================================================================


! ---> Combustion gaz USEBUC
!      Flamme de premelange modele EBU

!      Enthalpie du melange gazeux :
!         hors de la boucle pour eviter un appel par face.
!         Suppose que une entree est forcement IENTGF=1 ou IENTGB=1

if ( ippmod(icoebu) .eq. 1 .or.                                   &
     ippmod(icoebu) .eq. 3        ) then

  do ii = 1, nzfppp
    izone = ilzppp(ii)
!       Entree premelange froid ou dilution
    if ( ientgf(izone).eq.1 ) then
      tgazf    = tkent(izone)
      coefg(1) = fment(izone)
      coefg(2) = 1.d0 - fment(izone)
      coefg(3) = zero
      mode    = -1
      call cothht                                                 &
      !==========
       ( mode   , ngazg  , ngazgm , coefg  ,                      &
         npo    , npot   , th     , ehgazg ,                      &
         hgazf  , tgazf  )
      hgent(izone) = hgazf
!       Entree gaz brules (flamme pilote)
    elseif ( ientgb(izone).eq.1 ) then
      tgazb    = tkent(izone)
      coefg(1) = max(zero,(fment(izone)-fs(1))/(1.d0-fs(1)))
      coefg(3) = (fment(izone)-coefg(1))/fs(1)
      coefg(2) = 1.d0 - coefg(1) - coefg(3)
      mode    = -1
      call cothht                                                 &
      !==========
       ( mode   , ngazg , ngazgm  , coefg  ,                      &
         npo    , npot   , th     , ehgazg ,                      &
         hgazb , tgazb )
      hgent(izone) = hgazb
    endif
  enddo

endif


do ifac = 1, nfabor

  izone = izfppp(ifac)

!      ELEMENT ADJACENT A LA FACE DE BORD

  if ( itypfb(ifac).eq.ientre ) then

! ----  Traitement automatique des scalaires physiques particulieres

!       Entree premelange froid ou dilution

    if ( ientgf(izone).eq.1 ) then

!         - Fraction massique de gaz frais
      rcodcl(ifac,isca(iygfm),1) = 1.d0

!         - Fraction de melange
      if ( ippmod(icoebu) .eq. 2 .or.                             &
           ippmod(icoebu) .eq. 3      ) then
        rcodcl(ifac,isca(ifm),1) = fment(izone)
      endif

!          - Enthalpie du melange gazeux
      if ( ippmod(icoebu) .eq. 1 .or.                             &
           ippmod(icoebu) .eq. 3        ) then
        rcodcl(ifac,isca(iscalt),1) = hgent(izone)
      endif

    elseif ( ientgb(izone).eq.1 ) then

!       Entree gaz brules (flamme pilote)

!         - Fraction massique de gaz frais
        rcodcl(ifac,isca(iygfm),1) = zero

!         - Fraction de melange
      if ( ippmod(icoebu) .eq. 2 .or.                             &
           ippmod(icoebu) .eq. 3      ) then
        rcodcl(ifac,isca(ifm),1) = fment(izone)
      endif

!          - Enthalpie du melange gazeux
      if ( ippmod(icoebu) .eq. 1 .or.                             &
           ippmod(icoebu) .eq. 3        ) then
        rcodcl(ifac,isca(iscalt),1) = hgent(izone)
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
