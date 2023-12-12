!-------------------------------------------------------------------------------

! This file is part of code_saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2023 EDF S.A.
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

!===============================================================================
! Function :
! --------

!> \file lwctcl.f90
!> \brief Automatic boundary conditions for partially premixed flame combustion
!>        model (LWC)

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     nvar          total number of variables
!> \param[in]     itypfb        boundary face types
!> \param[out]    izfppp        boundary face zone number
!> \param[out]    rcodcl        boundary condition values:
!>                               - rcodcl(1) value of the dirichlet
!>                               - rcodcl(2) value of the exterior exchange
!>                                 coefficient (infinite if no exchange)
!>                               - rcodcl(3) value flux density
!>                                 (negative if gain) in w/m2 or roughness
!>                                 in m if icodcl=6
!>                                 -# for the velocity \f$ (\mu+\mu_T)
!>                                    \gradt \, \vect{u} \cdot \vect{n}  \f$
!>                                 -# for the pressure \f$ \Delta t
!>                                    \grad P \cdot \vect{n}  \f$
!>                                 -# for a scalar \f$ cp \left( K +
!>                                     \dfrac{K_T}{\sigma_T} \right)
!>                                     \grad T \cdot \vect{n} \f$
!_______________________________________________________________________________

subroutine lwctcl &
 ( itypfb , izfppp ,                                              &
   rcodcl )

!===============================================================================

!===============================================================================
! Module files
!===============================================================================

use paramx
use numvar
use optcal
use cstphy
use cstnum
use dimens, only : nvar
use entsor
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use mesh
use field
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

integer          itypfb(nfabor)
integer          izfppp(nfabor)

double precision rcodcl(nfabor,nvar,3)

! Local variables

integer          nbr
integer          igg, ifac, izone
integer          icke, ii, iel, iok
double precision qisqc, viscla, d2s3, uref2, rhomoy, dhy, xiturb
double precision hgazf , tgazf, hgazb, tgazb
double precision qcalc(nozppm), hgent(nozppm)
double precision coefg(ngazgm), rval(1)
double precision, dimension(:), pointer ::  brom
double precision, dimension(:), pointer :: viscl

!===============================================================================
! 1.  INITIALISATIONS
!===============================================================================


call field_get_val_s(ibrom, brom)
call field_get_val_s(iviscl, viscl)

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

if (irangp.ge.0) then
  call parrmx(nozapm,qimp  )
  call parrmx(nozapm,fment )
  call parrmx(nozapm,tkent )
  call parimx(nozapm,iqimp )
  call parimx(nozapm,ientgf)
  call parimx(nozapm,ientgb)
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
'@  Verifier uslwcc, et en particulier                        ',/,&
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

      if (icke.eq.1) then
        !   Calculation of turbulent inlet conditions using
        !     standard laws for a circular pipe
        !     (their initialization is not needed here but is good practice).
        call turbulence_bc_inlet_hyd_diam(ifac, uref2, dhy, rhomoy, viscla,  &
                                          rcodcl)
      else if (icke.eq.2) then

        ! Calculation of turbulent inlet conditions using
        !   the turbulence intensity and standard laws for a circular pipe
        !   (their initialization is not needed here but is good practice)

        call turbulence_bc_inlet_turb_intensity(ifac, uref2, xiturb, dhy,  &
                                                rcodcl)


      endif

    endif

  endif

enddo

!===============================================================================
! 3.  VERIFICATION DES DONNEES POUR LA FRACTION DE MELANGE
!                              ET LA TEMPERATURE DES GAZ FRAIS
!    (modele lwc)
!===============================================================================

! --- FRMEL et TGF (on n'en veut qu'un : on prend le max)
!     EBU nominal est a f homogene
!     On se limite pour l'instant a une temperature
!       des gaz frais identiques

!      FRMEL = 0.D0
!      TGF   = 0.D0
!      DO IFAC = 1, NFABOR
!        IF ( ITYPFB(IFAC).EQ.IENTRE ) THEN
!          IZONE = IZFPPP(IFAC)
!          IF ( IPPMOD(ICOEBU).EQ.0 .OR. IPPMOD(ICOEBU).EQ.1 ) THEN
!            FRMEL = MAX(FMENT(IZONE),FRMEL)
!          ENDIF
!          IF (IENTGF(IZONE).EQ.1) THEN
!            TGF = MAX(TKENT(IZONE),TGF)
!          ENDIF
!        ENDIF
!      ENDDO

!     IF(IRANGP    .GE.0) THEN
!       CALL PARMAX(FRMEL)
!       CALL PARMAX(TGF  )
!      ENDIF

! Attention, ici on modifie FMENT et TKENT sur les zones
!  presentes sur le proc local. Ca suffit pour le traitement qui suit.
!     DO IFAC = 1, NFABOR
!        IZONE = IZFPPP(IFAC)
!        IF ( ITYPFB(IFAC).EQ.IENTRE ) THEN
!          IF ( IPPMOD(ICOEBU).EQ.0 .OR. IPPMOD(ICOEBU).EQ.1 ) THEN
!            FMENT(IZONE) = FRMEL
!          ENDIF
!          IF (IENTGF(IZONE).EQ.1) THEN
!            TKENT(IZONE) = TGF
!          ENDIF
!        ENDIF
!      ENDDO

!===============================================================================
! 2.  REMPLISSAGE DU TABLEAU DES CONDITIONS LIMITES
!       ON BOUCLE SUR TOUTES LES FACES D'ENTREE
!                     =========================
!         ON DETERMINE LA FAMILLE ET SES PROPRIETES
!           ON IMPOSE LES CONDITIONS AUX LIMITES
!           POUR LES SCALAIRES
!    (modele ebu)
!===============================================================================


! ---> Combustion gaz USLWCC
!      Flamme de premelange modele LWC

!      Enthalpie du melange gazeux :
!         hors de la boucle pour eviter un appel par face.
!         Suppose que une entree est forcement IENTGF=1 ou IENTGB=1


do ii = 1, nzfppp
  izone = ilzppp(ii)
!       Entree 1
  if ( ientgf(izone).eq.1 ) then
    tgazf    = tkent(izone)
    coefg(1) = fment(izone)
    coefg(2) = 1.d0 - fment(izone)
    coefg(3) = zero
    hgazf = cs_gas_combustion_t_to_h(coefg, tgazf)
    hgent(izone) = hgazf
!       Entree 2
  elseif ( ientgb(izone).eq.1 ) then
    tgazb    = tkent(izone)
    coefg(1) = fment(izone)
    coefg(2) = 1.d0 - fment(izone)
    coefg(3) = zero
    hgazb = cs_gas_combustion_t_to_h(coefg, tgazb)
    hgent(izone) = hgazb
  endif
enddo

do ifac = 1, nfabor

  izone = izfppp(ifac)

!      ELEMENT ADJACENT A LA FACE DE BORD

  if ( itypfb(ifac).eq.ientre ) then

! ----  Traitement automatique des scalaires physiques particulieres

!       Entree gaz frais

    if ( ientgf(izone).eq.1 ) then

!         - Fraction massique de fuel
      rcodcl(ifac,isca(iyfm),1) = fment(izone)

!         - Variance de la fraction massique
      rcodcl(ifac,isca(iyfp2m),1) = zero

!         - Fraction de melange
        rcodcl(ifac,isca(ifm),1) = fment(izone)

!         - Variance de la fraction de melange
      rcodcl(ifac,isca(ifp2m),1) = zero

      if ( ippmod(icolwc).ge.2 ) then
        rcodcl(ifac,isca(icoyfp),1) = zero
      endif

!         - Enthalpie du melange gazeux
      if ( ippmod(icolwc) .eq. 1 .or.                             &
           ippmod(icolwc) .eq. 3 .or.                             &
           ippmod(icolwc) .eq. 5    ) then
        rcodcl(ifac,isca(iscalt),1) = hgent(izone)
      endif

    elseif ( ientgb(izone).eq.1 ) then

!       Entree gaz brule

!         - Fraction massique de fuel
      rcodcl(ifac,isca(iyfm),1) = zero

!         - Variance de la fraction massique
       rcodcl(ifac,isca(iyfp2m),1) = zero

!         - Fraction de melange
      rcodcl(ifac,isca(ifm),1) = fment(izone)

!         - Variance de la fraction de melange
      rcodcl(ifac,isca(ifp2m),1) = zero

      if ( ippmod(icolwc) .ge.2) then
        rcodcl(ifac,isca(icoyfp),1) = zero
      endif

!         - Enthalpie du melange gazeux
      if ( ippmod(icolwc) .eq. 1 .or.                             &
           ippmod(icolwc) .eq. 3 .or.                             &
           ippmod(icolwc) .eq. 5    ) then
        rcodcl(ifac,isca(iscalt),1) = hgent(izone)
      endif

    endif

  endif

enddo

! Calcul de FMIN/FMAX et HMIN/HMAX sur les entrees

fmin = 1.e+30
fmax =-1.e+30

do ifac = 1, nfabor

  izone = izfppp(ifac)

!      ELEMENT ADJACENT A LA FACE DE BORD

  if ( itypfb(ifac).eq.ientre ) then

    if ( fment(izone) .lt. fmin ) then
     fmin= fment(izone)
     hmin= hgent(izone)
    endif
    if ( fment(izone) .gt. fmax ) then
     fmax= fment(izone)
     hmax= hgent(izone)
    endif
  endif
enddo

if (irangp.ge.0) then

  nbr = 1

  rval(1) = hmax
  call parmxl(nbr,fmax,rval)
  hmax = rval(1)

  rval(1) = hmin
  call parmnl(nbr,fmin,rval)
  hmin = rval(1)

endif

!----
! FORMATS
!----

!----
! FIN
!----

return
end subroutine
