!-------------------------------------------------------------------------------

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

subroutine elphyv &
!================

 ( nvar   , nscal  ,                                              &
   mbrom  , izfppp ,                                              &
   dt     , rtp    , propce )

!===============================================================================
! FONCTION :
! --------

!   REMPLISSAGE DES VARIABLES PHYSIQUES : Version Electrique

!     ----> Effet Joule
!     ----> Arc Electrique
!     ----> Conduction Ionique

!      1) Masse Volumique
!      2) Viscosite moleculaire
!      3) Cp
!      4) Lambda/Cp moleculaire
!      4) Diffusivite moleculaire



! ATTENTION :
! =========


! Il est INTERDIT de modifier la viscosite turbulente VISCT ici
!        ========
!  (une routine specifique est dediee a cela : usvist)


!  Il FAUT AVOIR PRECISE ICP = 1
!     ==================
!    si on souhaite imposer une chaleur specifique
!    CP variable (sinon: ecrasement memoire).


!  Il FAUT AVOIR PRECISE IVISLS(Numero de scalaire) = 1
!     ==================
!     si on souhaite une diffusivite VISCLS variable
!     pour le scalaire considere (sinon: ecrasement memoire).




! Remarques :
! ---------

! Cette routine est appelee au debut de chaque pas de temps

!    Ainsi, AU PREMIER PAS DE TEMPS (calcul non suite), les seules
!    grandeurs initialisees avant appel sont celles donnees
!      - dans usipsu :
!             . la masse volumique (initialisee a RO0)
!             . la viscosite       (initialisee a VISCL0)
!      - dans usiniv :
!             . les variables de calcul  (initialisees a 0 par defaut
!             ou a la valeur donnee dans usiniv)

! On peut donner ici les lois de variation aux cellules
!     - de la masse volumique                      ROM    kg/m3
!         (et eventuellememt aux faces de bord     ROMB   kg/m3)
!     - de la viscosite moleculaire                VISCL  kg/(m s)
!     - de la chaleur specifique associee          CP     J/(kg degres)
!     - des "diffusivites" associees aux scalaires VISCLS kg/(m s)


! On dispose des types de faces de bord au pas de temps
!   precedent (sauf au premier pas de temps, ou les tableaux
!   ITYPFB et ITRIFB n'ont pas ete renseignes)




! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nvar             ! i  ! <-- ! total number of variables                      !
! nscal            ! i  ! <-- ! total number of scalars                        !
! mbrom            ! te ! <-- ! indicateur de remplissage de romb              !
!        !    !     !                                                !
! izfppp           ! te ! <-- ! numero de zone de la face de bord              !
! (nfabor)         !    !     !  pour le module phys. part.                    !
! dt(ncelet)       ! ra ! <-- ! time step (per cell)                           !
! rtp              ! ra ! <-- ! calculated variables at cell centers           !
!  (ncelet, *)     !    !     !  (at current time steps)                       !
! propce(ncelet, *)! ra ! <-- ! physical properties at cell centers            !
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
use optcal
use cstnum
use cstphy
use entsor
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

integer          mbrom
integer          izfppp(nfabor)

double precision dt(ncelet), rtp(ncelet,nflown:nvar)
double precision propce(ncelet,*)

! Local variables

integer          iel
integer          ipcray
integer          ipcvsl, ith   , iscal , ii
integer          iiii  , ipcsig, it
integer          iesp  , iesp1 , iesp2 , mode , isrrom

double precision tp    , delt  , somphi, val
double precision alpro , alpvis, alpcp , alpsig, alplab , alpkab
double precision rhonp1
double precision ym    (ngazgm),yvol  (ngazgm)
double precision coef(ngazgm,ngazgm)
double precision roesp (ngazgm),visesp(ngazgm),cpesp(ngazgm)
double precision sigesp(ngazgm),xlabes(ngazgm),xkabes(ngazgm)
double precision, dimension(:), pointer ::  crom
double precision, dimension(:), pointer :: viscl, cpro_cp

integer          ipass
data             ipass /0/
save             ipass

!===============================================================================
!===============================================================================
! 0 - INITIALISATIONS A CONSERVER
!===============================================================================

! Initialize variables to avoid compiler warnings

ipcvsl = 0
ipcsig = 0
ipcray = 0

! --- Initialisation memoire


ipass = ipass + 1

!     Sous relaxation de la masse volumique (pas au premier pas de temps)
if(ntcabs.gt.1.and.srrom.gt.0.d0) then
  isrrom = 1
else
  isrrom = 0
endif

!===============================================================================
! 1 - EFFET JOULE
!===============================================================================

!  -- Les lois doivent etre imposees par l'utilisateur
!       donc on ne fait rien.

!      IF ( IPPMOD(IELJOU).GE.1 ) THEN


!  -- Attention, dans les modules electriques, la chaleur massique, la
!       conductivite thermique et la conductivite electriques sont
!       toujours dans le tableau PROPCE
!       qu'elles soient physiquement variables ou non.

!       On n'utilisera donc PAS les variables
!          =====================
!                                CP0, VISLS0(ISCALT)
!                                VISLS0(IPOTR) et VISLS0(IPOTI)

!       Informatiquement, ceci se traduit par le fait que
!                                ICP>0, IVISLS(ISCALT)>0,
!                                IVISLS(IPOTR)>0 et IVISLS(IPOTI)>0

!       Les verifications ont ete faites dans elveri

!  -- Si la conductivite electrique est toujours la meme pour
!       le potentiel reel et le potentiel imaginaire, on pourrait
!       n'en avoir qu'une seule (modif dans varpos pour definir
!       IVISLS(IPOTI) = IVISLS(IPOTR)) et economiser NCEL reels .

!      call field_get_val_s(icrom, crom)
!      call field_get_val_s(iprpfl(iviscl), viscl)
!      call field_get_val_s(iprpfl(icp), cpro_cp)
!      IPCVSL = IPPROC(IVISLS(ISCALT))
!      IPCSIR = IPPROC(IVISLS(IPOTR))
!      IPCSII = IPPROC(IVISLS(IPOTI))

!      PROPCE(IEL,IPPROC(ITEMP)) =
!      crom(iel) =
!      viscl(iel) =
!      cpro_cp(iel) =
!      PROPCE(IEL,IPCVSL) =
!      PROPCE(IEL,IPCSIR) =
!      PROPCE(IEL,IPCSII) =

!      ENDIF

!===============================================================================
! 2 - ARC ELECTRIQUE
!===============================================================================

if ( ippmod(ielarc).ge.1 ) then

!      Un message une fois au moins pour dire
!                                    qu'on prend les valeurs sur fichier
  if(ipass.eq.1) then
    write(nfecra,1000)
  endif

!      Calcul de la temperature a partir de l'enthalpie

  mode = 1

  if ( ngazg .eq. 1 ) then
    ym(1) = 1.d0
    mode = 1
    do iel = 1, ncel
      call elthht(mode,ngazg,ym,rtp(iel,isca(iscalt)),               &
                                propce(iel,ipproc(itemp)))
    enddo
  else
    do iel = 1, ncel
      ym(ngazg) = 1.d0
      do iesp = 1, ngazg-1
        ym(iesp) = rtp(iel,isca(iycoel(iesp)))
        ym(ngazg) = ym(ngazg) - ym(iesp)
      enddo
      call elthht(mode,ngazg,ym,rtp(iel,isca(iscalt)),               &
                                propce(iel,ipproc(itemp)))
    enddo
  endif

!      Pointeurs pour les differentes variables

  call field_get_val_s(icrom, crom)
  call field_get_val_s(iprpfl(iviscl), viscl)
  if(icp.gt.0) call field_get_val_s(iprpfl(icp), cpro_cp)
  if(ivisls(iscalt).gt.0) then
    ipcvsl = ipproc(ivisls(iscalt))
  endif
  if ( ivisls(ipotr).gt.0 ) then
    ipcsig = ipproc(ivisls(ipotr))
  endif
  if ( ixkabe .gt. 0 ) then
     ipcray = ipproc(idrad)
  endif

!       Interpolation des donnees sur le fichier de donnees
!         en fonction de la temperature

  do iel = 1, ncel

!        Valeur de la temperature

    tp = propce(iel,ipproc(itemp))

!        On determine  le IT ou il faut interpoler

    it = 0
    if ( tp .le. th(1) ) then
      it = 1
    else if ( tp .ge. th(npo) ) then
      it = npo
    else
      do iiii = 1, npo-1
        if ( tp .gt. th(iiii) .and. tp .le. th(iiii+1) ) then
          it = iiii
        endif
      enddo
    endif
    if ( it .eq. 0 ) then
      write(nfecra,9900) tp
      call csexit(1)
    endif

!        Fraction massique

    ym(ngazg) = 1.d0
    do iesp = 1, ngazg-1
      ym(iesp)  = rtp(iel,isca(iycoel(iesp)))
      ym(ngazg) = ym(ngazg) - ym(iesp)
    enddo

!     Masse volumique, Viscosite, CP, Sigm et Lambda de chaque constituant

    if ( tp .le. th(1) ) then

!         Extrapolation : Valeur constante = 1ere valeur de la table

      do iesp = 1, ngazg
        roesp (iesp) = rhoel (iesp,1)
        visesp(iesp) = visel (iesp,1)
        cpesp (iesp) = cpel  (iesp,1)
        sigesp(iesp) = sigel (iesp,1)
        xlabes(iesp) = xlabel(iesp,1)
        if ( ixkabe .gt. 0 ) then
          xkabes(iesp) = xkabel(iesp,1)
        endif
      enddo

    else if ( tp .ge. th(npo) ) then

!         Extrapolation : valeur constante = derniere valeur de la table

      do iesp = 1, ngazg
        roesp (iesp) = rhoel (iesp,npo)
        visesp(iesp) = visel (iesp,npo)
        cpesp (iesp) = cpel  (iesp,npo)
        sigesp(iesp) = sigel (iesp,npo)
        xlabes(iesp) = xlabel(iesp,npo)
        if ( ixkabe .gt. 0 ) then
          xkabes(iesp) = xkabel(iesp,npo)
        endif
      enddo

    else

!         Interpolation

      delt = th(it+1) - th(it)
      do iesp = 1, ngazg

!          Masse volumique de chaque constituant

        alpro = (rhoel(iesp,it+1)-rhoel(iesp,it))/delt
        roesp(iesp)  = rhoel(iesp,it) + alpro*(tp-th(it))

!          Viscosite de chaque constituant

        alpvis = (visel(iesp,it+1)-visel(iesp,it))/delt
        visesp(iesp) = visel(iesp,it) + alpvis*(tp-th(it))

!          CP de chaque constituant

        alpcp = (cpel(iesp,it+1)-cpel(iesp,it))/delt
        cpesp(iesp) = cpel(iesp,it) + alpcp*(tp-th(it))

!          Conductivite electrique (Sigma) de chaque constituant

        alpsig = (sigel(iesp,it+1)-sigel(iesp,it))/delt
        sigesp(iesp) = sigel(iesp,it) + alpsig*(tp-th(it))

!          Conductivite thermique (Lambda) de chaque constituant

        alplab = (xlabel(iesp,it+1)-xlabel(iesp,it))/delt
        xlabes(iesp) = xlabel(iesp,it) + alplab*(tp-th(it))

!          Emission nette radiative ou Terme source radiatif
!          de chaque constituant

        if ( ixkabe .gt. 0 ) then
          alpkab = (xkabel(iesp,it+1)-xkabel(iesp,it))/delt
          xkabes(iesp) = xkabel(iesp,it) + alpkab*(tp-th(it))
        endif

      enddo

    endif

!       Masse volumique du melange (sous relaxee eventuellement)
!       ==========================

    rhonp1 = 0.d0
    do iesp = 1, ngazg
      rhonp1 = rhonp1+ym(iesp)/roesp(iesp)
    enddo
    rhonp1 = 1.d0/rhonp1
    if(isrrom.eq.1) then
      crom(iel) =                                        &
           srrom*crom(iel)+(1.d0-srrom)*rhonp1
    else
      crom(iel) = rhonp1
    endif

!        Fraction volumique de chaque constituant

    do iesp = 1, ngazg
      yvol(iesp) = ym(iesp)*roesp(iesp)/crom(iel)
      if ( yvol(iesp) .le. 0.d0 ) yvol(iesp) = epzero**2
    enddo

!       Viscosite moleculaire dynamique en kg/(m s)
!       ==========================================

    do iesp1 = 1, ngazg
      do iesp2 = 1, ngazg
        coef(iesp1,iesp2) = ( 1.d0                                &
               + sqrt(visesp(iesp1)/visesp(iesp2))                &
                *sqrt(sqrt(roesp(iesp2)/roesp(iesp1))))**2.       &
                / sqrt( 1.d0 + roesp(iesp1)/roesp(iesp2) )        &
            / sqrt(8.d0)
      enddo
    enddo

    viscl(iel) = 0.d0
    do iesp1=1,ngazg

      somphi = 0.d0
      do iesp2=1,ngazg
        if ( iesp1 .ne. iesp2 ) then
          somphi = somphi                                         &
                  +coef(iesp1,iesp2)*yvol(iesp2)/yvol(iesp1)
        endif
      enddo

      viscl(iel) = viscl(iel)                     &
                          +visesp(iesp1)/(1.d0+somphi)

    enddo

!       Chaleur specifique J/(kg degres)
!       ================================

    if(icp.gt.0) then

      cpro_cp(iel) = 0.d0
      do iesp = 1, ngazg
        cpro_cp(iel) = cpro_cp(iel) + ym(iesp)*cpesp(iesp)
      enddo

    endif

!       Lambda/Cp en kg/(m s)
!       ---------------------

    if(ivisls(iscalt).gt.0) then

      do iesp1=1,ngazg
        do iesp2=1,ngazg
          coef(iesp1,iesp2) = ( 1.d0                              &
                 + sqrt(xlabes(iesp1)/xlabes(iesp2))              &
                  *sqrt(sqrt(roesp(iesp2)/roesp(iesp1))))**2.d0   &
                  / sqrt( 1.d0 + roesp(iesp1)/roesp(iesp2) )      &
              / sqrt(8.d0)
        enddo
      enddo

!        On calcule d'abord juste Lambda

      propce(iel,ipcvsl) = 0.d0
      do iesp1=1,ngazg

        somphi = 0.d0
        do iesp2=1,ngazg
          if ( iesp1 .ne. iesp2 ) then
            somphi = somphi + coef(iesp1,iesp2)*yvol(iesp2)/yvol(iesp1)
          endif
        enddo

        propce(iel,ipcvsl) = propce(iel,ipcvsl)                   &
                            +xlabes(iesp1)/(1.d0+1.065*somphi)

      enddo

!        On divise par CP pour avoir Lambda/CP
!          On suppose Cp renseigne au prealable.

      if(icp.le.0) then

! --- Si CP est uniforme, on utilise CP0

        propce(iel,ipcvsl) = propce(iel,ipcvsl)/cp0

      else

! --- Si CP est non uniforme, on utilise le CP calcul au dessus
        propce(iel,ipcvsl) = propce(iel,ipcvsl)/cpro_cp(iel)

      endif
    endif

!       Conductivite electrique en S/m
!       ==============================

    if ( ivisls(ipotr).gt.0 ) then
      propce(iel,ipcsig) = 0.d0
      val = 0.d0
      do iesp=1,ngazg
        val = val + yvol(iesp)/sigesp(iesp)
      enddo

      propce(iel,ipcsig) = 1.d0/val
    endif

!       Emission nette radiative en W/m3
!       ================================

    if ( ixkabe .gt. 0 ) then
      propce(iel,ipcray) = 0.d0
      val = 0.d0
      do iesp=1,ngazg
        val = val + yvol(iesp)*xkabes(iesp)
      enddo

      propce(iel,ipcray) = val
    endif

  enddo

!       Diffusivite variable a l'exclusion de l'enthalpie et de IPOTR
!       -------------------------------------------------------------
!         Il n'y a pas d'autres scalaires, et la boucle ne fait donc rien


  do ii = 1, nscapp

! --- Numero du scalaire
    iscal = iscapp(ii)

! --- Si il s'agit de l'enthalpie son cas a deja ete traite plus haut
    ith = 0
    if (iscal.eq.iscalt) ith = 1

! --- Si il s'agit de Potentiel (IPOTR), son cas a deja ete traite
    if (iscal.eq.ipotr) ith = 1

! --- Si la variable est une fluctuation, sa diffusivite est
!       la meme que celle du scalaire auquel elle est rattachee :
!       il n'y a donc rien a faire ici : on passe directement
!       a la variable suivante sans renseigner PROPCE(IEL,IPCVSL).

    if ( ith.eq.0 .and. iscavr(iscal).le.0) then

! --- On ne traite ici que les variables non thermiques
!                                        et pas le potentiel (sigma)
!                                   et qui ne sont pas des fluctuations

      if(ivisls(iscal).gt.0) then

! --- Rang de Lambda du scalaire
!     dans PROPCE, prop. physiques au centre des elements       : IPCVSL

        ipcvsl = ipproc(ivisls(iscal))

        do iel = 1, ncel
          propce(iel,ipcvsl) = 1.d0
        enddo

      endif

    endif

  enddo

endif

!===============================================================================
! 3 - CONDUCTION IONIQUE
!===============================================================================

! POUR LE MOMENT CETTE OPTION N'EST PAS ACTIVEE

if ( ippmod(ielion).ge.1  ) then

!       Masse volumique
!       ---------------

  call field_get_val_s(icrom, crom)
  do iel = 1, ncel
    crom(iel) = 1.d0
  enddo

!       VISCOSITE
!       =========

  call field_get_val_s(iprpfl(iviscl), viscl)
  do iel = 1, ncel
    viscl(iel) = 1.d-2
  enddo

!       CHALEUR SPECIFIQUE VARIABLE J/(kg degres)
!       =========================================

  if(icp.gt.0) then
    call field_get_val_s(iprpfl(icp), cpro_cp)
    do iel = 1, ncel
      cpro_cp(iel) = 1000.d0
    enddo

  endif

!       Lambda/CP  VARIABLE en kg/(m s)
!       ===============================

  if (ivisls(iscalt).gt.0) then

    ipcvsl = ipproc(ivisls(iscalt))

    if(icp.le.0) then

! --- Si CP est uniforme, on utilise CP0

      do iel = 1, ncel
        propce(iel,ipcvsl) = 1.d0/cp0
      enddo

    else

! --- Si CP est non uniforme, on utilise PROPCE ci dessus
      do iel = 1, ncel
        propce(iel,ipcvsl) = 1.d0 /cpro_cp(iel)
      enddo

    endif

  endif

!       DIFFUSIVITE VARIABLE A L'EXCLUSION DE L'ENTHALPIE
!       ==================================================

  do ii = 1, nscapp

! --- Numero du scalaire
    iscal = iscapp(ii)

! --- Si il s'agit de l'enthqlpie son cas a deja ete traite plus haut
    ith = 0
    if (iscal.eq.iscalt) ith = 1

! --- Si la variable est une fluctuation, sa diffusivite est
!       la meme que celle du scalaire auquel elle est rattachee :
!       il n'y a donc rien a faire ici : on passe directement
!       a la variable suivante sans renseigner PROPCE(IEL,IPCVSL).

    if ( ith.eq.0 .and. iscavr(iscal).le.0) then

! --- On ne traite ici que les variables non thermiques
!                                   et qui ne sont pas des fluctuations

      if(ivisls(iscal).gt.0) then

! --- Rang de Lambda du scalaire
!     dans PROPCE, prop. physiques au centre des elements       : IPCVSL

        ipcvsl = ipproc(ivisls(iscal))

! --- Lambda en kg/(m s) au centre des cellules


        do iel = 1, ncel
          propce(iel,ipcvsl) = 1.d0
        enddo

      endif

    endif

  enddo

endif

!===============================================================================
! 4 - ON PASSE LA MAIN A L'UTILISATEUR (joule en particulier)
!===============================================================================

call uselph                                                       &
!==========
 ( nvar   , nscal  ,                                              &
   mbrom  , izfppp ,                                              &
   dt     )



! La masse volumique au bord est traitee dans phyvar (recopie de la valeur
!     de la cellule de bord).

!--------
! FORMATS
!--------

 1000 format(/,                                                   &
' Module electrique: proprietes physiques lues sur fichier',/)
 9900 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR DANS ELPHYV (MODULE ELECTRIQUE)      ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Tabulation echoue avec une temperature TP = ', E14.5       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

!----
! FIN
!----

return
end subroutine
