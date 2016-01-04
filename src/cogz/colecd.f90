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

subroutine colecd
!================

!===============================================================================
!  FONCTION  :
!  ---------

! LECTURE DU FICHIER DE DONNEES PHYSIQUE PARTICULIERE
!            RELATIF A LA COMBUSTION GAZ

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
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
use pointe
use entsor
use cstnum
use cstphy
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat

!===============================================================================

implicit none

! Arguments

! Local variables

character(len=150) :: chain1,chain2
character(len=12) :: nomcoe(ngazem)

integer          it, igg, ir, ige, iat, ios
integer          ncgm, nrgm
integer          inicoe , inicha
integer          idebch , ifinch , lonch
integer          ichai , ichcoe
integer          atgaze(ngazem,natom)
integer          igoxy(nrgazm), igfuel(nrgazm)
integer          ncoel
integer          iico2 , iih2o
integer          mode, icalck

double precision tmin , tmax
double precision kabse(ngazem)
double precision compog(ngazem,ngazgm)
double precision wmolce (ngazem)
double precision cpgaze(ngazem,npot)
double precision coefg(ngazgm), tgaz, efgaz(ngazgm)
double precision mfuel, mreac, epsi, nmolg, bilan
double precision moxyd

!===============================================================================
!===============================================================================
! -0. INITIALISATION et VERIFICATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

iico2 = 0
iih2o = 0

mfuel = 0.d0
moxyd = 0.d0

epsi = 1.d-9


if (indjon.ne.1.and.indjon.ne.0) then
  write(nfecra,9900) indjon
  call csexit (1)
  !==========
endif

!===============================================================================
! -1. UTILISATION DE JANAF
!===============================================================================

if (indjon.eq.1) then


! 1.1 LECTURE DU FICHIER DONNEES SPECIFIQUES
!===================================================

! --> Ouverture du fichier

  open ( unit=impfpp, file=ficfpp,                                &
          STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL',    &
                                          iostat=ios, err=99 )
  rewind ( unit=impfpp,err=99 )

! --> On se limite pour l'instant a de la combustion gaz
!     Plus exactement a une flamme de diffusion chimie 3 corps
!                  et a une flamme de premelange modele EBU
!     --------------------------------------------------------

! --> Lecture des donnees

! ---- Nb de constituants gazeux elementaires

  read ( impfpp,*,err=999,end=999 ) ngaze
  if ( ngaze.gt.ngazem ) then
    write(nfecra,9990) ngazem,ngaze
    call csexit (1)
    !==========
  endif

! ---- Nb de points de tabulation ENTH-TEMP

  read ( impfpp,*,err=999,end=999 ) npo
  if ( npo.gt.npot ) then
    write(nfecra,9991) npot, npo
    call csexit (1)
    !==========
  endif

! --- Temperature Min et Max

  read (impfpp,*,err=999,end=999) tmin
  read (impfpp,*,err=999,end=999) tmax

! ---- Lecture des noms des constituants elementaires gazeux

  do ige=1 , ngaze
    do inicoe=1,len(nomcoe(ige))
      NOMCOE(IGE)(INICOE:INICOE)=' '
    enddo
  enddo

  do inicha=1,len(chain1)
    CHAIN1(INICHA:INICHA)=' '
  enddo

  do inicha=1,len(chain2)
    CHAIN2(INICHA:INICHA)=' '
  enddo

  read ( impfpp,*,err=999,end=999 )
  read ( impfpp,1010,err=999,end=999 ) chain1
  call verlon (chain1 , idebch , ifinch , lonch)
  chain2(1:lonch)=chain1(idebch:ifinch)

  ige   =1
  ichcoe=0
  do ichai = 1, lonch
    IF (CHAIN2(ICHAI:ICHAI).NE.' ') THEN
      ichcoe=ichcoe+1
      nomcoe(ige)(ichcoe:ichcoe)=chain2(ichai:ichai)
    else
      if (ichcoe.ne.0) then
        ige=ige+1
        ichcoe=0
      endif
    endif
  enddo

 1010   format(a150)


! ---- Rayonnement
! Le rayonnement n'est autorise qu'avec des modeles permeatiques,
! car en adiabatique l'enthalpie est une grandeur algebrique qui
! ne prend pas en compte les pertes par rayonnement.

  if ( iirayo.gt.0 .and. ippmod(icod3p).ne.1                      &
       .and. ippmod(icoebu).ne.1 .and. ippmod(icoebu).ne.3        &
       .and. ippmod(icolwc).ne.1 .and. ippmod(icolwc).ne.3        &
       .and. ippmod(icolwc).ne.5 ) then
    write(nfecra,9993)                                            &
         iirayo,ippmod(icod3p),ippmod(icoebu),ippmod(icolwc)
    call csexit (1)
    !==========
  endif

! ---- Coefficient d'absorption des especes courantes

  read (impfpp,*,err=999,end=999 )                                &
                          ( kabse(ige),ige=1,ngaze )

! ---- Nb especes atomiques (C, H, O, N, ...)

  read (impfpp,*,err=999,end=999 ) nato
  if ( nato.gt.natom ) then
    write(nfecra,9994) natom,nato
    call csexit (1)
    !==========
  endif

! ---- Masse molaire especes atomiques
!      Composition des especes courantes en fonction des especes
!        elementaires

  do iat = 1, nato
    read (impfpp,*,err=999,end=999 ) wmolat(iat),                 &
                      ( atgaze(ige,iat),ige=1,ngaze )
  enddo

! ---- Nb especes globales

  read (impfpp,*,err=999,end=999 ) ngazg
!       On ne considere qu'UNE SEULE REACTION GLOBALE
!         NGAZG = NCGM = 3 par consequent (F, O, P)
  ncgm = 3
  if ( ngazg.ne.ncgm ) then
    write(nfecra,9995) ncgm,ngazg
    call csexit (1)
    !==========
  endif

! ---- Composition des especes globales en fonction des especes
!        courantes

  do igg = 1, ngazg
    read (impfpp,*,err=999,end=999 )                              &
                       ( compog(ige,igg),ige=1,ngaze )
  enddo

! ----- Nb de reactions globales

  read (impfpp,*,err=999,end=999 ) nrgaz
!       On ne considere qu'UNE SEULE REACTION GLOBALE
  nrgm = 1
  if ( nrgaz.ne.nrgm ) then
    write(nfecra,9996) nrgm,nrgaz
    call csexit (1)
    !==========
  endif

! ---- No des especes concernees par la rapport stoechio
!       Stoechio en especes globales des reactions

  do ir = 1, nrgaz
    read (impfpp,*,err=999,end=999 )                              &
        igfuel(ir),igoxy(ir),( stoeg(igg,ir),igg=1,ngazg )
  enddo

! --> Fermeture du fichier

  close(impfpp)

! 1.2 CALCULS DE DONNEES COMPLEMENTAIRES
!===============================================

! ---- Calcul des masses molaires des especes courantes

  do ige = 1, ngaze
    wmole(ige) = 0.d0
    do iat = 1, nato
      wmole(ige)= wmole(ige) + atgaze(ige,iat)*wmolat(iat)
    enddo
  enddo

! --- Discretisation de la temperature

  do it = 1, npo
    th(it) = dble(it-1)*(tmax-tmin)/dble(npo-1)+tmin
  enddo

! ---Calcul des enthalpies par appel a la subroutine PPTBHT

  ncoel  = ngaze
  do ige = 1, ngaze
    wmolce(ige) = wmole (ige)
  enddo

  call pptbht                                                     &
  !==========
 ( ncoel ,                                                        &
  nomcoe , ehgaze , cpgaze , wmolce )

! ---- Calcul des masses molaires des especes globales
!          de la tabulation temperature - enthalpie massique
!          et des coefficients d'absorption des especes globales
!             si RAYONNEMENT

  icalck = 0
  if ((ippmod(icod3p).eq.1.or.                           &
       ippmod(icoebu).eq.1.or.ippmod(icoebu).eq.3.or.    &
       ippmod(icolwc).eq.1.or.ippmod(icolwc).eq.3.or.    &
       ippmod(icolwc).eq.5 ).and.iirayo.ge.1) then
     icalck = 1
  endif

  do igg = 1 , ngazg
    wmolg(igg) = 0.d0
    nmolg      = 0.d0
    do ige = 1 , ngaze
      wmolg(igg) = wmolg(igg)+compog(ige,igg)*wmole(ige)
      nmolg      = nmolg     +compog(ige,igg)
    enddo
    do it = 1,npo
      ehgazg(igg,it) = 0.d0
      cpgazg(igg,it) = 0.d0
      if (icalck.eq.1) ckabsg(igg) = 0.d0
      do ige = 1 , ngaze
        ehgazg(igg,it) = ehgazg(igg,it)                           &
             + compog(ige,igg)*wmole(ige)*ehgaze(ige,it)
        cpgazg(igg,it) = cpgazg(igg,it)                           &
             + compog(ige,igg)*wmole(ige)*cpgaze(ige,it)
        if (icalck.eq.1) ckabsg(igg) = ckabsg(igg)                &
                         + compog(ige,igg)*kabse(ige)*wmole(ige)
      enddo
      ehgazg(igg,it) = ehgazg(igg,it)/wmolg(igg)
      cpgazg(igg,it) = cpgazg(igg,it)/wmolg(igg)
      if (icalck.eq.1) ckabsg(igg) = ckabsg(igg)/wmolg(igg)
    enddo
    wmolg(igg) = wmolg(igg) / nmolg
    do ige = 1 , ngaze
      compog(ige,igg) = compog(ige,igg) / nmolg
    enddo
  enddo

! ---- Calcul des coefficients molaires XCO2 , XH2O

  do ige = 1 , ngaze
    IF (NOMCOE(IGE)(1:3).EQ.'CO2') IICO2=IGE
    IF (NOMCOE(IGE)(1:3).EQ.'H2O') IIH2O=IGE
  enddo

  xco2 =  compog(iico2,3)
  xh2o =  compog(iih2o,3)

! ---- Calcul bilan pour verification
!        et taux de melange a la stochio pour chaque reaction

  do ir = 1, nrgaz
    do iat = 1, nato
      bilan = 0.d0
      do igg = 1, ngazg
        do ige = 1, ngaze
          bilan = bilan                                           &
              + stoeg(igg,ir)*compog(ige,igg)*atgaze(ige,iat)
        enddo
      enddo
      if ( abs(bilan) .gt. epsi ) then
        write(nfecra,9997) ir, iat
        call csexit (1)
        !==========
      endif
    enddo
    mfuel = stoeg(igfuel(ir),ir)*wmolg(igfuel(ir))
    moxyd = stoeg(igoxy(ir),ir)*wmolg(igoxy(ir))
    mreac = mfuel + moxyd
    fs(ir) = mfuel/mreac
  enddo

! ---- Calcul des coefficients de la fraction massique d'oxydant
!      du programme pdflwc

  coeff1 = compog(2,2)*wmole(2)/wmolg(2)

  coeff3 =   moxyd / mfuel

  coeff2 = coeff3*coeff1

  ! --- Conversion coefficients from global species to elementary species

  do igg = 1, ngazg
    do ige = 1, ngaze
     coefeg(ige,igg) = compog(ige,igg)*wmole(ige)/wmolg(igg)
    enddo
  enddo

  ! --- PCI calculation

  ! gas name storage
  namgas = nomcoe(1)

  pcigas = 0.d0

  do ir = 1, nrgaz

    do igg = 1, ngazg

      ! enthalpies of formation
      coefg(1)  = 0.d0
      coefg(2)  = 0.d0
      coefg(3)  = 0.d0
      coefg(igg) = 1.d0
      tgaz      = 300.d0

      mode = -1
      call cothht                                                   &
      !==========
        ( mode   , ngazg , ngazgm  , coefg  ,                     &
          npo    , npot   , th     , ehgazg ,                     &
          efgaz(igg)      , tgaz   )

      pcigas = pcigas + stoeg(igg,ir)*wmolg(igg)*efgaz(igg)

    enddo

    ! PCI dimension is J/kg of combustible
    pcigas = pcigas / (stoeg(1,ir)*wmolg(1))

  enddo

!===============================================================================
! -2. UTILISATION D'UNE TABULATION ENTHALPIE-TEMPERATURE
!===============================================================================

else

  open ( unit=impfpp, file=ficfpp,                                &
          STATUS='OLD', FORM='FORMATTED', ACCESS='SEQUENTIAL',    &
                                          iostat=ios, err=99 )
  rewind ( unit=impfpp,err=99 )

! ---- Nb de points de tabulation ENTH-TEMP

  read ( impfpp,*,err=999,end=999 ) npo
  if ( npo.gt.npot ) then
    write(nfecra,9991) npot, npo
    call csexit (1)
    !==========
  endif

! --- Tabulation ENTH-TEMP pour les especes globales

  do it = 1, npo
    read (impfpp,*,err=999,end=999) th(it),                       &
                 ehgazg(1,it),ehgazg(2,it),ehgazg(3,it)
  enddo

!       On ne considere qu'UNE SEULE REACTION GLOBALE
!         NGAZG = NCGM = 3 par consequent (F, O, P)
  ngazg = 3

!       On ne considere qu'UNE SEULE REACTION GLOBALE
  nrgaz = 1

! --- Masses molaires pour les especes globales

  read (impfpp,*,err=999,end=999) wmolg(1),wmolg(2),wmolg(3)

! --- Fraction de melange a la stoechiometrie

  read (impfpp,*,err=999,end=999) fs(1)

! --- Rayonnement

  if ( iirayo.gt.0 .and. ippmod(icod3p).ne.1                      &
       .and. ippmod(icoebu).ne.1 .and. ippmod(icoebu).ne.3        &
       .and. ippmod(icolwc).ne.1 .and. ippmod(icolwc).ne.3        &
       .and. ippmod(icolwc).ne.5 ) then
    write(nfecra,9993)                                            &
         iirayo,ippmod(icod3p),ippmod(icoebu),ippmod(icolwc)
    call csexit (1)
    !==========
  endif

! --- Coefficients d'absorption des especes globales

  read (impfpp,*,err=999,end=999) ckabsg(1),ckabsg(2),ckabsg(3)

! --- Coefficients molaires de CO2 et H2O dans les produits
!      (pour le rayonnement)

  read (impfpp,*,err=999,end=999) xco2, xh2o


! ---> Fermeture du fichier

  close (impfpp)


! ---- Calcul des coefficients de la fraction massique d oxydant
!      on considère que l'oxydant est un mélange d'O2 et N2

  coeff1 = ((wmolg(2)-0.028)/(0.032-0.028))* (0.032/wmolg(2))

  coeff3 = (1-fs(1))/fs(1)

  coeff2 = coeff3*coeff1

  ! --- Conversion coefficients from global species to elementary species

  do igg = 1, ngazg
    do ige = 1, ngaze
     coefeg(ige,igg) = compog(ige,igg)*wmole(ige)/wmolg(igg)
    enddo
  enddo

  ! --- PCI calculation

  ! gas name storage
  namgas = nomcoe(1)

  pcigas = 0.d0

  do ir = 1, nrgaz

    do igg = 1, ngazg

      ! enthalpies of formation
      coefg(1)  = 0.d0
      coefg(2)  = 0.d0
      coefg(3)  = 0.d0
      coefg(igg) = 1.d0
      tgaz      = 300.d0

      mode = -1
      call cothht                                                   &
      !==========
        ( mode   , ngazg , ngazgm  , coefg  ,                     &
          npo    , npot   , th     , ehgazg ,                     &
          efgaz(igg)      , tgaz   )

      pcigas = pcigas + stoeg(igg,ir)*wmolg(igg)*efgaz(igg)

    enddo

    ! dimension is J/kg of combustible
    pcigas = pcigas / (stoeg(1,ir)*wmolg(1))

  enddo

endif


return

!============================
! 3. SORTIE EN ERREUR
!============================

  99  continue
write ( nfecra,9998 )
call csexit (1)
!==========

  999 continue
write ( nfecra,9999 )
call csexit (1)
!==========


!--------
! FORMATS
!--------


 9900 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  L''indicateur INDJON doit etre un entier de valeur        ',/,&
'@    1 (utilisation de Janaf) ou 0 (tabulation utilisateur). ',/,&
'@                                                            ',/,&
'@  Il vaut ici ',I10                                          ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier INDJON dans usppmo.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9990 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  Le nombre d''especes courantes doit etre ',I10             ,/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9991 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  Le nombre de points de tabulation est limite a ',I10       ,/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9993 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  LE RAYONNEMENT NE PEUT ETRE ACTIVE QU''AVEC UN MODELE DE  ',/,&
'@   COMBUSTION EN CONDITIONS PERMEATIQUES.                   ',/,&
'@                                                            ',/,&
'@  Un mode de rayonnement a ete specifie                     ',/,&
'@   IIRAYO = ',I10                                            ,/,&
'@  Or dans usppmo on a :                                     ',/,&
'@   IPPMOD(ICOD3P) = ',I10                                    ,/,&
'@   IPPMOD(ICOEBU) = ',I10                                    ,/,&
'@   IPPMOD(ICOLWC) = ',I10                                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique et usppmo.               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9994 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  Le nombre d''especes elementaires est limite a ',I10       ,/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9995 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  Le nombre d''especes globales doit etre ',I10              ,/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9996 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  Le nombre de reactions globales doit etre ',I10            ,/,&
'@   Il vaut ',I10   ,' dans le fichier parametrique          ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9997 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  Probleme de conservation rencontre dans la                ',/,&
'@   reaction ',I10   ,' pour l''element ',I10                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9998 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  Erreur a l''ouverture du fichier parametrique.            ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
'@  Erreur a la lecture du fichier parametrique.              ',/,&
'@    Le fichier a ete ouvert mais est peut etre incomplet    ',/,&
'@    ou son format inadapte.                                 ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

end subroutine


