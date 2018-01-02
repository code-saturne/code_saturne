!-------------------------------------------------------------------------------

! This file is part of Code_Saturne, a general-purpose CFD tool.
!
! Copyright (C) 1998-2018 EDF S.A.
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

character(len=150) :: chain1,chain2,chain3
character(len=12) :: nomgaz

integer          it, igg, ir, ige, iat, ios, igf, igo, igp
integer          ncgm, nrgm
integer          inicoe , inicha
integer          idebch , ifinch , lonch
integer          ichai , ichcoe
integer          atgaze(ngazem,natom)
integer          iereac(ngazem)
integer          ncoel
integer          mode, icalck

double precision tmin , tmax
double precision kabse(ngazem)
double precision wmolce (ngazem)
double precision cpgaze(ngazem,npot)
double precision coefg(ngazgm), tgaz, efgaz(ngazgm)
double precision mfuel, mreac, epsi, nmolg, bilan
double precision moxyd

double precision, dimension(:,:), allocatable :: aa
double precision, dimension(:), allocatable :: bb, xx

!===============================================================================
!===============================================================================
! -0. INITIALISATION et VERIFICATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings

iio2  = 0
iico  = 0
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
    write(nfecra,9980) ngazem,ngaze
    call csexit (1)
    !==========
  endif

! ---- Nb de points de tabulation ENTH-TEMP

  read ( impfpp,*,err=999,end=999 ) npo
  if ( npo.gt.npot ) then
    write(nfecra,9981) npot, npo
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

  if (ige.lt.ngaze) then
    write(nfecra,9984) ngaze, ige
    call csexit(1)
  endif

  1010 format(a150)


! ---- Rayonnement

! FIXME ce test n'a plus sa place ici, iirayo n'est plus fixe
! dans le fichier de thermochimie

! Le rayonnement n'est autorise qu'avec des modeles permeatiques,
! car en adiabatique l'enthalpie est une grandeur algebrique qui
! ne prend pas en compte les pertes par rayonnement.

  if ( iirayo.gt.0 .and. ippmod(icod3p).ne.1                      &
       .and. ippmod(icoebu).ne.1 .and. ippmod(icoebu).ne.3        &
       .and. ippmod(icolwc).ne.1 .and. ippmod(icolwc).ne.3        &
       .and. ippmod(icolwc).ne.5 ) then
    write(nfecra,9982)                                            &
         iirayo,ippmod(icod3p),ippmod(icoebu),ippmod(icolwc)
    call csexit (1)
  endif

! ---- Coefficient d'absorption des especes courantes

  read (impfpp,*,err=999,end=999 ) ( kabse(ige),ige=1,ngaze )

! ---- Nb especes atomiques (C, H, O, N, ...)

  read (impfpp,*,err=999,end=999 ) nato
  if ( nato.gt.natom ) then
    write(nfecra,9983) natom,nato
    call csexit (1)
  endif

! ---- Masse molaire especes atomiques
!      Composition des especes courantes en fonction des especes
!        elementaires

  do iat = 1, nato
    read (impfpp,*,err=999,end=999 ) wmolat(iat),                 &
                      ( atgaze(ige,iat),ige=1,ngaze )
  enddo

! ---- Nb especes globales

  ! lecture du nombre d'especes globales dont la composition est connue
  ! souvent fuel et oxydant
  read (impfpp,*,err=999,end=999 ) ngazg
  ! On ne considere qu'UNE SEULE REACTION GLOBALE
  !   NGAZG = NCGM = 3 par consequent (F, O, P)
  ncgm = 3

  ! soit l'utilisateur equilibre la reaction, 3 especes globales sont alors lues
  ! soit on calcule l'equilibre et l'utilisateur doit donner les deux especes
  ! globales reactives
  if (ngazg.lt.2.or.ngazg.gt.ncgm) then
    write(nfecra,9985) ncgm,ngazg
    call csexit (1)
  endif

  ! -----------------------------------------------------------------------------
  ! a. Si les trois especes globales sont connues, l'utilisateur fournit aussi la
  !    stoechiometrique en nombre de mol d'especes globales
  ! -----------------------------------------------------------------------------

  if (ngazg.eq.3) then

    ! ---- Composition des especes globales en fonction des especes
    !      courantes
    do igg = 1, ngazg
      read (impfpp,*,err=999,end=999) (compog(ige,igg),ige=1,ngaze)
    enddo

    ! ----- Nb de reactions globales
    read (impfpp,*,err=999,end=999) nrgaz
    ! On ne considere qu'UNE SEULE REACTION GLOBALE
    nrgm = 1
    if (nrgaz.ne.nrgm) then
      write(nfecra,9990) nrgm,nrgaz
      call csexit(1)
    endif

    ! ---- No des especes concernees par la rapport stoechio
    !       Stoechio en especes globales des reactions
    ! igfuel(ir) : espece global fuel    de la reaction ir
    ! igoxy (ir) : espece global oxydant de la reaction ir
    ! igprod(ir) : espece global produit de la reaction ir
    do ir = 1, nrgaz
      read (impfpp,*,err=999,end=999 )                              &
          igfuel(ir),igoxy(ir),( stoeg(igg,ir),igg=1,ngazg )
    enddo
    ! Produit
    igprod(1) = 3

  ! -----------------------------------------------------------------------------
  ! b. si l'utilisateur ne fournit que les especes reactives, la reaction est
  !    equilibre automatiquement
  ! -----------------------------------------------------------------------------
  else

    ! --- Composition des especes globales en fonction des especes courantes
    ! igfuel(ir) : espece global fuel    de la reaction ir
    ! igoxy (ir) : espece global oxydant de la reaction ir
    ! igprod(ir) : espece global produit de la reaction ir
    ! iereac(igg): espece elementaire reactive de l'espece globale igg

    ! Fuel
    igfuel(1) = 1
    iereac(igfuel(1)) = 1
    read (impfpp,*,err=999,end=999) (compog(ige,igfuel(1)),ige=1,ngaze)
    if (compog(iereac(igfuel(1)),igfuel(1)).eq.0.d0) then
      write(nfecra,9989) iereac(igfuel(1)), igfuel(1), compog(iereac(igfuel(1)),igfuel(1))
      call csexit(1)
    endif

    ! Oxydant
    igoxy(1) = 2
    iereac(igoxy(1)) = 2
    read (impfpp,*,err=999,end=999) (compog(ige,igoxy(1)),ige=1,ngaze)
    if (compog(iereac(igoxy(1)),igoxy(1)).eq.0.d0) then
      write(nfecra,9989) iereac(igoxy(1)), igoxy(1), compog(iereac(igoxy(1)),igoxy(1))
      call csexit(1)
    endif

    ! Produit
    igprod(1) = 3
    iereac(igprod(1)) = 0
    ! La composition est calculee
    ! Si on a plus de 3 especes elementaires dans les produits, on a besoin de
    ! donnees supplementaires (part du produit en fonction du combustible)
    read (impfpp,*,err=999,end=999) (compog(ige,igprod(1)),ige=1,ngaze)

    ! Ces especes doivent etre positionnees a la fin pour eviter les pivots nuls
    ! lors de l'inversion de la matrice (on pourrait peut etre les deplacer...)
    do ige = 1, ngaze
      if (compog(ige,igprod(1)).gt.0.d0.and.ige.le.5) then
        write(nfecra,9987) trim(nomcoe(ige)), ige
        call csexit(1)
      endif
    enddo

    ! Nb de reactions globales : ON NE PEUT EQUILIBRER QU'UNE SEULE REACTION
    nrgaz = 1

  endif

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

  ! --- Calcul de la compostion des produits et les coefficients
  !     stoechiometriques

  igf = igfuel(1)
  igo = igoxy(1)
  igp = igprod(1)

  ! Automatic equilibrium of the reaction
  ! -------------------------------------

  if (ngazg.eq.2) then

    ! On rajoute l'espece globale produit
    ngazg = 3

    ! Composition des especes connues dans les produits, donnee par unite de
    ! masse combustible
    do ige = 1, ngaze
      compog(ige,igp) = compog(ige,igp) * wmole(igf)/wmole(ige)
    enddo

    ! Calcul des coefficients stoechiometriques
    allocate(aa(nato,ngaze),bb(nato),xx(ngaze))

    do iat = 1, nato
      do ige = 1, ngaze
        aa(iat,ige) = 0.d0
        xx(ige) = 0.d0
      enddo
      bb(iat) = 0.d0
    enddo

    do iat = 1, nato
      do ige = 1, ngaze
        ! fuel
        bb(iat) = bb(iat) + atgaze(ige,iat)*compog(ige,igf)
        ! oxydiser
        aa(iat,iereac(igo)) = aa(iat,iereac(igo)) - atgaze(ige,iat)*compog(ige,igo)
        ! products
        bb(iat) = bb(iat) - atgaze(ige,iat)*compog(ige,igp)
        if ((ige.ne.iereac(igf).and.ige.ne.iereac(igo)).and.   &
            compog(ige,igp).eq.0.d0                          ) then
          aa(iat,ige) = aa(iat,ige) + atgaze(ige,iat)
        endif
      enddo
    enddo

    call gauss(nato,nato,aa(1,2),xx(2),bb)

    ! maintenant, on connait les coefficients stoechio des especes globales
    nreact(igf) = -1.d0
    nreact(igo) = -xx(iereac(igo))
    nreact(igp) = 1.d0

    do igg = 1, ngazg
      ! les especes globales reactives sont deja connues
      if (igg.ne.igf.and.igg.ne.igo) then
        do ige = 1, ngaze
          ! les especes elementaires reactives aussi
          if (ige.ne.iereac(igf).and.ige.ne.iereac(igo).and. &
              compog(ige,igg).eq.0.d0                       ) then
            compog(ige,igg) = abs(xx(ige))
          endif
        enddo
      endif
    enddo

    ! Stoechio en especes globales des reactions
    do ir = 1, nrgaz
      do igg = 1, ngazg
        stoeg(igg,ir) = 0.d0
        do ige = 1, ngaze
          stoeg(igg,ir) = stoeg(igg,ir) + compog(ige,igg)*nreact(igg)
        enddo
      enddo
    enddo

    deallocate(aa,bb,xx)

  ! Equilibrium defined by the user
  ! -------------------------------
  else

    ! global stoechiometric coefficient, to write the reaction
    do igg = 1 , ngazg
      nmolg      = 0.d0
      do ige = 1 , ngaze
        nmolg      = nmolg + compog(ige,igg)
      enddo
      if (nmolg.eq.0.d0) then
        write(nfecra,9988) igg, nmolg
        call csexit(1)
      endif
      nreact(igg) = stoeg(igg,1)/nmolg
    enddo

  endif

  ! Reaction writing
  ! ----------------

  ! fuel
  chain1 = ""
  write(chain1,'(a)') trim(nomcoe(igf))
  nomcog(igf) = trim(chain1)

  ! oxydiser
  chain2 = ""
  do ige = 1, ngaze
    if (compog(ige,igo).eq.1.d0) then
      write(chain2,1) trim(chain2), trim(nomcoe(ige)),"+"
    elseif (compog(ige,igo).gt.0.d0) then
      write(chain2,2) trim(chain2), compog(ige,igo), trim(nomcoe(ige)),"+"
    endif
  enddo
  write(chain2,'(a)') trim(chain2(2:len(trim(chain2))-1))
  nomcog(igo) = trim(chain2)

  ! product
  chain3 = ""
  do ige = 1, ngaze
    if (compog(ige,igp).gt.0.d0) then
      write(chain3,2) trim(chain3), compog(ige,igp), trim(nomcoe(ige)), "+"
    endif
  enddo
  write(chain3,'(a)') trim(chain3(2:len(trim(chain3))-1))
  nomcog(igp) = trim(chain3)

  1 format(a,1x,a,1x,a)
  2 format(a,1x,f6.3,1X,a,1x,a)

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
    if (nmolg.eq.0.d0) then
      write(nfecra,9988) igg, nmolg
      call csexit(1)
    endif
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
    nomgaz = nomcoe(ige)
    IF (trim(nomgaz).EQ.'C(S)') IIC=IGE
    IF (trim(nomgaz).EQ.'CO'  ) IICO=IGE
    IF (trim(nomgaz).EQ.'O2'  ) IIO2=IGE
    IF (trim(nomgaz).EQ.'CO2' ) IICO2=IGE
    IF (trim(nomgaz).EQ.'H2O' ) IIH2O=IGE
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
        write(nfecra,9991) ir, iat
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
    write(nfecra,9981) npot, npo
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

  ! FIXME ce test n'a plus sa place ici, iirayo n'est plus fixe
  ! dans le fichier de thermochimie

  if ( iirayo.gt.0 .and. ippmod(icod3p).ne.1                      &
       .and. ippmod(icoebu).ne.1 .and. ippmod(icoebu).ne.3        &
       .and. ippmod(icolwc).ne.1 .and. ippmod(icolwc).ne.3        &
       .and. ippmod(icolwc).ne.5 ) then
    write(nfecra,9982)                                            &
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
write ( nfecra,9992 )
call csexit (1)
!==========

  999 continue
write ( nfecra,9993 )
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
 9980 format(                                                           &
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
 9981 format(                                                           &
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
 9982 format(                                                           &
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
 9983 format(                                                           &
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
 9984 format(                                                           &
"@                                                            ",/,&
"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@",/,&
"@                                                            ",/,&
"@ @@ ATTENTION : ARRET A L'ENTREE DES DONNEES (COLECD)       ",/,&
"@    =========                                               ",/,&
"@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ",/,&
"@                                                            ",/,&
"@  Le nombre d'especes courantes est fixe a ",I10             ,/,&
"@  or seulement ",I10   ," especes courantes sont listees    ",/,&
"@  dans le fichier parametrique.                             ",/,&
"@                                                            ",/,&
"@  Le calcul ne sera pas execute.                            ",/,&
"@                                                            ",/,&
"@  Verifier le fichier parametrique.                         ",/,&
"@                                                            ",/,&
"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@",/,&
"@                                                            ",/)
 9985 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
"@  Le nombre d'especes globales ne peut etre egal qu'a ",I10  ,/,&
"@  (equilibre automatique de la reaction) ou ",I10            ,/,&
"@  (composition des produits connue).                        ",/,&
"@  Il vaut ",I10," dans le fichier parametrique.             ",/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9987 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
"@  Les especes dont la part est specifiee doivent etre       ",/,&
"@  positionnees en dernier dans les produits.                ",/,&
"@  L'espece elementaire ",a6," est a la position ",i6         ,/,&
"@  dans la composition de l'espece globale produit dans le   ",/,&
"@  fichier parametrique.                                     ",/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9988 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
"@  Le nombre de moles dans l'espece globale ",i6,"           ",/,&
"@  vaut ",f10.4," .                                          ",/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9989 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES (COLECD)      ',/,&
'@    =========                                               ',/,&
'@      PHYSIQUE PARTICULIERE (COMBUSTION GAZ)                ',/,&
'@                                                            ',/,&
"@  L'espece elementaire reactive ", I10," n'est pas presente ",/,&
"@  dans l'espece globale ",I10,".                            ",/,&
"@  Sa proportion vaut ",f10.4," dans le fichier parametrique.",/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
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
'@  Le nombre de reactions globales doit etre ',I10            ,/,&
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
'@  Probleme de conservation rencontre dans la                ',/,&
'@   reaction ',I10   ,' pour l''element ',I10                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier le fichier parametrique.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 9992 format(                                                           &
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
 9993 format(                                                           &
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


