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
! Function:
! ---------

!> \file colecd.f90
!>
!> \brief Specific physic subroutine: gas combustion
!>
!> Read data file dp_C3P* for gas combustion
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
subroutine colecd

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
use cs_c_bindings

!===============================================================================

implicit none

! Arguments

! Local variables

character(len=150) :: chain1,chain2,chain3
character(len=12) :: nomgaz

integer          it, igg, ir, ige, iat, ios, igf, igo, igp, iehc
integer          ncgm, nrgm
integer          inicoe, inicha
integer          lonch, ichai, ichcoe
integer          iereac(ngazem)
integer          ncoel, icoel

double precision tmin, tmax
double precision kabse(ngazem)
double precision wmolce(ngazem)
double precision cpgaze(ngazem,npot)
double precision atgaze(ngazem, natom)
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
  write(nfecra, 9900) indjon
  call csexit(1)
endif

!===============================================================================
! -1. UTILISATION DE JANAF
!===============================================================================

if (indjon.eq.1) then

  ! 1.1 LECTURE DU FICHIER DONNEES SPECIFIQUES
  !===========================================

  ! --> Ouverture du fichier

  open(unit=impfpp, file=ficfpp, status='old', form='formatted',    &
       access='sequential', iostat=ios, err=99)
  rewind(unit=impfpp, err=99)

  ! --> On se limite pour l'instant a de la combustion gaz
  !     Plus exactement a une flamme de diffusion chimie 3 corps
  !                  et a une flamme de premelange modele EBU
  !     --------------------------------------------------------

  ! --> Lecture des donnees

  ! ---- Nb de constituants gazeux elementaires

  read(impfpp, *, err=999, end=999 ) ngaze
  if (ngaze.gt.ngazem) then
    write(nfecra,9980) ngazem, ngaze
    call csexit(1)
  endif

  ! ---- Nb de points de tabulation ENTH-TEMP

  read(impfpp, *, err=999, end=999) npo
  if (npo.gt.npot) then
    write(nfecra,9981) npot, npo
    call csexit(1)
  endif

  ! --- Temperature Min et Max

  read(impfpp,*,err=999,end=999) tmin
  read(impfpp,*,err=999,end=999) tmax

  ! ---- Lecture des noms des constituants elementaires gazeux

  do ige=1 , ngaze
    do inicoe=1,len(nomcoe(ige))
      nomcoe(ige)(inicoe:inicoe) = ' '
    enddo
  enddo

  do inicha=1,len(chain1)
    chain1(inicha:inicha) = ' '
  enddo

  do inicha=1,len(chain2)
    chain2(inicha:inicha) = ' '
  enddo

  read(impfpp, *, err=999, end=999)
  read(impfpp, 1010, err=999, end=999) chain1
  chain2 = adjustl(chain1)
  lonch = len(trim(chain2))

  ige = 1
  ichcoe = 0
  do ichai = 1, lonch
    if (chain2(ichai:ichai).ne.' ') then
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

  read(impfpp, *, err=999, end=999) (kabse(ige), ige = 1, ngaze)

  ! ---- Nb especes atomiques (C, H, O, N, ...)

  read(impfpp, *, err=999, end=999) nato
  if (nato.gt.natom) then
    write(nfecra,9983) natom, nato
    call csexit(1)
  endif

  ! ---- Masse molaire especes atomiques
  !      Composition des especes courantes en fonction des especes
  !        elementaires

  do iat = 1, nato
    read(impfpp,*,err=999,end=999 ) wmolat(iat),                 &
         (atgaze(ige, iat), ige = 1, ngaze)
  enddo

  ! ---- Nb especes globales

  ! lecture du nombre d'especes globales dont la composition est connue
  ! souvent fuel et oxydant
  read(impfpp,*,err=999,end=999 ) ngazg
  ! On ne considere qu'UNE SEULE REACTION GLOBALE
  !   NGAZG = NCGM = 3 par consequent (F, O, P)
  ncgm = 3

  ! soit l'utilisateur equilibre la reaction, 3 especes globales sont alors lues
  ! soit on calcule l'equilibre et l'utilisateur doit donner les deux especes
  ! globales reactives
  if (ngazg.lt.2.or.ngazg.gt.ncgm) then
    write(nfecra,9985) ncgm, ngazg
    call csexit (1)
  endif

  ! -----------------------------------------------------------------------------
  ! a. Si les trois especes globales sont connues, l'utilisateur fournit aussi la
  !    stoechiometrique en nombre de mol d'especes globales
  ! -----------------------------------------------------------------------------

  if (ngazg.eq.3) then

    ! Composition des especes globales en fonction des especes courantes
    do igg = 1, ngazg
      read (impfpp, *, err=999, end=999) (compog(ige,igg), ige = 1, ngaze)
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
      read(impfpp, *, err=999, end=999)                             &
           igfuel(ir), igoxy(ir), (stoeg(igg,ir), igg = 1, ngazg)
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

  ! --- Effective Heat of Combustion (J/kg)

  pcigas = 0.d0
  read(impfpp, '(a)', iostat=ios) chain1
  if (ios.eq.0) then
    iehc = index(chain1,"EHC")
    if (iehc.gt.0) read(chain1(iehc+3:), *) pcigas
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
    do igg = 1, ngazg
      nmolg = 0.d0
      do ige = 1, ngaze
        nmolg = nmolg + compog(ige,igg)
      enddo
      if (nmolg.eq.0.d0) then
        write(nfecra, 9988) igg, nmolg
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

  ! --- Enthalpy calculation of elementary species

  ! if user specifies a EHC, fuel enthalpy is not computed from PPTBHT
  ! but set to zero waiting its calculation from EHC
  if (pcigas.gt.0.d0) then
    ncoel = ngaze-1
    icoel = 2
    do it = 1, npo
      ehgaze(1,it) = 0.d0
    enddo
    ! fuel must be placed in first position in the elementary species
    if (compog(1, igfuel(1)).ne.1.d0) then
      write(nfecra, 9986) trim(nomcoe(1))
      call csexit(1)
    endif
  ! else, it is computed
  else
    ncoel = ngaze
    icoel = 1
  endif

  do ige = icoel, ngaze
    wmolce(ige) = wmole(ige)
  enddo

  call pptbht(ncoel, nomcoe(icoel), ehgaze(icoel,1), cpgaze(icoel,1), wmolce(icoel))

  ! --- Masses of global species (becoming molar masses below)

  do igg = 1, ngazg
    wmolg(igg) = 0.d0
    do ige = 1, ngaze
      wmolg(igg) = wmolg(igg)+compog(ige,igg)*wmole(ige)
    enddo
  enddo

  ! --- Absorption coefficients of global species in case of radiation calculation

  if ((ippmod(icod3p).eq.1.or.                           &
       ippmod(icoebu).eq.1.or.ippmod(icoebu).eq.3.or.    &
       ippmod(icolwc).eq.1.or.ippmod(icolwc).eq.3.or.    &
       ippmod(icolwc).eq.5 ).and.iirayo.ge.1) then
    do igg = 1, ngazg
      ckabsg(igg) = 0.d0
      do ige = 1, ngaze
        ckabsg(igg) = ckabsg(igg) + compog(ige,igg)*kabse(ige)*wmole(ige)
      enddo
      ckabsg(igg) = ckabsg(igg)/wmolg(igg)
    enddo
  endif

  ! --- Enthalpies and mass heat capacity of global species

  do igg = icoel , ngazg
    do it = 1, npo
      ehgazg(igg,it) = 0.d0
      cpgazg(igg,it) = 0.d0
      do ige = 1, ngaze
        ehgazg(igg,it) = ehgazg(igg,it)                           &
             + compog(ige,igg)*wmole(ige)*ehgaze(ige,it)
        cpgazg(igg,it) = cpgazg(igg,it)                           &
             + compog(ige,igg)*wmole(ige)*cpgaze(ige,it)
      enddo
      ehgazg(igg,it) = ehgazg(igg,it)/wmolg(igg)
      cpgazg(igg,it) = cpgazg(igg,it)/wmolg(igg)
    enddo
  enddo

  ! --- Molar masses of global species

  do igg = 1, ngazg
    nmolg = 0.d0
    do ige = 1, ngaze
      nmolg = nmolg + compog(ige,igg)
    enddo
    if (nmolg.eq.0.d0) then
      write(nfecra, 9988) igg, nmolg
      call csexit(1)
    endif
    wmolg(igg) = wmolg(igg) / nmolg
    do ige = 1, ngaze
      compog(ige,igg) = compog(ige,igg) / nmolg
    enddo
  enddo

  ! --- Estimation of fuel enthalpy in case of user EHC

  if (pcigas.gt.0.d0) then
    do it = 1, npo
      ehgazg(1,it) = 0.d0
      do igg = icoel, ngazg
        ehgazg(1,it) = ehgazg(1,it)                           &
             + stoeg(igg,1)*wmolg(igg)*ehgazg(igg,it)
      enddo
      ehgazg(1,it) = pcigas - ehgazg(1,it)/(wmolg(1)*stoeg(1,1))
    enddo
  endif

  ! --- Molar coefficients XCO2 , XH2O

  do ige = 1, ngaze
    nomgaz = nomcoe(ige)
    if (trim(nomgaz).EQ.'C(S)') IIC=IGE
    if (trim(nomgaz).EQ.'CO'  ) IICO=IGE
    if (trim(nomgaz).EQ.'O2'  ) IIO2=IGE
    if (trim(nomgaz).EQ.'CO2' ) IICO2=IGE
    if (trim(nomgaz).EQ.'H2O' ) IIH2O=IGE
  enddo

  xco2 = compog(iico2,3)
  xh2o = compog(iih2o,3)

  ! Balance verification and stoechiometric ratio for each reaction

  do ir = 1, nrgaz
    do iat = 1, nato
      bilan = 0.d0
      do igg = 1, ngazg
        do ige = 1, ngaze
          bilan = bilan + stoeg(igg,ir)*compog(ige,igg)*atgaze(ige,iat)
        enddo
      enddo
      if (abs(bilan) .gt. epsi) then
        write(nfecra,9991) ir, iat, bilan
        call csexit(1)
      endif
    enddo
    mfuel = stoeg(igfuel(ir),ir)*wmolg(igfuel(ir))
    moxyd = stoeg(igoxy(ir),ir)*wmolg(igoxy(ir))
    mreac = mfuel + moxyd
    fs(ir) = mfuel/mreac
  enddo

  ! Calcul des coefficients de la fraction massique d'oxydant
  ! du programme pdflwc

  coeff1 = compog(2,2)*wmole(2)/wmolg(2)
  coeff3 = moxyd / mfuel
  coeff2 = coeff3*coeff1

  ! Conversion coefficients from global species to elementary species

  do igg = 1, ngazg
    do ige = 1, ngaze
      coefeg(ige,igg) = compog(ige,igg)*wmole(ige)/wmolg(igg)
    enddo
  enddo

  ! PCI calculation

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

      efgaz(igg) = cs_gas_combustion_t_to_h(coefg, tgaz)

      pcigas = pcigas + stoeg(igg,ir)*wmolg(igg)*efgaz(igg)

    enddo

    ! PCI dimension is J/kg of combustible
    pcigas = pcigas / (stoeg(1,ir)*wmolg(1))

  enddo

!===============================================================================
! -2. UTILISATION D'UNE TABULATION ENTHALPIE-TEMPERATURE
!===============================================================================

else

  open (unit=impfpp, file=ficfpp, status='old', form='formatted',  &
        access='sequential', iostat=ios, err=99)
  rewind (unit=impfpp, err=99)

  ! Nb de points de tabulation ENTH-TEMP

  read(impfpp, *, err=999, end=999) npo
  if (npo.gt.npot) then
    write(nfecra, 9981) npot, npo
    call csexit (1)
  endif

  ! Tabulation ENTH-TEMP pour les especes globales

  do it = 1, npo
    read(impfpp, *, err=999, end=999) th(it),                       &
                 ehgazg(1,it),ehgazg(2,it),ehgazg(3,it)
  enddo

  ! On ne considere qu'UNE SEULE REACTION GLOBALE
  ! NGAZG = NCGM = 3 par consequent (F, O, P)
  ngazg = 3

  ! On ne considere qu'UNE SEULE REACTION GLOBALE
  nrgaz = 1

  ! Masses molaires pour les especes globales

  read (impfpp,*,err=999,end=999) wmolg(1),wmolg(2),wmolg(3)

  ! Fraction de melange a la stoechiometrie

  read (impfpp,*,err=999,end=999) fs(1)

  ! Rayonnement

  ! FIXME ce test n'a plus sa place ici, iirayo n'est plus fixe
  ! dans le fichier de thermochimie

  if (iirayo.gt.0 .and. ippmod(icod3p).ne.1                      &
      .and. ippmod(icoebu).ne.1 .and. ippmod(icoebu).ne.3        &
      .and. ippmod(icolwc).ne.1 .and. ippmod(icolwc).ne.3        &
      .and. ippmod(icolwc).ne.5) then
    write(nfecra,9982)                                           &
         iirayo,ippmod(icod3p),ippmod(icoebu),ippmod(icolwc)
    call csexit(1)
  endif

  ! Coefficients d'absorption des especes globales

  read (impfpp,*,err=999,end=999) ckabsg(1),ckabsg(2),ckabsg(3)

  ! Coefficients molaires de CO2 et H2O dans les produits
  ! (pour le rayonnement)

  read (impfpp,*,err=999,end=999) xco2, xh2o

  ! Fermeture du fichier

  close(impfpp)

  ! Calcul des coefficients de la fraction massique d oxydant
  ! on considère que l'oxydant est un mélange d'O2 et N2

  coeff1 = ((wmolg(2)-0.028)/(0.032-0.028))* (0.032/wmolg(2))

  coeff3 = (1-fs(1))/fs(1)

  coeff2 = coeff3*coeff1

  ! Conversion coefficients from global species to elementary species

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

      efgaz(igg) = cs_gas_combustion_t_to_h(coefg, tgaz)

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
write(nfecra, 9992)
call csexit(1)

  999 continue
write(nfecra, 9993)
call csexit(1)

!--------
! Formats
!--------

 9900 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The INDJON indicator must have value 1 (use Janaf)',        /,&
'@    or 0 (user tabulation).',                                 /,&
'@',                                                            /,&
'@  Its current value is', i10                                  /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9980 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The number of current species must be ', i10,               /,&
'@   Its value is ', i10, ' in the parameters file.',           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9981 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The number of tabulation points is limited to ', i10,       /,&
'@   Its value is ', i10, ' in the parameters file.',           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9982 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The radiation model can only be activated with a',          /,&
'@   combustion model in permeatic conditions.',                /,&
'@',                                                            /,&
'@  A radiation model was specified:',                          /,&
'@   iirayo = ', i10,                                           /,&
'@  But we have:',                                              /,&
'@   ippmod(icod3p) = ', i10,                                   /,&
'@   ippmod(icoebu) = ', i10,                                   /,&
'@   ippmod(icolwc) = ', i10,                                   /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9983 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The number of elementary species must be ', i10,            /,&
'@   Its value is ', i10, ' in the parameters file.',           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9984 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The number of current species is set to ', i10,             /,&
'@  but only ', i10, ' current species are listed',             /,&
'@  in the parameters file.',                                   /,&
'@',                                                            /,&
"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@",/,&
"@                                                            ",/)
 9985 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The number of global species can only be equal to ', i10,   /,&
'@  (automatic reaction balance) or ', i10,                     /,&
'@  (composition of known products).',                          /,&
'@   Its value is ', i10, ' in the parameters file.',           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9986 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  Fuel must be placed at first place in the list of'          /,&
'@  elementary species. First elementary specie is now ', a6,   /,&
'@  in the parameters file',                                    /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9987 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The species whose part is specified must be positionned',   /,&
'@  last among the products.',                                  /,&
'@  The elementary species ', a6, 'is at position ', i6,        /,&
'@  in the compositions of the global species product in the',  /,&
'@  parameters file.',                                          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9988 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The number of moles in the global species ', i6,            /,&
'@  is ', f10.4, '.',                                           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9989 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The reactive elementary species ', i10, 'is not present',   /,&
'@  in the global species ', i10, '.',                          /,&
'@  Its ratio is ', f10.4, ' in the parameters file.',          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9990 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  The number of global reactions must be ', i10,              /,&
'@   Its value is ', i10, ' in the parameters file.',           /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9991 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  Conservation problem encountered in reaction ', i10,        /,&
'@  for element ', i10,                                         /,&
'@  The molar balance gives ', f10.4, ' mol, whereas it should,',/&
'@  be zero.',                                                  /,&
'@',                                                            /,&
'@  Check the parameters file.',                                /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

 9992 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  Error opening the parameters file.',                        /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)
 9993 format(                                                     &
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /,&
'@ @@ ERROR:   STOP WHILE READING INPUT DATA (COLECD)',         /,&
'@    =====',                                                   /,&
'@             GAS COMBUSTION',                                 /,&
'@',                                                            /,&
'@  Error reading the parameters file.',                        /,&
'@    It may be incomplete or incorrectly formatted.',          /,&
'@',                                                            /,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@',                                                            /)

end subroutine colecd


