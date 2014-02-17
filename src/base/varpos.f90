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

!> \file varpos.f90
!> \brief Variables location initialization, according to calculation type
!> selected by the user.
!>
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Arguments
!------------------------------------------------------------------------------
!   mode          name          role
!------------------------------------------------------------------------------
!______________________________________________________________________________

subroutine varpos

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use numvar
use optcal
use cstphy
use cstnum
use entsor
use albase
use lagpar
use lagdim
use lagran
use parall
use ppppar
use ppthch
use coincl
use cpincl
use ppincl
use radiat
use ihmpre
use mesh
use field

!===============================================================================

implicit none

! Arguments

! Local variables

character     rubriq*64, cmoy4*4, cindfm*4
character     ficsui*32
character*80  f_label, f_name, s_label, s_name
integer       ivar  , iscal , irphas , iprop, id
integer       ii    , jj    , kk    , ll
integer       iok   , ippok , ipppst
integer       iest  , ivisph, ipropp
integer       imom  , jmom  , imold , jmold , jbmomt, jtcabs
integer       iiplus, iimoin, jmomok, idto  , jdto
integer       ierror, irtyp , itysup, nbval
integer       idffin, idfmji, imomr , ilecec
integer       ilsmom, ivers , inmfin
integer       impamx
integer       nfmtmo, nberro
integer       idtold(nbmomx)
integer       nprayc
integer       idttur

double precision gravn2

integer       ipass
data          ipass /0/
save          ipass

integer       nprmax
data          nprmax /0/
save          nprmax

integer       nppmax
data          nppmax /0/
save          nppmax

!===============================================================================
! 0. INITIALISATIONS
!===============================================================================

! Initialize variables to avoid compiler warnings
ippok = 0

ipass = ipass + 1

!     Nombre max pour les formats choisis
nfmtmo = 9999
!     Indefini a 4 caracteres
cindfm = 'YYYY'

!===============================================================================
! 1. PREMIER APPEL : VERIFICATIONS ET
!        POSITIONNEMENT DES VARIABLES  : IPR, IU ... ISCA, NVAR
!                       ET DES PROPRIETES PPALES
!===============================================================================

if (ipass.eq.1) then

! ---> 1.1 VERIFICATIONS
!      -----------------

  iok = 0

! --- ISCAVR(ISCAL) doit etre compris entre 0 et NSCAL.

  if (nscaus.gt.0) then
    do ii = 1, nscaus
      iscal = ii
      if (iscavr(iscal).gt.nscal.or.iscavr(iscal).lt.0) then
        write(nfecra,7030) iscal,ii,ii,iscavr(iscal),nscal
        iok = iok + 1
      endif
    enddo
  endif
  if (nscapp.gt.0) then
    do ii = 1, nscapp
      iscal = iscapp(ii)
      if (iscavr(iscal).gt.nscal.or.iscavr(iscal).lt.0) then
        write(nfecra,7031) iscal,ii,ii,iscavr(iscal),nscal
        iok = iok + 1
      endif
    enddo
  endif


! --- IVISLS(ISCAL) doit etre non initialise pour les variances
!     Il prend la valeur du scalaire associe
!     Tous les tests qui suivent sont utilises pour le message
!       d'erreur eventuel.

  if (nscaus.gt.0) then
    do jj = 1, nscaus
      ii    = jj
      iscal = iscavr(ii)
!       Si on a une variance avec ivisls initialise : erreur
      if    ( (iscal.gt.0.and.iscal.le.nscal).and.                &
               ivisls(ii).ne.-1                        ) then
        ll = 0
        do kk = 1, nscaus
          if (       kk .eq.iscal) ll = kk
        enddo
        do kk = 1, nscapp
          if (iscapp(kk).eq.iscal) ll = -kk
        enddo
        if (ll.gt.0) then
          write(nfecra,7040)                                      &
               ii,ii,jj,iscal,ll,jj,iscal,jj,ivisls(iscal)
        else
          write(nfecra,7041)                                      &
               ii,ii,jj,iscal,-ll,jj,iscal,jj,ivisls(iscal)
        endif
        iok = iok + 1
!     Si on n'a pas une variance std mais que ivisls est incorrect : erreur
      else if ( (iscal.le.0 .or.iscal.gt.nscal).and.                &
             (ivisls(ii).ne.0.and.ivisls(ii).ne.1) ) then
        write(nfecra,7050) ii,jj,jj,ivisls(ii)
        iok = iok + 1
      endif
    enddo
  endif

  if (nscapp.gt.0) then
    do jj = 1, nscapp
      ii    = iscapp(jj)
      iscal = iscavr(ii)
!       Si on a une variance avec ivisls initialise : erreur
      if    ( (iscal.gt.0.and.iscal.le.nscal).and.                &
               ivisls(ii).ne.-1                        ) then
        ll = 0
        do kk = 1, nscaus
          if (       kk .eq.iscal) ll = kk
        enddo
        do kk = 1, nscapp
          if (iscapp(kk).eq.iscal) ll = -kk
        enddo
        if (ll.gt.0) then
          write(nfecra,7042)                                      &
               ii,ii,jj,iscal,ll,jj,iscal,jj,ivisls(iscal)
        else
          write(nfecra,7043)                                      &
               ii,ii,jj,iscal,-ll,jj,iscal,jj,ivisls(iscal)
        endif
        iok = iok + 1
!       Si on n'a pas une variance std mais que ivisls est incorrect : erreur
      else if ( (iscal.le.0 .or.iscal.gt.nscal).and.                &
              (ivisls(ii).ne.0.and.ivisls(ii).ne.1) ) then
        write(nfecra,7051)ii,jj,jj,ivisls(ii)
        iok = iok + 1
      endif
    enddo
  endif

!       On initialise les ivisls des variances
  if (nscal.gt.0) then
    do ii = 1, nscal
      iscal = iscavr(ii)
      if (iscal.gt.0.and.iscal.le.nscal) then
        ivisls(ii) = ivisls(iscal)
        if (ivisls(ii).gt.0) then
          call field_set_key_int(ivarfl(isca(ii)), kivisl,    &
                                 iprpfl(ivisls(iscal)))
        endif
      endif
    enddo
  endif


! ---> VISCOSITE ALE
  if (iale.eq.1) then
    if (iortvm.ne.0 .and. iortvm.ne.1) then
      write(nfecra,7070) iortvm
      iok = iok + 1
    endif
  endif
! --- On s'arrete si quelque chose s'est mal passe

  if (iok.ne.0) then
    call csexit (1)
  endif


! ---> 1.2 POSITIONNEMENT DES PROPRIETES PRINCIPALES
!      --------------------------------

! --- Numerotation des proprietes presentes ici
!       Ceci depend du type de solveur branche derriere
!        (CP, Poly, Mono...)
!       Dans l'ideal, il y aurait donc plusieurs varpos.

!       Pour le moment, on fait les hypotheses suivantes :
!         Il y a toujours, pour toutes les phases,  rho, viscl, visct
!         Il y a toujours la pression totale (sauf en compressible)
!         Lorsqu'elles sont variables, on a les proprietes suivantes :
!           . cp    (par phase)
!           . visls (par scalaire)
!           . csmago (par phase) en LES dynamique
!         En ALE on a la viscosite de maillage
!         On a aussi les flux de masse porteurs :
!           . les variables u,v,w,p,turbulence sont portees par leur
!               phase (1 flux)
!           . les scalaires sont portes par leur phase (1 flux)
!           On suppose donc qu'il n'y a pas de scalaire en
!             taux de vide a convecter avec un flux particulier (ce
!             serait possible : ca rajoute un test, par exemple
!             if alpro...

!     ATTENTION : on s'arrange pour numeroter a la queue-leu-leu sans
!       trous les proprietes qui sont definies au centre des cellules
!       ceci permet ensuite de ne pas se fatiguer lors de la
!       construction de IPPPRO plus bas.
!      Cependant, pour les physiques particulieres, ce n'est pas le cas.

  ! Base properties, always present

  call add_property_field('density', 'Density', irom)
  icrom = iprpfl(irom)

  call add_property_field('molecular_viscosity', 'Laminar Viscosity', &
                          iviscl)

  call add_property_field('turbulent_viscosity', 'Turb Viscosity', &
                          ivisct)
  if (iturb.eq.0) call hide_property(ivisct)

  call add_property_field('courant_number', 'CFL', icour)
  call add_property_field('fourier_number', 'Fourier Number', ifour)

  ! Pression totale stockee dans IPRTOT, si on n'est pas en compressible
  ! (sinon Ptot=P* !)
  if (ippmod(icompf).lt.0) then
    call add_property_field('total_pressure', 'Total Pressure', &
                            iprtot)
  endif

  ! CP when variable
  if (icp.ne.0) then
    call add_property_field('specific_heat', 'Specific Heat', icp)
  endif

  ! Cs^2 si on est en LES dynamique
  if (iturb.eq.41) then
    call add_property_field('smagorinsky_constant^2', 'Csdyn2', ismago)
  else
    ismago = 0
  endif

  ! Viscosite de maillage en ALE
  if (iale.eq.1) then
    call add_property_field('mesh_viscosity_1', 'Mesh ViscX', ivisma(1))
    ! si la viscosite est isotrope, les trois composantes pointent
    !  au meme endroit
    if (iortvm.eq.0) then
      ivisma(2) = ivisma(1)
      ivisma(3) = ivisma(1)
    else
      call add_property_field('mesh_viscosity_2', 'Mesh ViscY', ivisma(2))
      call add_property_field('mesh_viscosity_3', 'Mesh ViscZ', ivisma(3))
    endif
  endif

  ! Estimateurs d'erreur
  do iest = 1, nestmx
    iestim(iest) = -1
  enddo

  if (iescal(iespre).gt.0) then
    write(f_name,  '(a14,i1)') 'est_error_pre_', iescal(iespre)
    write(f_label, '(a5,i1)') 'EsPre', iescal(iespre)
    call add_property_field(f_name, f_label, iestim(iespre))
  endif
  if (iescal(iesder).gt.0) then
    write(f_name,  '(a14,i1)') 'est_error_der_', iescal(iesder)
    write(f_label, '(a5,i1)') 'EsDer', iescal(iesder)
    call add_property_field(f_name, f_label, iestim(iesder))
  endif
  if (iescal(iescor).gt.0) then
    write(f_name,  '(a14,i1)') 'est_error_cor_', iescal(iescor)
    write(f_label, '(a5,i1)') 'EsCor', iescal(iescor)
    call add_property_field(f_name, f_label, iestim(iescor))
  endif
  if (iescal(iestot).gt.0) then
    write(f_name,  '(a14,i1)') 'est_error_tot_', iescal(iestot)
    write(f_label, '(a5,i1)') 'EsTot', iescal(iestot)
    call add_property_field(f_name, f_label, iestim(iestot))
  endif

!   Proprietes des scalaires : VISCLS si elle est variable
!     On utilisera IVISLS comme suit :
!       Pour le scalaire II
!         si IVISLS(II) = 0    : visc = VISLS0(II)
!         si IVISLS(II) .GT. 0 : visc = PROPCE(IEL ,IPPROC(IVISLS(II)))
!     Ceci permet de ne pas reserver des tableaux vides pour les
!       scalaires a viscosite constante

  ! Add a scalar diffusivity when defined as variable

  if (nscal.ge.1) then
    do ii = 1, nscal
      if (ivisls(ii).ne.0 .and. iscavr(ii).le.0) then
        ! Build name and label, using a general rule, with a
        ! fixed name for temperature or enthalpy
        id = ivarfl(isca(ii))
        call field_get_name(id, s_name)
        call field_get_label(id, s_label)
        if (ii.eq.iscalt) then
          s_name = 'thermal'
          s_label = 'Th'
        endif
        if (iscacp(ii).gt.0) then
          f_name  = trim(s_name) // '_conductivity'
          f_label = trim(s_label) // ' Cond'
        else
          f_name  = trim(s_name) // '_diffusivity'
          f_label = trim(s_label) // ' Diff'
        endif
        ! Special case for electric arcs: real and imaginary electric
        ! conductivity is the same (and ipotr < ipoti)
        if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
          if (ii.eq.ipotr) then
            f_name = 'elec_sigma'
            f_label = 'Sigma'
          else if (ii.eq.ipoti) then
            ivisls(ipoti) = ivisls(ipotr)
            cycle ! go to next scalar in loop, avoid creating property
          endif
        endif
        ! Now create matching property
        call add_property_field(f_name, f_label, ivisls(ii))
        call field_set_key_int(ivarfl(isca(ii)), kivisl, iprpfl(ivisls(ii)))
      endif
    enddo
  endif

  if (iscalt.gt.0) then
    if (ityturt(iscalt).gt.0.and.irovar.eq.1) then
      call add_property_field('thermal_expansion', 'Beta', ibeta)
    endif
  endif

!  Pour les fluctuations, le pointeur de la diffusivite
!    envoie directement sur la diffusivite du scalaire associe.

  do ii = 1, nscal
    if (ivisls(ii).gt.0) then
      if (iscavr(ii).gt.0) then
        ivisls(ii) = ivisls(iscavr(ii))
        if (ivisls(ii).gt.0) then
          call field_set_key_int(ivarfl(isca(ii)), kivisl,    &
                                 iprpfl(ivisls(iscavr(ii))))
        endif
      endif
    endif
  enddo

!     Numero max des proprietes ; ressert plus bas pour
!       ajouter celles relatives a la physique particuliere
  nprmax = nproce


! --- Positionnement dans les tableaux PROPCE

!   Au centre des cellules (tout sauf les flux de masse)
!     Pour les fluctuations, le pointeur de la diffusivite
!     envoie directement sur la diffusivite du scalaire associe.

!  On positionne en meme temps les pointeurs IPPPRO pour
!     le post traitement des proprietes physiques definies
!     aux cellules afin de ne pas en oublier.
!     Les pointeurs ont ete initialises a 1 (poubelle) dans iniini.
!     IPPPST commence a 2 car 1 est une poubelle.
!  Attention, IPPPST ressert plus bas.


!   - Repere pour la position des grandeurs posttraitees.

  ipppst = 1

!   - Pointeurs post pour les variables.

  do ivar = 1, nvar
    ipppst        = ipppst + 1
    ipprtp(ivar)  = ipppst
  enddo


!   - Numero des proprietes physiques definies au centre de cellules.
!     Et pointeur pour le post.
!      (attention, la viscosite turbulente n'est pas post traitable
!       pour les calculs laminaires, car elle n'existe pas)
!      (attention, en module electrique, on utilise la meme
!       conductivite electrique pour le potentiel reel et le potentiel
!       imaginaire


  iprop                 = irom
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = iviscl
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = ivisct
  if (iturb.eq.0) then
    ipppro(iprop)         = 1
  else
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
  endif
  iprop                 = icour
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  iprop                 = ifour
  ipppst                = ipppst + 1
  ipppro(iprop)         = ipppst

  if (ippmod(icompf).lt.0) then
    iprop                 = iprtot
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
  endif

  if (icp.gt.0) then
    iprop                 = icp
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
  endif

  if (ismago.gt.0) then
    iprop                 = ismago
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
  endif

  if (iale.eq.1) then
    iprop             = ivisma(1)
    ipppst            = ipppst + 1
    ipppro(iprop)     = ipppst
    if (iortvm.eq.1) then
      iprop             = ivisma(2)
      ipppst            = ipppst + 1
      ipppro(iprop)     = ipppst
      iprop             = ivisma(3)
      ipppst            = ipppst + 1
      ipppro(iprop)     = ipppst
    endif
  endif

  do iest = 1, nestmx
    if (iescal(iest).gt.0) then
      iprop                      = iestim(iest)
      ipppst                     = ipppst + 1
      ipppro(iprop)              = ipppst
    endif
  enddo

!     Conductivite electrique imaginaire :
!     La conductivite reelle et imaginaire sont dans le meme tableau.
!       En Joule, on ne reserve donc pas de propriete "viscosite"
!       pour le potentiel imaginaire.
!       Intervention 1/2

  if (ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
    do ii = 1, nscal
      if (ivisls(ii).gt.0.and.ii.ne.ipoti) then
        if (iscavr(ii).le.0) then
          iprop                 = ivisls(ii)
          ipppst                = ipppst + 1
          ipppro(iprop)         = ipppst
        endif
      endif
    enddo
  else
    do ii = 1, nscal
      if (ivisls(ii).gt.0) then
        if (iscavr(ii).le.0) then
          iprop                 = ivisls(ii)
          ipppst                = ipppst + 1
          ipppro(iprop)         = ipppst
        endif
      endif
    enddo
  endif
  if (iscalt.gt.0) then
    if (ityturt(iscalt).gt.0.and.irovar.eq.1) then!FIXME
      iprop                 = ibeta
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst
    endif
  endif

  nproce = iprop

! --- Modifications pour la physique particuliere
!      des entiers NPROCE

!      Sauvegarde pour la physique particuliere de IPROP
!      afin d'initialiser les positions des variables d'etat
!      Attention IPROPP est le dernier numero affecte pour les proprietes.
!                IPPPST est le rang de la derniere grandeur definie aux
!                       cellules pour le post traitement
  ipropp = nprmax
  call ppprop (ipropp,ipppst)
  !==========

!     On sauve IPROPP modifie pour les passages ulterieurs dans varpos
  nprmax = ipropp
!     On sauve IPPPST modifie pour les passages ulterieurs dans varpos
  nppmax = ipppst


! --- Verification de NPROCE

  if (nproce.gt.npromx) then
    write(nfecra,7200)nproce, npromx, nproce
    call csexit (1)
    !==========
  endif

  return

endif

!===============================================================================
! 2. DEUXIEME APPEL :
!        POSITIONNEMENT DES PROPRIETES POUR LE SCHEMA EN TEMPS,
!                                           LES MOMENTS ET FIN
!===============================================================================

!     Dans le cas ou on a un schema en temps d'ordre 2, il faut aussi
!       prevoir les proprietes au temps n-1. Ce sera fait au dernier appel

!     Dans le cas ou on calcule des moments, il faut en prevoir le nombre
!       et prevoir le nombre de tableaux necessaires pour stocker le
!       temps cumule de moyenne. On suppose que l'on peut s'aider
!       d'infos situees en tete de fichier suite (si on fait une
!       suite avec des moments non reinitialises).

if (ipass.eq.2) then

! ---> 2.1 PROPRIETES ADDITIONNELLES POUR LES ET SCHEMA EN TEMPS
!      ---------------------------------------------------------

! --- Initialisations par defaut eventuelles et verifications
!       des options utilisees ci-dessous pour decider si l'on
!       reserve des tableaux supplementaires pour des grandeurs
!       au pas de temps precedent

  iok = 0

!     Pression hydrostatique
  if (iphydr.eq.0.or.iphydr.eq.2) then
    icalhy = 0
  else if (iphydr.eq.1) then
    gravn2 = gx**2+gy**2+gz**2
    if (gravn2.lt.epzero**2) then
      icalhy = 0
    else
      icalhy = 1
    endif
  endif

!     Schemas en temps
!         en LES : Ordre 2 ; sinon Ordre 1
!         (en particulier, ordre 2 impossible en k-eps couple)
  if (ischtp.eq.-999) then
    if (itytur.eq.4) then
      ischtp = 2
    else
      ischtp = 1
    endif
  endif

!     Schemas en temps : variables deduites
!     Schema pour le Flux de masse
  if (istmpf.eq.-999) then
    if (ischtp.eq.1) then
      istmpf = 1
    else if (ischtp.eq.2) then
      istmpf = 2
    endif
  endif
  !     Masse volumique
  if (iroext.eq.-999) then
    if (ischtp.eq.1) then
      iroext = 0
    else if (ischtp.eq.2) then
      !       Pour le moment par defaut on ne prend pas l'ordre 2
      !              IROEXT = 1
      iroext = 0
    endif
  endif
  !     Viscosite
  if (iviext.eq.-999) then
    if (ischtp.eq.1) then
      iviext = 0
    else if (ischtp.eq.2) then
      !       Pour le moment par defaut on ne prend pas l'ordre 2
      !              IVIEXT = 1
      iviext = 0
    endif
  endif
  !     Chaleur massique
  if (icpext.eq.-999) then
    if (ischtp.eq.1) then
      icpext = 0
    else if (ischtp.eq.2) then
      !       Pour le moment par defaut on ne prend pas l'ordre 2
      !              ICPEXT = 1
      icpext = 0
    endif
  endif
  !     Termes sources NS,
  if (isno2t.eq.-999) then
    if (ischtp.eq.1) then
      isno2t = 0
      !            ELSE IF (ISCHTP.EQ.2.AND.IVISSE.EQ.1) THEN
    else if (ischtp.eq.2) then
      !       Pour le moment par defaut on prend l'ordre 2
      isno2t = 1
      !              ISNO2T = 0
    endif
  endif
  !     Termes sources turbulence (k-eps, Rij, v2f ou k-omega)
  !     On n'autorise de changer ISTO2T qu'en Rij (sinon avec
  !       le couplage k-eps/omega il y a pb)
  if (isto2t.eq.-999) then
    if (ischtp.eq.1) then
      isto2t = 0
    else if (ischtp.eq.2) then
      !       Pour le moment par defaut on ne prend pas l'ordre 2
      !              ISTO2T = 1
      isto2t = 0
    endif
  else if ( itytur.eq.2.or.iturb.eq.50             &
       .or.iturb.ne.60) then
    write(nfecra,8132) iturb,isto2t
    iok = iok + 1
  endif

  idttur = 0

  do iscal = 1, nscal
!     Termes sources Scalaires,
    if (isso2t(iscal).eq.-999) then
      if (ischtp.eq.1) then
        isso2t(iscal) = 0
      else if (ischtp.eq.2) then
!       Pour coherence avec Navier Stokes on prend l'ordre 2
!       mais de toute facon qui dit ordre 2 dit LES et donc
!       generalement pas de TS scalaire a interpoler.
        isso2t(iscal) = 1
!              ISSO2T(ISCAL) = 0
      endif
    endif
!     Diffusivite scalaires
    if (ivsext(iscal).eq.-999) then
      if (ischtp.eq.1) then
        ivsext(iscal) = 0
      else if (ischtp.eq.2) then
!       Pour le moment par defaut on ne prend pas l'ordre 2
!              IVSEXT(ISCAL) = 1
        ivsext(iscal) = 0
      endif
    endif

    ! ---> Model for turbulent fluxes u'T' (SGDH, GGDH, AFM, DFM)
    ityturt(iscal) = iturt(iscal)/10

    ! Index of the turbulent flux
    if (ityturt(iscal).eq.3) then
      idttur = idttur + 1
      ifltur(iscal) = idttur
    endif

  enddo

!     Pression hydrostatique
  if (iphydr.ne.0.and.iphydr.ne.1.and.iphydr.ne.2) then
    write(nfecra,8121) 'IPHYDR ',iphydr
    iok = iok + 1
  endif

!     Viscosite secondaire
  ivisph = ivisse
  if (ivisph.ne.0.and.ivisph.ne.1) then
    write(nfecra,8022) 'IVISSE ',ivisph
    iok = iok + 1
  endif

!     Schemas en temps

  !     Schema en temps global.
  if (ischtp.ne. 1.and.ischtp.ne.2) then
    write(nfecra,8101) 'ISCHTP',ischtp
    iok = iok + 1
  endif
  if (ischtp.eq. 2.and.idtvar.ne.0) then
    write(nfecra,8111) ischtp,idtvar
    iok = iok + 1
  endif
  if (ischtp.eq. 2.and.itytur.eq.2) then
    write(nfecra,8112) ischtp,iturb
    iok = iok + 1
  endif
  if (ischtp.eq.1.and.itytur.eq.4) then
    write(nfecra,8113) ischtp,iturb
  endif
  if (ischtp.eq. 2.and.iturb.eq.50) then
    write(nfecra,8114) ischtp,iturb
    iok = iok + 1
  endif
  if (ischtp.eq. 2.and.iturb.eq.51) then
    write(nfecra,8117) ischtp,iturb
    iok = iok + 1
  endif
  if (ischtp.eq. 2.and.iturb.eq.60) then
    write(nfecra,8115) ischtp,iturb
    iok = iok + 1
  endif
  if (ischtp.eq. 2.and.iturb.eq.70) then
    write(nfecra,8116) ischtp,iturb
    iok = iok + 1
  endif

  !     Schema en temps pour le flux de masse
  if (istmpf.ne. 2.and.istmpf.ne.0.and.            &
       istmpf.ne. 1) then
    write(nfecra,8121) 'ISTMPF',istmpf
    iok = iok + 1
  endif

  !     Schema en temps pour les termes sources de NS
  if (isno2t.ne.0.and.                                    &
       isno2t.ne. 1.and.isno2t.ne.2) then
    write(nfecra,8131) 'ISNO2T',isno2t
    iok = iok + 1
  endif
  !     Schema en temps pour les termes sources des grandeurs
  !     turbulentes
  if (isto2t.ne.0.and.                                    &
       isto2t.ne. 1.and.isto2t.ne.2) then
    write(nfecra,8131) 'ISTO2T',isto2t
    iok = iok + 1
  endif

  !     Schema en temps pour la masse volumique
  if (iroext.ne.0.and.                                    &
       iroext.ne. 1.and.iroext.ne.2) then
    write(nfecra,8131) 'IROEXT',iroext
    iok = iok + 1
  endif

  !     Schema en temps pour la viscosite
  if (iviext.ne.0.and.                                    &
       iviext.ne. 1.and.iviext.ne.2) then
    write(nfecra,8131) 'IVIEXT',iviext
    iok = iok + 1
  endif

  !     Schema en temps pour la chaleur specifique
  if (icpext.ne.0.and.                                    &
       icpext.ne. 1.and.icpext.ne.2) then
    write(nfecra,8131) 'ICPEXT',icpext
    iok = iok + 1
  endif

  do iscal = 1, nscal
!     Schema en temps pour les termes sources des scalaires
    if (isso2t(iscal).ne.0.and.                                    &
       isso2t(iscal).ne. 1.and.isso2t(iscal).ne.2) then
      write(nfecra,8141) iscal,'ISSO2T',isso2t(iscal)
      iok = iok + 1
    endif
!     Schema en temps pour la viscosite
    if (ivsext(iscal).ne.0.and.                                    &
       ivsext(iscal).ne. 1.and.ivsext(iscal).ne.2) then
      write(nfecra,8141) iscal,'IVSEXT',ivsext(iscal)
      iok = iok + 1
    endif
  enddo

!     Stop si probleme
  if (iok.gt.0) then
    call csexit(1)
  endif


! --- Reprise du dernier numero de propriete
  iprop  = nprmax

! --- Numeros de propriete

  !    Source term for weakly compressible algorithm (semi analytic scheme)
  if (idilat.eq.4) then
    do iscal = 1, nscal
      iprop         = iprop + 1
      iustdy(iscal) = iprop
    enddo
    itsrho = nscal + 1
    iprop          = iprop + 1
    iustdy(itsrho) = iprop
  endif

  ! The density at the previous time step is required if idilat>1 or if
  ! we perform a hydrostatic pressure correction (icalhy=1)
  if (iroext.gt.0.or.icalhy.eq.1.or.idilat.gt.1.or.ippmod(icompf).ge.0) then
    iprop         = iprop + 1
    iroma  = iprop
    call field_set_n_previous(iprpfl(irom), 1)
  endif
  !     Dans le cas d'une extrapolation de la viscosite totale
  if (iviext.gt.0) then
    iprop         = iprop + 1
    ivisla = iprop
    iprop         = iprop + 1
    ivista = iprop
    call field_set_n_previous(iprpfl(ivisct), 1)
  endif

  !     Proprietes des phases : CP s'il est variable
  if (icp.ne.0) then
    if (icpext.gt.0) then
      iprop         = iprop + 1
      icpa   = iprop
      call field_set_n_previous(iprpfl(icp), 1)
    endif
  endif

  !     On a besoin d'un tableau pour les termes sources de Navier Stokes
  !       a extrapoler. Ce tableau est NDIM
  if (isno2t.gt.0) then
    iprop         = iprop + 1
    itsnsa = iprop
  endif
  if (isto2t.gt.0) then
    iprop         = iprop + 1
    itstua = iprop
  endif

!     Proprietes des scalaires : termes sources pour theta schema
!       et VISCLS si elle est variable
  if (nscal.ge.1) then
    do iscal = 1, nscal
      if (isso2t(iscal).gt.0) then
        iprop         = iprop + 1
        itssca(iscal) = iprop
      endif
      if (ivisls(iscal).ne.0) then
        if (ivsext(iscal).gt.0) then
          iprop         = iprop + 1
          ivissa(iscal) = iprop
        endif
      endif
    enddo
  endif

! --- Sauvegarde du dernier numero de propriete
  nprmax = iprop


! --- Reprise du dernier NPROCE et du dernier NPPMAX
  iprop                 = nproce
  ipppst                = nppmax

! --- Positionnement des PROPCE

  !    Source term for weakly compressible algorithm (semi analytic scheme)
  if (idilat.eq.4) then
    do iscal = 1, nscal
      iprop                 = iprop + 1
      ipproc(iustdy(iscal)) = iprop
    enddo
    iprop                  = iprop + 1
    ipproc(iustdy(itsrho)) = iprop
  endif

  ! Variables schema en temps
  if (iroext.gt.0.or.icalhy.eq.1.or.idilat.gt.1.or.ippmod(icompf).ge.0) then
    iprop                 = iprop  + 1
    ipproc(iroma ) = iprop
  endif
  if (iviext.gt.0) then
    iprop                 = iprop  + 1
    ipproc(ivisla) = iprop
  endif
  if (iviext.gt.0) then
    iprop                 = iprop  + 1
    ipproc(ivista) = iprop
  endif
  if (icpext.gt.0) then
    iprop                 = iprop + 1
    ipproc(icpa  ) = iprop
  endif
  if (isno2t.gt.0) then
    iprop                 = iprop + 1
    ipproc(itsnsa) = iprop
    !     Ce tableau est NDIM :
    iprop                 = iprop + ndim-1
  endif
  if (isto2t.gt.0) then
    iprop                 = iprop + 1
    ipproc(itstua) = iprop
    !     Ce tableau est 2, 7 ou 4 selon le modele de turbulence :
    if    (itytur.eq.2) then
      iprop                 = iprop + 2-1
    else if (itytur.eq.3) then
      iprop                 = iprop + 7-1
      if (iturb.eq.32) then
        iprop                 = iprop + 1
      endif
    else if (iturb.eq.50) then
      iprop                 = iprop + 4-1
    else if (iturb.eq.70) then
      iprop                 = iprop + 1-1
    endif
  endif

  do ii = 1, nscal
! Termes source des scalaires pour theta schema
    if (isso2t(ii).gt.0) then
      iprop                 = iprop + 1
      ipproc(itssca(ii))    = iprop
    endif

    if (ivisls(ii).gt.0) then
      if (iscavr(ii).le.0) then
        if (ivsext(ii).gt.0) then
          iprop                 = iprop + 1
          ipproc(ivissa(ii)  )  = iprop
        endif
      endif
    endif
  enddo
  do ii = 1, nscal
    if (ivisls(ii).gt.0) then
      if (iscavr(ii).gt.0.and.iscavr(ii).le.nscal) then
        if (ivsext(ii).gt.0) then
          ipproc(ivissa(ii)  )  = ipproc(ivissa(iscavr(ii))  )
        endif
      endif
    endif
  enddo

! --- Sauvegarde du dernier NPROCE et du dernier NPPMAX
  nproce = iprop
  nppmax = ipppst


! ---> 2.2 CALCUL DE LA TAILLE DU TABLEAU DES TEMPS CUMULES POUR LES MOMENTS
!      ---------------------------------------------------------------------

!     Pour verification des definitions de moments
  iok = 0

! --- Calcul du nombre de moments definis (et verif qu'il n y a pas de trous)
  nbmomt = 0
  inmfin = 0
  do imom = 1, nbmomx
!     Si on n'est pas a la fin de la liste
    if (inmfin.eq.0) then
!       Si il y en a, ca en fait en plus
      if (idfmom(1,imom).ne.0) then
        nbmomt = nbmomt + 1
!       Si il n'y en a pas, c'est la fin de la liste
      else
        inmfin = 1
      endif
!     Si on est a la fin de la liste, il n'en faut plus
    else
      if (idfmom(1,imom).ne.0) then
        iok = iok + 1
      endif
    endif
  enddo

  if (iok.ne.0) then
    write(nfecra,8200)nbmomt+1,idfmom(1,nbmomt+1),nbmomt,nbmomt
    do imom = 1, nbmomx
      write(nfecra,8201)imom,idfmom(1,imom)
    enddo
    write(nfecra,8202)
  endif

! --- Verification de IDFMOM
  iok = 0
  do imom = 1, nbmomx
    idffin = 0
    do jj = 1, ndgmox
      idfmji = idfmom(jj,imom)
      if (idffin.eq.0) then
        if (idfmji.lt.-nprmax) then
          iok = iok + 1
          write(nfecra,8210)jj,imom,idfmji,nprmax
        else if (idfmji.gt.nvar) then
          iok = iok + 1
          write(nfecra,8211)jj,imom,idfmji,nvar
        else if (idfmji.lt.0) then
          if (ipproc(-idfmji).le.0) then
            iok = iok + 1
            write(nfecra,8212)jj,imom,idfmji,-idfmji,             &
                 ipproc(-idfmji)
          endif
        else if (idfmji.eq.0) then
          idffin = 1
        endif
      else
        if (idfmji.ne.0) then
          iok = iok + 1
          write(nfecra,8213)imom,jj,idfmji
        endif
      endif
    enddo
  enddo

! --- Verification de NTDMOM (>0)
  do imom = 1, nbmomt
    if (ntdmom(imom).lt.0) then
      iok = iok + 1
      write(nfecra,8214)imom,ntdmom(imom)
    endif
  enddo

  if (iok.ne.0) then
    call csexit(1)
  endif


! --- Indicateur pour juger de l'utilite de consulter le fichier suite
!     Il n'est pas indispensable de consulter avec succes le fichier suite
!       lorsqu'il n'y a aucun moment a relire ou qu'on ne fait pas de suite
!       (si ILSMOM=0)
!       ie ILSMOM = 0 == (ISUITE=0 OR IMOOLD(IMOM)=-1 pour tout IMOM)
!       ie ILSMOM = 1 == (ISUITE=1 AND IMOOLD(IMOM) different de -1 pour un IMOM)
  ilsmom = 0
  if (isuite.eq.1.and.nbmomt.gt.0) then
    do imom = 1, nbmomt
      if (imoold(imom).ne.-1) then
        ilsmom = 1
      endif
    enddo
  endif

! --- Lecture du fichier suite (debut : info sur les moments et sur le ntpabs)
  if (ilsmom.eq.1) then

!     Ouverture
!        (ILECEC=1:lecture)
    ilecec = 1
    ficsui = 'auxiliary'
    call opnsui(ficsui,len(ficsui),ilecec,impamx,ierror)
    !==========
    if (ierror.ne.0) then
      write(nfecra,8300) ficsui
      call csexit (1)
    endif

!     Tests : type de fichier et dimension des supports (puisqu'on y
!       relira des infos, autant s'arreter au plus vite
!       s'il n'est pas correct)

!     Type (fichier auxiliaire)

    itysup = 0
    nbval  = 1
    irtyp  = 1
    rubriq = 'version_fichier_suite_auxiliaire'
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                ivers,ierror)

    if (ierror.ne.0) then
      write(nfecra,8301)ficsui
      call csexit (1)
    endif

    nberro=0

!     Nb de moments
!       Si le nombre de moments n'a pas pu etre relu, on devra s'arreter
    itysup = 0
    nbval  = 1
    irtyp  = 1
    rubriq = 'nombre_moyennes_temps'
    jbmomt = 0
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                jbmomt,ierror)
    nberro=nberro+ierror
    if (ierror.ne.0) then
      write(nfecra,8311)
    endif

!     Nombre de pas de temps final du calcul precedent
    rubriq = 'nbre_pas_de_temps'
    itysup = 0
    nbval  = 1
    irtyp  = 1
    jtcabs = 0
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                jtcabs,ierror)
    nberro=nberro+ierror
    if (ierror.ne.0) then
      write(nfecra,8312)
    endif

    if (nberro.ne.0) then
      call csexit (1)
    endif

  endif


! --- Avec JTCABS et JBMOMT
!     on complete et on verifie la correspondance nouveaux -> anciens moments

!       Si suite et il existe un IMOOLD different de -1 :
!                      si IMOOLD mal renseigne erreur
!                         (valeur non admissible ou demande de relecture
!                             avec NTDBMO pas coherent)
!                      si IMOOLD pas renseigne pour NTDMOM > JTCABS : initialise
!                                              pour IMOM < JBMOMT+1 : on relit
!                                              pour IMOM > JBMOMT : initialise
!                      si on pointe deux fois sur le meme : erreur
!                                  (simplifie la lecture lecamo)
  iok = 0

  if (ilsmom.eq.1) then
    do imom = 1, nbmomt
      if (imoold(imom).gt.jbmomt.or.imoold(imom).eq.0              &
                             .or.imoold(imom).lt.-2) then
        write(nfecra,8400) imom,imoold(imom),jbmomt
        iok = iok + 1
      else if (imoold(imom).le.jbmomt.and.imoold(imom).gt.0        &
                             .and.ntdmom(imom).gt.jtcabs) then
        write(nfecra,8401) imom,imoold(imom),ntdmom(imom),jtcabs
        iok = iok + 1
      else if (imoold(imom).eq.-2.and.ntdmom(imom).gt.jtcabs) then
        imoold(imom)=-1
      else if (imoold(imom).eq.-2.and.imom.le.jbmomt) then
        imoold(imom)=imom
      else if (imoold(imom).eq.-2.and.imom.gt.jbmomt) then
        imoold(imom)=-1
      endif
    enddo
    do imom = 1, nbmomt
      if (imoold(imom).gt.0) then
        do jmom = 1, imom-1
          if (imoold(imom).eq.imoold(jmom)) then
            write(nfecra,8402) imom,jmom,imoold(imom)
            iok = iok + 1
          endif
        enddo
      endif
    enddo

!       Si pas suite : si IMOOLD pas renseigne, on initialise
!                      sinon erreur
  else if (isuite.eq.0) then
    do imom = 1, nbmomt
      if (imoold(imom).ne.-2) then
        write(nfecra,8403) isuite,imom,imoold(imom)
        iok = iok + 1
      else
        imoold(imom)=-1
      endif
    enddo
  endif

  if (iok.ne.0) then
    call csexit (1)
  endif

! --- D'autres informations sont lues dans le fichier suite pour
!     l'initialisation du numero dans le fichier suite du cumul temporel
!       associe a chaque moyenne (IDTOLD)
!     IDTOLD(IMOM) donne pour la moyenne (du calcul courant) IMOM
!       le numero dans le fichier suite du cumul temporel associe
!       a la moyenne du calcul precedent a laquelle correspond IMOM (ouf          !)

!     Initialisation
  do imom = 1, nbmomx
    idtold(imom) = 0
  enddo

!     Si on ne lit rien dans le fichier suite, IDTOLD restera nul
  if (ilsmom.eq.1) then

    nberro = 0

!     On ne lit des choses que s'il y a des moments stockes
    if ( jbmomt.gt.0 ) then

!        On lit les infos pour tous les moments courants
      do imom = 1, nbmomt
!           On s'interesse uniquement aux moments qui suivnt d'anciens moments
        if (imoold(imom).gt.0) then
!            Si le numero de l'ancien moment est incompatible avec le format
!              on cherche numero_cumul_temps_momentYYYY qui n existe pas
!              ca genere une erreur
          if (imom.le.nfmtmo) then
            write(cmoy4,'(i4.4)')imoold(imom)
          else
            cmoy4 = cindfm
          endif
          itysup = 0
          nbval  = 1
          irtyp  = 1
          rubriq = 'numero_cumul_temps_moment'//cmoy4
          call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,     &
                      irtyp,idtold(imom),ierror)
          nberro=nberro+ierror
          if (idtold(imom).eq.0.or.ierror.ne.0) then
            write(nfecra,8313)imom,imoold(imom)
            iok = iok + 1
          endif
        endif
      enddo
    endif

    if (nberro.ne.0) then
      call csexit (1)
    endif

!     Fermeture du fichier suite auxiliaire
    call clssui(impamx,ierror)

    if (ierror.ne.0) then
      write(nfecra,8390) ficsui
    endif

  endif


! --- Remplissage de IDTMOM (pointeur dans DTCMOM pour chaque moment)
!       > 0 : pointe sur un tableau NCEL dans PROPCE (DT cumule non uniforme)
!             remplissage sans trous de 1 a n
!       < 0 : pointe sur une case d'un tableau NPROMX  (DT cumule uniforme)
!             remplissage sans trous de -1 a -p

  iiplus = 0
  iimoin = 0

  do imom = 1, nbmomx
    idtmom(imom) = 0
  enddo

!       Pour les moments reinitialises (suite) ou initialises (non suite)
!         un nouveau tableau, sauf si leur calcul commence au meme instant
  do imom = 1, nbmomt
    imold = imoold(imom)
!        Si on (re)initialise IMOM
    if (imold.eq.-1) then
!          On cherche si on en a deja vu (re)initialisees au meme moment
      imomr = 0
      do jmom = 1, imom-1
        jmold = imoold(jmom)
        if (jmold.eq.-1.and.ntdmom(jmom).eq.ntdmom(imom)) then
          imomr = jmom
        endif
      enddo
!          Si oui : on utilise le meme tableau
      if (imomr.gt.0) then
        idtmom(imom) = idtmom(imomr)
!          Si non : on a besoin d'un nouveau tableau ou reel
      else
        if (idtvar.eq.2.or.idtvar.eq.-1) then
          iiplus = iiplus + 1
          idtmom(imom) = iiplus
        else
          iimoin = iimoin - 1
          idtmom(imom) = iimoin
        endif
      endif
    endif
  enddo

!       Pour les moments IMOM relus dans IMOLD
!         (ie suite + IMOOLD non egal a -1)
  if (isuite.eq.1) then
    do imom = 1, nbmomt
      imold = imoold(imom)
      if (imold.gt.0) then
!         On regarde si le DTcumule du IMOLD a deja ete vu
!            (indicateur JMOMOK)
        idto = idtold(imom)
        jmomok = 0
        do jmom = 1, imom-1
          jmold = imoold(jmom)
          if (jmold.gt.0) then
            jdto = idtold(jmom)
            if (jdto.eq.idto) then
              jmomok = jmom
            endif
          endif
        enddo
!         Si on le voit pour la premiere fois, ca en fait un de plus
        if (jmomok.eq.0) then
!           Si DT non uniforme dans le present calcul ou que
!              DT cumule etait non uniforme dans l'ancien, on a un
!              DT cumule non uniforme (sinon, il est uniforme)
          if (idtvar.eq.2.or.idtvar.eq.-1.or.idto.gt.0) then
            iiplus = iiplus + 1
            idtmom(imom) = iiplus
          else
            iimoin = iimoin - 1
            idtmom(imom) = iimoin
          endif
!         Si on l'a deja rencontre, on pointe au meme endroit
        else
          idtmom(imom) = idtmom(jmomok)
        endif
      endif
    enddo
  endif

!       Verification de IDTMOM : normalement, jamais ca plante ici.
  iok = 0
  do imom = 1, nbmomt
    if (idtmom(imom).eq.0) then
      iok = iok + 1
      write(nfecra,8410)imom, idtmom(imom)
    endif
  enddo
  if (iok.ne.0) then
    call csexit (1)
  endif

! --- Calcul du nombre de tableaux NCEL "temps cumule"
  nbdtcm = 0
  do imom = 1, nbmomt
    nbdtcm=max(idtmom(imom),nbdtcm)
  enddo

! ---> 2.3 POSITIONNEMENT DANS PROPCE DES MOMENTS ET DU TEMPS CUMULE
!      -------------------------------------------------------------

! --- Reprise du dernier numero de propriete
  iprop  = nprmax

! --- Numeros de propriete
  do imom = 1, nbmomt
    iprop  = iprop + 1
    icmome(imom) = iprop
  enddo
  do ii = 1, nbdtcm
    iprop  = iprop + 1
    icdtmo(ii) = iprop
  enddo

! --- Sauvegarde du dernier numero de propriete
  nprmax = iprop

! --- Reprise des derniers NPROCE et NPPMAX (PROPCE et POST-TRAITEMENT)
  iprop                 = nproce
  ipppst                = nppmax

! --- Positionnement
  do imom = 1, nbmomt
    iprop                = iprop + 1
    ipproc(icmome(imom)) = iprop
    ipppst               = ipppst + 1
    ipppro(iprop)        = ipppst
  enddo
  do ii = 1, nbdtcm
    iprop                = iprop + 1
    ipproc(icdtmo(ii))   = iprop
    ipppst               = ipppst + 1
    ipppro(iprop)        = ipppst
  enddo

! --- Sauvegarde du dernier NPROCE et NPPMAX
  nproce = iprop
  nppmax = ipppst

!        En compressible, avec l'algorithme en pression, celle-ci est
!        une grandeur instationnaire => istat = 1
  if (ippmod(icompf).ge.0) then
    istat(ipr) = 1
  endif

! ---> 2.4 POINTEURS POST-PROCESSING / LISTING / HISTORIQUES / CHRONOS
!      ---------------------------------------------------------------------

! --- Les pointeurs ont ete initialises a 1 (poubelle) dans iniini.

!     On posttraitera les variables localisees au centre des cellules.

!      IPPRTP(IVAR)  pour RTP a ete complete plus haut.

!      IPPPRO(IPPROC(II)) pour PROPCE a ete complete plus haut
!        au fur et a mesure (voir ppprop en particulier)

!      Le rang de la derniere propriete pour le post est IPPPST.


  if (idtvar.gt.0) then
    ipppst = ipppst + 1
    ippdt  = ipppst
    nppmax = ipppst
  endif

  if (ipucou.eq.1) then
    ipppst = ipppst + 1
    ipptx  = ipppst
    ipppst = ipppst + 1
    ippty  = ipppst
    ipppst = ipppst + 1
    ipptz  = ipppst
    nppmax = ipppst
  endif


! Verification de la limite sur IPPPST

  if (ipppst.gt.nvppmx) then
    write(nfecra,8900)ipppst,nvppmx
    call csexit (1)
    !==========
  endif

  return

endif

!===============================================================================
! 3. TROISIEME APPEL :
!        RESERVATION D'UNE PLACE DANS PROPCE SI RAYONNEMENT
!        ET LAGRANGIEN AVEC THERMIQUE DES PARTICULES
!===============================================================================

if (ipass.eq.3) then

  if ( iirayo.gt.0 ) then

! --- Reprise du dernier numero de propriete
    iprop  = nprmax

! --- Numeros de propriete
    iprop        = iprop + 1
    ilumin       = iprop
    iprop        = iprop + 1
    iqx          = iprop
    iprop        = iprop + 1
    iqy          = iprop
    iprop        = iprop + 1
    iqz          = iprop


    do irphas = 1, nrphas

      iprop                 = iprop + 1
      itsre(irphas)         = iprop
      iprop                 = iprop + 1
      itsri(irphas)         = iprop
      iprop                 = iprop + 1
      iabs(irphas)          = iprop
      iprop                 = iprop + 1
      iemi(irphas)          = iprop
      iprop                 = iprop + 1
      icak(irphas)          = iprop

    enddo

    nprayc = iprop - nprmax

! --- Sauvegarde du dernier numero de propriete
    nprmax = iprop


! --- Reprise des derniers NPROCE et NPPMAX (PROPCE et POST-TRAITEMENT)
    iprop         = nproce
    ipppst        = nppmax

! --- Positionnement
    iprop          = iprop + 1
    ipproc(ilumin) = iprop
    ipppst         = ipppst + 1
    ipppro(iprop)  = ipppst

    iprop         = iprop + 1
    ipproc(iqx)   = iprop
    ipppst        = ipppst + 1
    ipppro(iprop) = ipppst

    iprop         = iprop + 1
    ipproc(iqy)   = iprop
    ipppst        = ipppst + 1
    ipppro(iprop) = ipppst

    iprop         = iprop + 1
    ipproc(iqz)   = iprop
    ipppst        = ipppst + 1
    ipppro(iprop) = ipppst

! Positionnement de ITSRE, ITSRI, ICAK, IABS et IEME
! Leur dimensionnement n'est pas le meme si on est en charbon ou non


    do irphas = 1, nrphas
!
      iprop                = iprop + 1
      ipproc(itsre(irphas)) = iprop
      ipppst               = ipppst + 1
      ipppro(iprop)        = ipppst

      iprop                = iprop + 1
      ipproc(itsri(irphas)) = iprop
      ipppst               = ipppst + 1
      ipppro(iprop)        = ipppst
!
      iprop                = iprop + 1
      ipproc(iabs(irphas))  = iprop
      ipppst               = ipppst + 1
      ipppro(iprop)        = ipppst

      iprop                = iprop + 1
      ipproc(iemi(irphas))  = iprop
      ipppst               = ipppst + 1
      ipppro(iprop)        = ipppst

      iprop                = iprop + 1
      ipproc(icak(irphas)) = iprop
      ipppst               = ipppst + 1
      ipppro(iprop)        = ipppst

    enddo


! --- Sauvegarde du dernier NPROCE et NPPMAX
    nproce = iprop
    nppmax = ipppst


! --- Reprise du dernier numero de propriete
    iprop = nprmax


! --- Numeros de propriete
    iprop        = iprop + 1
    ixlam        = iprop
    iprop        = iprop + 1
    iepa         = iprop
    iprop        = iprop + 1
    ieps         = iprop
    iprop        = iprop + 1
    ifnet        = iprop

! --- Sauvegarde du dernier numero de propriete
    nprmax = iprop

    if (iihmpr.eq.1) then

      call uirapr &
      !==========
    ( nprayc, nrphas, ipppro, ipproc,                   &
      ilumin, iqx, iqy, iqz,                            &
      itsre, itsri, iabs, iemi, icak)

    endif

!
! --- Verification de NPROCE

    if (nproce.gt.npromx) then
      write(nfecra,7200) nproce, npromx, nproce
      call csexit (1)
      !==========
    endif

  endif

  return

endif


!===============================================================================
! 4. Formats
!===============================================================================

#if defined(_CS_LANG_FR)

 7030 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ASSOCIE A UNE VARIANCE INCORRECT               ',/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire utilisateur           ',I10   ,') represente  ',/,&
'@     une variance puisque                                   ',/,&
'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,&
'@  Les valeurs de ISCAVR doivent cependant etre              ',/,&
'@    superieures ou egales a                  0              ',/,&
'@    inferieures ou egales a NSCAL = ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier ISCAVR.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7031 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ASSOCIE A UNE VARIANCE INCORRECT               ',/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire physique particuliere ',I10   ,') represente  ',/,&
'@     une variance puisque    ',/,                         &
'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,&
'@  Les valeurs de ISCAVR doivent cependant etre              ',/,&
'@    superieures ou egales a                  0              ',/,&
'@    inferieures ou egales a NSCAL = ',I10                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier ISCAVR.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7040 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10                                           ,/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire utilisateur           ',I10   ,') represente  ',/,&
'@    la variance des fluctuations du scalaire ',I10           ,/,&
'@    (scalaire utilisateur           ',I10   ,') puisque     ',/,&
'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,&
'@                                                            ',/,&
'@  L''indicateur de diffusivite IVISLS(',I10  ,')            ',/,&
'@    ne doit pas etre renseigne.                             ',/,&
'@  Il sera pris automatiquement egal a celui du scalaire     ',/,&
'@    associe, soit ',I10                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7041 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10                                           ,/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire utilisateur           ',I10   ,') represente  ',/,&
'@    la variance des fluctuations du scalaire ',I10           ,/,&
'@    (scalaire physique particuliere ',I10   ,') puisque     ',/,&
'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,&
'@                                                            ',/,&
'@  L''indicateur de diffusivite IVISLS(',I10  ,')            ',/,&
'@    ne doit pas etre renseigne.                             ',/,&
'@  Il sera pris automatiquement egal a celui du scalaire     ',/,&
'@    associe, soit ',I10                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7042 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10                                           ,/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire physique particuliere ',I10   ,') represente  ',/,&
'@    la variance des fluctuations du scalaire ',I10           ,/,&
'@    (scalaire utilisateur           ',I10   ,') puisque     ',/,&
'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,&
'@                                                            ',/,&
'@  L''indicateur de diffusivite IVISLS(ISCAPP(',I10  ,'))    ',/,&
'@    ne doit pas etre renseigne.                             ',/,&
'@  Il sera pris automatiquement egal a celui du scalaire     ',/,&
'@    associe, soit ',I10                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7043 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10                                           ,/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire physique particuliere ',I10   ,') represente  ',/,&
'@    la variance des fluctuations du scalaire ',I10           ,/,&
'@    (scalaire physique particuliere ',I10   ,') puisque     ',/,&
'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,&
'@                                                            ',/,&
'@  L''indicateur de diffusivite IVISLS(ISCAPP(',I10  ,'))    ',/,&
'@    ne doit pas etre renseigne.                             ',/,&
'@  Il sera pris automatiquement egal a celui du scalaire     ',/,&
'@    associe, soit ',I10                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7050 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10                                           ,/,&
'@                                                            ',/,&
'@  L''indicateur de diffusivite variable                     ',/,&
'@    du scalaire utilisateur           ',I10                  ,/,&
'@    IVISLS(',I10   ,')                                      ',/,&
'@    doit etre un entier egal a 0 ou 1.                      ',/,&
'@    Il vaut ici ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7051 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10                                           ,/,&
'@                                                            ',/,&
'@  L''indicateur de diffusivite variable                     ',/,&
'@    du scalaire physique particuliere ',I10                  ,/,&
'@    IVISLS(ISCAPP(',I10   ,'))                              ',/,&
'@    doit etre un entier egal a 0 ou 1.                      ',/,&
'@    Il vaut ici ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7070 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    L''INDICATEUR IORTVM NE PEUT PRENDRE QUE LES VALEURS    ',/,&
'@      0 OU 1.                                               ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usalin.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE PROPRIETES TROP GRAND                        ',/,&
'@                                                            ',/,&
'@  Le type de calcul defini                                  ',/,&
'@    correspond aux nombres de proprietes suivants           ',/,&
'@      au centre des cellules       : NPROCE = ',I10          ,/,&
'@  Le nombre de proprietes maximal prevu                     ',/,&
'@                      dans paramx.h est NPROMX = ',I10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@  NPROMX doit valoir au moins ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8022 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8101 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 1 ou 2                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8111 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    AVEC UN SCHEMA EN TEMPS D ORDRE 2 : ISCHTP = ', I10      ,/,&
'@    IL FAUT UTILISER UN PAS DE TEMPS CONSTANT ET UNIFORME   ',/,&
'@    OR IDTVAR = ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8112 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN K-EPSILON (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l ordre 2 avec le    ',/,&
'@   couplage des termes sources du k-epsilon.                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8113 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 1 (ISCHTP = ',I10   ,/,&
'@    EN LES (ITURB = ',I10,' )'                               ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres.              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8114 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN PHI_FBAR (ITURB = ',I10,' )'                          ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources du k-epsilon.                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8117 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN BL-V2/K  (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources du k-epsilon.                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8115 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN K-OMEGA   (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources du k-omega.                  ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8116 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN SPALART   (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources de Spalart-Allmaras.         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8121 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A  0, 1 OU 2            ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8131 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2             ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8132 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Avec le modele de turbulence choisi, ITURB = ',I10         ,/,&
'@    la valeur de ISTO2T (extrapolation des termes sources   ',/,&
'@    pour les variables turbulentes) ne doit pas etre modifie',/,&
'@    or ISTO2T a ete force a ',I10                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8141 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES SCALAIRE ',I10 ,/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2             ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@      SUR LA LISTE DES MOYENNES TEMPORELLES                 ',/,&
'@                                                            ',/,&
'@    La valeur de IDFMOM(1,',I10   ,' est ',I10               ,/,&
'@      ceci indique que l''on a defini ',I10   ,' moyennes   ',/,&
'@      temporelles en renseignant les IMOM = ',I10            ,/,&
'@      premieres cases du tableau IDFMOM(.,IMOM).            ',/,&
'@    Les cases suivantes devraient etre nulles.              ',/,&
'@                                                            ',/,&
'@    Ce n''est cependant pas le cas :                        ',/,&
'@                                                            ',/,&
'@        IMOM    IDFMOM(1,IMOM)                              ',/,&
'@  ----------------------------                              '  )
 8201 format(                                                     &
'@  ',I10   ,'        ',     I10                                 )
 8202 format(                                                     &
'@  ----------------------------                              ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier les parametres.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@      SUR LES VARIABLES COMPOSANT LES MOYENNES TEMPORELLES  ',/,&
'@                                                            ',/,&
'@    IDFMOM(',I10   ,',',I10   ,') = ',I10                    ,/,&
'@      Les valeurs negatives renvoient a des proprietes      ',/,&
'@      physiques, or il n y en a que NPRMAX = ', I10          ,/,&
'@    La valeur de IDFMOM est donc erronee.                   ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier les parametres.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8211 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@      SUR LES VARIABLES COMPOSANT LES MOYENNES TEMPORELLES  ',/,&
'@                                                            ',/,&
'@    IDFMOM(',I10   ,',',I10   ,') = ',I10                    ,/,&
'@      Les valeurs positives renvoient a des variables de    ',/,&
'@      calcul, or il n y en a que NVAR   = ', I10             ,/,&
'@    La valeur de IDFMOM est donc erronee.                   ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier les parametres.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8212 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@      SUR LES VARIABLES COMPOSANT LES MOYENNES TEMPORELLES  ',/,&
'@                                                            ',/,&
'@    La valeur                                               ',/,&
'@      IDFMOM(',I10   ,',',I10   ,') = ',I10                  ,/,&
'@      n''est pas une propriete associee aux cellules        ',/,&
'@      (IPPROC(',I10   ,') = ',I10   ,')                     ',/,&
'@    La valeur de IDFMOM est donc erronee.                   ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier les parametres.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8213 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@      SUR LES VARIABLES COMPOSANT LES MOYENNES TEMPORELLES  ',/,&
'@                                                            ',/,&
'@    Le tableau IDFMOM(JJ,IMOM) pour IMOM = ',I10             ,/,&
'@      doit etre renseigne continuement. Or ici,             ',/,&
'@      IDFMOM(',I10,',IMOM) est non nul (=',I10   ,')        ',/,&
'@      alors qu il existe II < JJ pour lequel                ',/,&
'@      IDFMOM(II,IMOM) est nul.                              ',/,&
'@    La valeur de IDFMOM est donc erronee.                   ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier les parametres.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8214 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@  SUR L''INSTANT DE DEBUT DE CALCUL DES MOYENNES TEMPORELLES',/,&
'@                                                            ',/,&
'@    La variable NTDMOM(IMOM)   pour IMOM = ',I10             ,/,&
'@      doit etre renseignee pour indiquer                    ',/,&
'@      a partir de quel pas de temps (absolu) doit etre      ',/,&
'@      calculee la moyenne temporelle IMOM correspondante.   ',/,&
'@      NTDMOM(IMOM) vaur ici ',I10                            ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier les parametres.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8300 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@      ERREUR A L''OUVERTURE DU FICHIER SUITE AUXILIAIRE     ',/,&
'@                                                            ',/,&
'@    Pour permettre de realiser une suite de calcul en       ',/,&
'@      prenant en compte les moyennes temporelles,           ',/,&
'@      on cherche a relire le fichier suite auxiliaire.      ',/,&
'@    Une erreur se produit a son ouverture.                  ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier l''existence et le nom (',A13,') du            ',/,&
'@        fichier suite dans le repertoire de travail.        ',/,&
'@    Il est possible de s''affranchir du fichier suite       ',/,&
'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,&
'@      ou en ne calculant pas de moyennes.                   ',/,&
'@      Voir alors usipsu.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8301 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========                                     AUXILIAIRE',/,&
'@      TYPE DE FICHIER INCORRECT                             ',/,&
'@                                                            ',/,&
'@    Le fichier ',A13      ,' ne semble pas etre un fichier  ',/,&
'@      suite auxiliaire.                                     ',/,&
'@                                                            ',/,&
'@    Pour permettre de realiser une suite de calcul en       ',/,&
'@      prenant en compte les moyennes temporelles,           ',/,&
'@      on cherche a relire le fichier suite auxiliaire.      ',/,&
'@    Une erreur se produit a sa lecture.                     ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier que le fichier suite utilise correspond bien   ',/,&
'@        a un fichier suite auxiliaire.                      ',/,&
'@    Il est possible de s''affranchir du fichier suite       ',/,&
'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,&
'@      ou en ne calculant pas de moyennes.                   ',/,&
'@      Voir alors usipsu.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8311 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========  AUXILIAIRE (varpos)                          ',/,&
'@                                                            ',/,&
'@    Pour permettre de realiser une suite de calcul en       ',/,&
'@      prenant en compte les moyennes temporelles,           ',/,&
'@      on cherche a relire le fichier suite auxiliaire.      ',/,&
'@    Erreur a la lecture du nombre de moyennes               ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire.                   ',/,&
'@    Il est possible de s''affranchir du fichier suite       ',/,&
'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,&
'@      ou en ne calculant pas de moyennes.                   ',/,&
'@      Voir alors usipsu.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8312 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========  AUXILIAIRE (varpos)                          ',/,&
'@                                                            ',/,&
'@    Pour permettre de realiser une suite de calcul en       ',/,&
'@      prenant en compte les moyennes temporelles,           ',/,&
'@      on cherche a relire le fichier suite auxiliaire.      ',/,&
'@    Erreur a la lecture du numero du pas de temps precedent ',/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire.                   ',/,&
'@    Il est possible de s''affranchir du fichier suite       ',/,&
'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,&
'@      ou en ne calculant pas de moyennes.                   ',/,&
'@      Voir alors usipsu.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8313 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE         ',/,&
'@    =========  AUXILIAIRE (varpos)                          ',/,&
'@                                                            ',/,&
'@    Pour permettre de realiser une suite de calcul en       ',/,&
'@      prenant en compte les moyennes temporelles,           ',/,&
'@      on cherche a relire le fichier suite auxiliaire.      ',/,&
'@    Erreur a la lecture du numero du cumul temporel         ',/,&
'@      de la moyenne                  ', I10                  ,/,&
'@      associee a l''ancienne moyenne ', I10                  ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut pas etre execute.                     ',/,&
'@                                                            ',/,&
'@    Verifier le fichier suite auxiliaire.                   ',/,&
'@    Il est possible de s''affranchir du fichier suite       ',/,&
'@      auxiliaire en reinitialisant les moyennes (IMOOLD)    ',/,&
'@      ou en ne calculant pas de moyennes.                   ',/,&
'@      Voir alors usipsu.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8390 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ERREUR A LA FERMETURE DU FICHIER SUITE      ',/,&
'@    =========   AUXILIAIRE (varpos)                         ',/,&
'@                                                            ',/,&
'@    Probleme sur le fichier de nom (',A13,')                ',/,&
'@                                                            ',/,&
'@    Le calcul se poursuit...                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER SUITE (VARPOS)',/,&
'@    =========                                     AUXILIAIRE',/,&
'@      ERREUR A LA LECTURE DU FICHIER SUITE                  ',/,&
'@                                                            ',/,&
'@    On souhaite faire correspondre la moyenne   ',I10        ,/,&
'@            du present calcul avec la moyenne   ',I10        ,/,&
'@            du calcul precedent, or, le numero des anciennes',/,&
'@            moyennes  doit etre strictement positif et      ',/,&
'@            inferieur ou egal a                 ',I10        ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier les valeurs de IMOOLD.                         ',/,&
'@    Verifier que le fichier suite est le bon.               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8401 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@      SUR LA CORRESPONDANCE AVEC LES ANCIENNES MOYENNES     ',/,&
'@                                                            ',/,&
'@    On souhaite initialiser la nouvelle moyenne IMOM = ',I10 ,/,&
'@      en relisant la moyenne IMOOLD(IMOM) = ',I10            ,/,&
'@      dans le fichier suite.                                ',/,&
'@    Or on a specifie que le pas de temps initial pour le    ',/,&
'@      calcul de la moyenne IMOM etait NTDMOM(IMOM) = ',I10   ,/,&
'@      et le calcul stocke dans le fichier suite correspond  ',/,&
'@      au pas de temps = ',I10                                ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier la valeur de  NTDMOM                           ',/,&
'@      verifier que le fichier suite est le bon.             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8402 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@      SUR LA CORRESPONDANCE AVEC LES ANCIENNES MOYENNES     ',/,&
'@                                                            ',/,&
'@    Deux moyennes distinctes ',I10,' et ', I10               ,/,&
'@      sont initialises avec la meme ',I10                    ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Verifier IMOOLD dans usipsu.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8403 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@      SUR LA CORRESPONDANCE AVEC LES ANCIENNES MOYENNES     ',/,&
'@                                                            ',/,&
'@    On ne souhaite pas faire un calcul suite puisque        ',/,&
'@      ISUITE = ',I10                                         ,/,&
'@      mais on a renseigne le tableau IMOOLD de              ',/,&
'@      correspondance nouvelles -> anciennes moyennes :      ',/,&
'@      IMOOLD(',I10,') = ',I10                                ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Ne pas modifier IMOOLD.                                 ',/,&
'@      ou realiser un calcul suite.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8410 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA CONSTRUCTION DE IDTMOM  (varpos) ',/,&
'@    =========                                               ',/,&
'@      CALCUL DES MOYENNES TEMPORELLES                       ',/,&
'@                                                            ',/,&
'@    IDTMOM(',I10,') = ', I10                                 ,/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8900 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE VARIABLES A SUIVRE TROP GRAND                ',/,&
'@                                                            ',/,&
'@  Le type de calcul defini                                  ',/,&
'@    correspond a un nombre de variables a suivre dans       ',/,&
'@    le listing et le post-processing egal a      ',I10       ,/,&
'@  Le nombre de variables a suivre maximal prevu             ',/,&
'@                      dans paramx.h est NVPPMX = ',I10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres.                                  ',/,&
'@  Contacter l assistance.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 7030 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    INCORRECT SCALAR ASSOCIATED TO A VARIANCE               ',/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (user scalar                    ',I10   ,') represents  ',/,&
'@     one variance as                                        ',/,&
'@    ISCAVR(',I10   ,') has a value',I10   ,' (non-zero )    ',/,&
'@  However, the values of ISCAVR must be                     ',/,&
'@    larger or equal to                0                     ',/,&
'@    lower or equal to       NSCAL = ',I10                    ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   ISCAVR.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7031 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    INCORRECT SCALAR ASSOCIATED TO A VARIANCE               ',/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (scalar of  a particualr physics',I10   ,') represents  ',/,&
'@     one variance as         ',/,                         &
'@  ISCAVR(ISCAPP(',I10   ,')) has a value',I10   ,'(non-zero)',/,&
'@  However, the values of ISCAVR must be                     ',/,&
'@    larger or equal to                       0              ',/,&
'@    lower or equal to       NSCAL = ',I10                    ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   ISCAVR.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7040 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10                                           ,/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (user scalar                    ',I10   ,') represents  ',/,&
'@    the variance of fluctuations of a scalar ',I10           ,/,&
'@    (user scalar                    ',I10   ,') since       ',/,&
'@    ISCAVR(',I10   ,') has a value',I10   ,'(non-zero)      ',/,&
'@                                                            ',/,&
'@  The diffusivity index        IVISLS(',I10  ,')            ',/,&
'@    has not been set.                                       ',/,&
'@  It will be automatically set equal to that of the         ',/,&
'@    associated scalar ',I10                                  ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7041 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10                                           ,/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (user scalar                    ',I10   ,') represents  ',/,&
'@    the variance of fluctuations of a scalar ',I10           ,/,&
'@    (scalar of a specific physics ',I10   ,') since         ',/,&
'@    ISCAVR(',I10   ,') has a value ',I10   ,' (non-zero)    ',/,&
'@                                                            ',/,&
'@  The diffusivity index        IVISLS(',I10  ,')            ',/,&
'@    has not been set.                                       ',/,&
'@  It will be automatically set equal to that of the         ',/,&
'@    associated scalar ',I10                                  ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7042 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10                                           ,/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (scalar of a specific physics ',I10   ,') represents    ',/,&
'@    the variance of fluctuations of a scalar ',I10           ,/,&
'@    (user scalar                    ',I10   ,') since       ',/,&
'@ ISCAVR(ISCAPP(',I10   ,')) has a value ',I10   ,'(non-zero)',/,&
'@                                                            ',/,&
'@  The diffusivity index        IVISLS(ISCAPP(',I10  ,'))    ',/,&
'@    has not been set.                                       ',/,&
'@  It will be automatically set equal to that of the         ',/,&
'@    associated scalar ',I10                                  ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7043 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10                                           ,/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (scalar of a specific physics ',I10   ,') represents    ',/,&
'@    the variance of fluctuations of a scalar ',I10           ,/,&
'@    (scalar of a specific physics ',I10   ,') since         ',/,&
'@ ISCAVR(ISCAPP(',I10   ,')) has a value',I10   ,' (non-zero)',/,&
'@                                                            ',/,&
'@  The diffusivity index        IVISLS(ISCAPP(',I10  ,'))    ',/,&
'@    has not been set.                                       ',/,&
'@  It will be automatically set equal to that of the         ',/,&
'@    associated scalar ',I10                                  ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7050 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10                                           ,/,&
'@                                                            ',/,&
'@  The variable diffusivity index of the                     ',/,&
'@    user scalar                       ',I10                  ,/,&
'@    IVISLS(',I10   ,')                                      ',/,&
'@    must be set equal to 0 or 1.                            ',/,&
'@    Here it is  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7051 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10                                           ,/,&
'@                                                            ',/,&
'@  The variable diffusivity index of the                     ',/,&
'@    scalar of a specific physics    ',I10                    ,/,&
'@    IVISLS(ISCAPP(',I10   ,'))                              ',/,&
'@    must be set equal to 0 or 1.                            ',/,&
'@    Here it is  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IVISLS.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7070 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    THE INDEX  IORTVM    CANNOT HAVE VALUES OTHER THAN      ',/,&
'@      0 OR 1.                                               ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   usalin.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF VARIABLES TOO LARGE                          ',/,&
'@                                                            ',/,&
'@  The type of calculation defined                           ',/,&
'@    corresponds  to the following number of properties      ',/,&
'@      at the cell centres          : NPROCE = ',I10          ,/,&
'@  The maximum number of properties allowed                  ',/,&
'@                      in   paramx   is  NPROMX = ',I10       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@  NPROMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8022 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0 OR 1               ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8101 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 1 OR 2               ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8111 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    WITH A SECOND ORDER SCHEME IN TIME: ISCHTP = ', I10      ,/,&
'@    IT IS NECESSARY TO USE A CONSTANT AND UNIFORM TIME STEP ',/,&
'@    BUT IDTVAR = ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8112 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2ND ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    WITH K-EPSILON (ITURB = ',I10,' )'                       ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8113 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   :      AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 1st ORDER SCHEME HAS BEEN IMPOSSED   (ISCHTP = ',I10   ,/,&
'@    FOR LES (ITURB = ',I10,' )'                              ,/,&
'@                                                            ',/,&
'@  The calculation will   be executed                        ',/,&
'@                                                            ',/,&
'@  It is recommended to verify  parameters.                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8114 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR PHI_FBAR (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8117 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA                    ',/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR BL-V2/K  (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8115 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR K-OMEGA   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-omega.                 ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8116 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR SPALART   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of Spalart-Allmaras.        ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8121 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8131 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8132 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA'                    ,/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  With the chosen turbulence model   , ITURB = ',I10         ,/,&
'@    the value of ISTO2T (extrapolation of the source terms  ',/,&
'@    for the turbulent variables) cannot be modified         ',/,&
'@    yet ISTO2T has been forced to ',I10                      ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8141 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITAL DATA FOR SCALARS    ',I10 ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8200 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@      ON THE LIST OF TEMPORAL AVERAGES                      ',/,&
'@                                                            ',/,&
'@    The value of IDFMOM(1,',I10   ,' is  ',I10               ,/,&
'@      this indicates that ',I10   ,' temporal averages have ',/,&
'@      been defined to find out the   IMOM = ',I10            ,/,&
'@      first locations of the array IDFMOM(.,IMOM).          ',/,&
'@    The follwing locations should be zero.                  ',/,&
'@                                                            ',/,&
'@    This however, is not the case  :                        ',/,&
'@                                                            ',/,&
'@        IMOM    IDFMOM(1,IMOM)                              ',/,&
'@  ----------------------------                              '  )
 8201 format(                                                     &
'@  ',I10   ,'        ',     I10                                 )
 8202 format(                                                     &
'@  ----------------------------                              ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify   parameters.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8210 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@    ON THE VARIABLES THAT CONSTITUTE THE TEMPORAL AVERAGES  ',/,&
'@                                                            ',/,&
'@    IDFMOM(',I10   ,',',I10   ,') = ',I10                    ,/,&
'@      The negative value reflect physical properties        ',/,&
'@      but there is none in          NPRMAX = ', I10          ,/,&
'@    The value of IDFMOM is wrongly set.                     ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify   parameters.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8211 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@    ON THE VARIABLES THAT CONSTITUTE THE TEMPORAL AVERAGES  ',/,&
'@                                                            ',/,&
'@    IDFMOM(',I10   ,',',I10   ,') = ',I10                    ,/,&
'@      The positive values reflect variables of the          ',/,&
'@      calculation, yet there none in NVAR   = ', I10         ,/,&
'@    The value of IDFMOM is wrongly set.                     ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify   parameters.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8212 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@    ON THE VARIABLES THAT CONSTITUTE THE TEMPORAL AVERAGES  ',/,&
'@                                                            ',/,&
'@    The value                                               ',/,&
'@      IDFMOM(',I10   ,',',I10   ,') = ',I10                  ,/,&
'@      is not a property associated with the cells           ',/,&
'@      (IPPROC(',I10   ,') = ',I10   ,')                     ',/,&
'@    The value of IDFMOM is wrongly set.                     ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify   parameters.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8213 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@    ON THE VARIABLES THAT CONSTITUTE THE TEMPORAL AVERAGES  ',/,&
'@                                                            ',/,&
'@    The array  IDFMOM(JJ,IMOM) for  IMOM = ',I10             ,/,&
'@      must be assigned continuously.     Yet here,          ',/,&
'@      IDFMOM(',I10,',IMOM) is not zero (=',I10   ,')        ',/,&
'@      while it exists    II < JJ for which                  ',/,&
'@      IDFMOM(II,IMOM) is zero.                              ',/,&
'@    The value of IDFMOM is wrongly set.                     ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify   parameters.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8214 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@  ON THE INSTANT OF THE START OF THE CALCULATION OF THE     ',/,&
'@  TEMPORAL AVERAGES                                         ',/,&
'@                                                            ',/,&
'@    The variable NTDMOM(IMOM)   for IMOM = ',I10             ,/,&
'@      must be set to indicate from which                    ',/,&
'@      time step (absolute) the calculation ot the           ',/,&
'@      corresponding average IMOM must start.                ',/,&
'@      NTDMOM(IMOM) here is  ',I10                            ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify   parameters.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8300 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP WHILE READING THE AUXILIARY RESTARTING ',/,&
'@    =========                                           FILE',/,&
'@      ERROR OPENING THE RESTARTING AUXILIARY FILE           ',/,&
'@                                                            ',/,&
'@    In order to restart a calculation taking into account   ',/,&
'@      the temporal averages, it is necessary to             ',/,&
'@      read the auxiliary restarting file.                   ',/,&
'@    An error has occur during while opening it.             ',/,&
'@                                                            ',/,&
'@    The calculation cannot be executed                      ',/,&
'@                                                            ',/,&
'@    Verify the existence and the name (',A13,') of          ',/,&
'@        the restarting file on the working directory.       ',/,&
'@    It is possible to liberate the auxiliary restarting     ',/,&
'@      file and reinitialise the averages        (IMOOLD)    ',/,&
'@      or not to compute the averages.                       ',/,&
'@      Check parameters.                                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8301 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP WHILE READING THE AUXILIARY RESTARTING ',/,&
'@    =========                                           FILE',/,&
'@      INCORRECT FILE TYPE                                   ',/,&
'@                                                            ',/,&
'@    The file   ',A13      ,' does not look like a auxiliary ',/,&
'@      restarting file.                                      ',/,&
'@                                                            ',/,&
'@    In order to restart a calculation taking into account   ',/,&
'@      the temporal averages, it is necessary to             ',/,&
'@      read the auxiliary restarting file.                   ',/,&
'@    An error has occur during while reading it.             ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify  that the restarting file corresponds to an      ',/,&
'@       auxiliary restarting file.                           ',/,&
'@    It is possible to liberate the auxiliary restarting     ',/,&
'@      file and reinitialise teh averages        (IMOOLD)    ',/,&
'@      or not to compute the averages.                       ',/,&
'@      Check parameters.                                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8311 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP WHILE READING THE AUXILIARY RESTARTING ',/,&
'@    =========  FILE       (varpos)                          ',/,&
'@                                                            ',/,&
'@    In order to restart a calculation taking into account   ',/,&
'@      the temporal averages, it is necessary to             ',/,&
'@      read the auxiliary restarting file.                   ',/,&
'@    Error while reading the number of averages              ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify  the auxiliary restarting file.                  ',/,&
'@    It is possible to liberate the auxiliary restarting     ',/,&
'@      file and reinitialise teh averages        (IMOOLD)    ',/,&
'@      or not to compute the averages.                       ',/,&
'@      Check parameters.                                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8312 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP WHILE READING THE AUXILIARY RESTARTING ',/,&
'@    =========  FILE       (varpos)                          ',/,&
'@                                                            ',/,&
'@    In order to restart a calculation taking into account   ',/,&
'@      the temporal averages, it is necessary to             ',/,&
'@      read the auxiliary restarting file.                   ',/,&
'@    Error while reading the number of averages              ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify  the auxiliary restarting file.                  ',/,&
'@    It is possible to liberate the auxiliary restarting     ',/,&
'@      file and reinitialise teh averages        (IMOOLD)    ',/,&
'@      or not to compute the averages.                       ',/,&
'@      Check parameters.                                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8313 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP WHILE READING THE AUXILIARY RESTARTING ',/,&
'@    =========  FILE       (varpos)                          ',/,&
'@                                                            ',/,&
'@    In order to restart a calculation taking into account   ',/,&
'@      the temporal averages, it is necessary to             ',/,&
'@      read the auxiliary restarting file.                   ',/,&
'@    Error while reading the number of temporal cumulative   ',/,&
'@      of the average                 ', I10                  ,/,&
'@      associated with the previous   ', I10                  ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify the auxiliary restarting file.                   ',/,&
'@    Verify  the auxiliary restarting file.                  ',/,&
'@    It is possible to liberate the auxiliary restarting     ',/,&
'@      file and reinitialise teh averages        (IMOOLD)    ',/,&
'@      or not to compute the averages.                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8390 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : ERROR WHILE CLOSING THE AUXILIARY RESTARTING',/,&
'@    =========   FILE       (varpos)                         ',/,&
'@                                                            ',/,&
'@    Problem with the file named    (',A13,')                ',/,&
'@                                                            ',/,&
'@    The calculation continues ...                           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8400 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP WHILE READING THE RESTARTING AUXILIARY ',/,&
'@    =========                            FILE  (VARPOS)     ',/,&
'@      ERROR WHILE READING THE RESTARTING FILE               ',/,&
'@                                                            ',/,&
'@    Trying to match the temporal average        ',I10        ,/,&
'@    of the present calculation with the average ',I10        ,/,&
'@    of the previous calculation, yet the number of previous ',/,&
'@    averages must be strictly positive and lower or         ',/,&
'@    equal to                                    ',I10        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify the values of    IMOOLD.                         ',/,&
'@    Verify  that the restarting file is correct.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8401 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@      ON THE CORRESPONDANCE WITH PREVIOUS AVERAGES          ',/,&
'@                                                            ',/,&
'@    Trying to initialise the new average        IMOM = ',I10 ,/,&
'@      while reading the average IMOOLD(IMOM) = ',I10         ,/,&
'@      in the restarting file                                ',/,&
'@    Yet it has been specified that initial time step for the',/,&
'@    calculation of the average IMOM is NTDMOM(IMOM) = ',I10  ,/,&
'@    and the calculation stored in the restarting file       ',/,&
'@    corresponds to a time step  = ',I10                      ,/,&
'@                                                            ',/,&
'@    The calculation cannot be executed                      ',/,&
'@                                                            ',/,&
'@    Modify the value of    NTDMOM.                          ',/,&
'@      Verify  that the restarting file is correct.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8402 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@      ON THE CORRESPONDANCE WITH PREVIOUS AVERAGES          ',/,&
'@                                                            ',/,&
'@    Two diferent averages    ',I10,' and ', I10              ,/,&
'@      are initialised with the same ',I10                    ,/,&
'@                                                            ',/,&
'@    The calculation cannot be executed                      ',/,&
'@                                                            ',/,&
'@    Verify   IMOOLD.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8403 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@      ON THE CORRESPONDANCE WITH PREVIOUS AVERAGES          ',/,&
'@                                                            ',/,&
'@    A restarting calculation cannot be executed since       ',/,&
'@      ISUITE = ',I10                                         ,/,&
'@      but the array IMOOLD matching the old -> new         ',/, &
'@      averages has been filled                              ',/,&
'@      IMOOLD(',I10,') = ',I10                                ,/,&
'@                                                            ',/,&
'@    The calculation cannot be executed                      ',/,&
'@                                                            ',/,&
'@    Do not modify   IMOOLD.                                 ',/,&
'@      or run a restart calculation.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8410 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE CONSTRUCTION OF IDTMOM  (varpos)',/,&
'@    =========                                               ',/,&
'@      COMPUTATION OF THE TEMPORAL AVERAGES                  ',/,&
'@                                                            ',/,&
'@    IDTMOM(',I10,') = ', I10                                 ,/,&
'@                                                            ',/,&
'@    The calculation cannot be executed                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8900 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NAME OF THE VARIABLE TO BE CONTINUED TOO LARGE         ',/,&
'@                                                            ',/,&
'@  The type of calcultion defined                            ',/,&
'@    corresponds to a number of variables to continue in     ',/,&
'@    the listing and  post-processing equal to    ',I10       ,/,&
'@  The maximum number of variables to continue in            ',/,&
'@                           paramx.h is  NVPPMX = ',I10       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   parameters.                                      ',/,&
'@  Contact help.                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!===============================================================================
! 5. End
!===============================================================================

return
end subroutine varpos

!===============================================================================

!> \function add_property_field_nd
!
!> \brief add field defining a property field defined on cells,
!>        with default options
!
!> It is recommended not to define property names of more than 16
!> characters, to get a clear execution listing (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     name          field name
!> \param[in]     label         field default label, or empty
!> \param[in]     dim           field dimension
!> \param[out]    iprop         matching field property id
!_______________________________________________________________________________

subroutine add_property_field_nd &
 ( name, label, dim, iprop )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(in)          :: dim
integer, intent(out)         :: iprop

! Local variables

integer  id, type_flag, location_id
logical  has_previous, interleaved

!===============================================================================

type_flag = FIELD_INTENSIVE + FIELD_PROPERTY
location_id = 1 ! variables defined on cells
has_previous = .false.
interleaved = .true.

! Test if the field has already been defined
call field_get_id_try(trim(name), id)
if (id .ge. 0) then
  write(nfecra,1000) trim(name)
  call csexit (1)
endif

! Create field

call field_create(name, type_flag, location_id, dim, interleaved, has_previous, &
                  id)

call field_set_key_int(id, keyvis, 1)
call field_set_key_int(id, keylog, 1)

if (len(trim(label)).gt.0) then
  call field_set_key_str(id, keylbl, trim(label))
endif

! Property number and mapping to field

iprop = nproce + 1
nproce = nproce + dim

call varpos_check_nproce

iprpfl(iprop) = id
ipproc(iprop) = iprop
nomprp(iprop) = label

return

!---
! Formats
!---

#if defined(_CS_LANG_FR)
 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERREUR :    ARRET A L''ENTREE DES DONNEES               ',/,&
'@    ========                                                ',/,&
'@     LE CHAMP : ', a, 'EST DEJA DEFINI.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#else
 1000 format(                                                     &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP              ',/,&
'@    ======                                                  ',/,&
'@     FIELD: ', a, 'HAS ALREADY BEEN DEFINED.                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
#endif

end subroutine add_property_field_nd

!===============================================================================

!> \function add_property_field
!
!> \brief add field defining a property field defined on cells,
!>        with default options
!
!> It is recommended not to define property names of more than 16
!> characters, to get a clear execution listing (some advanced writing
!> levels take into account only the first 16 characters).
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     name          field name
!> \param[in]     label         field default label, or empty
!> \param[out]    iprop         matching field property id
!_______________________________________________________________________________

subroutine add_property_field &
 ( name, label, iprop )

!===============================================================================
! Module files
!===============================================================================

use field

!===============================================================================

implicit none

! Arguments

character(len=*), intent(in) :: name, label
integer, intent(out)         :: iprop

!===============================================================================

call add_property_field_nd(name, label, 1, iprop)

return

end subroutine add_property_field

!===============================================================================

!> \function hide_property
!
!> \brief disable logging and postprocessing for an unused property
!
!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!> \param[in]     iprop         property id
!_______________________________________________________________________________

subroutine hide_property &
 ( iprop )

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar
use field

!===============================================================================

implicit none

! Arguments

integer, intent(in) :: iprop

! Local variables

integer  id

!===============================================================================

ipppro(iprop) = 1
id = iprpfl(iprop)
call field_set_key_int(id, keyvis, 0)
call field_set_key_int(id, keylog, 0)

return

end subroutine hide_property

!===============================================================================

!> \function varpos_check_nproce

!> \brief check npromx is sufficient for the required number of properties.

!-------------------------------------------------------------------------------
! Arguments
!______________________________________________________________________________.
!  mode           name          role                                           !
!______________________________________________________________________________!
!_______________________________________________________________________________

subroutine varpos_check_nproce

!===============================================================================
! Module files
!===============================================================================

use paramx
use dimens
use entsor
use numvar

!===============================================================================

implicit none

! Arguments

! Local variables

if (nproce .gt. npromx) then
  write(nfecra,1000) nproce, npromx
  call csexit (1)
endif

return

!---
! Formats
!---

#if defined(_CS_LANG_FR)

 1000 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ERREUR :    ARRET A L''ENTREE DES DONNEES'               ,/,&
'@    ========'                                                ,/,&
'@     NOMBRE DE PROPRIETES TROP GRAND'                        ,/,&
'@'                                                            ,/,&
'@  Le type de calcul defini'                                  ,/,&
'@    correspond a un nombre de proprietes NPROCE >= ', i10    ,/,&
'@  Le nombre de proprietes maximal prevu'                     ,/,&
'@                      dans paramx    est NPROMX  = ', i10    ,/,&
'@'                                                            ,/,&
'@  Le calcul ne sera pas execute.'                            ,/,&
'@'                                                            ,/,&
'@  Verifier les parametres'                                   ,/,&
'@'                                                            ,/,&
'@  Si NPROMX est augmente, le code doit etre reinstalle.'     ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#else

 1000 format(                                                     &
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/,&
'@ @@ ERROR:      STOP AT THE INITIAL DATA SETUP'              ,/,&
'@    ======'                                                  ,/,&
'@     NUMBER OF PROPERTIES TOO LARGE'                         ,/,&
'@'                                                            ,/,&
'@  The type of calculation defined'                           ,/,&
'@    corresponds to a number of properties NPROCE >= ', i10   ,/,&
'@  The maximum number of properties allowed'                  ,/,&
'@                      in   paramx     is  NPROMX  = ', i10   ,/,&
'@'                                                            ,/,&
'@  The calculation cannot be executed'                        ,/,&
'@'                                                            ,/,&
'@  Verify   parameters.'                                      ,/,&
'@'                                                            ,/,&
'@  If NVARMX is increased, the code must be reinstalled.'     ,/,&
'@'                                                            ,/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@'                                                            ,/)

#endif

end subroutine varpos_check_nproce

!===============================================================================
