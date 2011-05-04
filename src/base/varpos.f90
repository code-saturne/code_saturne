!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2009 EDF S.A., France

!     contact: saturne-support@edf.fr

!     The Code_Saturne Kernel is free software; you can redistribute it
!     and/or modify it under the terms of the GNU General Public License
!     as published by the Free Software Foundation; either version 2 of
!     the License, or (at your option) any later version.

!     The Code_Saturne Kernel is distributed in the hope that it will be
!     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with the Code_Saturne Kernel; if not, write to the
!     Free Software Foundation, Inc.,
!     51 Franklin St, Fifth Floor,
!     Boston, MA  02110-1301  USA

!-------------------------------------------------------------------------------

subroutine varpos &
!================

(nmodpp )

!===============================================================================
!  FONCTION  :
!  ---------

! INIT DES POSITIONS DES VARIABLES SELON
!   LE TYPE DE CALCUL INDIQUE PAR L'UTILISATEUR

!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! nmodpp           ! e  ! --> ! nombre de modeles phys.part. actives           !
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

!===============================================================================

implicit none

! Arguments

integer       nmodpp

! Local variables


character     rubriq*64, cmoy4*4, cindfm*4
character     ficsui*32
integer       ivar  , ipp   , iscal , iphas , irphas , iprop , iprofl
integer       ii    , jj    , kk    , ll
integer       iok   , ippok , ipppst
integer       iflum , icondl
integer       iest  , iiflaa, ivisph, iprofa, ipropp
integer       imom  , jmom  , imold , jmold , jbmomt, jtcabs
integer       iiplus, iimoin, jmomok, idto  , jdto
integer       ierror, irtyp , itysup, nbval
integer       idffin, idfmji, imomr , ilecec
integer       ilsmom, ivers , inmfin
integer       impamx
integer       nfmtmo, nberro
integer       idtold(nbmomx)
integer       nprayc, nprayb

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
! 1. PREMIER APPEL : VERIFICATION DU NOMBRE DE PHASES
!                    RETURN
!===============================================================================

if(ipass.eq.1) then

  iok = 0

  if(nphas.le.0) then
    write(nfecra,5000) nphas
    iok = iok + 1
  endif
  if(nphas.gt.nphsmx) then
    write(nfecra,5001) nphas, nphsmx, nphas
    iok = iok + 1
  endif

  if(iok.ne.0) then
    call csexit (1)
    !==========
  endif

  return

endif

!===============================================================================
! 2. SECOND APPEL : CALCUL DE NSCAPP
!                   VERIFICATION DU NOMBRE DE SCALAIRES
!                   CONSTRUCTION DE ISCAPP
!                   CALCUL DE NSCAL
!                   RETURN

!  C'est juste avant ce second appel que les modeles de combustion
!    auront ete renseignes. C'est dans la section ci-dessous qu'on en
!    en deduira NSCAPP (avant meme les verifications).
!  A la sortie de cette section, NSCAL, NSCAUS et NSCAPP sont connus.
!  On renseignera egalement ici les valeurs de IPHSCA, ISCAVR, IVISLS
!    pour les scalaires physiques particulieres en question.
!  On en profite aussi pour remplir ITYTUR puisque ITURB vient d'etre
!    defini.
!===============================================================================

if(ipass.eq.2) then

! ---> Remplissage de ITYTUR
  do iphas = 1, nphas
    itytur(iphas) = iturb(iphas)/10
  enddo

! ---> Coherence modele
!     Rq : ATTENTION il faudrait renforcer le blindage

  iok   = 0
  nmodpp = 0
  do ipp = 2, nmodmx
    if ( ippmod(ipp).ne.-1 ) then
      nmodpp = nmodpp+1
      ippok = ipp
    endif
  enddo
  if ( nmodpp.gt.1 ) then
    write(nfecra,6000)
    iok = iok + 1
  endif

  if ( nmodpp.eq.1 ) then
    if ( ippmod(ippok).lt.0 .or. ippmod(ippok).gt.5 ) then
      write(nfecra,6001)
      iok = iok + 1
    endif
  endif

  if(iok.ne.0) then
    call csexit (1)
    !==========
  endif

! ---> On positionne l'indicateur global IPPMOD(IPHPAR)
!         0 : pas de physique particuliere
!         1 : physique particuliere enclenchee
!         2 : physique particuliere avec definition du coefficient
!             d'absorption par fichier parametrique pour le rayonnement
  ippmod(iphpar) = 0
  if (nmodpp.gt.0) then
    ippmod(iphpar) = 1
    if (ippmod(icompf).eq.-1 .and. ippmod(iatmos).eq.-1           &
                             .and. ippmod(iaeros).eq.-1)          &
         ippmod(iphpar) = 2
  endif

! ---> Lecture donnees thermochimie

  call pplecd
  !==========

! ---> Calcul de NSCAPP

  call ppcsca
  !==========


! ---> Verifications

  iok = 0

  if(nscaus.lt.0) then
    write(nfecra,6010) nscaus
    iok = iok + 1
  endif

  if(nscaus.gt.0 .or. nscapp.gt.0) then
    if((nscaus+nscapp).gt.nscamx) then

      if(nscapp.le.0) then
        write(nfecra,6011)                                        &
             nscaus,       nscamx,nscamx       ,nscaus
      else
        write(nfecra,6012)                                        &
             nscaus,nscapp,nscamx,nscamx-nscapp,nscaus+nscapp
      endif
      iok = iok + 1
    endif
  endif

  if(iok.ne.0) then
    call csexit (1)
    !==========
  endif

! ---> Calcul de ISCAPP et NSCAL
!      On prefere que l'identite porte sur les scalaires utilisateurs,
!        ca minimisera peut etre des erreurs utilisateur

  iscal = max(0,nscaus)
  if (nscapp.gt.0) then
    do ii = 1, nscapp
      iscal = iscal + 1
      iscapp(ii) = iscal
    enddo

    call ppvarp
    !==========

  endif

  nscal = iscal

  return

endif


!===============================================================================
! 3. TROISIEME APPEL : VERIFICATIONS ET
!        POSITIONNEMENT DES VARIABLES  : IPR, IU ... ISCA, NVAR
!                       ET DES PROPRIETES PPALES
!===============================================================================

if(ipass.eq.3) then


! ---> 3.1 VERIFICATIONS
!      -----------------

  iok = 0

! ---  NPHAS et NSCAL ont deja ete verifies, mais on ne sait jamais.

  if(nphas.le.0) then
    write(nfecra,5000) nphas
    iok = iok + 1
  endif
  if(nphas.gt.nphsmx) then
    write(nfecra,5001) nphas, nphsmx, nphas
    iok = iok + 1
  endif

  if(nscal.lt.0) then
    write(nfecra,7010) nscal, nscaus, nscapp
    iok = iok + 1
  endif
  if(nscal.gt.nscamx) then
    write(nfecra,7011) nscal, nscamx, nscaus, nscapp, nscal
    iok = iok + 1
  endif

! --- IPHSCA(ISCAL) doit etre compris entre 0 et NPHAS.

  if(nscaus.gt.0) then
    do ii = 1, nscaus
      iscal = ii
      if(iphsca(iscal).gt.nphas.or.iphsca(iscal).lt.0) then
        write(nfecra,7020) iscal, iscal, ii,ii,                   &
             iphsca(iscal), nphas
        iok = iok + 1
      endif
    enddo
  endif
  if(nscapp.gt.0) then
    do ii = 1, nscapp
      iscal = iscapp(ii)
      if(iphsca(iscal).gt.nphas.or.iphsca(iscal).lt.0) then
        write(nfecra,7021) iscal, iscal, ii,ii,                   &
             iphsca(iscal), nphas
        iok = iok + 1
      endif
    enddo
  endif

! --- ISCAVR(ISCAL) doit etre compris entre 0 et NSCAL.

  if(nscaus.gt.0) then
    do ii = 1, nscaus
      iscal = ii
      if(iscavr(iscal).gt.nscal.or.iscavr(iscal).lt.0) then
        write(nfecra,7030) iscal,ii,ii,iscavr(iscal),nscal
        iok = iok + 1
      endif
    enddo
  endif
  if(nscapp.gt.0) then
    do ii = 1, nscapp
      iscal = iscapp(ii)
      if(iscavr(iscal).gt.nscal.or.iscavr(iscal).lt.0) then
        write(nfecra,7031) iscal,ii,ii,iscavr(iscal),nscal
        iok = iok + 1
      endif
    enddo
  endif


! --- IVISLS(ISCAL) doit etre non initialise pour les variances
!     Il prend la valeur du scalaire associe
!     Tous les tests qui suivent sont utilises pour le message
!       d'erreur eventuel.

  if(nscaus.gt.0) then
    do jj = 1, nscaus
      ii    = jj
      iscal = iscavr(ii)
!       Si on a une variance avec ivisls initialise : erreur
      if    ( (iscal.gt.0.and.iscal.le.nscal).and.                &
               ivisls(ii).ne.-1                        ) then
        ll = 0
        do kk = 1, nscaus
          if(       kk .eq.iscal) ll = kk
        enddo
        do kk = 1, nscapp
          if(iscapp(kk).eq.iscal) ll = -kk
        enddo
        if(ll.gt.0) then
          write(nfecra,7040)                                      &
               ii,ii,jj,iscal,ll,jj,iscal,jj,ivisls(iscal)
        else
          write(nfecra,7041)                                      &
               ii,ii,jj,iscal,-ll,jj,iscal,jj,ivisls(iscal)
        endif
        iok = iok + 1
!     Si on n'a pas une variance std mais que ivisls est incorrect : erreur
      elseif( (iscal.le.0 .or.iscal.gt.nscal).and.                &
             (ivisls(ii).ne.0.and.ivisls(ii).ne.1) ) then
        write(nfecra,7050) ii,jj,jj,ivisls(ii)
        iok = iok + 1
      endif
    enddo
  endif

  if(nscapp.gt.0) then
    do jj = 1, nscapp
      ii    = iscapp(jj)
      iscal = iscavr(ii)
!       Si on a une variance avec ivisls initialise : erreur
      if    ( (iscal.gt.0.and.iscal.le.nscal).and.                &
               ivisls(ii).ne.-1                        ) then
        ll = 0
        do kk = 1, nscaus
          if(       kk .eq.iscal) ll = kk
        enddo
        do kk = 1, nscapp
          if(iscapp(kk).eq.iscal) ll = -kk
        enddo
        if(ll.gt.0) then
          write(nfecra,7042)                                      &
               ii,ii,jj,iscal,ll,jj,iscal,jj,ivisls(iscal)
        else
          write(nfecra,7043)                                      &
               ii,ii,jj,iscal,-ll,jj,iscal,jj,ivisls(iscal)
        endif
        iok = iok + 1
!       Si on n'a pas une variance std mais que ivisls est incorrect : erreur
      elseif( (iscal.le.0 .or.iscal.gt.nscal).and.                &
              (ivisls(ii).ne.0.and.ivisls(ii).ne.1) ) then
        write(nfecra,7051)ii,jj,jj,ivisls(ii)
        iok = iok + 1
      endif
    enddo
  endif

!       On initialise les ivisls des variances
  if(nscal.gt.0) then
    do ii = 1, nscal
      iscal = iscavr(ii)
      if(iscal.gt.0.and.iscal.le.nscal) then
        ivisls(ii) = ivisls(iscal)
      endif
    enddo
  endif

! ---> IPHSCA
!      Pour les variances de fluctuations, les valeurs de IPHSCA
!        ne doivent pas avoir ete modifiees par l'utilisateur
!        Elles sont prises egales aux valeurs correspondantes
!        pour le scalaire associe.

  if(nscaus.gt.0) then
    do jj = 1, nscaus
      ii    = jj
      iscal = iscavr(ii)
      if(iscal.gt.0.and.iscal.le.nscal)then
        if(iphsca(ii).eq.0) then
          iphsca(ii) = iphsca(iscal)
        else
          ll = 0
          do kk = 1, nscaus
            if(       kk .eq.iscal) ll = kk
          enddo
          do kk = 1, nscapp
            if(iscapp(kk).eq.iscal) ll = -kk
          enddo
          if(ll.gt.0) then
            write(nfecra,7060)ii,                                 &
                 ii,jj,iscal,ll,jj,iscal,                         &
                 jj,iphsca(iscal)
          else
            write(nfecra,7061)ii,                                 &
                 ii,jj,iscal,-ll,jj,iscal,                        &
                 jj,iphsca(iscal)
          endif
          iok = iok + 1
        endif
      endif
    enddo
  endif

  if(nscapp.gt.0) then
    do jj = 1, nscapp
      ii    = iscapp(jj)
      iscal = iscavr(ii)
      if(iscal.gt.0.and.iscal.le.nscal)then
        if(iphsca(ii).eq.0) then
          iphsca(ii) = iphsca(iscal)
        else
          ll = 0
          do kk = 1, nscaus
            if(       kk .eq.iscal) ll = kk
          enddo
          do kk = 1, nscapp
            if(iscapp(kk).eq.iscal) ll = -kk
          enddo
          if(ll.gt.0) then
            write(nfecra,7062)ii,                                 &
                 ii,jj,iscal,ll,jj,iscal,                         &
                 jj,iphsca(iscal)
          else
            write(nfecra,7063)ii,                                 &
                 ii,jj,iscal,-ll,jj,iscal,                        &
                 jj,iphsca(iscal)
          endif
          iok = iok + 1
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

  if(iok.ne.0) then
    call csexit (1)
  endif


! ---> 3.2 POSITIONNEMENT DES VARIABLES  : IPR, IU ... ISCA, NVAR
!      --------------------------------

  ivar = 0

  do iphas = 1, nphas

! --- Pression : supposons ici qu'il n'y a qu'une pression
!       quelque soit le nombre de phases.
!     Le reste du code devrait se preter a l'extension, du fait
!       qu'il n'y a qu'ici qu'on fait cette hypothese (evidemment, il
!       faudrait ecrire un schema a deux pressions et remplacer navsto,
!       mais c'est un "detail"...relativement a ce qui nous occupe dans
!       le present sous-programme)

    if(iphas.eq.1) then
      ivar          = ivar + 1
      ipr   (iphas) = ivar
    else
      ipr   (iphas) = ipr(1)
    endif

! --- Vitesse
    ivar          = ivar + 1
    iu    (iphas) = ivar
    ivar          = ivar + 1
    iv    (iphas) = ivar
    ivar          = ivar + 1
    iw    (iphas) = ivar

! --- Turbulence
    if (itytur(iphas).eq.2) then
      ivar          = ivar + 1
      ik    (iphas) = ivar
      ivar          = ivar + 1
      iep   (iphas) = ivar
    elseif(itytur(iphas).eq.3) then
      ivar          = ivar + 1
      ir11  (iphas) = ivar
      ivar          = ivar + 1
      ir22  (iphas) = ivar
      ivar          = ivar + 1
      ir33  (iphas) = ivar
      ivar          = ivar + 1
      ir12  (iphas) = ivar
      ivar          = ivar + 1
      ir13  (iphas) = ivar
      ivar          = ivar + 1
      ir23  (iphas) = ivar
      ivar          = ivar + 1
      iep   (iphas) = ivar
    elseif(iturb(iphas).eq.50) then
      ivar          = ivar + 1
      ik    (iphas) = ivar
      ivar          = ivar + 1
      iep   (iphas) = ivar
      ivar          = ivar + 1
      iphi  (iphas) = ivar
      ivar          = ivar + 1
      ifb   (iphas) = ivar
    elseif(iturb(iphas).eq.60) then
      ivar          = ivar + 1
      ik    (iphas) = ivar
      ivar          = ivar + 1
      iomg  (iphas) = ivar
    elseif (iturb(iphas).eq.70) then
      ivar          = ivar + 1
      inusa (iphas) = ivar
    endif

  enddo

! --- Scalaires
  if(nscapp.ge.1) then
    do jj = 1, nscapp
      ii       = iscapp(jj)
      ivar     = ivar + 1
      isca(ii) = ivar
    enddo
  endif
  if(nscaus.ge.1) then
    do jj = 1, nscaus
      ii       = jj
      ivar     = ivar + 1
      isca(ii) = ivar
    enddo
  endif

! --- Vitesse de maillage en ALE
  if (iale.eq.1) then
    ivar = ivar + 1
    iuma = ivar
    ivar = ivar + 1
    ivma = ivar
    ivar = ivar + 1
    iwma = ivar
  endif

! --- Nombre total de variables
  nvar = ivar

! --- Verification de NVAR

  if(nvar.gt.nvarmx) then
    write(nfecra,7100)nvar,nvarmx,nvar
    call csexit (1)
  endif


! --- Maintenant on peut faire ceci :

  do iphas = 1, nphas
    istat (ipr(iphas)) = 0
    iconv (ipr(iphas)) = 0
    if (iturb(iphas).eq.50) then
      istat(ifb(iphas))  = 0
      iconv(ifb(iphas))  = 0
!     Pour fb, on sait qu'on a un terme diagonal, meme si ISTAT=0,
!       donc on ne decalera pas la diagonale
      idircl(ifb(iphas)) = 0
    endif
  enddo
  if (iale.eq.1) then
    istat(iuma) = 0
    iconv(iuma) = 0
    imgr (iuma) = 1
    istat(ivma) = 0
    iconv(ivma) = 0
    imgr (ivma) = 1
    istat(iwma) = 0
    iconv(iwma) = 0
    imgr (iwma) = 1
  endif


! ---> 3.3 POSITIONNEMENT DES PROPRIETES PRINCIPALES
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


  iprop = 0

!   Proprietes des phases : proprietes toujours presentes
  do iphas = 1, nphas
    iprop         = iprop + 1
    irom  (iphas) = iprop
    iprop         = iprop + 1
    iviscl(iphas) = iprop
    iprop         = iprop + 1
    ivisct(iphas) = iprop
    iprop         = iprop + 1
    icour (iphas) = iprop
    iprop         = iprop + 1
    ifour (iphas) = iprop
  enddo

!  Pression totale stockee dans IPRTOT, si on n'est pas en compressible
!  (sinon Ptot=P* !)
  if (ippmod(icompf).lt.0) then
    do iphas = 1, nphas
      if (iphas.eq.1) then
        iprop         = iprop + 1
        iprtot(1)     = iprop
      else
        iprtot(iphas) = iprtot(1)
      endif
    enddo
  endif

!  Proprietes des phases : CP s'il est variable
  do iphas = 1, nphas
    if(icp(iphas).ne.0) then
      iprop         = iprop + 1
      icp   (iphas) = iprop
    endif
  enddo

!  Proprietes des phases : Cs^2 si on est en LES dynamique
  do iphas = 1, nphas
    if(iturb(iphas).eq.41) then
      iprop         = iprop + 1
      ismago(iphas) = iprop
    else
      ismago(iphas) = -1
    endif
  enddo

!  Viscosite de maillage en ALE
  if (iale.eq.1) then
    iprop     = iprop + 1
    ivisma(1) = iprop
!     si la viscosite est isotrope, les trois composantes pointent
!       au meme endroit
    if (iortvm.eq.0) then
      ivisma(2) = iprop
      ivisma(3) = iprop
    else
      iprop     = iprop + 1
      ivisma(2) = iprop
      iprop     = iprop + 1
      ivisma(3) = iprop
    endif
  endif

!   Proprietes des phases : estimateurs d'erreur
  do iphas = 1, nphas
    do iest = 1, nestmx
      iprop              = iprop + 1
      iestim(iest,iphas) = iprop
    enddo
  enddo


!   Proprietes des scalaires : VISCLS si elle est variable
!     On utilisera IVISLS comme suit :
!       Pour le scalaire II
!         si IVISLS(II) = 0    : visc = VISLS0(II)
!         si IVISLS(II) .GT. 0 : visc = PROPCE(IEL ,IPPROC(IVISLS(II)))
!     Ceci permet de ne pas reserver des tableaux vides pour les
!       scalaires a viscosite constante
  if(nscal.ge.1) then
    do ii = 1, nscal
      if(ivisls(ii).ne.0) then
        iprop      = iprop + 1
        ivisls(ii) = iprop
      endif
    enddo
  endif


!   Proprietes des variables : flux de masse porteur

  iprofl = iprop
  do iphas = 1, nphas
    iprop               = iprop + 1
    if(iphas.eq.1) then
      ifluma(ipr (iphas)) = iprop
    endif
    ifluma(iu  (iphas)) = iprop
    ifluma(iv  (iphas)) = iprop
    ifluma(iw  (iphas)) = iprop
    if(itytur(iphas).eq.2) then
      ifluma(ik  (iphas)) = iprop
      ifluma(iep (iphas)) = iprop
    elseif(itytur(iphas).eq.3) then
      ifluma(ir11(iphas)) = iprop
      ifluma(ir22(iphas)) = iprop
      ifluma(ir33(iphas)) = iprop
      ifluma(ir12(iphas)) = iprop
      ifluma(ir13(iphas)) = iprop
      ifluma(ir23(iphas)) = iprop
      ifluma(iep (iphas)) = iprop
    elseif(iturb(iphas).eq.50) then
      ifluma(ik  (iphas)) = iprop
      ifluma(iep (iphas)) = iprop
      ifluma(iphi(iphas)) = iprop
      ifluma(ifb (iphas)) = iprop
    elseif(iturb(iphas).eq.60) then
      ifluma(ik  (iphas)) = iprop
      ifluma(iomg(iphas)) = iprop
    elseif(iturb(iphas).eq.70) then
      ifluma(inusa(iphas))= iprop
    endif
  enddo
  do iscal = 1, nscal
    ifluma(isca(iscal)) = ifluma(iu(iphsca(iscal)))
  enddo
  if (iale.eq.1) then
    ifluma(iuma) = ifluma(ipr(1))
    ifluma(ivma) = ifluma(ipr(1))
    ifluma(iwma) = ifluma(ipr(1))
  endif
!     Nombre total de flux de masse
!       IPROFL ressert plus bas.
  nfluma = iprop - iprofl

!     Numero max des proprietes ; ressert plus bas pour
!       ajouter celles relatives a la physique particuliere
  nprmax = iprop


! --- Positionnement dans les tableaux PROPCE, PROFA, PROFB

!   Au centre des cellules (tout sauf les flux de masse)
!     Pour les fluctuations, le pointeur de la diffusivite
!     envoie directement sur la diffusivite du scalaire associe.

!  On positionne en meme temps les pointeurs IPPPRO pour
!     le post traitement des proprietes physiques definies
!     aux cellules afin de ne pas en oublier.
!     Les pointeurs ont ete initialises a 1 (poubelle).
!       dans iniini (avant usini1).
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

  iprop = 0
  do iphas = 1, nphas

    iprop                 = iprop  + 1
    ipproc(irom  (iphas)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
    iprop                 = iprop  + 1
    ipproc(iviscl(iphas)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
    iprop                 = iprop  + 1
    ipproc(ivisct(iphas)) = iprop
    if (iturb(iphas).eq.0) then
      ipppro(iprop)         = 1
    else
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst
    endif
    iprop                 = iprop  + 1
    ipproc(icour (iphas)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst
    iprop                 = iprop  + 1
    ipproc(ifour (iphas)) = iprop
    ipppst                = ipppst + 1
    ipppro(iprop)         = ipppst

    if (ippmod(icompf).lt.0) then
      iprop                 = iprop  + 1
      ipproc(iprtot(iphas)) = iprop
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst
    endif

    if(icp(iphas).gt.0) then
      iprop                 = iprop + 1
      ipproc(icp   (iphas)) = iprop
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst
    endif

    if(ismago(iphas).ne.-1) then
      iprop                 = iprop  + 1
      ipproc(ismago(iphas)) = iprop
      ipppst                = ipppst + 1
      ipppro(iprop)         = ipppst
    endif

    if (iale.eq.1) then
      iprop             = iprop + 1
      ipproc(ivisma(1)) = iprop
      ipppst            = ipppst + 1
      ipppro(iprop)     = ipppst
      if (iortvm.eq.1) then
        iprop             = iprop + 1
        ipproc(ivisma(2)) = iprop
        ipppst            = ipppst + 1
        ipppro(iprop)     = ipppst
        iprop             = iprop + 1
        ipproc(ivisma(3)) = iprop
        ipppst            = ipppst + 1
        ipppro(iprop)     = ipppst
      endif
    endif

    do iest = 1, nestmx
      if(iescal(iest,iphas).gt.0) then
        iprop                      = iprop + 1
        ipproc(iestim(iest,iphas)) = iprop
        ipppst                     = ipppst + 1
        ipppro(iprop)              = ipppst
      endif
    enddo

  enddo

!     Conductivite electrique imaginaire :
!     La conductivite reelle et imaginaire sont dans le meme tableau.
!       En Joule, on ne reserve donc pas de propriete "viscosite"
!       pour le potentiel imaginaire.
!       Intervention 1/2

  if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
    do ii = 1, nscal
      if(ivisls(ii).gt.0.and.ii.ne.ipoti) then
        if(iscavr(ii).le.0) then
          iprop                 = iprop + 1
          ipproc(ivisls(ii)  )  = iprop
          ipppst                = ipppst + 1
          ipppro(iprop)         = ipppst
        endif
      endif
    enddo
  else

    do ii = 1, nscal
      if(ivisls(ii).gt.0) then
        if(iscavr(ii).le.0) then
          iprop                 = iprop + 1
          ipproc(ivisls(ii)  )  = iprop
          ipppst                = ipppst + 1
          ipppro(iprop)         = ipppst
        endif
      endif
    enddo

 endif

  do ii = 1, nscal
    if(ivisls(ii).gt.0) then
      if(iscavr(ii).gt.0.and.iscavr(ii).le.nscal) then
        ipproc(ivisls(ii)  )  = ipproc(ivisls(iscavr(ii))  )
      endif
    endif
  enddo

!     Conductivite electrique imaginaire :
!     La conductivite reelle et imaginaire sont dans le meme tableau.
!       En Joule, le pointeur sur la  "viscosite" du potentiel imaginaire
!         renvoie sur celle du potentiel reel.
!       Intervention 2/2

  if(ippmod(ieljou).eq.2 .or. ippmod(ieljou).eq.4) then
    ipproc(ivisls(ipoti)  )  = ipproc(ivisls(ipotr)  )
  endif

  nproce = iprop

!   Au centre des faces de bord (rho et flux de masse)

  iprop = 0
  do iphas = 1, nphas
    iprop                 = iprop + 1
    ipprob(irom  (iphas)) = iprop
  enddo
  do iflum = 1, nfluma
    iprop                 = iprop + 1
    ipprob(iprofl+iflum)  = iprop
  enddo
  nprofb = iprop

!   Au centre des faces internes (flux de masse)

  iprop = 0
  do iflum = 1, nfluma
    iprop                 = iprop + 1
    ipprof(iprofl+iflum)  = iprop
  enddo
  nprofa = iprop


! --- Modifications pour la physique particuliere
!      des entiers NPROCE, NPROFA, NPROFB

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


! --- Verification de NPROCE, NPROFA, NPROFB

  if(nproce.gt.npromx.or.                                         &
     nprofa.gt.npromx.or.nprofb.gt.npromx) then
    write(nfecra,7200)nproce, nprofa, nprofb, npromx,             &
         max(max(nproce,nprofa),nprofb)
    call csexit (1)
    !==========
  endif

  return

endif

!===============================================================================
! 4. QUATRIEME APPEL :
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

if(ipass.eq.4) then


! ---> 4.1 PROPRIETES ADDITIONNELLES POUR LES ET SCHEMA EN TEMPS
!      ---------------------------------------------------------

! --- Initialisations par defaut eventuelles et verifications
!       des options utilisees ci-dessous pour decider si l'on
!       reserve des tableaux supplementaires pour des grandeurs
!       au pas de temps precedent

  iok = 0

!     Pression hydrostatique
  if(iphydr.eq.0) then
    icalhy = 0
  elseif(iphydr.eq.1) then
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
  do iphas = 1, nphas
    if(ischtp(iphas).eq.-999) then
      if(itytur(iphas).eq.4) then
        ischtp(iphas) = 2
      else
        ischtp(iphas) = 1
      endif
    endif
  enddo

!     Schemas en temps : variables deduites
  do iphas = 1, nphas
!     Schema pour le Flux de masse
    if(istmpf(iphas).eq.-999) then
      if(ischtp(iphas).eq.1) then
        istmpf(iphas) = 1
      elseif(ischtp(iphas).eq.2) then
        istmpf(iphas) = 2
      endif
    endif
!     Masse volumique
    if(iroext(iphas).eq.-999) then
      if(ischtp(iphas).eq.1) then
        iroext(iphas) = 0
      elseif(ischtp(iphas).eq.2) then
!       Pour le moment par defaut on ne prend pas l'ordre 2
!              IROEXT(IPHAS) = 1
        iroext(iphas) = 0
      endif
    endif
!     Viscosite
    if(iviext(iphas).eq.-999) then
      if(ischtp(iphas).eq.1) then
        iviext(iphas) = 0
      elseif(ischtp(iphas).eq.2) then
!       Pour le moment par defaut on ne prend pas l'ordre 2
!              IVIEXT(IPHAS) = 1
        iviext(iphas) = 0
      endif
    endif
!     Chaleur massique
    if(icpext(iphas).eq.-999) then
      if(ischtp(iphas).eq.1) then
        icpext(iphas) = 0
      elseif(ischtp(iphas).eq.2) then
!       Pour le moment par defaut on ne prend pas l'ordre 2
!              ICPEXT(IPHAS) = 1
        icpext(iphas) = 0
      endif
    endif
!     Termes sources NS,
    if(isno2t(iphas).eq.-999) then
      if(ischtp(iphas).eq.1) then
        isno2t(iphas) = 0
!            ELSEIF(ISCHTP(IPHAS).EQ.2.AND.IVISSE(IPHAS).EQ.1) THEN
      elseif(ischtp(iphas).eq.2) then
!       Pour le moment par defaut on prend l'ordre 2
        isno2t(iphas) = 1
!              ISNO2T(IPHAS) = 0
      endif
    endif
!     Termes sources turbulence (k-eps, Rij, v2f ou k-omega)
!     On n'autorise de changer ISTO2T qu'en Rij (sinon avec
!       le couplage k-eps/omega il y a pb)
    if(isto2t(iphas).eq.-999) then
      if(ischtp(iphas).eq.1) then
        isto2t(iphas) = 0
      elseif(ischtp(iphas).eq.2) then
!       Pour le moment par defaut on ne prend pas l'ordre 2
!              ISTO2T(IPHAS) = 1
        isto2t(iphas) = 0
      endif
    else if( itytur(iphas).eq.2.or.iturb(iphas).eq.50             &
         .or.iturb(iphas).ne.60) then
      write(nfecra,8132) iphas,iturb(iphas),isto2t(iphas)
      iok = iok + 1
    endif
  enddo

  do iscal = 1, nscal
!     Termes sources Scalaires,
    iphas = iphsca(iscal)
    if(isso2t(iscal).eq.-999) then
      if(ischtp(iphas).eq.1) then
        isso2t(iscal) = 0
      elseif(ischtp(iphas).eq.2) then
!       Pour coherence avec Navier Stokes on prend l'ordre 2
!       mais de toute facon qui dit ordre 2 dit LES et donc
!       generalement pas de TS scalaire a interpoler.
        isso2t(iscal) = 1
!              ISSO2T(ISCAL) = 0
      endif
    endif
!     Diffusivite scalaires
    if(ivsext(iscal).eq.-999) then
      iphas = iphsca(iscal)
      if(ischtp(iphas).eq.1) then
        ivsext(iscal) = 0
      elseif(ischtp(iphas).eq.2) then
!       Pour le moment par defaut on ne prend pas l'ordre 2
!              IVSEXT(ISCAL) = 1
        ivsext(iscal) = 0
      endif
    endif
  enddo


!     Pression hydrostatique
  if (iphydr.ne.0.and.iphydr.ne.1) then
    WRITE(NFECRA,8021) 'IPHYDR ',IPHYDR
    iok = iok + 1
  endif

!     Viscosite secondaire
  do iphas = 1, nphas
    ivisph = ivisse(iphas)
    if (ivisph.ne.0.and.ivisph.ne.1) then
      WRITE(NFECRA,8022) IPHAS,'IVISSE ',IVISPH
      iok = iok + 1
    endif
  enddo

!     Schemas en temps
  do iphas = 1, nphas

!     Schema en temps global.
    if(ischtp(iphas).ne. 1.and.ischtp(iphas).ne.2) then
      WRITE(NFECRA,8101) IPHAS,'ISCHTP',ISCHTP(IPHAS)
      iok = iok + 1
    endif
    if(ischtp(iphas).eq. 2.and.idtvar.ne.0) then
      write(nfecra,8111) iphas,ischtp(iphas),idtvar
      iok = iok + 1
    endif
    if(ischtp(iphas).eq. 2.and.itytur(iphas).eq.2) then
      write(nfecra,8112) iphas,ischtp(iphas),iturb(iphas)
      iok = iok + 1
    endif
    if(ischtp(iphas).eq.1.and.itytur(iphas).eq.4) then
      write(nfecra,8113) iphas,ischtp(iphas),iturb(iphas)
    endif
    if(ischtp(iphas).eq. 2.and.iturb(iphas).eq.50) then
      write(nfecra,8114) iphas,ischtp(iphas),iturb(iphas)
      iok = iok + 1
    endif
    if(ischtp(iphas).eq. 2.and.iturb(iphas).eq.60) then
      write(nfecra,8115) iphas,ischtp(iphas),iturb(iphas)
      iok = iok + 1
    endif
    if(ischtp(iphas).eq. 2.and.iturb(iphas).eq.70) then
      write(nfecra,8116) iphas,ischtp(iphas),iturb(iphas)
      iok = iok + 1
    endif

!     Schema en temps pour le flux de masse
    if(istmpf(iphas).ne. 2.and.istmpf(iphas).ne.0.and.            &
       istmpf(iphas).ne. 1) then
      WRITE(NFECRA,8121) IPHAS,'ISTMPF',ISTMPF(IPHAS)
      iok = iok + 1
    endif

!     Schema en temps pour les termes sources de NS
    if(isno2t(iphas).ne.0.and.                                    &
       isno2t(iphas).ne. 1.and.isno2t(iphas).ne.2) then
      WRITE(NFECRA,8131) IPHAS,'ISNO2T',ISNO2T(IPHAS)
      iok = iok + 1
    endif
!     Schema en temps pour les termes sources des grandeurs
!     turbulentes
    if(isto2t(iphas).ne.0.and.                                    &
       isto2t(iphas).ne. 1.and.isto2t(iphas).ne.2) then
      WRITE(NFECRA,8131) IPHAS,'ISTO2T',ISTO2T(IPHAS)
      iok = iok + 1
    endif

!     Schema en temps pour la masse volumique
    if(iroext(iphas).ne.0.and.                                    &
       iroext(iphas).ne. 1.and.iroext(iphas).ne.2) then
      WRITE(NFECRA,8131) IPHAS,'IROEXT',IROEXT(IPHAS)
      iok = iok + 1
    endif

!     Schema en temps pour la viscosite
    if(iviext(iphas).ne.0.and.                                    &
       iviext(iphas).ne. 1.and.iviext(iphas).ne.2) then
      WRITE(NFECRA,8131) IPHAS,'IVIEXT',IVIEXT(IPHAS)
      iok = iok + 1
    endif

!     Schema en temps pour la chaleur specifique
    if(icpext(iphas).ne.0.and.                                    &
       icpext(iphas).ne. 1.and.icpext(iphas).ne.2) then
      WRITE(NFECRA,8131) IPHAS,'ICPEXT',ICPEXT(IPHAS)
      iok = iok + 1
    endif

  enddo

  do iscal = 1, nscal
!     Schema en temps pour les termes sources des scalaires
    if(isso2t(iscal).ne.0.and.                                    &
       isso2t(iscal).ne. 1.and.isso2t(iscal).ne.2) then
      WRITE(NFECRA,8141) ISCAL,'ISSO2T',ISSO2T(ISCAL)
      iok = iok + 1
    endif
!     Schema en temps pour la viscosite
    if(ivsext(iscal).ne.0.and.                                    &
       ivsext(iscal).ne. 1.and.ivsext(iscal).ne.2) then
      WRITE(NFECRA,8141) ISCAL,'IVSEXT',IVSEXT(ISCAL)
      iok = iok + 1
    endif
  enddo

!     Stop si probleme
  if(iok.gt.0) then
    call csexit(1)
  endif


! --- Reprise du dernier numero de propriete
  iprop  = nprmax

! --- Numeros de propriete
  do iphas = 1, nphas
!     On a besoin de la masse volumique si on l'extrapole ou si ICALHY
    if(iroext(iphas).gt.0.or.icalhy.eq.1) then
      iprop         = iprop + 1
      iroma (iphas) = iprop
    endif
!     Dans le cas d'une extrapolation de la viscosite totale
    if(iviext(iphas).gt.0) then
      iprop         = iprop + 1
      ivisla(iphas) = iprop
      iprop         = iprop + 1
      ivista(iphas) = iprop
    endif
  enddo
!     Proprietes des phases : CP s'il est variable
  do iphas = 1, nphas
    if(icp(iphas).ne.0) then
      if(icpext(iphas).gt.0) then
        iprop         = iprop + 1
        icpa  (iphas) = iprop
      endif
    endif
  enddo
!     On a besoin d'un tableau pour les termes sources de Navier Stokes
!       a extrapoler. Ce tableau est NDIM
  do iphas = 1, nphas
    if(isno2t(iphas).gt.0) then
      iprop         = iprop + 1
      itsnsa(iphas) = iprop
    endif
    if(isto2t(iphas).gt.0) then
      iprop         = iprop + 1
      itstua(iphas) = iprop
    endif
  enddo
!     Proprietes des scalaires : termes sources pour theta schema
!       et VISCLS si elle est variable
  if(nscal.ge.1) then
    do iscal = 1, nscal
      if(isso2t(iscal).gt.0) then
        iprop         = iprop + 1
        itssca(iscal) = iprop
      endif
      if(ivisls(iscal).ne.0) then
        if(ivsext(iscal).gt.0) then
          iprop         = iprop + 1
          ivissa(iscal) = iprop
        endif
      endif
    enddo
  endif
!     Proprietes des variables : flux de masse porteur
!       On en ajoute un (et un seul) par phase s'il existe une phase
!         qui en a besoin.
!       On est donc dans l'hypothese implicite qu'il n'y a qu'un seul
!         flux de masse pour toutes les variables
!         et ceci n'est pas optimal si il y a plusieurs phases traitees
!         differemment.
!       On suppose que si une phase en a besoin, on en ajoute donc
!         autant qu'il y a de flux de masse, cad NFLUMA

!     On les initialise a -1 pour iniva0
  do ivar = 1, nvarmx
    ifluaa(ivar) = -1
  enddo
!     On regarde s'il y en a besoin
  iiflaa = 0
  do iphas = 1, nphas
    if(istmpf(iphas).ne.1) iiflaa = 1
  enddo
!     On les affecte
  iprofa = iprop
  if(iiflaa.eq.1) then
    do iphas = 1, nphas
      if(iphas.eq.1) then
        iprop               = iprop + 1
        ifluaa(ipr (iphas)) = iprop
      endif
      ifluaa(iu  (iphas)) = iprop
      ifluaa(iv  (iphas)) = iprop
      ifluaa(iw  (iphas)) = iprop
      if(itytur(iphas).eq.2) then
        ifluaa(ik  (iphas)) = iprop
        ifluaa(iep (iphas)) = iprop
      elseif(itytur(iphas).eq.3) then
        ifluaa(ir11(iphas)) = iprop
        ifluaa(ir22(iphas)) = iprop
        ifluaa(ir33(iphas)) = iprop
        ifluaa(ir12(iphas)) = iprop
        ifluaa(ir13(iphas)) = iprop
        ifluaa(ir23(iphas)) = iprop
        ifluaa(iep (iphas)) = iprop
      elseif(iturb(iphas).eq.50) then
        ifluaa(ik  (iphas)) = iprop
        ifluaa(iep (iphas)) = iprop
        ifluaa(iphi(iphas)) = iprop
        ifluaa(ifb (iphas)) = iprop
      elseif(iturb(iphas).eq.60) then
        ifluaa(ik  (iphas)) = iprop
        ifluaa(iomg(iphas)) = iprop
      elseif (iturb(iphas).eq.70) then
        ifluaa(inusa(iphas))= iprop
      endif
    enddo
    do iscal = 1, nscal
      ifluaa(isca(iscal)) = ifluaa(iu(iphsca(iscal)))
    enddo
  endif

! --- Sauvegarde du dernier numero de propriete
  nprmax = iprop


! --- Reprise du dernier NPROCE et du dernier NPPMAX
  iprop                 = nproce
  ipppst                = nppmax

! --- Positionnement des PROPCE
  do iphas = 1, nphas

!     Variables schema en temps
    if(iroext(iphas).gt.0.or.icalhy.eq.1) then
      iprop                 = iprop  + 1
      ipproc(iroma (iphas)) = iprop
    endif
    if(iviext(iphas).gt.0) then
      iprop                 = iprop  + 1
      ipproc(ivisla(iphas)) = iprop
    endif
    if(iviext(iphas).gt.0) then
      iprop                 = iprop  + 1
      ipproc(ivista(iphas)) = iprop
    endif
    if(icpext(iphas).gt.0) then
      iprop                 = iprop + 1
      ipproc(icpa  (iphas)) = iprop
    endif
    if(isno2t(iphas).gt.0) then
      iprop                 = iprop + 1
      ipproc(itsnsa(iphas)) = iprop
!     Ce tableau est NDIM :
      iprop                 = iprop + ndim-1
    endif
    if(isto2t(iphas).gt.0) then
      iprop                 = iprop + 1
      ipproc(itstua(iphas)) = iprop
!     Ce tableau est 2, 7 ou 4 selon le modele de turbulence :
      if    (itytur(iphas).eq.2) then
        iprop                 = iprop + 2-1
      elseif(itytur(iphas).eq.3) then
        iprop                 = iprop + 7-1
      elseif(iturb(iphas).eq.50) then
        iprop                 = iprop + 4-1
      elseif(iturb(iphas).eq.70) then
        iprop                 = iprop + 1-1
      endif
    endif
  enddo

  do ii = 1, nscal
! Termes source des scalaires pour theta schema
    if(isso2t(ii).gt.0) then
      iprop                 = iprop + 1
      ipproc(itssca(ii))    = iprop
    endif

    if(ivisls(ii).gt.0) then
      if(iscavr(ii).le.0) then
        if(ivsext(ii).gt.0) then
          iprop                 = iprop + 1
          ipproc(ivissa(ii)  )  = iprop
        endif
      endif
    endif
  enddo
  do ii = 1, nscal
    if(ivisls(ii).gt.0) then
      if(iscavr(ii).gt.0.and.iscavr(ii).le.nscal) then
        if(ivsext(ii).gt.0) then
          ipproc(ivissa(ii)  )  = ipproc(ivissa(iscavr(ii))  )
        endif
      endif
    endif
  enddo

! --- Sauvegarde du dernier NPROCE et du dernier NPPMAX
  nproce = iprop
  nppmax = ipppst


! --- Reprise du dernier NPROFB
  iprop                 = nprofb

! --- Positionnement des PROPFB
  do iphas = 1, nphas
!     Variables schema en temps : rhoa (pas pour icalhy)
    if(iroext(iphas).gt.0) then
      iprop                 = iprop  + 1
      ipprob(iroma (iphas)) = iprop
    endif
  enddo
!     Variables schema en temps : flux de masse A
  if(iiflaa.eq.1) then
    do iflum = 1, nfluma
      iprop                 = iprop + 1
      ipprob(iprofa+iflum)  = iprop
    enddo
  endif

! --- Sauvegarde du dernier NPROFB
  nprofb = iprop


! --- Reprise du dernier NPROFA
  iprop                 = nprofa

! --- Positionnement des PROPFA
  if(iiflaa.eq.1) then
    do iflum = 1, nfluma
      iprop                 = iprop + 1
      ipprof(iprofa+iflum)  = iprop
    enddo
  endif

! --- Sauvegarde du dernier NPROFA
  nprofa = iprop


! ---> 4.2 CALCUL DE LA TAILLE DU TABLEAU DES TEMPS CUMULES POUR LES MOMENTS
!      ---------------------------------------------------------------------

!     Pour verification des definitions de moments
  iok = 0

! --- Calcul du nombre de moments definis (et verif qu'il n y a pas de trous)
  nbmomt = 0
  inmfin = 0
  do imom = 1, nbmomx
!     Si on n'est pas a la fin de la liste
    if(inmfin.eq.0) then
!       Si il y en a, ca en fait en plus
      if(idfmom(1,imom).ne.0) then
        nbmomt = nbmomt + 1
!       Si il n'y en a pas, c'est la fin de la liste
      else
        inmfin = 1
      endif
!     Si on est a la fin de la liste, il n'en faut plus
    else
      if(idfmom(1,imom).ne.0) then
        iok = iok + 1
      endif
    endif
  enddo

  if(iok.ne.0) then
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
      if(idffin.eq.0) then
        if(idfmji.lt.-nprmax) then
          iok = iok + 1
          write(nfecra,8210)jj,imom,idfmji,nprmax
        elseif(idfmji.gt.nvar) then
          iok = iok + 1
          write(nfecra,8211)jj,imom,idfmji,nvar
        elseif(idfmji.lt.0) then
          if (ipproc(-idfmji).le.0) then
            iok = iok + 1
            write(nfecra,8212)jj,imom,idfmji,-idfmji,             &
                 ipproc(-idfmji)
          endif
        elseif(idfmji.eq.0) then
          idffin = 1
        endif
      else
        if(idfmji.ne.0) then
          iok = iok + 1
          write(nfecra,8213)imom,jj,idfmji
        endif
      endif
    enddo
  enddo

! --- Verification de NTDMOM (>0)
  do imom = 1, nbmomt
    if(ntdmom(imom).lt.0) then
      iok = iok + 1
      write(nfecra,8214)imom,ntdmom(imom)
    endif
  enddo

  if(iok.ne.0) then
    call csexit(1)
  endif


! --- Indicateur pour juger de l'utilite de consulter le fichier suite
!     Il n'est pas indispensable de consulter avec succes le fichier suite
!       lorsqu'il n'y a aucun moment a relire ou qu'on ne fait pas de suite
!       (si ILSMOM=0)
!       ie ILSMOM = 0 == (ISUITE=0 OR IMOOLD(IMOM)=-1 pour tout IMOM)
!       ie ILSMOM = 1 == (ISUITE=1 AND IMOOLD(IMOM) different de -1 pour un IMOM)
  ilsmom = 0
  if(isuite.eq.1.and.nbmomt.gt.0) then
    do imom = 1, nbmomt
      if(imoold(imom).ne.-1) then
        ilsmom = 1
      endif
    enddo
  endif

! --- Lecture du fichier suite (debut : info sur les moments et sur le ntpabs)
  if(ilsmom.eq.1) then

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
    RUBRIQ = 'version_fichier_suite_auxiliaire'
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
    RUBRIQ = 'nombre_moyennes_temps'
    jbmomt = 0
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                jbmomt,ierror)
    nberro=nberro+ierror
    if(ierror.ne.0) then
      write(nfecra,8311)
    endif

!     Nombre de pas de temps final du calcul precedent
    RUBRIQ = 'nbre_pas_de_temps'
    itysup = 0
    nbval  = 1
    irtyp  = 1
    jtcabs = 0
    call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,irtyp,     &
                jtcabs,ierror)
    nberro=nberro+ierror
    if(ierror.ne.0) then
      write(nfecra,8312)
    endif

    if(nberro.ne.0) then
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

  if(ilsmom.eq.1) then
    do imom = 1, nbmomt
      if(imoold(imom).gt.jbmomt.or.imoold(imom).eq.0              &
                             .or.imoold(imom).lt.-2) then
        write(nfecra,8400) imom,imoold(imom),jbmomt
        iok = iok + 1
      elseif(imoold(imom).le.jbmomt.and.imoold(imom).gt.0         &
                             .and.ntdmom(imom).gt.jtcabs) then
        write(nfecra,8401) imom,imoold(imom),ntdmom(imom),jtcabs
        iok = iok + 1
      elseif(imoold(imom).eq.-2.and.ntdmom(imom).gt.jtcabs) then
        imoold(imom)=-1
      elseif(imoold(imom).eq.-2.and.imom.le.jbmomt) then
        imoold(imom)=imom
      elseif(imoold(imom).eq.-2.and.imom.gt.jbmomt) then
        imoold(imom)=-1
      endif
    enddo
    do imom = 1, nbmomt
      if(imoold(imom).gt.0) then
        do jmom = 1, imom-1
          if(imoold(imom).eq.imoold(jmom)) then
            write(nfecra,8402) imom,jmom,imoold(imom)
            iok = iok + 1
          endif
        enddo
      endif
    enddo

!       Si pas suite : si IMOOLD pas renseigne, on initialise
!                      sinon erreur
  elseif(isuite.eq.0) then
    do imom = 1, nbmomt
      if(imoold(imom).ne.-2) then
        write(nfecra,8403) isuite,imom,imoold(imom)
        iok = iok + 1
      else
        imoold(imom)=-1
      endif
    enddo
  endif

  if(iok.ne.0) then
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
  if(ilsmom.eq.1) then

    nberro = 0

!     On ne lit des choses que s'il y a des moments stockes
    if ( jbmomt.gt.0 ) then

!        On lit les infos pour tous les moments courants
      do imom = 1, nbmomt
!           On s'interesse uniquement aux moments qui suivnt d'anciens moments
        if(imoold(imom).gt.0) then
!            Si le numero de l'ancien moment est incompatible avec le format
!              on cherche numero_cumul_temps_momentYYYY qui n existe pas
!              ca genere une erreur
          if(imom.le.nfmtmo) then
            WRITE(CMOY4,'(I4.4)')IMOOLD(IMOM)
          else
            cmoy4 = cindfm
          endif
          itysup = 0
          nbval  = 1
          irtyp  = 1
          RUBRIQ = 'numero_cumul_temps_moment'//CMOY4
          call lecsui(impamx,rubriq,len(rubriq),itysup,nbval,     &
                      irtyp,idtold(imom),ierror)
          nberro=nberro+ierror
          if(idtold(imom).eq.0.or.ierror.ne.0) then
            write(nfecra,8313)imom,imoold(imom)
            iok = iok + 1
          endif
        endif
      enddo
    endif

    if(nberro.ne.0) then
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
    if(imold.eq.-1) then
!          On cherche si on en a deja vu (re)initialisees au meme moment
      imomr = 0
      do jmom = 1, imom-1
        jmold = imoold(jmom)
        if(jmold.eq.-1.and.ntdmom(jmom).eq.ntdmom(imom)) then
          imomr = jmom
        endif
      enddo
!          Si oui : on utilise le meme tableau
      if(imomr.gt.0) then
        idtmom(imom) = idtmom(imomr)
!          Si non : on a besoin d'un nouveau tableau ou reel
      else
        if(idtvar.eq.2.or.idtvar.eq.-1) then
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
  if(isuite.eq.1) then
    do imom = 1, nbmomt
      imold = imoold(imom)
      if(imold.gt.0) then
!         On regarde si le DTcumule du IMOLD a deja ete vu
!            (indicateur JMOMOK)
        idto = idtold(imom)
        jmomok = 0
        do jmom = 1, imom-1
          jmold = imoold(jmom)
          if(jmold.gt.0) then
            jdto = idtold(jmom)
            if(jdto.eq.idto) then
              jmomok = jmom
            endif
          endif
        enddo
!         Si on le voit pour la premiere fois, ca en fait un de plus
        if(jmomok.eq.0) then
!           Si DT non uniforme dans le present calcul ou que
!              DT cumule etait non uniforme dans l'ancien, on a un
!              DT cumule non uniforme (sinon, il est uniforme)
          if(idtvar.eq.2.or.idtvar.eq.-1.or.idto.gt.0) then
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
    if(idtmom(imom).eq.0) then
      iok = iok + 1
      write(nfecra,8410)imom, idtmom(imom)
    endif
  enddo
  if(iok.ne.0) then
    call csexit (1)
  endif

! --- Calcul du nombre de tableaux NCEL "temps cumule"
  nbdtcm = 0
  do imom = 1, nbmomt
    nbdtcm=max(idtmom(imom),nbdtcm)
  enddo

! ---> 4.3 POSITIONNEMENT DANS PROPCE DES MOMENTS ET DU TEMPS CUMULE
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



! ---> 4.4  POSITIONNEMENT DES CONDITIONS AUX LIMITES
!      ---------------------------------------------------------------------

! --- Numerotation des tableaux NFABOR de type COEFA/COEFB presents ici
!       On suppose pour le moment que seules les variables de calcul
!         disposent de conditions aux limites. Ceci pourra ensuite etre
!         etendu aux variables pysiques en declarant un pointeur du type
!         de ICLRTP (par exemple ICLPRO)
!       On suppose que les seules variables qui ont 2 types de cl sont
!         les variables vitesse en k-epsilon, k-omega et en LES, et la
!         pression si IPHYDR=1

  icondl = 0
  do ivar = 1, nvar
    icondl = icondl + 1
    iclrtp(ivar,icoef ) = icondl
    iclrtp(ivar,icoeff) = icondl
  enddo
  do iphas = 1, nphas
    if( itytur(iphas).eq.2 .or. itytur(iphas).eq.4                &
         .or. iturb(iphas).eq.60 .or. iturb(iphas).eq.70          &
         ) then
      ivar = iu(iphas)
      icondl = icondl + 1
      iclrtp(ivar,icoeff) = icondl
      ivar = iv(iphas)
      icondl = icondl + 1
      iclrtp(ivar,icoeff) = icondl
      ivar = iw(iphas)
      icondl = icondl + 1
      iclrtp(ivar,icoeff) = icondl
    endif
  enddo
  if (iphydr.eq.1) then
    do iphas = 1, nphas
      icondl = icondl + 1
      iclrtp(ipr(iphas),icoeff) = icondl
    enddo
  endif

! Compressible
  if (ippmod(icompf).ge.0) then
    do iphas = 1, nphas
      icondl = icondl + 1
      iclrtp(isca(ienerg(iphas)),icoeff) = icondl
    enddo
  endif

  ncofab = icondl


! ---> 4.5 POINTEURS POST-PROCESSING / LISTING / HISTORIQUES / CHRONOS
!      ---------------------------------------------------------------------

! --- Les pointeurs ont ete initialises a 1 (poubelle).
!       dans iniini (avant usini1).

!     On posttraitera les variables localisees au centre des cellules.

!      IPPRTP(IVAR)  pour RTP a ete complete plus haut.

!      IPPPRO(IPPROC(II)) pour PROPCE a ete complete plus haut
!        au fur et a mesure (voir ppprop en particulier)

!      Le rang de la derniere propriete pour le post est IPPPST.


  if(idtvar.gt.0) then
    ipppst = ipppst + 1
    ippdt  = ipppst
    nppmax = ipppst
  endif

  if(ipucou.eq.1) then
    ipppst = ipppst + 1
    ipptx  = ipppst
    ipppst = ipppst + 1
    ippty  = ipppst
    ipppst = ipppst + 1
    ipptz  = ipppst
    nppmax = ipppst
  endif


! Verification de la limite sur IPPPST

  if(ipppst.gt.nvppmx) then
    write(nfecra,8900)ipppst,nvppmx
    call csexit (1)
    !==========
  endif

  return

endif

!===============================================================================
! 5. CINQUIEME APPEL :
!        RESERVATION D'UNE PLACE DANS PROPCE SI RAYONNEMENT
!        ET LAGRANGIEN AVEC THERMIQUE DES PARTICULES
!===============================================================================

if (ipass.eq.5) then

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
      ipproc(icak(irphas))  = iprop
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
    itparo       = iprop
    iprop        = iprop + 1
    iqinci       = iprop
    iprop        = iprop + 1
    ixlam        = iprop
    iprop        = iprop + 1
    iepa         = iprop
    iprop        = iprop + 1
    ieps         = iprop
    iprop        = iprop + 1
    ifnet        = iprop
    iprop        = iprop + 1
    ifconv       = iprop
    iprop        = iprop + 1
    ihconv       = iprop

! --- Sauvegarde du dernier numero de propriete
    nprmax = iprop

! --- Reprise du dernier NPROFB (PROPFB)
    iprop   = nprofb

! --- Positionnement
    iprop          = iprop + 1
    ipprob(itparo) = iprop

    iprop          = iprop + 1
    ipprob(iqinci) = iprop

    iprop          = iprop + 1
    ipprob(ixlam)  = iprop

    iprop          = iprop + 1
    ipprob(iepa)   = iprop

    iprop          = iprop + 1
    ipprob(ieps)   = iprop

    iprop          = iprop + 1
    ipprob(ifnet)  = iprop

    iprop          = iprop + 1
    ipprob(ifconv) = iprop

    iprop          = iprop + 1
    ipprob(ihconv) = iprop

    nprayb = iprop - nprofb

    if (iihmpr.eq.1) then

      call uirapr &
      !==========
    ( nprayc, nprayb, nrphas, ipppro, ipproc,           &
      ilumin, iqx, iqy, iqz,                            &
      itsre, itsri, iabs, iemi, icak)

    endif

!
! --- Sauvegarde du dernier NPROFB
    nprofb = iprop

! --- Verification de NPROCE, NPROFA, NPROFB

    if(nproce.gt.npromx.or.                                       &
       nprofa.gt.npromx.or.nprofb.gt.npromx) then
      write(nfecra,7200)nproce, nprofa, nprofb, npromx,           &
           max(max(nproce,nprofa),nprofb)
      call csexit (1)
      !==========
    endif

  endif

  return

endif


!===============================================================================
! 6. FORMATS
!===============================================================================

#if defined(_CS_LANG_FR)
 5000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE PHASES ERRONE                                ',/,&
'@                                                            ',/,&
'@  Le nombre de phases doit etre un entier strictement       ',/,&
'@    positif. Il vaut ici                NPHAS  = ',I10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE PHASES TROP GRAND                            ',/,&
'@                                                            ',/,&
'@  Le nombre de phases                                       ',/,&
'@    - demande          dans usini1   est NPHAS  = ',I10      ,/,&
'@    - maximal autorise dans paramx.h est NPHSMX = ',I10      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@  NPHSMX doit valoir au moins ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     PLUSIEURS MODELES PHYSIQUES PARTICULIERES ACTIVES      ',/,&
'@                                                            ',/,&
'@  Un seul modele physique particuliere peut etre active a la',/,&
'@    fois.                                                   ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les indicateurs de IPPMOD dans usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     SELECTION INCORRECTE DU MODELE PHYSIQUE PARTICULIERE   ',/,&
'@                                                            ',/,&
'@  Les valeurs des indicateurs du tableau IPPMOD ne sont pas ',/,&
'@    admissibles                                             ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier les indicateurs de IPPMOD dans usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE SCALAIRES ERRONE                             ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires utilisateur doit etre un entier    ',/,&
'@    positif ou nul. Il vaut ici   NSCAUS  = ',I10            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE SCALAIRES TROP GRAND                         ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires utilisateurs                       ',/,&
'@    demande    dans usini1           est  NSCAUS = ',I10     ,/,&
'@  Le nombre de scalaires total                              ',/,&
'@    autorise   dans paramx.h         est  NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  La valeur maximale autorisee de NSCAUS                    ',/,&
'@                          est donc  NSCAMX        = ',I10    ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier NSCAUS dans usini1.                              ',/,&
'@                                                            ',/,&
'@  NSCAMX doit valoir au moins ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6012 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE SCALAIRES TROP GRAND                         ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires utilisateurs                       ',/,&
'@    demande    dans usini1           est  NSCAUS = ',I10     ,/,&
'@  Le nombre de scalaires pour les physiques particulieres   ',/,&
'@    necessaire avec le modele choisi est  NSCAPP = ',I10     ,/,&
'@  Le nombre de scalaires total                              ',/,&
'@    autorise   dans paramx.h         est  NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  La valeur maximale autorisee de NSCAUS                    ',/,&
'@    avec le modele choisi est donc NSCAMX-NSCAPP = ',I10     ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier NSCAUS dans usini1.                              ',/,&
'@                                                            ',/,&
'@  NSCAMX doit valoir au moins ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE SCALAIRES ERRONE                             ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires doit etre un entier                ',/,&
'@    positif ou nul. Il vaut ici   NSCAL   = ',I10            ,/,&
'@    Remarque : NSCAUS = ',I10                                ,/,&
'@               NSCAPP = ',I10                                ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE SCALAIRES TROP GRAND                         ',/,&
'@                                                            ',/,&
'@  Le nombre de scalaires necessaire  est  NSCAL  = ',I10     ,/,&
'@  Le nombre de scalaires total                              ',/,&
'@    autorise   dans paramx.h         est  NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@    Remarque : NSCAUS = ',I10                                ,/,&
'@               NSCAPP = ',I10                                ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX doit valoir au moins ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHASE PORTEUSE INCORRECTE POUR LE SCALAIRE ',I10         ,/,&
'@                                                            ',/,&
'@  La phase porteuse du scalaire ',I10                        ,/,&
'@          (scalaire utilisateur ',I10   ,') indiquee dans   ',/,&
'@    usini1 est IPHSCA(',I10   ,') = ',I10                    ,/,&
'@  Elle devrait etre comprise entre 0 et NPHAS = ',I10        ,/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IPHSCA dans usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7021 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHASE PORTEUSE INCORRECTE POUR LE SCALAIRE ',I10         ,/,&
'@                                                            ',/,&
'@  La phase porteuse du scalaire ',I10                        ,/,&
'@          (scalaire physique particuliere ',I10   ,')       ',/,&
'@           est IPHSCA(ISCAPP(',I10   ,')) = ',I10            ,/,&
'@  Elle devrait etre comprise entre 0 et NPHAS = ',I10        ,/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IPHSCA.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7030 format(                                                           &
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
'@  Verifier ISCAVR dans usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7031 format(                                                           &
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
'@  Verifier ISCAVR dans usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7040 format(                                                           &
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
'@  Verifier IVISLS dans usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7041 format(                                                           &
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
'@  Verifier IVISLS dans usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7042 format(                                                           &
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
 7043 format(                                                           &
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
 7050 format(                                                           &
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
'@  Verifier IVISLS dans usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7051 format(                                                           &
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
 7060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA PHASE PORTEUSE   ',/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire utilisateur           ',I10   ,') represente  ',/,&
'@    la variance des fluctuations du scalaire ',I10           ,/,&
'@    (scalaire utilisateur           ',I10   ,') puisque     ',/,&
'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,&
'@                                                            ',/,&
'@  La phase porteuse IPHSCA(',I10   ,') du scalaire          ',/,&
'@    ne doit pas etre renseignee.                            ',/,&
'@  Elle sera automatiquement prise identique a la phase      ',/,&
'@    porteuse du scalaire associe, soit ',I10                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IPHSCA dans usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7061 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA PHASE PORTEUSE   ',/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire utilisateur           ',I10   ,') represente  ',/,&
'@    la variance des fluctuations du scalaire ',I10           ,/,&
'@    (scalaire physique particuliere ',I10   ,') puisque     ',/,&
'@    ISCAVR(',I10   ,') vaut ',I10   ,' (non nul)            ',/,&
'@                                                            ',/,&
'@  La phase porteuse IPHSCA(',I10   ,') du scalaire          ',/,&
'@    ne doit pas etre renseignee.                            ',/,&
'@  Elle sera automatiquement prise identique a la phase      ',/,&
'@    porteuse du scalaire associe, soit ',I10                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IPHSCA dans usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7062 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA PHASE PORTEUSE   ',/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire physique particuliere ',I10   ,') represente  ',/,&
'@    la variance des fluctuations du scalaire ',I10           ,/,&
'@    (scalaire utilisateur           ',I10   ,') puisque     ',/,&
'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,&
'@                                                            ',/,&
'@  La phase porteuse IPHSCA(ISCAPP(',I10   ,')) du scalaire  ',/,&
'@    ne doit pas etre renseignee.                            ',/,&
'@  Elle sera automatiquement prise identique a la phase      ',/,&
'@    porteuse du scalaire associe, soit ',I10                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IPHSCA.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7063 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',I10   ,' NE PAS MODIFIER LA PHASE PORTEUSE   ',/,&
'@                                                            ',/,&
'@  Le scalaire ',I10                                          ,/,&
'@    (scalaire physique particuliere ',I10   ,') represente  ',/,&
'@    la variance des fluctuations du scalaire ',I10           ,/,&
'@    (scalaire physique particuliere ',I10   ,') puisque     ',/,&
'@    ISCAVR(ISCAPP(',I10   ,')) vaut ',I10   ,' (non nul)    ',/,&
'@                                                            ',/,&
'@  La phase porteuse IPHSCA(ISCAPP(',I10   ,')) du scalaire  ',/,&
'@    ne doit pas etre renseignee.                            ',/,&
'@  Elle sera automatiquement prise identique a la phase      ',/,&
'@    porteuse du scalaire associe, soit ',I10                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier IPHSCA.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7070 format(                                                           &
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
 7100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE VARIABLES TROP GRAND                         ',/,&
'@                                                            ',/,&
'@  Le type de calcul defini dans usini1                      ',/,&
'@    correspond a un nombre de variables NVAR   = ',I10       ,/,&
'@  Le nombre de variables maximal prevu                      ',/,&
'@                      dans paramx.h est NVARMX = ',I10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@  NVARMX doit valoir au moins ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE PROPRIETES TROP GRAND                        ',/,&
'@                                                            ',/,&
'@  Le type de calcul defini dans usini1                      ',/,&
'@    correspond aux nombres de proprietes suivants           ',/,&
'@      au centre des cellules       : NPROCE = ',I10          ,/,&
'@      au centre des faces internes : NPROFA = ',I10          ,/,&
'@      au centre des faces de bord  : NPROFB = ',I10          ,/,&
'@  Le nombre de proprietes maximal prevu                     ',/,&
'@                      dans paramx.h est NPROMX = ',I10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@  NPROMX doit valoir au moins ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8021 format(                                                           &
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
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8022 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHASE ',I10                                              ,/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0 OU 1                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8101 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PHASE ',I10                                              ,/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 1 ou 2                ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8111 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    AVEC UN SCHEMA EN TEMPS D ORDRE 2 : ISCHTP = ', I10      ,/,&
'@    IL FAUT UTILISER UN PAS DE TEMPS CONSTANT ET UNIFORME   ',/,&
'@    OR IDTVAR = ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8112 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN K-EPSILON (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l ordre 2 avec le    ',/,&
'@   couplage des termes sources du k-epsilon.                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8113 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 1 (ISCHTP = ',I10   ,/,&
'@    EN LES (ITURB = ',I10,' )'                               ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier usini1.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8114 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN V2F       (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources du k-epsilon.                ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8115 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN K-OMEGA   (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources du k-omega.                  ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8116 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    ON IMPOSE UN SCHEMA EN TEMPS D ORDRE 2 (ISCHTP = ',I10   ,/,&
'@    EN SPALART   (ITURB = ',I10,' )'                         ,/,&
'@                                                            ',/,&
'@   La version courante ne supporte pas l''ordre 2 avec le   ',/,&
'@   couplage des termes sources de Spalart-Allmaras.         ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Modifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8121 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A  0, 1 OU 2            ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8131 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2             ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8132 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  Avec le modele de turbulence choisi, ITURB = ',I10         ,/,&
'@    la valeur de ISTO2T (extrapolation des termes sources   ',/,&
'@    pour les variables turbulentes) ne doit pas etre modifie',/,&
'@    or ISTO2T a ete force a ',I10                            ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8141 format(                                                           &
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
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8200 format(                                                           &
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
 8201 format(                                                           &
'@  ',I10   ,'        ',     I10                                 )
 8202 format(                                                           &
'@  ----------------------------                              ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas  execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8210 format(                                                           &
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
'@    Verifier usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8211 format(                                                           &
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
'@    Verifier usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8212 format(                                                           &
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
'@    Verifier usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8213 format(                                                           &
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
'@    Verifier usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8214 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA VERIFICATION DES DONNEES         ',/,&
'@    =========                                               ',/,&
'@  SUR L''INSTANT DE DEBUT DE CALCUL DES MOYENNES TEMPORELLES',/,&
'@                                                            ',/,&
'@    La variable NTDMOM(IMOM)   pour IMOM = ',I10             ,/,&
'@      doit etre renseignee dans usini1 pour indiquer        ',/,&
'@      a partir de quel pas de temps (absolu) doit etre      ',/,&
'@      calculee la moyenne temporelle IMOM correspondante.   ',/,&
'@      NTDMOM(IMOM) vaur ici ',I10                            ,/,&
'@                                                            ',/,&
'@    Le calcul ne peut etre execute.                         ',/,&
'@                                                            ',/,&
'@    Verifier usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8300 format(                                                           &
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
'@      Voir alors usini1.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8301 format(                                                           &
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
'@      Voir alors usini1.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8311 format(                                                           &
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
'@      Voir alors usini1.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8312 format(                                                           &
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
'@      Voir alors usini1.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8313 format(                                                           &
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
'@      Voir alors usini1.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8390 format(                                                           &
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
 8400 format(                                                           &
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
'@    Verifier les valeurs de IMOOLD dans usini1.             ',/,&
'@    Verifier que le fichier suite est le bon.               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8401 format(                                                           &
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
'@    Modifier la valeur de  NTDMOM dans usini1               ',/,&
'@      verifier que le fichier suite est le bon.             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8402 format(                                                           &
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
'@    Verifier IMOOLD dans usini1      .                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8403 format(                                                           &
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
'@    Ne pas modifier IMOOLD dans usini1.                     ',/,&
'@      ou realiser un calcul suite.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8410 format(                                                           &
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
 8900 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     NOMBRE DE VARIABLES A SUIVRE TROP GRAND                ',/,&
'@                                                            ',/,&
'@  Le type de calcul defini dans usini1                      ',/,&
'@    correspond a un nombre de variables a suivre dans       ',/,&
'@    le listing et le post-processing egal a      ',I10       ,/,&
'@  Le nombre de variables a suivre maximal prevu             ',/,&
'@                      dans paramx.h est NVPPMX = ',I10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@  Contacter l assistance.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 5000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA                    ',/,&
'@    =========                                               ',/,&
'@     WRTONG    NUMBER OF PHASES                             ',/,&
'@                                                            ',/,&
'@  The number of phases must be an integer strictly          ',/,&
'@    positive. Here it has a value of    NPHAS  = ',I10       ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF PHASES TOO LARGE                             ',/,&
'@                                                            ',/,&
'@  The number of phases                                      ',/,&
'@    - requested        in   usini1   is  NPHAS  = ',I10      ,/,&
'@    - maximmum authorised in paramx.h is NPHSMX = ',I10      ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@  NPHSMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     TOO MANY SPECIFIC PHYSICS MODULES ACTIVATED            ',/,&
'@                                                            ',/,&
'@  Only one specific physics module can be active for one    ',/,&
'@    given calculation.                                      ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Modify the indices of       IPPMOD in   usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     WRONG SELLECTION OF THE MODEL FOR SPECIFIC PHYSICS     ',/,&
'@                                                            ',/,&
'@  The values of the indices of the array IPPMOD are not     ',/,&
'@    admissible                                              ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Modify the indices of       IPPMOD in   usppmo.           ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     ERRONEOUS NUMBER OF SCALARS                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars must be an integer either     ',/,&
'@   positive or zero. Here is      NSCAUS  = ',I10            ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@  requested in      usini1           is   NSCAUS = ',I10     ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx.h         is   NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  The maximmum value allowed of   NSCAUS                    ',/,&
'@                          is in   NSCAMX        = ',I10      ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS in   usini1.                              ',/,&
'@                                                            ',/,&
'@  NSCAMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6012 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The number of users scalars                               ',/,&
'@     requested in   usini1           is   NSCAUS = ',I10     ,/,&
'@  The number of scalars necessary for the specific physics'  ,/,&
'@    with the chosen model is              NSCAPP = ',I10     ,/,&
'@  The total number of scalars                               ',/,&
'@    allowed    in   paramx.h         est  NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@  The maximum value allowed for  NSCAUS                     ',/,&
'@    with the chosen model is       NSCAMX-NSCAPP = ',I10     ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   NSCAUS in   usini1.                              ',/,&
'@                                                            ',/,&
'@  NSCAMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     WRONG NUMBER OF SCALARS                                ',/,&
'@                                                            ',/,&
'@  The number of scalars must be an integer either           ',/,&
'@    positive or zero. Here it is  NSCAL   = ',I10            ,/,&
'@    Note     : NSCAUS = ',I10                                ,/,&
'@               NSCAPP = ',I10                                ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF SCALARS TOO LARGE                            ',/,&
'@                                                            ',/,&
'@  The necessary number of scalars is      NSCAL  = ',I10     ,/,&
'@  The number of scalars available in                        ',/,&
'@                    paramx.h           is NSCAMX = ',I10     ,/,&
'@                                                            ',/,&
'@    Note     : NSCAUS = ',I10                                ,/,&
'@               NSCAPP = ',I10                                ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@  NSCAMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    CARRIER PHASE IS INCORRECT FOR THE SCALAR  ',I10         ,/,&
'@                                                            ',/,&
'@  The carrier phase of the scalar ' ,I10                     ,/,&
'@    (user scalar          ',I10   ,') indicated in          ',/,&
'@    usini1 is  IPHSCA(',I10   ,') = ',I10                    ,/,&
'@    It should be between zero and       NPHAS = ',I10        ,/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IPHSCA in   usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7021 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    CARRIER PHASE IS INCORRECT FOR THE SCALAR  ',I10         ,/,&
'@                                                            ',/,&
'@  The carrier phase of the scalar ' ,I10                     ,/,&
'@          (scalar in paricular phisics    ',I10   ,')       ',/,&
'@           is  IPHSCA(ISCAPP(',I10   ,')) = ',I10            ,/,&
'@  It should be between zero and    NPHAS = ',I10             ,/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IPHSCA.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7030 format(                                                           &
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
'@  Verify   ISCAVR in   usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7031 format(                                                           &
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
'@  Verify   ISCAVR in   usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7040 format(                                                           &
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
'@  Verify   IVISLS in   usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7041 format(                                                           &
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
'@  Verify   IVISLS in   usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7042 format(                                                           &
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
 7043 format(                                                           &
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
 7050 format(                                                           &
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
'@  Verify   IVISLS in   usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7051 format(                                                           &
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
 7060 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10   ,' DO NOT MODIFY THE CARRIER PHASE     ',/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (user scalar                    ',I10   ,') represents  ',/,&
'@    the variance of fluctuations of the scalar',I10          ,/,&
'@    (user scalar                    ',I10   ,') since       ',/,&
'@    ISCAVR(',I10   ,') has a value ',I10   ,' (non-zero)    ',/,&
'@                                                            ',/,&
'@  The carrier phase IPHSCA(',I10   ,') of the scalar        ',/,&
'@    must not be set.                                        ',/,&
'@  It will automatically be set equal to the carrier phase   ',/,&
'@    of the associated scalar ',I10                           ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IPHSCA in   usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7061 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10   ,' DO NOT MODIFY THE CARRIER PHASE     ',/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (user scalar                    ',I10   ,') represents  ',/,&
'@    the variance of fluctuations of the scalar',I10          ,/,&
'@    (scalar of specific physics   ',I10   ,') since         ',/,&
'@    ISCAVR(',I10   ,')has a value ',I10   ,' (non-zero)     ',/,&
'@                                                            ',/,&
'@  The carrier phase IPHSCA(',I10   ,') of the scalar        ',/,&
'@    must not be set.                                        ',/,&
'@  It will automatically be set equal to the carrier phase   ',/,&
'@    of the associated scalar ',I10                           ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IPHSCA in   usini1.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7062 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10   ,' DO NOT MODIFY THE CARRIER PHASE     ',/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (scalar of specific physics   ',I10   ,') represents    ',/,&
'@    the variance of fluctuations of the scalar',I10          ,/,&
'@    (user scalar                    ',I10   ,') since       ',/,&
'@ ISCAVR(ISCAPP(',I10   ,'))has a value ',I10   ,' (non-zero)',/,&
'@                                                            ',/,&
'@  The carrier phase IPHSCA(ISCAPP(',I10   ,')) of the scalar',/,&
'@    must not be set.                                        ',/,&
'@  It will automatically be set equal to the carrier phase   ',/,&
'@    of the associated scalar ',I10                           ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IPHSCA.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7063 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',I10   ,' DO NOT MODIFY THE CARRIER PHASE     ',/,&
'@                                                            ',/,&
'@  The scalar  ',I10                                          ,/,&
'@    (scalar of specific physics   ',I10   ,') represents    ',/,&
'@    the variance of fluctuations of the scalar',I10          ,/,&
'@    (scalar of specific physics   ',I10   ,') since         ',/,&
'@ ISCAVR(ISCAPP(',I10   ,'))has a value ',I10   ,' (non-zero)',/,&
'@                                                            ',/,&
'@  The carrier phase IPHSCA(ISCAPP(',I10   ,')) of the scalar',/,&
'@    must not be set.                                        ',/,&
'@  It will automatically be set equal to the carrier phase   ',/,&
'@    of the associated scalar ',I10                           ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   IPHSCA.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7070 format(                                                           &
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
 7100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF VARIABLES TOO LARGE                          ',/,&
'@                                                            ',/,&
'@  The type of calculation defined in usini1                 ',/,&
'@    corresponds to a number of variables NVAR  = ',I10       ,/,&
'@  The maximum number of variables allowed                   ',/,&
'@                      in   paramx.h is  NVARMX = ',I10       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@  NVARMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NUMBER OF VARIABLES TOO LARGE                          ',/,&
'@                                                            ',/,&
'@  The type of calculation defined in usini1                 ',/,&
'@    corresponds  to the following number of properties      ',/,&
'@      at the cell centres          : NPROCE = ',I10          ,/,&
'@      at the internal face centres : NPROFA = ',I10          ,/,&
'@      at the boundary face centres : NPROFB = ',I10          ,/,&
'@  The maxumum number of properties allowed                  ',/,&
'@                      in   paramx.h is  NPROMX = ',I10       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@  NPROMX must be at least     ',I10                          ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8021 format(                                                           &
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
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8022 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    PHASE ',I10                                              ,/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0 OR 1               ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8101 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@    PHASE ',I10                                              ,/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 1 OR 2               ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8111 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA FOR PHASE'  ,I10    ,/,&
'@    =========                                               ',/,&
'@    WITH A SECOND ORDER SCHEME IN TIME: ISCHTP = ', I10      ,/,&
'@    IT IS NECESSARY TO USE A CONSTANT AND UNIFORM TIME STEP ',/,&
'@    BUT IDTVAR = ',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verifier usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8112 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA FOR PHASE'  ,I10    ,/,&
'@    =========                                               ',/,&
'@    A 2ND ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    WITH K-EPSILON (ITURB = ',I10,' )'                       ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8113 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   :     AT THE INITIAL DATA FOR PHASE   ',I10    ,/,&
'@    =========                                               ',/,&
'@    A 1st ORDER SCHEME HAS BEEN IMPOSSED   (ISCHTP = ',I10   ,/,&
'@    FOR LES (ITURB = ',I10,' )'                              ,/,&
'@                                                            ',/,&
'@  The calculation will   be executed                        ',/,&
'@                                                            ',/,&
'@  It is recommended to verify  usini1.                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8114 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA FOR  PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR V2F       (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-epsilon.               ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8115 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA FOR  PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR K-OMEGA   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of k-omega.                 ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8116 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA FOR  PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    A 2nd ORDER SCHEME HAS BEEN IMPOSED    (ISCHTP = ',I10   ,/,&
'@    FOR SPALART   (ITURB = ',I10,' )'                        ,/,&
'@                                                            ',/,&
'@   The current version does not support the 2nd order with  ',/,&
'@   coupling of the source terms of Spalart-Allmaras.        ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Modify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8121 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA FOR  PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8131 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA FOR  PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1 OR 2            ',/,&
'@    HERE IT IS  ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8132 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA FOR  PHASE ',I10    ,/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  With the chosen turbulence model   , ITURB = ',I10         ,/,&
'@    the value of ISTO2T (extrapolation of the source terms  ',/,&
'@    for the turbulent variables) cannot be modified         ',/,&
'@    yet ISTO2T has been forced to ',I10                      ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8141 format(                                                           &
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
'@  Verify   usini1.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8200 format(                                                           &
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
 8201 format(                                                           &
'@  ',I10   ,'        ',     I10                                 )
 8202 format(                                                           &
'@  ----------------------------                              ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify   usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8210 format(                                                           &
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
'@    Verify   usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8211 format(                                                           &
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
'@    Verify   usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8212 format(                                                           &
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
'@    Verify   usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8213 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@    ON THE VARIABLES THAT CONSTITUTE THE TEMPORAL AVERAGES  ',/,&
'@                                                            ',/,&
'@    The array  IDFMOM(JJ,IMOM) for  IMOM = ',I10             ,/,&
'@      must be assigned continuously.   . Yet here,          ',/,&
'@      IDFMOM(',I10,',IMOM) is not zero (=',I10   ,')        ',/,&
'@      while it exists    II < JJ for which                  ',/,&
'@      IDFMOM(II,IMOM) is zero.                              ',/,&
'@    The value of IDFMOM is wrongly set.                     ',/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify   usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8214 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE VERIFICATION OF DATA            ',/,&
'@    =========                                               ',/,&
'@  ON THE INSTANT OF THE START OF THE CALCULATION OF THE     ',/,&
'@  TEMPORAL AVERAGES                                         ',/,&
'@                                                            ',/,&
'@    The variable NTDMOM(IMOM)   for IMOM = ',I10             ,/,&
'@      must be assigned in usini1 to indicate from which     ',/,&
'@      time step (absolute) the calculation ot the           ',/,&
'@      corresponding average IMOM must start.                ',/,&
'@      NTDMOM(IMOM) here is  ',I10                            ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@    Verify   usini1.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8300 format(                                                           &
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
'@      Look at    usini1.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8301 format(                                                           &
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
'@      Look at    usini1.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8311 format(                                                           &
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
'@      Look at    usini1.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8312 format(                                                           &
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
'@      Look at    usini1.                                    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8313 format(                                                           &
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
 8390 format(                                                           &
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
 8400 format(                                                           &
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
'@    Verify the values of    IMOOLD in   usini1.             ',/,&
'@    Verify  that the restarting file is correct.            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8401 format(                                                           &
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
'@    Modify the value of    NTDMOM in   usini1               ',/,&
'@      Verify  that the restarting file is correct.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8402 format(                                                           &
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
'@    Verify   IMOOLD in   usini1      .                      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8403 format(                                                           &
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
'@    Do not modify   IMOOLD in   usini1.                     ',/,&
'@      or run a restart calculation.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8410 format(                                                           &
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
 8900 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING   : STOP AT THE INITIAL DATA VERIFICATION       ',/,&
'@    =========                                               ',/,&
'@     NAME OF THE VARIABLE TO BE CONTINUED TOO LARGE         ',/,&
'@                                                            ',/,&
'@  The type of calcultion defined in usini1                  ',/,&
'@    corresponds to a number of variables to continue in     ',/,&
'@    the listing and  post-processing equal to    ',I10       ,/,&
'@  The maximum number of variables to continue in            ',/,&
'@                           paramx.h is  NVPPMX = ',I10       ,/,&
'@                                                            ',/,&
'@  The calculation cannot be executed                        ',/,&
'@                                                            ',/,&
'@  Verify   usini1.                                          ',/,&
'@  Contact help.                                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif



!===============================================================================
! 5. FIN
!===============================================================================

return
end subroutine
