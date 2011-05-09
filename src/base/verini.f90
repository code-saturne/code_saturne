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

subroutine verini &
!================

 ( iok    )

!===============================================================================
!  FONCTION  :
!  ---------

! VERIFICATION DES PARAMETRES DE CALCUL
!   APRES INTERVENTION UTILISATEUR
!   (COMMONS)
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
use cstnum
use dimens
use numvar
use optcal
use mltgrd
use cstphy
use entsor
use albase
use alstru
use parall
use period
use ppthch
use ppppar
use ppincl
use lagpar
use lagran
use radiat
use mesh

!===============================================================================

implicit none

! Arguments

integer          iok

! Local variables

character        chaine*80, chain2*80
integer          ii    , iis   , jj    , iisct , iphas
integer          iscal , iest  , iiesca, ivar
integer          iuiph , iviph , iwiph, ikiph, ieiph
integer          itrbph, nbsccp
integer          ipp   , imgrok, nbccou
integer          iokpre, indest, itests, iiidef, istop
integer          iresop, ipolop
double precision testth
double precision arakfr

!===============================================================================

! Initialize variables to avoid warnings

jj = 0

!===============================================================================
! 1. ENTREES SORTIES entsor : formats 1000
!===============================================================================

! --- Dimension

if (ndim.ne.3) then
  write(nfecra,1100)ndim
  iok = iok + 1
endif

! --- Suite, Chrono, Historiques, Listing

if(ntchr.ne.-1.and.ntchr.le.0) then
  WRITE(NFECRA,1210) 'NTCHR  (Periode   Sortie Chrono.)',NTCHR
  iok = iok + 1
endif

do ipp = 2, nvppmx
  if(ichrvr(ipp).ne.1.and.ichrvr(ipp).ne.0) then
    chaine=nomvar(ipp)
    write(nfecra,1220)chaine(1:8),ipp,ichrvr(ipp)
    iok = iok + 1
  endif
enddo

if(ncapt.lt.0.or.ncapt.gt.ncaptm) then
  write(nfecra,1230)ncaptm,ncapt
  iok = iok + 1
endif

if(nthist.le.0.and.nthist.ne.-1) then
  WRITE(NFECRA,1210) 'NTHIST (Periode   Sortie Histo. )',NTHIST
  iok = iok + 1
endif

if(nthsav.lt.-1) then
  WRITE(NFECRA,1200) 'NTHSAV (Periode   Save   Histo. )',NTHSAV
  iok = iok + 1
endif

do ipp = 2, nvppmx
  if( ihisvr(ipp,1).gt.ncapt.or.                                  &
     (ihisvr(ipp,1).lt.0.and.ihisvr(ipp,1).ne.-1) ) then
    chaine=nomvar(ipp)
    write(nfecra,1240)chaine(1:8),ipp,ncapt,ihisvr(ipp,1)
    iok = iok + 1
  endif
enddo

do ipp = 2, nvppmx
  if( (ihisvr(ipp,1).gt.0.and.ihisvr(ipp,1).lt.ncapt)) then
    do jj = 1, ihisvr(ipp,1)
      if(ihisvr(ipp,jj+1).le.0.or.ihisvr(ipp,jj+1).gt.ncapt) then
        chaine=nomvar(ipp)
        write(nfecra,1250)                                        &
          chaine(1:8),ipp,jj+1,ncapt,ihisvr(ipp,jj+1)
        iok = iok + 1
      endif
    enddo
  endif
enddo

do ipp = 2, nvppmx
  if( ilisvr(ipp).ne.0.and.ilisvr(ipp).ne.1) then
    chaine=nomvar(ipp)
    write(nfecra,1260)chaine(1:8),ipp,ilisvr(ipp)
    iok = iok + 1
  endif
enddo

do ipp = 2, nvppmx
  if(itrsvr(ipp).ne.0) then
    if(itrsvr(ipp).le.0.or.itrsvr(ipp).gt.nvar) then
      chaine=nomvar(ipp)
      write(nfecra,1270)chaine(1:8),ipp,nvar,itrsvr(ipp)
      iok = iok + 1
    endif
  endif
enddo

if(ntlist.ne.-1.and.ntlist.le.0) then
  WRITE(NFECRA,1210) 'NTLIST (Periode   Sortie Listing)',NTLIST
  iok = iok + 1
endif

! --- Post traitement automatique (bord)

if(ipstdv.ne.1.and.                                               &
   mod(ipstdv,ipstyp).ne.0.and.                                   &
   mod(ipstdv,ipstcl).ne.0.and.                                   &
   mod(ipstdv,ipstft).ne.0) then
  WRITE(NFECRA,1300) 'IPSTDV',                                    &
                     'IPSTYP',IPSTYP,                             &
                     'IPSTCL',IPSTCL,                             &
                     'IPSTFT',IPSTFT,                             &
                     'IPSTFO',IPSTFO,                             &
                     'IPSTDV',IPSTDV
  iok = iok + 1
endif


!===============================================================================
! 2. OPTIONS DU CALCUL : TABLEAUX DE optcal : formats 2000
!===============================================================================

! --- Dimensions

if(nphas.lt.0.or.nphas.gt.nphsmx) then
  WRITE(NFECRA,2000)'NPHAS ',NPHSMX,NPHAS
  iok = iok + 1
endif
if(nscal.lt.0.or.nscal.gt.nscamx) then
  WRITE(NFECRA,2000)'NSCAL ',NSCAMX,NSCAL
  iok = iok + 1
endif
if(nscaus.lt.0.or.nscaus.gt.nscamx) then
  WRITE(NFECRA,2000)'NSCAUS',NSCAMX,NSCAUS
  iok = iok + 1
endif
if(nscapp.lt.0.or.nscapp.gt.nscamx) then
  WRITE(NFECRA,2000)'NSCAPP',NSCAMX,NSCAPP
  iok = iok + 1
endif
if(nvar.lt.0.or.nvar.gt.nvarmx) then
  WRITE(NFECRA,2000)'NVAR  ',NVARMX,NVAR
  iok = iok + 1
endif

! --- Rho et visc constants ou variables

do iphas = 1, nphas
  if(irovar(iphas).ne.0.and.irovar(iphas).ne.1) then
    WRITE(NFECRA,2201)'IROVAR',IROVAR(IPHAS)
    iok = iok + 1
  endif
  if(ivivar(iphas).ne.0.and.ivivar(iphas).ne.1) then
    WRITE(NFECRA,2201)'IVIVAR',IVIVAR(IPHAS)
    iok = iok + 1
  endif
enddo

! --- Definition des equations, schema en temps, schema convectif

do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if( (iconv (ii).ne.0.and.iconv (ii).ne.1).or.                 &
        (istat (ii).ne.0.and.istat (ii).ne.1).or.                 &
        (idiff (ii).ne.0.and.idiff (ii).ne.1).or.                 &
        (idifft(ii).ne.0.and.idifft(ii).ne.1).or.                 &
        (thetav(ii).gt.1.d0.or.thetav(ii).lt.0.d0).or.            &
        (blencv(ii).gt.1.d0.or.blencv(ii).lt.0.d0).or.            &
        (ischcv(ii).ne.0.and.ischcv(ii).ne.1).or.                 &
        (isstpc(ii).ne.0.and.isstpc(ii).ne.1)   ) then
      chaine=nomvar(ipp)
      write(nfecra,2100) chaine(1:8),                             &
                         istat(ii),iconv(ii),                     &
                         idiff(ii),idifft(ii),                    &
                         thetav(ii),blencv(ii),ischcv(ii),        &
                         isstpc(ii)
      iok = iok + 1
    endif
  endif
enddo


!     Extrap de rho : necessairement rho variable
do iphas = 1, nphas
  if(irovar(iphas).eq.0.and.iroext(iphas).gt.0) then
    write(nfecra,2005)iroext(iphas),irovar(iphas)
    iok = iok + 1
  endif
enddo

!     Coherence des vitesses
!       Pour le moment, theta est fixe automatiquement dans modini
!       On conserve quand meme le test pour plus tard puisqu'il est ecrit.
do iphas = 1, nphas
  if(abs(thetav(iv(iphas))-thetav(iu(iphas))).gt.epzero.or.       &
     abs(thetav(iw(iphas))-thetav(iu(iphas))).gt.epzero) then
    write(nfecra,2111) thetav(iu(iphas)),thetav(iv(iphas)), &
                       thetav(iw(iphas))
    iok = iok + 1
  endif
enddo


!     Theta pression : vaut 1
!       Pour le moment, theta est fixe automatiquement dans modini
!         (on ne devrait donc jamais voir cet affichage)
!       On conserve quand meme le test pour plus tard puisqu'il est ecrit.
do iphas = 1, nphas
  jj = ipr(iphas)
  if(abs(thetav(jj)-1.0d0).gt.epzero) then
    ipp    = ipprtp(jj)
    chaine=nomvar(ipp)
    write(nfecra,2112) thetav(jj)
    iok = iok + 1
  endif
enddo

!     Flux de masse et proprietes
!       Pour le moment, theta est fixe automatiquement dans modini
!       Donc pas de test sur sa valeur ici

!     On verifie que toutes les phases sont traitees pareil
testth = thetfl(1)
itests = istmpf(1)
do iphas = 1, nphas

  if(abs(testth-thetfl(iphas)).gt.epzero) then
    write(nfecra,2113) testth,thetfl(iphas)
    iok = iok + 1
  endif
  if(itests.ne.istmpf(iphas)) then
    write(nfecra,2114) itests,istmpf(iphas)
    iok = iok + 1
  endif
enddo

!     En LES il y a des verification de coherence supplementaires.
!       (simple avertissement si on s'ecarte des choix std)
!        mais stop si on fait plus de 5% d'upwind
!     Schema centre sans/avec test de pente, nwsrsm
do iphas = 1, nphas
  if(itytur(iphas).eq.4) then
    do ii = 1,3
      if(ii.eq.1) jj = iu(iphas)
      if(ii.eq.2) jj = iv(iphas)
      if(ii.eq.3) jj = iw(iphas)
      ipp    = ipprtp(jj)
      chaine=nomvar(ipp)
      if(abs(thetav(jj)-0.5d0).gt.epzero) then
        write(nfecra,2121) chaine(1:8),thetav(jj)
      endif
      if (blencv(jj).lt.0.95d0) then
        write(nfecra,2127) chaine(1:8),blencv(jj)
        iok = iok + 1
      elseif(abs(blencv(jj)-1.d0).gt.epzero) then
        write(nfecra,2122) chaine(1:8),blencv(jj)
      endif
      if(isstpc(jj).eq.0) then
        write(nfecra,2123) chaine(1:8),isstpc(jj)
      endif
    enddo
  endif
  if(itytur(iphas).eq.4.or.ischtp(iphas).eq.2) then
    do ii = 1,3
      if(ii.eq.1) jj = iu(iphas)
      if(ii.eq.2) jj = iv(iphas)
      if(ii.eq.3) jj = iw(iphas)
      ipp    = ipprtp(jj)
      chaine=nomvar(ipp)
      iiidef = 10
      if(nswrsm(jj).ne.iiidef) then
        write(nfecra,2125) chaine(1:8),iiidef,nswrsm(jj)
      endif
    enddo
    jj = ipr(iphas)
    ipp    = ipprtp(jj)
    chaine=nomvar(ipp)
    iiidef = 5
    if(nswrsm(jj).ne.iiidef) then
      write(nfecra,2125) chaine(1:8),iiidef,nswrsm(jj)
    endif
  endif
enddo
do ii = 1, nscal
  iphas = 1
  if(itytur(iphas).eq.4) then
    jj    = isca(ii)
    ipp   = ipprtp(jj)
    chaine=nomvar(ipp)
    if(abs(thetav(jj)-0.5d0).gt.epzero) then
      write(nfecra,2121) chaine(1:8),thetav(jj)
    endif
    if (blencv(jj).lt.0.95d0) then
      write(nfecra,2127) chaine(1:8),blencv(jj)
      iok = iok + 1
    elseif(abs(blencv(jj)-1.d0).gt.epzero) then
      write(nfecra,2122) chaine(1:8),blencv(jj)
    endif
    if(isstpc(jj).eq.1) then
      write(nfecra,2124) chaine(1:8),isstpc(jj)
    endif
  endif
  if(itytur(iphas).eq.4.or.ischtp(iphas).eq.2) then
    iiidef = 10
    if(nswrsm(jj).ne.iiidef) then
      write(nfecra,2125) chaine(1:8),iiidef,nswrsm(jj)
    endif
  endif
enddo

!     Test du theta de la viscosite secondaire, du flux de masse et
!     de la viscosite par rapport a celui de la vitesse
do iphas = 1, nphas
  jj = iu(iphas)
  if( abs(thetav(jj)-1.d0).lt.epzero.and.                         &
       (istmpf(iphas).eq.2.or.                                    &
        isno2t(iphas).ne.0.or.                                    &
        isto2t(iphas).ne.0.or.                                    &
        iroext(iphas).ne.0.or.                                    &
        iviext(iphas).ne.0.or.                                    &
        icpext(iphas).ne.0   ) ) then
    write(nfecra,2131) thetav(jj),                          &
         istmpf(iphas),isno2t(iphas),isto2t(iphas),               &
         iroext(iphas),iviext(iphas),icpext(iphas)
  endif
  if( abs(thetav(jj)-0.5d0).lt.epzero.and.                        &
       (istmpf(iphas).ne.2.or.                                    &
        isno2t(iphas).ne.1.or.                                    &
        isto2t(iphas).ne.1.or.                                    &
        iroext(iphas).ne.1.or.                                    &
        iviext(iphas).ne.1.or.                                    &
        icpext(iphas).ne.1   ) ) then
    write(nfecra,2132) thetav(jj),                          &
         istmpf(iphas),isno2t(iphas),isto2t(iphas),               &
         iroext(iphas),iviext(iphas),icpext(iphas)
  endif
enddo
do iscal = 1, nscal
  iphas = 1
  if(isso2t(iscal).ne.isno2t(iphas))then
    write(nfecra,2133) iscal,isso2t(iscal),isno2t(iphas)
  endif
  if(ivsext(iscal).ne.iviext(iphas))then
    write(nfecra,2134) iscal,ivsext(iscal),iviext(iphas)
  endif
enddo

!     Test du theta de la diffusivite des scalaires et de Cp : ils doivent etre
!       variables en (en espace) si on les extrapole (en temps) (...)
do iphas = 1, nphas
  if( icpext(iphas).gt.0 .and. icp(iphas).le.0 ) then
    write(nfecra,2135) icpext(iphas), icp(iphas)
    iok = iok + 1
  endif
enddo
do iscal = 1, nscal
  if( ivsext(iscal).gt.0 .and. ivisls(iscal).le.0 ) then
    write(nfecra,2136) iscal, ivsext(iscal), ivisls(iscal)
    iok = iok + 1
  endif
enddo



!     Pour les tests suivants : Utilise-t-on un estimateur d'erreur ?
indest = 0
do iphas = 1, nphas
  do iest = 1, nestmx
    iiesca = iescal(iest,iphas)
    if(iiesca.gt.0) then
      indest = 1
    endif
  enddo
enddo

!     Estimateurs incompatibles avec calcul a champ de vitesse
!       fige (on ne fait rien, sauf ecrire des betises dans le listing)
 if(indest.eq.1.and.iccvfg.eq.1) then
   write(nfecra,2137)
   iok = iok + 1
 endif

!     A priori, pour le moment, l'ordre 2 en temps
!       (rho, visc, termes sources N.S, theta vitesse)
!     est incompatible avec
!       - estimateurs
!       - ipucou
!       - iphydr et icalhy
!       - dt variable en espace ou en temps et stationnaire
!     Ici on s'arrete si on n'est pas dans le cas du schema std
do iphas = 1, nphas
  iuiph = iu(iphas)
  iviph = iv(iphas)
  iwiph = iw(iphas)
  if( (abs(thetav(iuiph)-1.0d0).gt.1.d-3).or.                     &
      (abs(thetav(iviph)-1.0d0).gt.1.d-3).or.                     &
      (abs(thetav(iwiph)-1.0d0).gt.1.d-3).or.                     &
      (    thetsn(iphas)       .gt.0.d0 ).or.                     &
      (    isno2t(iphas)       .gt.0    ).or.                     &
      (    thetro(iphas)       .gt.0.d0 ).or.                     &
      (    iroext(iphas)       .gt.0    ).or.                     &
      (    thetvi(iphas)       .gt.0.d0 ).or.                     &
      (    iviext(iphas)       .gt.0    )    ) then
     if(indest.eq.1.or.ipucou.eq.1.or.                            &
        iphydr.eq.1.or.icalhy.eq.1.or.                            &
        idtvar.eq.1.or.idtvar.eq.2.or.idtvar.lt.0) then
       write(nfecra,2140)                                         &
            thetav(iuiph),thetav(iviph),thetav(iwiph),            &
            isno2t(iphas),thetsn(iphas),                          &
            iroext(iphas),thetro(iphas),                          &
            iviext(iphas),thetvi(iphas)
       iok = iok + 1
     endif
   endif
 enddo


!     Iterations sur navsto
!     Doit etre un entier superieur ou egal a 1
!     Pour le moment, on interdit NTERUP > 1 et
!       estimateurs, matrices poids, pression hydrostatique
!       et algo stationnaire

if(nterup.le.0) then
  WRITE(NFECRA,3100) 'NTERUP', NTERUP
  iok = iok + 1
endif

if(nterup.gt.1) then

  if(ipucou.eq.1.or.indest.eq.1.or.                               &
       ippmod(icompf).ge.0.or.iccvfg.eq.1.or.                     &
       iphydr.eq.1.or.icalhy.eq.1.or.idtvar.eq.-1) then
    write(nfecra,2141) nterup
    iok = iok + 1
  endif

endif

!     A priori, pour le moment, l'ordre 2 en temps
!       n'est pas pris en compte en k-eps, v2f ou k-omega couple : on s'arrete
do iphas = 1, nphas
  if (itytur(iphas).eq.2 .and.ikecou(iphas).eq.1) then
    if((    thetst(iphas)       .gt.0.d0 ).or.                    &
       (    isto2t(iphas)       .gt.0    ).or.                    &
       (abs(thetav(ik (iphas))-1.0d0).gt.epzero).or.              &
       (abs(thetav(iep(iphas))-1.0d0).gt.epzero) ) then
      write(nfecra,2142)iturb(iphas),ikecou(iphas),         &
           thetst(iphas),isto2t(iphas),                           &
           thetav(ik (iphas)),thetav(iep(iphas))
      iok = iok + 1
    endif
  endif
  if (iturb(iphas).eq.50.and.ikecou(iphas).eq.1) then
    if((    thetst(iphas)       .gt.0.d0 ).or.                    &
       (    isto2t(iphas)       .gt.0    ).or.                    &
       (abs(thetav(ik  (iphas))-1.0d0).gt.epzero).or.             &
       (abs(thetav(iep (iphas))-1.0d0).gt.epzero).or.             &
       (abs(thetav(iphi(iphas))-1.0d0).gt.epzero).or.             &
       (abs(thetav(ifb (iphas))-1.0d0).gt.epzero) ) then
      write(nfecra,2143)iturb(iphas),ikecou(iphas),         &
           thetst(iphas),isto2t(iphas),                           &
           thetav(ik  (iphas)),thetav(iep (iphas)),               &
           thetav(iphi(iphas)),thetav(ifb (iphas))
      iok = iok + 1
    endif
  endif
  if (iturb(iphas).eq.60.and.ikecou(iphas).eq.1) then
    if((    thetst(iphas)       .gt.0.d0 ).or.                    &
       (    isto2t(iphas)       .gt.0    ).or.                    &
       (abs(thetav(ik  (iphas))-1.0d0).gt.epzero).or.             &
       (abs(thetav(iomg(iphas))-1.0d0).gt.epzero) ) then
      write(nfecra,2144)iturb(iphas),ikecou(iphas),         &
           thetst(iphas),isto2t(iphas),                           &
           thetav(ik  (iphas)),thetav(iomg(iphas))
      iok = iok + 1
    endif
  endif
  if (iturb(iphas).eq.70) then
    if((    thetst(iphas)       .gt.0.d0 ).or.                    &
       (    isto2t(iphas)       .gt.0    ).or.                    &
       (abs(thetav(inusa(iphas))-1.0d0).gt.epzero) ) then
      write(nfecra,2145)iturb(iphas),                       &
           thetst(iphas),isto2t(iphas),                           &
           thetav(inusa  (iphas))
      iok = iok + 1
    endif
  endif
enddo

!     A priori, pour le moment, l'ordre 2 en temps
!       (rho, visc, cp, termes sources N.S., Turb., Scal., theta)
!     est incompatible avec les physiques particulieres
!     Ici on s'arrete si on n'est pas dans le cas du schema std

if(ippmod(iphpar).ge.1) then
  istop = 0
  do ivar = 1, nvar
    if( (abs(thetav(ivar)-1.0d0).gt.1.d-3) ) istop = 1
  enddo
  do iphas = 1, nphas
    if((    thetsn(iphas)       .gt.0.d0 ).or.                    &
       (    isno2t(iphas)       .gt.0    ).or.                    &
       (    thetro(iphas)       .gt.0.d0 ).or.                    &
       (    iroext(iphas)       .gt.0    ).or.                    &
       (    thetvi(iphas)       .gt.0.d0 ).or.                    &
       (    iviext(iphas)       .gt.0    ).or.                    &
       (    thetcp(iphas)       .gt.0.d0 ).or.                    &
       (    icpext(iphas)       .gt.0    )    ) istop = 1
  enddo
  do iscal = 1, nscal
    if((    thetss(iscal)       .gt.0.d0 ).or.                    &
       (    isso2t(iscal)       .gt.0    ).or.                    &
       (    thetvs(iphas)       .gt.0.d0 ).or.                    &
       (    ivsext(iphas)       .gt.0    )    ) istop = 1
  enddo

  if(istop.ne.0) then
    write(nfecra,2146)
    iok = iok + 1
  endif
endif

!     A priori, pour le moment, l'ordre 2 en temps
!       n'est pas pris en compte pour les termes issus du Lagrangien.
!       On pourrait le signaler et continuer : on s'arrete.
if(iilagr .eq. 2) then
  do iphas = 1, nphas
    if((    thetsn(iphas)       .gt.0.d0 ).or.                    &
      (    isno2t(iphas)       .gt.0    ).or.                     &
       (    thetst(iphas)       .gt.0.d0 ).or.                    &
       (    isto2t(iphas)       .gt.0    ) ) then
      write(nfecra,2147)thetsn(iphas),isno2t(iphas),        &
                        thetst(iphas),isto2t(iphas)
      iok = iok + 1
    endif
  enddo
  do iscal = 1, nscal
    if (iscsth(iscal).eq.1.or.iscsth(iscal).eq.2) then
      if((    thetss(iscal)       .gt.0.d0 ).or.                  &
         (    isso2t(iscal)       .gt.0    )) then
        write(nfecra,2148)                                        &
    'lagrangien ',ISCAL,THETSS(ISCAL),ISSO2T(ISCAL),'uslag1'
        iok = iok + 1
      endif
    endif
  enddo
endif

!     A priori, pour le moment, l'ordre 2 en temps
!       n'est pas pris en compte pour les termes issus du rayonnement.
!       On pourrait le signaler et continuer : on s'arrete.
if (iirayo.gt.0) then
  do iphas = 1, nphas
    do iscal = 1, nscal
      if (iscal.eq.iscalt(iphas)) then
        if((    thetss(iscal)       .gt.0.d0 ).or.                &
           (    isso2t(iscal)       .gt.0    )) then
          write(nfecra,2148)                                      &
         'rayonnement',ISCAL,THETSS(ISCAL),ISSO2T(ISCAL),'usray1'
          iok = iok + 1
        endif
      endif
    enddo
  enddo
endif

! --- Algorithme stationnaire
if (idtvar.lt.0) then
  do ipp = 2, nvppmx
    ii = itrsvr(ipp)
    if (ii.ge.1) then
      if (relaxv(ii).gt.1d0.or.relaxv(ii).lt.0d0) then
        chaine=nomvar(ipp)
        write(nfecra,2149) chaine(1:8),relaxv(ii)
        iok = iok + 1
      endif
    endif
  enddo
  do iphas = 1, nphas
    if((relaxv(iv(iphas)).ne.relaxv(iu(iphas)))                   &
   .or.(relaxv(iw(iphas)).ne.relaxv(iu(iphas))) ) then
      write(nfecra,2150) relaxv(iu(iphas)),relaxv(iv(iphas)),     &
           relaxv(iw(iphas))
      iok = iok + 1
    endif
  enddo
!       L'algorithme stationnaire n'est pas compatible avec le module Lagrangien
  if (iilagr.ne.0) then
    write(nfecra,2151) iilagr
    iok = iok + 1
  endif
!       L'algorithme stationnaire n'est pas compatible avec la LES
  do iphas = 1, nphas
    if (itytur(iphas).eq.4) then
      write(nfecra,2152) iturb(iphas)
      iok = iok + 1
    endif
  enddo
endif

! --- Reconstruction des gradients
if(imrgra.ne.0.and.imrgra.ne.1.and.                               &
   imrgra.ne.2.and.imrgra.ne.3.and.                               &
   imrgra.ne.4 ) then
  WRITE(NFECRA,2205) 'IMRGRA',IMRGRA
  iok = iok + 1
endif

! On verifie l'angle de non orthogonalite de selection du
!   voisinage etendu dans le cas du moindre carre qui l'utilise

if(imrgra.eq.3) then
  if(anomax.gt.pi*0.5d0.or.anomax.lt.0.d0) then
    write(nfecra,2206) anomax, imrgra
  endif
endif

! Extrapolation : indetermination possible par mc,
!     necessitant un traitement particulier dans gradmc,
!     pour lequel on fait certaines hypotheses
if(imrgra.eq.1.or.imrgra.eq.2.or.imrgra.eq.3) then
  if((nphas.gt.1).or.                                             &
     ((abs(extrag(ipr(1))-1.d0).gt.epzero).and.                   &
      (abs(extrag(ipr(1))     ).gt.epzero)  )) then
    write(nfecra,2207) imrgra, nphas, extrag(ipr(1))
    iok = iok + 1
  endif
endif

! Les nombres de sweeps n'ont pas a etre verifies :
!  ce sont simplement des entiers (negatifs si on veut etre sur de ne
!  *jamais* entrer dans les boucles)

do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if(imligr(ii).gt.1) then
      chaine=nomvar(ipp)
      write(nfecra,2300) chaine(1:8),ii,imligr(ii)
      iok = iok + 1
    endif
  endif
enddo

do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if(ircflu(ii).ne.1.and.ircflu(ii).ne.0) then
      chaine=nomvar(ipp)
      write(nfecra,2310) chaine(1:8),ii,ircflu(ii)
      iok = iok + 1
    endif
  endif
enddo
! Non reconstruction des flux en SOLU n'a pas de sens pour la convection
do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if(ircflu(ii).eq.0.and.ischcv(ii).eq.0.and.                   &
                           blencv(ii).ne.0.d0) then
      chaine=nomvar(ipp)
      write(nfecra,2311) chaine(1:8),ii,ircflu(ii),ii,ischcv(ii)
      iok = iok + 1
    endif
  endif
enddo

! Il n'y a pas besoin de test sur les epsilons
!   Ce sont simplement des reels
!   Une valeur negative indique qu'on veut atteindre
!   le nombre d'iterations maximal

do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if(climgr(ii).lt.1.d0) then
      chaine=nomvar(ipp)
      write(nfecra,2320) chaine(1:8),ii,climgr(ii)
      iok = iok + 1
    endif
  endif
enddo


! EXTRAG non nul permis uniquement pour la pression.
!        et dans ce cas egal a 1
do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if(abs(extrag(ii)     ).ge.epzero) then
      iokpre = 0
      do iphas = 1, nphas
        if (ii.eq.ipr(iphas)) then
          iokpre = 1
          if(abs(extrag(ii)-1.d0).ge.epzero) then
            chaine=nomvar(ipp)
            write(nfecra,2330) chaine(1:8),ii,extrag(ii)
            iok = iok + 1
          endif
        endif
      enddo
      if(iokpre.eq.0) then
        chaine=nomvar(ipp)
        write(nfecra,2331) chaine(1:8),ii,extrag(ii)
        iok = iok + 1
      endif
    endif
  endif
enddo


! --- Solveurs iteratifs

! Il n'y a pas besoin de test sur les epsilons
!   Ce sont simplement des reels
!   Une valeur negative indique qu'on veut atteindre
!   le nombre d'iterations maximal
! Il n'y a pas besoin de test sur le nombre d'iterations
!   Ce sont simplement des entiers
!   Une valeur negative indique au'on veut sortir de suite

do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if(iresol(ii).ne.-1) then
      iresop = mod(iresol(ii),1000)
      ipolop = (iresol(ii)-iresop)/1000
      if ((iresop.lt.0.or.iresop.gt.3).or.                        &
          (iresop.eq.1.and.ipolop.ne.0)) then
        chaine=nomvar(ipp)
        write(nfecra,2400) chaine(1:8),ii,iresol(ii)
        iok = iok + 1
      endif
    endif
    if (idircl(ii).ne.0.and.idircl(ii).ne.1) then
      chaine=nomvar(ipp)
      write(nfecra,2401) chaine(1:8),ii,idircl(ii)
      iok = iok + 1
    endif
  endif
enddo

do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if ((iresol(ii).eq.0.and.iconv(ii).eq.1).or.                  &
        (iresol(ii).eq.1.and.iconv(ii).eq.0)) then
      chaine=nomvar(ipp)
      write(nfecra,2410) chaine(1:8),ii,iresol(ii),iconv(ii)
    endif
  endif
enddo

! --- Le multigrille sera verifie d'un seul bloc plus bas.


! --- Suite de calcul

if(isuite.ne.0.and.isuite.ne.1) then
  WRITE(NFECRA,2200) 'ISUITE',ISUITE
  iok = iok + 1
endif
if(ileaux.ne.0.and.ileaux.ne.1) then
  WRITE(NFECRA,2200) 'ILEAUX',ILEAUX
  iok = iok + 1
endif
if(iecaux.ne.0.and.iecaux.ne.1) then
  WRITE(NFECRA,2200) 'IECAUX',IECAUX
  iok = iok + 1
endif
! En LES, on previent que ce n'est pas malin de ne pas relire le fichier
!   auxiliaire
if(iecaux.eq.0.or.ileaux.eq.0) then
  do iphas = 1, nphas
    if(itytur(iphas).eq.4) then
      write(nfecra,2420) iturb(iphas),ileaux,iecaux
    endif
  enddo
endif

! --- Reperage du temps et marche en temps

! Le nombre de pas de temps pourrait bien etre negatif : pas de test.

if(inpdt0.ne.0.and.inpdt0.ne.1) then
  WRITE(NFECRA,2200) 'INPDT0',INPDT0
  iok = iok + 1
endif

if(idtvar.lt.-1.or.idtvar.gt.2) then
  WRITE(NFECRA,2500) 'IDTVAR',IDTVAR
  iok = iok + 1
endif

if(idtvar.gt.0.and.varrdt.lt.0.d0) then
  WRITE(NFECRA,2510)'VARRDT', VARRDT
  iok = iok + 1
endif

if(dtref .lt.0.d0) then
  WRITE(NFECRA,2510)'DTREF ', DTREF
  iok = iok + 1
endif

if(dtmin.le.0.d0   .or. dtmax.le.0.d0 .or.                        &
   dtmin.gt.dtmax                       ) then
  write(nfecra,2520) dtmin, dtmax
  if(idtvar.gt.0) then
    iok = iok + 1
  endif
endif

do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if(cdtvar(ii).le.0.d0) then
      chaine=nomvar(ipp)
      write(nfecra,2530) chaine(1:8),ii,cdtvar(ii)
      iok = iok + 1
    endif
  endif
enddo

if(iptlro.ne.0.and.iptlro.ne.1) then
  WRITE(NFECRA,2200) 'IPTLRO',IPTLRO
  iok = iok + 1
endif

!     Si on est en pas de temps constant, on ne touche pas le pas de temps,
!     mais on indique juste les depassements de critere local.
if(idtvar.eq.0.and.iptlro.eq.1) then
  write(nfecra,2540) iptlro,idtvar
endif
if(idtvar.eq.-1.and.iptlro.eq.1) then
  write(nfecra,2541) iptlro,idtvar
endif

! --- Turbulence

!    Modele

do iphas = 1, nphas
  itrbph = iturb(iphas)
  if ( itrbph.ne. 0.and.itrbph.ne.10.and.itrbph.ne.20.and.        &
       itrbph.ne.21.and.itrbph.ne.30.and.itrbph.ne.31.and.        &
       itrbph.ne.40.and.itrbph.ne.41.and.itrbph.ne.42.and.        &
       itrbph.ne.50.and.itrbph.ne.60.and.itrbph.ne.70  ) then
    WRITE(NFECRA,2600) 'ITURB  ',ITRBPH
    iok = iok + 1
  endif

  ! In lagrangian with two-way coupling, k-omega SST is forbidden (not
  ! properly implemented)
  if (itrbph.eq.60 .and. iilagr.eq.2) then
     write(nfecra,2601) iilagr
     iok = iok + 1
  endif

enddo

!     Methode des vortex pour la LES

if (ivrtex.ne.0 .and.ivrtex.ne.1) then
  WRITE(NFECRA,2200) 'IVRTEX ',IVRTEX
  iok = iok + 1
endif
! On impose qu'il n'y ait qu'une seule phase
if (ivrtex.eq.1 .and. nphas.gt.1) then
  write(nfecra,2605)nphas,ivrtex
  iok = iok + 1
endif
iphas = 1
if(ivrtex.eq.1.and.itytur(iphas).ne.4) then
  write(nfecra,2606)itytur(iphas),ivrtex
  ivrtex = 0
endif

!    Nb de variables

if(nscal.ge.1) then
  do iphas = 1, nphas
    if(iscalt(iphas).gt.nscal) then
      write(nfecra,2610)                                          &
                 'NUMERO DU SCALAIRE TEMPERATURE ',ISCALT(IPHAS), &
                 'NOMBRE DE SCALAIRES            ',NSCAL
      iok = iok + 1
    endif
    if(  (nvar.lt. 4+nscal               ) .or.                   &
         (nvar.lt. 6+nscal.and.itytur(iphas).eq.2).or.            &
         (nvar.lt.11+nscal.and.itytur(iphas).eq.3).or.            &
         (nvar.lt. 8+nscal.and.iturb(iphas).eq.50).or.            &
         (nvar.lt. 6+nscal.and.iturb(iphas).eq.60).or.            &
         (nvar.lt. 5+nscal.and.iturb(iphas).eq.70)      ) then
      write(nfecra,2610)                                          &
                 'NOMBRE DE VARIABLES            ',NVAR,          &
                 'NOMBRE DE SCALAIRES            ',NSCAL
      iok = iok + 1
    endif
  enddo
endif

do iphas = 1, nphas

  if(ideuch(iphas).lt.0.or.ideuch(iphas).gt.2) then
    WRITE(NFECRA,2211)'IDEUCH',IDEUCH(IPHAS)
    iok = iok + 1
  endif
  if (ideuch(iphas).ne.0 .and.                                    &
       (iturb(iphas).eq.0 .or. iturb(iphas).eq.10 .or.            &
       itytur(iphas).eq.4 .or. iturb(iphas).eq.7 )) then
     write(nfecra,2209)iturb(iphas),ideuch(iphas)
     iok = iok + 1
  endif
  if(ilogpo(iphas).ne.0.and.ilogpo(iphas).ne.1) then
    WRITE(NFECRA,2201)'ILOGPO',ILOGPO(IPHAS)
    iok = iok + 1
  endif

!      Specifique k-epsilon, v2f et k-omega

 if(itytur(iphas).eq.2 .or. iturb(iphas).eq.50                    &
       .or. iturb(iphas).eq.60 ) then
    if( (nvar.le.5.and.itytur(iphas).eq.2) .or.                   &
        (nvar.le.7.and.iturb(iphas).eq.50) .or.                   &
        (nvar.le.5.and.iturb(iphas).eq.60)     ) then
      write(nfecra,2610)                                          &
                 'NOMBRE DE VARIABLES            ',NVAR,          &
                 'OPTION POUR LA TURBULENCE      ',ITURB(IPHAS)
      iok = iok + 1
    endif
!     Le choix de ICLKEP n'est possible qu'en k-eps ou v2f
    if (iturb(iphas).ne.60) then
      if(iclkep(iphas).ne.0.and.iclkep(iphas).ne.1) then
        WRITE(NFECRA,2201)'ICLKEP',ICLKEP(IPHAS)
        iok = iok + 1
      endif
    endif
    if(ikecou(iphas).ne.0.and.ikecou(iphas).ne.1) then
      WRITE(NFECRA,2201)'IKECOU',IKECOU(IPHAS)
      iok = iok + 1
    endif
!     En k-eps a prod lin et en v2f on force IKECOU a 0
    if (ikecou(iphas).eq.1 .and.                                  &
         (iturb(iphas).eq.21 .or. iturb(iphas).eq.50)) then
      write(nfecra,2208)iturb(iphas),ikecou(iphas)
      iok = iok + 1
    endif
!     En stationnaire on force IKECOU a 0
    if (ikecou(iphas).ne.0.and.idtvar.lt.0) then
      write(nfecra,2210)ikecou(iphas)
      iok = iok + 1
    endif

    if(igrhok(iphas).ne.0.and.igrhok(iphas).ne.1) then
      WRITE(NFECRA,2201)'IGRHOK',IGRHOK(IPHAS)
      iok = iok + 1
    endif
    if(igrake(iphas).ne.0.and.igrake(iphas).ne.1) then
      WRITE(NFECRA,2201)'IGRAKE',IGRAKE(IPHAS)
      iok = iok + 1
    endif
!        IF( IGRAKE.EQ.1.AND.(GX**2+GY**2+GZ**2).LE.EPZERO**2 ) THEN
!          WRITE(NFECRA,2620)'IGRAKE',IGRAKE,GX,GY,GZ
!          IOK = IOK + 1
!        ENDIF
    if(nscal.gt.0) then
      if(iscalt(iphas).le.0.and.                                  &
           (gx**2+gy**2+gz**2).ge.epzero**2) then
        write(nfecra,2621) gx,gy,gz,iscalt(iphas)
        if(igrake(iphas).eq.1) then
          WRITE(NFECRA,2622)'IGRAKE',IGRAKE(IPHAS)
        endif
!MO            IOK = IOK + 1
      endif
    endif

!     Si RELAXV(IK) a ete modifie par l'utilisateur mais que IKECOU n'est
!     pas egal a 0, on previent que ce sera sans effet
!     Sinon on verifie que RELAXV(IK) est bien compris entre 0 et 1
!     (en stationnaire cela a deja ete fait plus haut)
    ikiph = ik(iphas)
    if (itytur(iphas).eq.6) then
      ieiph = iomg(iphas)
    else
      ieiph = iep(iphas)
    endif
    if ( (abs(relaxv(ikiph)+999.d0).gt.epzero .or.                &
          abs(relaxv(ieiph)+999.d0).gt.epzero ) .and.             &
          ikecou(iphas).ne.0) write(nfecra,2623)                  &
            relaxv(ikiph),relaxv(ieiph)
    if (ikecou(iphas).eq.0 .and. idtvar.ge.0) then
      if(relaxv(ikiph).gt.1.d0.or.relaxv(ikiph).lt.0.d0 .or.      &
         relaxv(ieiph).gt.1.d0.or.relaxv(ieiph).lt.0.d0) then
        write(nfecra,2624) relaxv(ikiph),relaxv(ieiph)
        iok = iok + 1
      endif
    endif

  endif

!     Specifique Rij-epsilon

  if(itytur(iphas).eq.3) then
    if(nvar.le.10) then
      write(nfecra,2610)                                          &
                 'NOMBRE DE VARIABLES            ',NVAR,          &
                 'OPTION POUR LA TURBULENCE      ',ITURB(IPHAS)
      iok = iok + 1
    endif
    if(irijnu(iphas).ne.0.and.irijnu(iphas).ne.1) then
      WRITE(NFECRA,2201)'IRIJNU',IRIJNU(IPHAS)
      iok = iok + 1
    endif
    if(irijrb(iphas).ne.0.and.irijrb(iphas).ne.1) then
      WRITE(NFECRA,2201)'IRIJRB',IRIJRB(IPHAS)
      iok = iok + 1
    endif
    if (iturb(iphas).eq.30) then
!     echo de paroi et implicitation speciale de la diffusion de epsilon
!     seulement en Rij standard
      if(irijec(iphas).ne.0.and.irijec(iphas).ne.1) then
        WRITE(NFECRA,2201)'IRIJEC',IRIJEC(IPHAS)
        iok = iok + 1
      endif
      if(idifre(iphas).ne.0.and.idifre(iphas).ne.1) then
        WRITE(NFECRA,2201)'IDIFRE',IDIFRE(IPHAS)
        iok = iok + 1
      endif
    endif
    if(igrari(iphas).ne.0.and.igrari(iphas).ne.1) then
      WRITE(NFECRA,2201)'IGRARI',IGRARI(IPHAS)
      iok = iok + 1
    endif
!        IF( IGRARI.EQ.1.AND.(GX**2+GY**2+GZ**2).LE.EPZERO**2 ) THEN
!          WRITE(NFECRA,2620)'IGRARI',IGRARI,GX,GY,GZ
!          IOK = IOK + 1
!        ENDIF
    if(iclsyr(iphas).ne.0.and.iclsyr(iphas).ne.1) then
      WRITE(NFECRA,2201)'ICLSYR',ICLSYR(IPHAS)
      iok = iok + 1
    endif
    if(iclptr(iphas).ne.0.and.iclptr(iphas).ne.1) then
      WRITE(NFECRA,2201)'ICLPTR',ICLPTR(IPHAS)
      iok = iok + 1
    endif
    if(nscal.gt.0) then
      if(iscalt(iphas).le.0.and.                                  &
           (gx**2+gy**2+gz**2).ge.epzero**2) then
        write(nfecra,2621)gx,gy,gz,iscalt(iphas)
        if(igrari(iphas).eq.1) then
          WRITE(NFECRA,2622)'IGRARI',IGRARI(IPHAS)
        endif
!MO            IOK = IOK + 1
      endif
    endif
  endif

!     Specifique LES

  if(itytur(iphas).eq.4) then
    if(idries(iphas).ne.1.and.idries(iphas).ne.0) then
      WRITE(NFECRA,2201)'IDRIES',IDRIES(IPHAS)
      iok = iok + 1
    endif
    if(idries(iphas).ne.0.and.(iturb(iphas).eq.41.or.iturb(iphas).eq.42)) then
      write(nfecra,2630) idries(iphas),iturb(iphas)
      iok = iok + 1
    endif
!         La reduction du voisinage etendu peut degrader
!         les resultats du modele dynamique en LES
    if(iturb(iphas).eq.41.and.imrgra.eq.3) then
      write(nfecra,2607) iturb(iphas), imrgra
    endif
  endif

enddo

! --- Stokes


if(iprco .ne.0.and.iprco .ne.1) then
  WRITE(NFECRA,2200) 'IPRCO ',IPRCO
  iok = iok + 1
endif
if(iprco.eq.1) then
  do iphas = 1, nphas
    if(irevmc(iphas).ne.0.and.irevmc(iphas).ne.1.and.             &
                              irevmc(iphas).ne.2) then
      WRITE(NFECRA,2211) 'IREVMC',IREVMC(IPHAS)
      iok = iok + 1
    endif
    arakfr = arak(iphas)
    if (idtvar.lt.0) arakfr=arakfr*relaxv(iu(iphas))
    if(arakfr.gt.1.d0 .or. arakfr.lt.0.d0) then
      WRITE(NFECRA,2640) 'ARAK  ',ARAKFR
      iok = iok + 1
    endif
    if( relaxv(ipr(iphas)).gt.1d0 .or.                            &
        relaxv(ipr(iphas)).lt.0d0     ) then
      write(nfecra,2625) relaxv(ipr(iphas))
      iok = iok + 1
    endif
  enddo
endif

! --- Couplage U-P

if (ipucou.ne.0.and.ipucou.ne.1) then
  WRITE(NFECRA,2200) 'IPUCOU ',IPUCOU
  iok = iok + 1
endif

! Incompatibilite pour le moment des matrices poids avec un theta schema
! Si theta n'est pas egal a 1 pour la vitesse (incompatibilite du
! pas de temps variable aussi)

do iphas = 1, nphas
  jj = iu(iphas)
  if((abs(thetav(jj)-1.0d0).gt.epzero).and.                       &
    ((idtvar.ne.0).or.(ipucou.eq.1))) then
    write(nfecra,2204) thetav(jj),idtvar,ipucou
  endif
enddo

! --- Prise en compte de la pression hydrostatique

!     IPHYDR et ICALHY sont testes dans varpos
!       (valeur 0 ou 1)

if (iphydr.eq.1.and.nphas.gt.1) then
  write(nfecra,2202) iphydr,nphas
  iok = iok +1
endif

! --- Champ de vitesse fige

if (iccvfg.ne.0.and.iccvfg.ne.1) then
  WRITE(NFECRA,2200) 'ICCVFG ',ICCVFG
  iok = iok + 1
endif

! --- Interpolation face des viscosites

if(imvisf.ne.0.and.imvisf.ne.1) then
  WRITE(NFECRA,2200) 'IMVISF',IMVISF
  iok = iok + 1
endif

! --- Traitement de la temperature pour couplage SYRTHES
!     Verification du nombre de couplages

if(itbrrb.ne.0 .and. itbrrb.ne.1       ) then
  WRITE(NFECRA,2200) 'ITBRRB',ITBRRB
  iok = iok + 1
endif

!     On regarde si ICPSYR a des valeurs realistes
if(nscal.gt.0) then
  do iscal = 1, nscal
    if(icpsyr(iscal).ne.0.and.icpsyr(iscal).ne.1) then
      chaine=nomvar(ipprtp(isca(iscal)))
      WRITE(NFECRA,2650)CHAINE(1:8),'ICPSYR',ISCAL,ICPSYR(ISCAL)
      iok = iok + 1
    endif
  enddo
endif

!     On compte le nombre de scalaires couples
nbsccp = 0
if(nscal.gt.0) then
  do iscal = 1, nscal
    nbsccp = nbsccp+icpsyr(iscal)
  enddo
endif

!     On regarde s'il y a du couplage

call nbcsyr (nbccou)
!==========

!     S'il n'y a pas de couplage
if(nbccou.eq.0) then

!       et qu'il n'y a pas zero scalaire couple, on s'arrete
  if(nbsccp.ne.0) then
    write(nfecra,2660)nbsccp,nscal
    iok = iok + 1
  endif

!     Sinon, s'il y a du couplage
else

!       et qu'il n'y a pas un et un seul scalaire couple, on s'arrete
  if(nbsccp.ne.1) then
    write(nfecra,2661)nscal,nbsccp
    iok = iok + 1
  endif

!          que le scalaire couple n'est pas la temperature, on s'arrete
!          attention : tout est pret, mais il n'y a pas eu de valid.
!          en outre, ca permet de bien verifier que l'utilisateur
!          ne s'est pas trompe dans 99% des cas d'utilisation
!         (en compressible, on couple l'energie)
  do iscal = 1, nscal
    if(icpsyr(iscal).eq.1) then
      if(ippmod(icompf).lt.0) then
        if(abs(iscsth(iscal)).ne.1) then
          write(nfecra,2662)iscal,iscsth(iscal)
          iok = iok + 1
        endif
      else
        if(iscsth(iscal).ne.3) then
          write(nfecra,2663)iscal,iscal,iscsth(iscal)
          iok = iok + 1
        endif
      endif
    endif
  enddo

endif


! --- Estimateurs  d'erreur pour Navier-Stokes

do iphas = 1, nphas
  do iest = 1, nestmx
    iiesca = iescal(iest,iphas)
    if (iiesca.ne.0.and.iiesca.ne.1.and.iiesca.ne.2) then
      write(nfecra,2664) iest,iest,iiesca,            &
                         iespre,iesder,iescor,iestot
      iok = iok + 1
    endif
  enddo
enddo


! --- Distance a la paroi

if(ineedy.eq.1) then

  if(abs(icdpar).ne.1.and.abs(icdpar).ne.2) then
    write(nfecra,2700) icdpar
    iok = iok + 1
  endif
  if(nitmay.lt.1) then
    WRITE(NFECRA,3100) 'NITMAY',NITMAY
    iok = iok + 1
  endif
  if(imligy.gt.1) then
    WRITE(NFECRA,2750) 'IMLIGY',IMLIGY
    iok = iok + 1
  endif
  if(ircfly.ne.1.and.ircfly.ne.0) then
    WRITE(NFECRA,2200) 'IRCFLY',IRCFLY
    iok = iok + 1
  endif
  if(ischcy.ne.1.and.ischcy.ne.0) then
    WRITE(NFECRA,2200) 'ISCHCY',ISCHCY
    iok = iok + 1
  endif
  if(isstpy.ne.1.and.isstpy.ne.0) then
    WRITE(NFECRA,2200) 'ISSTPY',ISSTPY
    iok = iok + 1
  endif
  if(imgrpy.ne.1.and.imgrpy.ne.0) then
    WRITE(NFECRA,2200) 'IMGRPY',IMGRPY
    iok = iok + 1
  endif
  if(ntcmxy.lt.1) then
    WRITE(NFECRA,3100) 'NTCMXY',NTCMXY
    iok = iok + 1
  endif

  if(blency.gt.1.d0.or.blency.lt.0.d0) then
    WRITE(NFECRA,2710) 'BLENCY',BLENCY
    iok = iok + 1
  endif
  if(climgy.lt.1.d0) then
    WRITE(NFECRA,2720) 'CLIMGY',CLIMGY
    iok = iok + 1
  endif
  if(abs(extray-1.d0).gt.epzero.and.abs(extray).gt.epzero) then
    WRITE(NFECRA,2730) 'EXTRAY',EXTRAY
    iok = iok + 1
  endif
  if(coumxy.le.0.d0) then
    WRITE(NFECRA,2740) 'COUMXY',COUMXY
    iok = iok + 1
  endif
  if(yplmxy.le.0.d0) then
    WRITE(NFECRA,2740) 'YPLMXY',YPLMXY
    iok = iok + 1
  endif

endif


!===============================================================================
! 2. MULTIGRILLE : TABLEAUX DU MULTIGRILLE : formats 3000
!===============================================================================

! --- Options generales

do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if(imgr(ii).ne.0.and.imgr(ii).ne.1) then
      chaine=nomvar(ipp)
      WRITE(NFECRA,3000) CHAINE(1:8),'IMGR  ',II,IMGR(II)
      iok = iok + 1
    endif
    do iphas = 1, nphas
      if(imgr(ii).eq.1.and.iconv(ii).eq.1) then
        write(nfecra,3001)
        iok = iok + 1
      endif
    enddo
    if(imgr(ii).eq.1) then
      if(ncymax(ii).le.0) then
        WRITE(NFECRA,3010) CHAINE(1:8),'NCYMAX',II,NCYMAX(II)
        iok = iok + 1
      endif
      if(nitmgf(ii).le.0) then
        WRITE(NFECRA,3010) CHAINE(1:8),'NITMGF',II,NITMGF(II)
        iok = iok + 1
      endif
    endif
  endif
enddo
imgrok = 0
do ipp = 2, nvppmx
  ii = itrsvr(ipp)
  if(ii.ge.1) then
    if(imgr(ii).eq.1) then
      imgrok = 1
    endif
  endif
enddo

! --- Options specifiques de niveau plus eleve

if(imgrok.eq.1.and.ncegrm.le.0)then
  WRITE(NFECRA,3100)'NCEGRM',NCEGRM
  iok = iok + 1
endif
if(imgrok.eq.1.and.ngrmax.le.0)then
  WRITE(NFECRA,3100)'NGRMAX',NGRMAX
  iok = iok + 1
endif

!===============================================================================
! 3. TABLEAUX DE cstphy : formats 4000
!===============================================================================

do iphas = 1, nphas

! --- Constantes physiques de chaque phase

  if(ro0(iphas)   .lt.0d0) then
    WRITE(NFECRA,2511)'RO0   ', RO0(IPHAS)
    iok = iok + 1
  endif
  if(viscl0(iphas).lt.0d0) then
    WRITE(NFECRA,2511)'VISCL0', VISCL0(IPHAS)
    iok = iok + 1
  endif

! --- Turbulence

!    On a besoin de UREF si on initialise la turbulence (debut de calcul
!      en turbulence ou suite laminaire->turbulent)
!    Ici on met juste un avertissement, car sans UREF l'utilisateur peut
!      initialiser la turbulence a la main. Un test complementaire sera fait
!      dans inivar.
  if(itytur(iphas).eq.2.or.itytur(iphas).eq.3                     &
       .or.iturb(iphas).eq.50.or.iturb(iphas).eq.60               &
       .or.iturb(iphas).eq.70) then
    if(uref(iphas)  .lt.0.d0) then
      write(nfecra,4100) uref(iphas)
    endif
  endif

  if(iturb(iphas).eq.10) then
    if(xlomlg(iphas).le.0.d0) then
      WRITE(NFECRA,2511)'XLOMLG', XLOMLG(IPHAS)
      iok = iok + 1
    endif
  endif

!     LES
  if(itytur(iphas).eq.4) then
    if(xlesfl(iphas).lt.0.d0) then
      WRITE(NFECRA,2511) 'XLESFL', XLESFL(IPHAS)
      iok = iok + 1
    endif
    if(ales  (iphas).lt.0.d0) then
      WRITE(NFECRA,2511) 'ALES  ', ALES(IPHAS)
      iok = iok + 1
    endif
    if(bles  (iphas).lt.0.d0) then
      WRITE(NFECRA,2511) 'BLES  ', BLES(IPHAS)
      iok = iok + 1
    endif
    if(csmago(iphas).lt.0.d0) then
      WRITE(NFECRA,2511) 'CSMAGO', CSMAGO(IPHAS)
      iok = iok + 1
    endif
    if(cwale(iphas).lt.0.d0) then
      WRITE(NFECRA,2511) 'CWALE', CWALE(IPHAS)
      iok = iok + 1
    endif
    if(idries(iphas).eq.1.and.cdries(iphas).lt.0) then
      WRITE(NFECRA,2511) 'CDRIES', CDRIES(IPHAS)
      iok = iok + 1
    endif
    if(iturb(iphas).eq.41) then
      if(xlesfd(iphas).lt.0.d0) then
        WRITE(NFECRA,2511) 'XLESFD', XLESFD(IPHAS)
        iok = iok + 1
      endif
      if(smagmx(iphas).lt.0.d0) then
        WRITE(NFECRA,2511) 'SMAGMX', SMAGMX(IPHAS)
        iok = iok + 1
      endif
    endif
  endif

enddo

! --- Scalaires

if(nscal.gt.0) then

!     Scalaire passif, temperature, enthalpie, energie
  do ii = 1, nscal
    if (iscsth(ii).lt.-1.or.iscsth(ii).gt.3) then
      chaine=nomvar(ipprtp(isca(ii)))
      write(nfecra,4300)chaine(1:8),ii,iscsth(ii)
      iok = iok + 1
    endif
  enddo

!     Scalaire associe dans le cas des variances
  do ii = 1, nscal
    if (iscavr(ii).gt.nscal.or.iscavr(ii).lt.0) then
      chaine=nomvar(ipprtp(isca(ii)))
      write(nfecra,4320)chaine(1:8),ii,nscal,iscavr(ii)
      iok = iok + 1
    endif
  enddo

!     On verifie que le scalaire associe a une fluctuation n'est
!     pas une fluctuation
  do ii = 1, nscal
    if (iscavr(ii).gt.0) then
      if (iscavr(iscavr(ii)).gt.0) then
        chaine=nomvar(ipprtp(isca(ii)))
        chain2=nomvar(ipprtp(isca(iscavr(ii))))
        write(nfecra,4321)chaine(1:8),chain2(1:8),ii,iscavr(ii),  &
             iscavr(ii),iscavr(iscavr(ii))
        iok = iok + 1
      endif
    endif
  enddo

!     Mode de clipping dans le cas des variances
!       pour les scalaires non variance, on n'utilise pas ICLVFL.
!       On demande donc a l'utilisateur de ne pas y toucher
!       (ca permet d'etre sur qu'il sait ce qu'il fait)
  do ii = 1, nscal
    if (iscavr(ii).le.nscal.and.iscavr(ii).gt.0) then
      if(iclvfl(ii).ne.0.and.                                     &
         iclvfl(ii).ne.1.and.iclvfl(ii).ne.2) then
        chaine=nomvar(ipprtp(isca(ii)))
        write(nfecra,4330)chaine(1:8),ii,iclvfl(ii)
        iok = iok + 1
      endif
    elseif (iscavr(ii).eq.0) then
      if(iclvfl(ii).ne.-1) then
        chaine=nomvar(ipprtp(isca(ii)))
        write(nfecra,4331)chaine(1:8),ii,iclvfl(ii)
        iok = iok + 1
      endif
    endif
  enddo

!     Valeur de la diffusivite positive (si cste)
!        si pas cste, on verifiera apres usphyv
  do ii = 1, nscal
    if (ivisls(ii).le.0.and.visls0(ii).lt.0d0) then
      chaine=nomvar(ipprtp(isca(ii)))
      write(nfecra,4340)chaine(1:8),ii,ii,ivisls(ii),visls0(ii)
      iok = iok + 1
    endif
  enddo

!     Valeur du sigma positif
  do ii = 1, nscal
    if (sigmas(ii).le.0d0) then
      chaine=nomvar(ipprtp(isca(ii)))
      write(nfecra,4350)chaine(1:8),ii,sigmas(ii)
      iok = iok + 1
    endif
  enddo

!     Si on n'utilise pas la borne inf de clipping, on demande
!      a l'utilisateur de ne pas y toucher (ca permet d'etre sur
!      qu'il sait ce qu'il fait)
  do ii = 1, nscal
    if(iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.               &
       iclvfl(ii).ne.2.and.scamin(ii).ne.-grand) then
      chaine=nomvar(ipprtp(isca(ii)))
      write(nfecra,4360)chaine(1:8),ii,scamin(ii),ii,iclvfl(ii)
      iok = iok + 1
    endif
  enddo

!     Si on n'utilise pas la borne sup de clipping, on demande
!      a l'utilisateur de ne pas y toucher (ca permet d'etre sur
!      qu'il sait ce qu'il fait)
  do ii = 1, nscal
    if(iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.               &
       iclvfl(ii).ne.2.and.scamax(ii).ne.grand) then
      chaine=nomvar(ipprtp(isca(ii)))
      write(nfecra,4361)chaine(1:8),ii,scamax(ii),ii,iclvfl(ii)
      iok = iok + 1
    endif
  enddo

!     Valeur de la borne sup de clipping si on l'utilise
  do ii = 1, nscal
    if(iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.               &
       iclvfl(ii).eq.2.and.scamax(ii).le.0.d0) then
      chaine=nomvar(ipprtp(isca(ii)))
      write(nfecra,4370)chaine(1:8),ii,scamax(ii)
      iok = iok + 1
    endif
  enddo

!     Warning sur Rvarfl < 0
  do ii = 1, nscal
    if(iscavr(ii).gt.0.and.iscavr(ii).le.nscal.and.               &
                           rvarfl(ii).le.0.d0) then
      chaine=nomvar(ipprtp(isca(ii)))
      write(nfecra,4380)chaine(1:8),ii,rvarfl(ii)
      iok = iok + 1
    endif
  enddo


!  Rien a verifier sur SCAMIN SCAMAX

!     Si CP0 est utilise (resolution d'un scalaire en temperature
!       et CP constant), il doit etre positif
  do iphas = 1, nphas
    if(icp(iphas).eq.0) then
      if (cp0(iphas).lt.0.d0) then
        iisct = 0
        do iis = 1, nscal
          if (abs(iscsth(iis)).eq.1) then
            iisct = 1
          endif
        enddo
        if (iisct.eq.1) then
          WRITE(NFECRA,2511)'CP0   ',CP0(IPHAS)
          iok = iok + 1
        endif
      endif
    endif
  enddo

endif

!===============================================================================
! 4. TABLEAUX DE period : formats 5000
!===============================================================================

! --- periodicite de rotation incompatible avec couplage
!       renforce vitesse pression et l'ALE
if(iperot.gt.0.and. (ipucou.ne.0.or.iale.ne.0)) then
  write(nfecra,5003)iperio,ipucou,iale
  iok = iok + 1
endif

! --- periodicite incompatible avec le mode de calcul
!       direct de la distance a la paroi
if(iperio.eq.1.and.ineedy.eq.1.and.abs(icdpar).eq.2) then
  write(nfecra,5005)iperio, icdpar
  iok = iok + 1
endif


! --- periodicite de rotation incompatible avec le rayonnement DOM
!if(iperio.gt.0.and.iirayo.gt.0) then ! de translation aussi ?
if(iperot.gt.0.and.iirayo.gt.0) then
  if(iirayo.eq.1) then
    write(nfecra,5008) iperio , iirayo
    iok = iok + 1
  endif
endif

! --- periodicite de rotation douteuse avec rij
!      (et donc a fortiori avec ordre 2 sur Rij)
if(iperot.gt.0) then
  do iphas = 1, nphas
    if(itytur(iphas).eq.3) then
      write(nfecra,5009)iperio,iturb(iphas)
!            IOK = IOK + 1
    endif
  enddo
endif

! --- periodicite de rotation douteuse avec ordre 2 vitesse
if(iperot.gt.0) then
  do iphas = 1, nphas
    iuiph = iu(iphas)
    iviph = iv(iphas)
    iwiph = iw(iphas)
    if( (abs(thetav(iuiph)-0.5d0).lt.1.d-3).or.                   &
        (abs(thetav(iviph)-0.5d0).lt.1.d-3).or.                   &
        (abs(thetav(iwiph)-0.5d0).lt.1.d-3)) then
      write(nfecra,5010)iperio,                                   &
        thetav(iuiph),thetav(iviph),thetav(iwiph)
!            IOK = IOK + 1
    endif
  enddo
endif


!===============================================================================
! 5. TABLEAUX DE parall : formats 6000 (limitations)
!===============================================================================

! --- parallelisme incompatible avec lagrangien
if(irangp.ge.0.and.iilagr.ne.0) then
  write(nfecra,6002)irangp,iilagr
  iok = iok + 1
endif

! --- parallelisme incompatible avec le mode de calcul
!       direct de la distance a la paroi
if(irangp.ge.0.and.ineedy.eq.1.and.abs(icdpar).eq.2) then
  write(nfecra,6005)irangp, icdpar
  iok = iok + 1
endif

!===============================================================================
! 6. METHODE ALE (albase, alstru) : formats 7000
!===============================================================================

if (iale.ne.0 .and. iale.ne.1) then
  write(nfecra,7000)iale
  iok = iok + 1
endif

if (iale.eq.1) then

  if (nalinf.lt.0) then
    write(nfecra,7010)nalinf
    iok = iok + 1
  endif

  if (alpnmk.lt.0.d0 .or. alpnmk.gt.1.d0  .or.                    &
       gamnmk.lt.0.d0 .or. gamnmk.gt.1.d0  .or.                   &
       betnmk.lt.0.d0 .or. betnmk.gt.0.5d0 ) then
    write(nfecra,7020)alpnmk,betnmk,gamnmk
    iok = iok + 1
  endif

  if (nalimx.le.0) then
    write(nfecra,7030)nalimx
    iok = iok + 1
  endif
  if (epalim.le.0.d0) then
    write(nfecra,7040)epalim
    iok = iok + 1
  endif

  if (italin.ne.-999 .and. italin.ne.0 .and. italin.ne.1) then
    write(nfecra,7050)italin
    iok = iok + 1
  endif

endif

!===============================================================================
! 7. COMPRESSIBLE : formats 8000
!===============================================================================

if(ippmod(icompf).ge.0) then
  do iphas = 1, nphas
    if(t0(iphas).le.0.d0.or.p0(iphas).le.0.d0) then
      write(nfecra,8000)t0(iphas),p0(iphas)
      iok = iok + 1
    endif
  enddo
endif

!===============================================================================
! 8. FORMATS VERIFICATION
!===============================================================================

#if defined(_CS_LANG_FR)

 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    NDIM   DOIT ETRE UN ENTIER EGAL A 3                     ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute                            ',/,&
'@                                                            ',/,&
'@  La dimension de l espace est 3, meme pour les calculs     ',/,&
'@    physiquement bidimensionnels.                           ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@  Verifier le fichier maillage.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A33,                          ' DOIT ETRE UN ENTIER   ',/,&
'@    SUPERIEUR OU EGAL A -1                                  ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1210 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A33,                          ' DOIT ETRE UN ENTIER   ',/,&
'@    STRICTEMENT POSITIF OU EGAL A -1                        ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1220 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ICHRVR(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  ICHRVR indique si la variable doit etre incluse dans les  ',/,&
'@    fichiers de post-traitement                             ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1230 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    NCAPT  DOIT ETRE UN ENTIER INFERIEUR OU EGAL A ',I10     ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  NCAPT  est le nombre de sondes utilisees pour les         ',/,&
'@    historiques.                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1240 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IHISVR(',I10   ,',1) DOIT ETRE UN ENTIER EGAL A -1      ',/,&
'@      POSITIF OU NUL ET INFERIEUR OU EGAL A NCAPT =',I10     ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  IHISVR(I,1) indique les sondes a utiliser pour la variable',/,&
'@    I (-1 signifiant que toutes sont utilisees)             ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1250 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IHISVR(',I10   ,',',I10   ,') DOIT ETRE UN ENTIER       ',/,&
'@      STRICTEMENT POSITIF ET                                ',/,&
'@      INFERIEUR OU EGAL A NCAPT = ',I10                      ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute                            ',/,&
'@                                                            ',/,&
'@  IHISVR(I,j+1) indique le numero de la jieme sonde a       ',/,&
'@    utiliser pour la variable a post-traiter numero I       ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1260 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ILISVR(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute                            ',/,&
'@                                                            ',/,&
'@  ILISVR(I) indique si la variable I sera suivie lors des   ',/,&
'@    impressions dans le listing                             ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1270 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ITRSVR(',I10   ,') DOIT ETRE UN ENTIER COMPRIS ENTRE    ',/,&
'@      0 ET NVAR=',I10                                        ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE PARAMETRE ',A6,' DOIT ETRE UN ENTIER MULTIPLE DE     ',/,&
'@      DE CERTAINS DES ENTIERS SUIVANTS :                    ',/,&
'@      ',A6,' = ',I10                                         ,/,&
'@      ',A6,' = ',I10                                         ,/,&
'@      ',A6,' = ',I10                                         ,/,&
'@      ',A6,' = ',I10                                         ,/,&
'@    IL VAUT ICI ',A6,' = ', I10                              ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute                            ',/,&
'@                                                            ',/,&
'@  Ce parametre precise les variables supplementaires        ',/,&
'@    a post-traiter.                                         ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER                              ',/,&
'@      STRICTEMENT POSITIF ET INFERIEUR OU EGAL A ',I10       ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2005 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ON DEMANDE UNE EXTRAPOLATION TEMPORELLE DE RHO AVEC     ',/,&
'@      IROEXT(IPHAS) = ',I10                                  ,/,&
'@    CECI EST INCOMPATIBLE AVEC RHO CONSTANT                 ',/,&
'@      IROVAR(IPHAS) = ',I10                                  ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     PARAMETRES DU SCHEMA NUMERIQUE POUR LA VARIABLE ',A8    ,/,&
'@                                                            ',/,&
'@ Parametre               ISTAT     ICONV                    ',/,&
'@ Valeurs acceptees     0 ou  1   0 ou  1                    ',/,&
'@ Valeurs entrees ici',I10     ,I10                           ,/,&
'@                                                            ',/,&
'@ Parametre                         IDIFF     IDIFFT         ',/,&
'@ Valeurs acceptees               0 ou  1   0 ou  1          ',/,&
'@ Valeurs entrees ici',10X     ,I10      ,I10                 ,/,&
'@                                                            ',/,&
'@ Parametre               THETAV   BLENCV    ISCHCV    ISSTPC',/,&
'@ Valeurs acceptees     [0.; 1.] [0.; 1.]   0 ou  1   0 ou  1',/,&
'@ Valeurs entrees ici      ',E14.5     ,E14.5,I10  ,I10       ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2111 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@   INCOMPATIBILITE POUR LE SCHEMA EN TEMPS                  ',/,&
'@                                                            ',/,&
'@   Schema en temps pour la vitesse phase                    ',/,&
'@      THETA n''a pas la meme valeur pour les 3 composantes  ',/,&
'@                                                            ',/,&
'@ Parametre THETAV              U          V          W      ',/,&
'@ Valeurs entrees ici ',E10.2,E10.2,E10.2                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2112 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS          ',/,&
'@                                                            ',/,&
'@  LE PARAMETRE THETAV POUR LA PRESSION DOIT VALOIR 1        ',/,&
'@                                                            ',/,&
'@  Il vaut ici ',E14.5                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2113 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS          ',/,&
'@                                                            ',/,&
'@   EN L.E.S.                                                ',/,&
'@   LE PARAMETRE THETFL DU SCHEMA EN TEMPS POUR LE FLUX DE   ',/,&
'@     MASSE EST DIFFERENT DE CELUI DE LA PHASE 1 ',E14.5      ,/,&
'@     THETFL A ETE IMPOSE ICI A ',E14.5                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute                             ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2114 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS          ',/,&
'@                                                            ',/,&
'@   LE PARAMETRE ISTMPF DU SCHEMA EN TEMPS POUR LE FLUX DE   ',/,&
'@     MASSE EST DIFFERENT DE CELUI DE LA PHASE 1 ',I10        ,/,&
'@     ISTMPF A ETE IMPOSE ICI A ',I10                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute                             ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2121 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,&
'@                                                            ',/,&
'@   EN L.E.S.                                                ',/,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE THETAV DU SCHEMA ',/,&
'@     EN TEMPS DE LA VARIABLE ',A8  ,' EST 0.5               ',/,&
'@     THETAV A ETE IMPOSE ICI A ',E14.5                       ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute                                    ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2122 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,&
'@                                                            ',/,&
'@   EN L.E.S.                                                ',/,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE BLENCV DU SCHEMA ',/,&
'@     CONVECTIF DE LA VARIABLE ',A8  ,' EST 1.0              ',/,&
'@     BLENCV A ETE IMPOSE ICI A ',E14.5                       ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute                                    ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2123 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,&
'@                                                            ',/,&
'@   EN L.E.S.                                                ',/,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE ISSTPC DU SCHEMA ',/,&
'@     CONVECTIF DE LA VARIABLE ',A8  ,' EST 1                ',/,&
'@     ISSTPC A ETE IMPOSE ICI A ',I10                         ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute                                    ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2124 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,&
'@                                                            ',/,&
'@   EN L.E.S.                                                ',/,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE ISSTPC DU SCHEMA ',/,&
'@     CONVECTIF DE LA VARIABLE ',A8  ,' EST 0                ',/,&
'@     ISSTPC A ETE IMPOSE ICI A ',I10                         ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute                                    ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2125 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,&
'@                                                            ',/,&
'@   ORDRE 2 EN TEMPS OU LES                                  ',/,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE NSWRSM POUR      ',/,&
'@     LA VARIABLE ',A8  ,' EST  ',I10                        ,/, &
'@     NSWRSM A ETE IMPOSE ICI A ',I10                         ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute                                    ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2127 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,&
'@                                                            ',/,&
'@   EN L.E.S.                                                ',/,&
'@   LA VALEUR RECOMMANDEE POUR LE PARAMETRE BLENCV DU SCHEMA ',/,&
'@     CONVECTIF DE LA VARIABLE ',A8  ,' EST 1.0              ',/,&
'@     BLENCV A ETE IMPOSE ICI A ',E14.5                       ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute                             ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2131 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX DU SCHEMA EN TEMPS                                 ',/,&
'@                                                            ',/,&
'@     LE SCHEMA EN TEMPS POUR LA VITESSE EST D ORDRE 1       ',/,&
'@       (THETAV = ',E10.2 ,')                                ',/,&
'@     CERTAINS TERMES SONT CEPENDANT PRIS A L''ORDRE 2 AVEC  ',/,&
'@       LES CHOIX SUIVANTS :                                 ',/,&
'@                                                            ',/,&
'@ Parametres       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT ',/,&
'@ Valeurs entrees ',6I7                                       ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2132 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX DU SCHEMA EN TEMPS                                 ',/,&
'@                                                            ',/,&
'@     LE SCHEMA EN TEMPS POUR LA VITESSE EST D ORDRE 2       ',/,&
'@       (THETAV = ',E10.2 ,')                                ',/,&
'@     CERTAINS TERMES SONT CEPENDANT PRIS A L''ORDRE 1 AVEC  ',/,&
'@       LES CHOIX SUIVANTS :                                 ',/,&
'@                                                            ',/,&
'@ Parametres       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT ',/,&
'@ Valeurs entrees ',6I7                                       ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute.                                   ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2133 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,&
'@                                                            ',/,&
'@   SCALAIRE ',I10,' ISSO2T = ',I10                           ,/,&
'@     EST DIFFERENT DE ISNO2T                                ',/,&
'@     ISNO2T(IPHAS) = ',I10                                   ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute                                    ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2134 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX NON STANDARD DU SCHEMA EN TEMPS                    ',/,&
'@                                                            ',/,&
'@   SCALAIRE ',I10,' IVSEXT = ',I10                           ,/,&
'@     EST DIFFERENT DE IVIEXT                                ',/,&
'@     IVIEXT(IPHAS) = ',I10                                   ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute                                    ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2135 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX INCOMPATIBLE POUR LE SCHEMA EN TEMPS               ',/,&
'@                                                            ',/,&
'@     La  chaleur massique est extrapolee en temps avec      ',/,&
'@       ICPEXT(IPHAS) = ',I10                                 ,/,&
'@     Pour cela, elle doit etre variable, or                 ',/,&
'@       ICP(IPHAS)    = ',I10                                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute                             ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1 ',/,&
'@    - desactiver le choix d''extrapolation de Cp en temps   ',/,&
'@      ou                                                    ',/,&
'@    - imposer Cp variable                                   ',/,&
'@         (et le renseigner alors dans usphyv)               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2136 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX INCOMPATIBLE POUR LE SCHEMA EN TEMPS               ',/,&
'@                                                            ',/,&
'@   Scalaire ISCAL = ',I10                                    ,/,&
'@     La  diffusivite      est extrapolee en temps avec      ',/,&
'@       IVSEXT(ISCAL) = ',I10                                 ,/,&
'@     Pour cela, elle doit etre variable, or                 ',/,&
'@       IVISLS(ISCAL) = ',I10                                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute                             ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1 ',/,&
'@    - desactiver le choix d''extrapolation en temps         ',/,&
'@                                     de la diffusivite      ',/,&
'@      ou                                                    ',/,&
'@    - imposer la diffusivite variable                       ',/,&
'@         (et la renseigner alors dans usphyv)               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2137 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :      A L''ENTREE DES DONNEES                ',/,&
'@    =========                                               ',/,&
'@   CHOIX INCOMPATIBLE POUR LES ESTIMATEURS D''ERREUR        ',/,&
'@                                                            ',/,&
'@  On a active un ou plusieurs estimateurs d''erreur pour    ',/,&
'@    Navier-Stokes dans un calcul a champ de vitesse fige.   ',/,&
'@    Le ou les estimateurs ne seront pas calcules.           ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute                             ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou        ',/,&
'@    usini1 :                                                ',/,&
'@      desactiver les estimateurs d erreur ou                ',/,&
'@               le calcul a champ de vitesse fige (ICCVFG)   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2140 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  On souhaite utiliser un schema en temps d''ordre 2 :      ',/,&
'@      U,V,W : THETA = ',3E12.4                               ,/,&
'@      Termes sources Navier-Stokes: ISNO2T = ',I10           ,/,&
'@                                    THETSN = ',E12.4         ,/,&
'@      Masse volumique             : IROEXT = ',I10           ,/,&
'@                                    THETRO = ',E12.4         ,/,&
'@      Viscosite                   : IVIEXT = ',I10           ,/,&
'@                                    THETVI = ',E12.4         ,/,&
'@  La version actuelle ne le permet pas lorsque l''une des   ',/,&
'@    options suivantes a ete activee (c''est le cas ici) :   ',/,&
'@    - utilisation d''un estimateur d''erreur (IESCAL)       ',/,&
'@    - couplage instationnaire (IPUCOU)                      ',/,&
'@    - prise en compte specifique de la pression             ',/,&
'@      hydrostatique (IPHYDR et ICALHY)                      ',/,&
'@    - pas de temps variable en temps ou en espace  ou       ',/,&
'@      algorithme stationnaire (IDTVAR)                      ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2141 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,&
'@                                                            ',/,&
'@  On souhaite utiliser un couplage                          ',/,&
'@    vitesse-pression par point fixe (NTERUP = ',I10          ,/,&
'@  La version actuelle ne le permet pas lorsque l''une des   ',/,&
'@    options suivantes a ete activee (c''est le cas ici) :   ',/,&
'@    - utilisation d''un estimateur d''erreur (IESCAL)       ',/,&
'@    - couplage instationnaire (IPUCOU)                      ',/,&
'@    - prise en compte specifique de la pression             ',/,&
'@      hydrostatique (IPHYDR et ICALHY)                      ',/,&
'@    - algorithme stationnaire (IDTVAR=-1)                   ',/,&
'@    - module compressible (IPPMOD(ICOMPF)>=0)               ',/,&
'@    - champ de vitesse fige (ICCVFG=1)                      ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2142 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Avec le modele de turbulence                              ',/,&
'@    k-epsilon (ITURB = ',I10 ,') couple (IKECOU = ',I10,') :',/,&
'@    la version courante ne permet pas de traiter les        ',/,&
'@    equations du modele k-epsilon a l''ordre 2 en temps avec',/,&
'@    couplage.                                               ',/,&
'@    Une ou plusieurs valeurs parmi les suivantes ne sont    ',/,&
'@    donc pas permises :                                     ',/,&
'@                                                            ',/,&
'@       THETST    ISTO2T     THETA K   THETA EPS             ',/,&
'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2143 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Avec le modele de turbulence                              ',/,&
'@    v2f (ITURB = ',I10   ,') couple (IKECOU = ',I10,    ') :',/,&
'@    la version courante ne permet pas de traiter les        ',/,&
'@    equations du modele v2f a l''ordre 2 en temps avec      ',/,&
'@    couplage.                                               ',/,&
'@    Une ou plusieurs valeurs parmi les suivantes ne sont    ',/,&
'@    donc pas permises :                                     ',/,&
'@                                                            ',/,&
'@       THETST    ISTO2T     THETA K   THETA EPS             ',/,&
'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,&
'@     THETA PHI    THETA FB                                  ',/,&
'@ ',      E12.4,      E12.4                                   ,/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2144 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Avec le modele de turbulence                              ',/,&
'@    k-omega (ITURB = ',I10   ,') couple (IKECOU = ',I10,') :',/,&
'@    la version courante ne permet pas de traiter les        ',/,&
'@    equations du modele k-omega a l''ordre 2 en temps avec  ',/,&
'@    couplage.                                               ',/,&
'@    Une ou plusieurs valeurs parmi les suivantes ne sont    ',/,&
'@    donc pas permises :                                     ',/,&
'@                                                            ',/,&
'@       THETST    ISTO2T     THETA K   THETA OMEGA           ',/,&
'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2145 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Avec le modele de turbulence                              ',/,&
'@    Spallart-Allmaras (ITURB = ',I10   ,')                  ',/,&
'@    la version courante ne permet pas de traiter les        ',/,&
'@    l''ordre 2 en temps.                                    ',/,&
'@    couplage.                                               ',/,&
'@    Une ou plusieurs valeurs parmi les suivantes ne sont    ',/,&
'@    donc pas permises :                                     ',/,&
'@                                                            ',/,&
'@       THETST    ISTO2T     THETA K   THETA OMEGA           ',/,&
'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2146 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  La version actuelle ne permet pas de modifier le schema en',/,&
'@    temps lorsqu''une physique particuliere est activee     ',/,&
'@    (combustion, charbon, electrique).                      ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface, usini1   ',/,&
'@    et usppmo.                                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2147 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Les termes sources provenant du module                    ',/,&
'@    Lagrangien ne sont pas traites a l''ordre 2 en temps    ',/,&
'@    dans la version courante malgre le choix utilisateur    ',/,&
'@    suivant :                                               ',/,&
'@                                                            ',/,&
'@     Pour Navier-Stokes    Pour la turbulence               ',/,&
'@       THETSN    ISNO2T      THETST    ISTO2T               ',/,&
'@ ',     E12.4,      I10,      E12.4,      I10                ,/,&
'@                                                            ',/,&
'@  (Les autres termes sources pourraient etre traites a      ',/,&
'@   l''ordre 2)                                              ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface, usini1   ',/,&
'@  et uslag1.                                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2148 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS DE CALCUL INCOMPATIBLES AVEC LE SCHEMA EN TEMPS ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Les termes sources provenant du module ',A11               ,/,&
'@    ne sont pas traites a l''ordre 2 en temps               ',/,&
'@    dans la version courante malgre le choix utilisateur    ',/,&
'@    suivant :                                               ',/,&
'@                                                            ',/,&
'@       Pour le scalaire ',I10                                ,/,&
'@       THETSS    ISSO2T                                     ',/,&
'@ ',     E12.4,      I10                                      ,/,&
'@                                                            ',/,&
'@  (Les autres termes sources pourraient etre traites a      ',/,&
'@   l''ordre 2)                                              ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface, usini1   ',/,&
'@  et ',A6                                                    ,/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2149 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     ALGORITHME STATIONNAIRE                                ',/,&
'@     COEFFICIENT DE RELAXATION POUR LA VARIABLE ',A8         ,/,&
'@                                                            ',/,&
'@ RELAXV doit etre un reel compris entre 0 et 1              ',/,&
'@ Il vaut ici ',E14.5                                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2150  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@   INCOMPATIBILITE ALGORITHME STATIONNAIRE                  ',/,&
'@                                                            ',/,&
'@   Coefficient de relaxation vitesse phase ',I10             ,/,&
'@      RELAXV n''a pas la meme valeur pour les 3 composantes ',/,&
'@                                                            ',/,&
'@ Parametre RELAXV              U          V          W      ',/,&
'@ Valeurs entrees ici ',E10.2,E10.2,E10.2                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2151  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@   L''ALGORITHME STATIONNAIRE N''EST PAS COMPATIBLE AVEC    ',/,&
'@    LE MODULE LAGRANGIEN QUI EST UNE APPROCHE INSTATIONNAIRE',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  L''indicateur IILAGR a ete positionne a ',I10              ,/,&
'@    dans uslag1 (module lagrangien active pour IILAGR>0).   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2152  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@   L''ALGORITHME STATIONNAIRE N''EST PAS COMPATIBLE AVEC    ',/,&
'@    LA L.E.S. QUI EST UNE MODELISATION INSTATIONNAIRE DE    ',/,&
'@    LA TURBULENCE                                           ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  L''indicateur ITURB a ete positionne a ',I10               ,/,&
'@    dans usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2200 format(                                                           &
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
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2201 format(                                                           &
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
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2202 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LA PRISE EN COMPTE EXPLICITE DE LA PRESSION             ',/,&
'@    HYDROSTATIQUE N''EST IMPLANTEE QUE POUR LE MONOPHASIQUE ',/,&
'@    IPHYDR VAUT ICI ',I10                                    ,/,&
'@    ET NPHAS        ',I10                                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
! 2203 format(
!     &'@                                                            ',/,
!     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!     &'@                                                            ',/,
!     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
!     &'@    =========                                               ',/,
!     &'@  ON DEMANDE LA PRISE EN COMPTE EXPLICITE DE LA PRESSION    ',/,
!     &'@    HYDROSTATIQUE POUR LES CONDITIONS DE SORTIE EN PRESSION ',/,
!     &'@    AVEC ACCELERATION DE LA PESANTEUR NULLE                 ',/,
!     &'@                                                            ',/,
!     &'@  ICALHY VAUT ICI ',I10                                      ,/,
!     &'@  G      VAUT ICI ',3E14.5                                   ,/,
!     &'@                                                            ',/,
!     &'@  Le calcul ne sera pas execute.                            ',/,
!     &'@                                                            ',/,
!     &'@  Verifier les parametres donnes via l''interface ou usini1.',/,
!     &'@                                                            ',/,
!     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!     &'@                                                            ',/)
 2204 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@   DONNEES NON ADMISSIBLES POUR LE SCHEMA EN TEMPS          ',/,&
'@                                                            ',/,&
'@  ON DEMANDE LA PRISE UN SCHEMA EN TEMPS POUR LA VITESSE    ',/,&
'@    D''ORDRE 2 AVEC UN PAS DE TEMPS NON CONSTANT, UN        ',/,&
'@    ALGORITHME STATIONNAIRE OU LES MATRICES POIDS           ',/,&
'@                                                            ',/,&
'@  THETAV VAUT ICI ',E14.5,' POUR LA VITESSE                 ',/,&
'@  ALORS QUE IDTVAR VAUT ',I10                                ,/,&
'@  ET IPUCOU VAUT        ',I10                                ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2205 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1, 2, 3 OU 4       ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2206 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ANOMAX DOIT ETRE UN REEL POSITIF OU NUL ET              ',/,&
'@                             INFERIEUR OU EGAL A PI/2       ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  On demande la reconstruction des gradients par moindres   ',/,&
'@    carres sur voisinage etendu reduit (IMRGRA = ',I10  ,').',/,&
'@    Le critere est base sur l''angle de non orthogonalite   ',/,&
'@    des faces ANOMAX qui doit etre fourni en radians et     ',/,&
'@    compris dans les bornes indiquees ci-dessus.            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2207 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    L''UTILISATION DE LA METHODE DE CALCUL DE GRADIENT PAR  ',/,&
'@      MOINDRES CARRES EST IMPOSSIBLE AVEC                   ',/,&
'@        NPHAS SUPERIEUR A 1   OU                            ',/,&
'@        EXTRAG(IPR(1)) DIFFERENT DE 0 ET 1                  ',/,&
'@                                                            ',/,&
'@    ON A ICI                                                ',/,&
'@        IMRGRA         = ',I10                               ,/,&
'@        NPHAS          = ',I10                               ,/,&
'@        EXTRAG(IPR(1)) = ',E14.5                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2208 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    EN K-EPS PROD LIN (ITURB=21) ET EN V2F (ITURB=50)       ',/,&
'@    IKECOU DOIT ETRE EGAL A 0                               ',/,&
'@    ITURB  VAUT ICI ',I10                                    ,/,&
'@    IKECOU VAUT ICI ',I10                                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2209 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE MODELE DE PAROI A DEUX ECHELLES (IDEUCH=1 OU 2)      ',/,&
'@    EST INCOMPATIBLE AVEC UN CALCUL EN LAMINAIRE, EN        ',/,&
'@    LONGUEUR DE MELANGE, EN SPALART-ALLMARAS OU EN L.E.S.   ',/,&
'@    ON A ICI ITURB=',I10                                     ,/,&
'@         ET IDEUCH=',I10                                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2210 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    L''ALGORITHME STATIONNAIRE EST INCOMPATIBLE AVEC LE     ',/,&
'@    COUPLAGE DES TERMES SOURCES EN K-EPS, V2F OU K-OMEGA    ',/,&
'@                                                            ',/,&
'@    ON A ICI IKECOU=',I10                                    ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2211 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2             ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IMLIGR(',I10   ,') DOIT ETRE UN ENTIER                  ',/,&
'@      INFERIEUR OU EGAL A 1                                 ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  IMLIGR(I) indique le mode de limitation des gradients     ',/,&
'@    pour la variable I                                      ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2310 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IRCFLU(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  IRCFLU(I) indique si les flux sont reconstruits           ',/,&
'@    pour la variable I                                      ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2311 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IRCFLU(',I10   ,') = ',I10   ,' EST INCOMPATIBLE AVEC   ',/,&
'@    ISCHCV(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  IRCFLU(I) = 0 indique que les flux ne sont pas            ',/,&
'@    reconstruits pour la variable I.                        ',/,&
'@  ISCHCV(I) = 0 (schema SOLU) demande une reconstruction    ',/,&
'@    pour les flux convectifs.                               ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2320 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    CLIMGR(',I10   ,') DOIT ETRE UN REEL                    ',/,&
'@      SUPERIEUR OU EGAL A 1                                 ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  CLIMGR(I) est le coefficient de limitation des gradients  ',/,&
'@    pour la variable I                                      ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2330 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    EXTRAG(',I10   ,') DOIT ETRE UN REEL EGAL A 0 OU 1      ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2331 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    EXTRAG(',I10   ,') DOIT ETRE NUL                        ',/,&
'@      (Des valeurs non nulles sont autorisees pour la       ',/,&
'@       pression uniquement)                                 ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2400 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IRESOL(',I10   ,') DOIT ETRE UN ENTIER EGAL             ',/,&
'@                                     A -1 OU A IPOL*1000+ J ',/,&
'@      AVEC IPOL LE DEGRE DU POLYNOME DE PRECONDITIONNEMENT  ',/,&
'@        ET J    = 0 POUR GRADIENT CONJUGUE                  ',/,&
'@                = 1 POUR JACOBI   (IPOL = 0 DANS CE CAS)    ',/,&
'@                = 2 POUR BI-CGSTAB                          ',/,&
'@                = 3 POUR GMRES                              ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  IRESOL(I) indique le solveur lineaire a utiliser          ',/,&
'@    pour la variable I                                      ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2401 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IDIRCL(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  IDIRCL(I) indique si le code doit decaler la diagonale de ',/,&
'@    la matrice de la variable I en l''absence de Dirichlet  ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2410 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    RISQUE D ECHEC A LA RESOLUTION DU SYSTEME LINEAIRE      ',/,&
'@    IRESOL(',I10   ,') = ',I10                               ,/,&
'@      ET LA VARIABLE EST CONVECTEE (ICONV = ',I10,')        ',/,&
'@                                                            ',/,&
'@  Le calcul sera engage.                                    ',/,&
'@                                                            ',/,&
'@  Le solveur iteratif choisi peut ne pas converger sur le   ',/,&
'@    systeme lineaire resultant du type de probleme considere',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2420 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    RISQUE DE PERTE D''INFORMATION EN CALCUL SUITE          ',/,&
'@                                                            ',/,&
'@  Le calcul sera engage.                                    ',/,&
'@                                                            ',/,&
'@  Un modele de LES a ete active par ITURB(IPHAS) = ',I10     ,/,&
'@    mais on a desactive l''ecriture ou la lecture du fichier',/,&
'@    suite auxiliaire :                                      ',/,&
'@    ILEAUX = ',I10   ,'    IECAUX = ',I10                    ,/,&
'@  Bien que ce fichier ne soit pas necessaire a la poursuite ',/,&
'@    d''un calcul, il contient neanmoins des informations    ',/,&
'@    qui permettent d''eviter les perturbations numeriques   ',/,&
'@    au moment des suites.                                   ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2500 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A -1, 0, 1 OU 2         ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2510 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN REEL POSITIF                        ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2511 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN REEL POSITIF                        ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2520 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@     DTMIN DOIT ETRE SUPERIEUR OU EGAL A 0. ET              ',/,&
'@                     INFERIEUR OU EGAL A DTMAX              ',/,&
'@      ICI DTMIN = ',E14.5     ,' ET DTMAX = ',E14.5          ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Quand le pas de temps n est pas uniforme et constant,     ',/,&
'@    les reels DTMIN et DTMAX bornent ses variations.        ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2530 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    CDTVAR(',I10   ,') DOIT ETRE UN REEL STRICTEMENT POSITIF',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  CDTVAR(I) est le coefficient multiplicatif applique au pas',/,&
'@    de temps pour la resolution de la variable I.           ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2540 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@  ON DEMANDE UNE LIMITATION DU PAS DE TEMPS LIEE AUX EFFETS ',/,&
'@    DE DENSITE (IPTLRO = ',I10   ,') AVEC UNE OPTION DE     ',/,&
'@    PAS DE TEMPS FIXE (IDTVAR = ',I10   ,')                 ',/,&
'@                                                            ',/,&
'@  Le calcul sera engage, mais le pas de temps ne sera pas   ',/,&
'@    clippe. Le code indiquera juste les eventuels           ',/,&
'@    depassements de critere local.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2541 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@  ON DEMANDE UNE LIMITATION DU PAS DE TEMPS LIEE AUX EFFETS ',/,&
'@    DE DENSITE (IPTLRO = ',I10   ,') AVEC UN ALGORITHME     ',/,&
'@    STATIONNAIRE (IDTVAR = ',I10   ,')                      ',/,&
'@                                                            ',/,&
'@  Le calcul sera engage, l''option IPTLRO ignoree.          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2600 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER EGAL A 0, 10, 20, 21, 30, 31,',/,&
'@    40, 41, 42, 50, OU 60'                                   ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2601 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE MODELE DE TURBULENCE K-OMEGA SST N''EST PAS          ',/,&
'@     COMPATIBLE AVEC LE LAGRANGIEN EN COUPLAGE INVERSE      ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2605 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LA METHODE DES VORTEX NE PEUT ETRE ACTIVEE QUE POUR UN  ',/,&
'@    CALCUL MONOPHASIQUE.                                    ',/,&
'@    ON A ICI                                                ',/,&
'@    NPHAS =',I10                                             ,/,&
'@    IVRETX=',I10                                             ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2606 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LA METHODE DES VORTEX NE CONCERNE QUE LES CALCULS EN LES',/,&
'@    ON A ICI                                                ',/,&
'@    ITURB(1)=',I10                                           ,/,&
'@    IVRETX  =',I10                                           ,/,&
'@                                                            ',/,&
'@  Le calcul sera execute en ignorant le mot clef IVRTEX.    ',/,&
'@    (il est repositione a 0)                                ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2607 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ON DEMANDE LA REDUCTION DU VOISINAGE ETENDU POUR LE     ',/,&
'@    CALCUL DES GRADIENTS PAR MOINDRE CARRES, POUR UN        ',/,&
'@    CALCUL EN L.E.S AVEC LE MODELE DYNAMIQUE.               ',/,&
'@    MODELE DYNAMIQUE.                                       ',/,&
'@      ITURB =',I10                                           ,/,&
'@      IMRGRA=',I10                                           ,/,&
'@                                                            ',/,&
'@  Le calcul sera engage.                                    ',/,&
'@                                                            ',/,&
'@  Le calcul de la moyenne locale de la constante de         ',/,&
'@    Smagorinsky dynamique peut etre degrade.                ',/,&
'@                                                            ',/,&
'@  Il est conseille :                                        ',/,&
'@    - soit d''utiliser une methode de calculs des gradients ',/,&
'@      par moindre carres sur voisinage etendu complet       ',/,&
'@      (IMRGRA = 2)                                          ',/,&
'@    - soit de calculer sa propre moyenne de la constante    ',/,&
'@      dynamique dans la subroutine USSMAG.                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2610 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    DONNEES INCOHERENTES                                    ',/,&
'@    ',A31,I10                                                ,/,&
'@    ',A31,I10                                                ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
! 2620 FORMAT
!     & (' ATTENTION : ON DEMANDE ',A6,' = ',I10,/,
!     &  '                  AVEC GRAVITE = ',3E14.5)
 2621 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@  ON DEMANDE LA PRISE EN COMPTE DE L''ACCELERATION DE LA    ',/,&
'@    PESANTEUR ',3E14.5                                       ,/,&
'@    SANS RESOLUTION D''UNE VARIABLE TEMPERATURE OU ENERGIE  ',/,&
'@    (ISCALT = ',I10   ,')                                   ',/,&
'@                                                            ',/,&
'@  Le calcul sera engage.                                    ',/,&
'@                                                            ',/,&
'@  Il n''y a pas d''incompatibilite a prendre en compte      ',/,&
'@    l''acceleration de la pesanteur sans effets thermiques, ',/,&
'@    mais, souvent, la masse volumique depend essentiellement',/,&
'@    de la temperature et la combinaison des options du      ',/,&
'@    present calcul est caracteristique d''un oubli.         ',/,&
'@                                                            ',/,&
'@  Il est conseille de verifier les parametres donnes via    ',/,&
'@  l''interface ou usini1.                                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2622 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@  ON DEMANDE LA PRISE EN COMPTE DES TERMES DE GRAVITE DANS  ',/,&
'@    LES EQUATIONS DE LA TURBULENCE (',A6,' = ',I10   ,')    ',/,&
'@    SANS RESOLUTION D''UNE VARIABLE TEMPERATURE OU ENERGIE  ',/,&
'@                                                            ',/,&
'@  Le calcul sera engage.                                    ',/,&
'@                                                            ',/,&
'@  Si des effets de gravite sont recherches, il convient de  ',/,&
'@    s assurer que la masse volumique est variable.          ',/,&
'@  Le nombre de Prandtl turbulent sera pris egal a 1.        ',/,&
'@  Elle peut varier en fonction d autres grandeurs que       ',/,&
'@    la temperature ou l enthalpie ; si c est le cas, ce     ',/,&
'@    message pourra etre ignore ; sinon, verifier usini1     ',/,&
'@    ou imposer une variation de masse volumique dans usphyv.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2623 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@ Le coefficient RELAXV des variables de la turbulence a ete ',/,&
'@ modifie alors que IKECOU ne vaut pas 0. Il vaut ici        ',/,&
'@ - pour k                  : ',E12.4                         ,/,&
'@ - pour epsilon (ou omega) : ',E12.4                         ,/,&
'@                                                            ',/,&
'@ La modification sera sans effet (RELAXV n''est utile que   ',/,&
'@ si IKECOU=0)                                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2624 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@ Le coefficient RELAXV des variables de la turbulence doit  ',/,&
'@ etre un reel compris entre 0 et 1. Il vaut ici :           ',/,&
'@ - pour k                  : ',E12.4                         ,/,&
'@ - pour epsilon (ou omega) : ',E12.4                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2625 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@ Le coefficient RELAXV de la pression doit etre un reel     ',/,&
'@ compris entre 0 et 1. Il vaut ici : ',E12.4                 ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2630 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@  On demande la prise en compte de l''amortissement de      ',/,&
'@    Van Driest (IDRIES(IPHAS) = ', I10,')'                   ,/,&
'@    avec un modele LES incompatible (ITURB(IPHAS) = ',I10,')',/,&
'@    (modele dynamique et modele WALE)                       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas realise.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2640 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN REEL INCLUS DANS L''INTERVALLE [0;1]',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2650 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    ',A6,'(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  ICPSYR(I) est l indicateur de couplage du scalaire I avec ',/,&
'@    SYRTHES.                                                ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2660 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Aucun couplage avec SYRTHES n''a ete defini.              ',/,&
'@  Le nombre de scalaires couples est cependant ',I10         ,/,&
'@    (Le nombre de scalaires total est ici ',I10   ,')       ',/,&
'@  Verifier le couplage SYRTHES-Noyau.                       ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2661 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Un couplage avec SYRTHES a ete defini.                    ',/,&
'@  Le couplage avec SYRTHES necessite de disposer d''un      ',/,&
'@    scalaire couple (et d''un seul).                        ',/,&
'@  Le nombre de scalaires total   est ici ',I10               ,/,&
'@  Le nombre de scalaires couples est ici ',I10               ,/,&
'@  Verifier le couplage SYRTHES-Noyau.                       ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2662 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Un couplage avec SYRTHES a ete defini.                    ',/,&
'@  Si le code est couple a SYRTHES, le scalaire couple doit  ',/,&
'@    etre la temperature.                                    ',/,&
'@  Le scalaire couple est ici le scalaire ',I10               ,/,&
'@    Ce n''est pas une temperature car                       ',/,&
'@                    ISCSTH(',I10   ,') = ',I10               ,/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@  Pour coupler un scalaire qui n''est pas la temperature,   ',/,&
'@    contacter l''equipe de developpement.                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2663 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Un couplage avec SYRTHES a ete defini.                    ',/,&
'@  En compressible, si le code est couple a SYRTHES, le      ',/,&
'@    scalaire couple doit etre l''energie.                   ',/,&
'@  Le scalaire couple est ici le scalaire ',I10               ,/,&
'@    Ce n''est pas l''energie car                            ',/,&
'@                    ISCSTH(',I10   ,') = ',I10               ,/,&
'@                                                            ',/,&
'@  Contacter l''equipe de developpement.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2664 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    CHOIX DU CALCUL DES ESTIMATEURS D''ERREUR               ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  L''indicateur IESCAL relatif a                            ',/,&
'@    l''estimateur d''erreur numero IEST = ',I10   ,' pour   ',/,&
'@    Navier-Stokes doit etre un entier egal a 0, 1 ou 2.     ',/,&
'@  Il vaut ici : IESCAL(',I10  ,') = ',I10                    ,/,&
'@                                                            ',/,&
'@  Rq. : les valeurs possibles de IEST sont :                ',/,&
'@        IESPRE = ',I10                                       ,/,&
'@        IESDER = ',I10                                       ,/,&
'@        IESCOR = ',I10                                       ,/,&
'@        IESTOT = ',I10                                       ,/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2700 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    CHOIX DU MODE DE CALCUL DE LA DISTANCE A LA PAROI       ',/,&
'@                                                            ',/,&
'@  ICDPAR DOIT ETRE UN ENTIER EGAL A -2, -1, 1 ou 2          ',/,&
'@  IL VAUT ICI ', I10                                         ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2710 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN REEL INCLUS DANS LINTERVALLE [0.;1.]',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2720 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN REEL SUPERIEUR OU EGAL A 1.         ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2730 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN REEL EGAL A 0. OU A 1.              ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2740 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN REEL STRICTEMENT POSITIF.           ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2750 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER INFERIEUR OU EGAL A 1        ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ',A6,'(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0 OU 1    ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@  L algorithme de resolution multigrille algebrique         ',/,&
'@  n est pas compatible avec les variables convectees.       ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ',A6,'(',I10   ,') DOIT ETRE UN ENTIER                  ',/,&
'@      STRICTEMENT POSITIF                                   ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE UN ENTIER STRICTEMENT POSITIF          ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ENTREE DES DONNEES                          ',/,&
'@    =========                                               ',/,&
'@    LA VITESSE DE REFERENCE UREF N''A PAS ETE INITIALISEE   ',/,&
'@    OU A ETE MAL INITIALISEE (VALEUR NEGATIVE).             ',/,&
'@    ELLE VAUT ICI ',E14.5                                    ,/,&
'@                                                            ',/,&
'@  Le calcul pourra etre execute si la turbulence est        ',/,&
'@  initialisee a partir d''un fichier suite de calcul ou par ',/,&
'@  la routine usiniv.                                        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    ISCSTH(',I10   ,') DOIT ETRE UN ENTIER EGAL A -1,0,1,2,3',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4320 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    ISCAVR(',I10   ,') DOIT ETRE UN ENTIER                  ',/,&
'@      POSITIF OU NUL ET                                     ',/,&
'@      INFERIEUR OU EGAL A NSCAL = ',I10                      ,/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Si ISCAVR(I) est nul, le scalaire I n est pas une variance',/,&
'@  Si ISCAVR(I) est positif, le scalaire I est une variance :',/,&
'@    il s agit de la variance des fluctuations du scalaire J ',/,&
'@    dont le numero est ISCAVR(I)                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4321 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE SCALAIRE ',A8  ,'EST DEFINI COMME LA FLUCTUATION     ',/,&
'@    DU SCALAIRE ',A8                                         ,/,&
'@    (ISCAVR(',I10   ,') = ',I10   ,'),                      ',/,&
'@    QUI EST LUI-MEME DEFINI COMME UNE FLUCTUATION           ',/,&
'@    (ISCAVR(',I10   ,') = ',I10   ,').                      ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Si ISCAVR(I) est positif, le scalaire I est une variance :',/,&
'@    il s agit de la variance des fluctuations du scalaire J ',/,&
'@    dont le numero est ISCAVR(I) et on a donc forcement     ',/,&
'@    ISCAVR(J) = 0                                           ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4330 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    ICLVFL(',I10   ,') DOIT ETRE UN ENTIER EGAL A 0, 1 OU 2 ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  ICLVFL(I) indique le mode de clipping du scalaire I       ',/,&
'@    lorsqu il s agit d une variance de fluctuations.        ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4331 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    ICLVFL(',I10   ,') N EST UTILISE QUE POUR LES VARIANCES ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@    ALORS QUE LE SCALAIRE N EST PAS UNE VARIANCE.           ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  ICLVFL(I) indique le mode de clipping du scalaire I       ',/,&
'@    lorsqu il s agit d une variance de fluctuations.        ',/,&
'@    Il n est pas utilise pour les autres scalaires.         ',/,&
'@  L utilisateur est invite a ne pas modifier ICLVFL pour    ',/,&
'@    les scalaires qui ne sont pas des variances.            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4340 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    VISLS0(',I10   ,') DOIT ETRE UN REEL POSITIF            ',/,&
'@      QUAND IVISLS(',I10   ,') = ',I10                       ,/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  VISLS0(I) est le coefficient de diffusion moleculaire du  ',/,&
'@    scalaire et doit etre positif quand ivisls est different',/,&
'@    de 1 (dans ce cas, un coefficient variable est donne    ',/,&
'@    dans usphyv).                                           ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4350 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    SIGMAS(',I10   ,') DOIT ETRE UN REEL POSITIF            ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  SIFMAS(I) est le nombre de Prandtl turbulent associe      ',/,&
'@    au scalaire I.                                          ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4360 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    SCAMIN(',I10   ,') VAUT ICI ',E14.5                      ,/,&
'@      AVEC ICLVFL(',I10   ,') = ',I10                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  SCAMIN(I) est la valeur minimale acceptee pour le         ',/,&
'@    scalaire I. Lorsque le scalaire est une variance        ',/,&
'@    (ISCAVR(I) > 0) la valeur de SCAMIN n est prise en      ',/,&
'@    compte que si ICLVFL(I) = 2                             ',/,&
'@  Si l utilisateur souhaite effectivement que le            ',/,&
'@    scalaire I (en fait, une variance) soit limite a SCAMIN ',/,&
'@    (positif) il faut imposer ICLVFL = 2 dans usini1.       ',/,&
'@  Si l utilisateur souhaite utiliser l option ICLVFL = 1    ',/,&
'@    il est invite a ne pas modifier SCAMIN dans usini1.     ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4361 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    SCAMAX(',I10   ,') VAUT ICI ',E14.5                      ,/,&
'@      AVEC ICLVFL(',I10   ,') = ',I10                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  SCAMAX(I) est la valeur maximale acceptee pour le         ',/,&
'@    scalaire I. Lorsque le scalaire est une variance        ',/,&
'@    (ISCAVR(I) > 0) la valeur de SCAMAX n est prise en      ',/,&
'@    compte que si ICLVFL(I) = 2                             ',/,&
'@  Si l utilisateur souhaite effectivement que le            ',/,&
'@    scalaire I (en fait, une variance) soit limite a SCAMAX ',/,&
'@    (positif) il faut imposer ICLVFL = 2 dans usini1.       ',/,&
'@  Si l utilisateur souhaite utiliser l option ICLVFL = 1    ',/,&
'@    il est invite a ne pas modifier SCAMAX dans usini1.     ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4370 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    SCAMAX(',I10   ,') DOIT ETRE UN REEL STRICTEMENT POSITIF',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  SCAMAX(I) est la valeur maximale acceptee pour le         ',/,&
'@    scalaire I, ici une variance                            ',/,&
'@  Avec ICLVFL(I) = 2, la valeur de SCAMAX doit donc etre    ',/,&
'@   strictment positive.                                     ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4380 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    RVARFL(',I10   ,') DOIT ETRE UN REEL STRICTEMENT POSITIF',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  RVARFL(I) est le coefficient R pour le scalaire I (qui est',/,&
'@    une variance) intervenant dans le terme de dissipation :',/,&
'@    - (1/R) rho scalaire epsilon/k                          ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LA PERIODICITE DE ROTATION N''EST PAS COMPATIBLE AVEC LE',/,&
'@      COUPLAGE VITESSE PRESSION RENFORCE  OU LA METHODE     ',/,&
'@      ALE DANS LA VERSION COURANTE                          ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,&
'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,&
'@    par IPERIO = ',I10,   ')                                ',/,&
'@    et certaines periodicites sont de rotation.             ',/,&
'@  L''indicateur IPUCOU a ete positionne a ',I10              ,/,&
'@    dans l''interface ou usini1 (couplage renforce pour     ',/,&
'@    IPUCOU=1).                                              ',/,&
'@  L''indicateur IALE a ete positionne a ',I10                ,/,&
'@    dans l''interface ou usini1 (methode activee si IALE=1) ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5005 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    MODE DE CALCUL DE LA DISTANCE A LA PAROI INCOMPATIBLE   ',/,&
'@      AVEC LA PERIODICITE                                   ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,&
'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,&
'@    par IPERIO = ',I10,   ').                               ',/,&
'@  Les parametres de calcul specifies necessitent le calcul  ',/,&
'@    la distance a la paroi (Rij-eps LRR avec echo de paroi, ',/,&
'@    LES avec amortissement de van Driest ou k-omega SST).   ',/,&
'@  Le mode de calcul de la distance a la paroi defini par    ',/,&
'@    ICDPAR = ',I10,   ' ne permet pas de prendre en compte  ',/,&
'@    la periodicite.                                         ',/,&
'@                                                            ',/,&
'@  Utiliser ICDPAR = 1 ou -1.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
! 5007 format(
!     &'@                                                            ',/,
!     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!     &'@                                                            ',/,
!     &'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,
!     &'@    =========                                               ',/,
!     &'@    LA PERIODICITE             N''EST PAS COMPATIBLE AVEC LE',/,
!     &'@      RAYONNEMENT SEMI-TRANSPARENT  (ORDONNEES DISCRETES)   ',/,
!     &'@      DANS LA VERSION COURANTE                              ',/,
!     &'@                                                            ',/,
!     &'@  Le calcul ne peut etre execute.                           ',/,
!     &'@                                                            ',/,
!     &'@  L''indicateur IPERIO a ete positionne a ',I10              ,/,
!     &'@    dans usini1 (periodicite activee pour IPERIO=1).        ',/,
!     &'@  L''indicateur IIRAYO(',I10,') a ete positionne a ',I10     ,/,
!     &'@    dans usray1.                                            ',/,
!     &'@                                                            ',/,
!     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!     &'@                                                            ',/)
 5008 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LA PERIODICITE DE ROTATION N''EST PAS COMPATIBLE AVEC LE',/,&
'@      RAYONNEMENT SEMI-TRANSPARENT  (ORDONNEES DISCRETES)   ',/,&
'@      DANS LA VERSION COURANTE                              ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,&
'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,&
'@    par IPERIO = ',I10,   ')                                ',/,&
'@    et certaines periodicites sont de rotation.             ',/,&
'@  L''indicateur IIRAYO a ete positionne a ',I10              ,/,&
'@    dans l''interface ou usray1.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5009 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    DES DEFAUTS PEUVENT SE RENCONTRER LORS DE L''UTILISATION',/,&
'@      DE LA PERIODICITE DE ROTATION EN RIJ-EPSILON.         ',/,&
'@                                                            ',/,&
'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,&
'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,&
'@    par IPERIO = ',I10,   ')                                ',/,&
'@    et certaines periodicites sont de rotation.             ',/,&
'@  L''indicateur ITURB(IPHAS) a ete positionne a ',I10        ,/,&
'@                                                            ',/,&
'@  Le calcul peut etre execute.                              ',/,&
'@    Les defauts eventuels evoques proviennent de la prise en',/,&
'@    compte de la rotation du tenseur de viscosite orthotrope',/,&
'@    Il a cependant en general une influence faible de sorte ',/,&
'@    que les tests menes jusqu''a present n''ont pas fait    ',/,&
'@    apparaitre de probleme.                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION :       A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    DES DEFAUTS PEUVENT SE RENCONTRER LORS DE L''UTILISATION',/,&
'@      DE LA PERIODICITE DE ROTATION AVEC UN SCHEMA D''ORDRE ',/,&
'@      DEUX EN TEMPS POUR LA VITESSE.                        ',/,&
'@                                                            ',/,&
'@  La variable COMMANDE_PERIO a ete renseignee dans le       ',/,&
'@    lanceur (la periodicite a ete activee, ce qui se traduit',/,&
'@    par IPERIO = ',I10,   ')                                ',/,&
'@    et certaines periodicites sont de rotation.             ',/,&
'@  Les indicateurs THETAV des trois composantes Ux, Uy, Uz   ',/,&
'@    de la vitesse                                           ',/,&
'@    ont ete positionnes (dans usini1 ou par defaut suite aux',/,&
'@    options de calcul selectionnees) aux valeurs suivantes :',/,&
'@    THETAV(IU(IPHAS))  THETAV(IV(IPHAS))  THETAV(IW(IPHAS)) ',/,&
'@        ',E14.5     ,'     ',E14.5     ,'     ',E14.5        ,/,&
'@                                                            ',/,&
'@  Le calcul peut etre execute.                              ',/,&
'@    Les defauts eventuels evoques proviennent de la prise en',/,&
'@    compte de la rotation du tenseur gradient de vitesse.   ',/,&
'@    En effet, on met en oeuvre une methode explicite qui est',/,&
'@    susceptible de faire chuter l''ordre du schema en temps.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 6002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    LE PARALLELISME N''EST PAS COMPATIBLE AVEC LE MODULE    ',/,&
'@      LAGRANGIEN DANS LA VERSION COURANTE                   ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Le processeur courant est de rang ',I10                    ,/,&
'@  L''indicateur IILAGR a ete positionne a ',I10              ,/,&
'@    dans uslag1 (module lagrangien active pour IILAGR>0).   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6005 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    MODE DE CALCUL DE LA DISTANCE A LA PAROI INCOMPATIBLE   ',/,&
'@      AVEC LE PARALLELISME                                  ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Le processeur courant est de rang ',I10                    ,/,&
'@  Les parametres de calcul specifies necessitent le calcul  ',/,&
'@    la distance a la paroi (Rij-eps LRR avec echo de paroi, ',/,&
'@    LES avec amortissement de van Driest ou k-omega SST).   ',/,&
'@  Le mode de calcul de la distance a la paroi defini par    ',/,&
'@    ICDPAR = ',I10,   ' ne permet pas de prendre en compte  ',/,&
'@    le parallelisme.                                        ',/,&
'@                                                            ',/,&
'@  Utiliser ICDPAR = 1 ou -1.                                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    INDICATEUR DE METHODE ALE                               ',/,&
'@                                                            ',/,&
'@  IALE DOIT VALOIR 0 OU 1                                   ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usalin.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    NOMBRE D''ITERATIONS D''INITIALISATION DU FLUIDE EN ALE ',/,&
'@                                                            ',/,&
'@  NALINF DOIT ETRE UN ENTIER POSITIF                        ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usalin.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    COEFFICIENTS DE LA METHODE DE NEWMARK NON VALIDES       ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  ALPNMK doit etre compris entre 0 et 1                     ',/,&
'@  BETNMK doit etre compris entre 0 et 1/2                   ',/,&
'@  GAMNMK doit etre compris entre 0 et 1                     ',/,&
'@  On a ici :                                                ',/,&
'@                                                            ',/,&
'@       ALPNMK      BETNMK      GAMNMK                       ',/,&
'@ ',     E12.4,      E12.4,      E12.4                        ,/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface, usini1   ',/,&
'@    ou usalin.                                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    NOMBRE D''ITERATIONS MAX DE COUPLAGE IMPLICITE EN ALE   ',/,&
'@                                                            ',/,&
'@  NALIMX DOIT ETRE UN ENTIER POSITIF                        ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usalin.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7040 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    PRECISION DU COUPLAGE IMPLICITE EN ALE                  ',/,&
'@                                                            ',/,&
'@  EPALIM DOIT ETRE UN REEL STRICTEMENT POSITIF              ',/,&
'@    IL VAUT ICI ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usalin.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7050 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========                                               ',/,&
'@    ITERATION D''INITIALISATION DE L''ALE                   ',/,&
'@                                                            ',/,&
'@  ITALIN DOIT VALOIR 0 OU 1                                 ',/,&
'@    IL VAUT ICI ',I10                                        ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes dans usalin.               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A L''ENTREE DES DONNEES               ',/,&
'@    =========   MODULE COMPRESSIBLE                         ',/,&
'@                                                            ',/,&
'@    T0 ET P0 DOIVENT ETRE DES REELS STRICTEMENT POSITIFS    ',/,&
'@    ILS VALENT ICI :                                        ',/,&
'@                   T0 = ',E14.5                              ,/,&
'@                   P0 = ',E14.5                              ,/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@  Verifier les parametres donnes via l''interface ou usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    NDIM   must have value equal to 3                       ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@  space dimension can only be 3 even for 2D simulations     ',/,&
'@                                                            ',/,&
'@  Check the input data given via User Interface or in usini1',/,&
'@  Check the mesh file                                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A33,                          ' MUST BE AN INTEGER   ',/, &
'@    SUPERIEUR or EGAL A -1                                  ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1210 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A33,                          ' MUST BE AN INTEGER    ',/,&
'@    LARGER OR EQUAL TO  1 or EQUAL TO  -1                   ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1220 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ICHRVR(',I10   ,') MUST BE AN INTEGER EQUAL  0  OR 1    ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  ICHRVR defines whether the variable should be included in ',/,&
'@    post-processing files                                   ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1230 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    NCAPT  MUST BE AN INTEGER LESS THAN or EGAL A ',I10     ,/, &
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The computation CANNOT  start.                           ',/,&
'@                                                            ',/,&
'@  NCAPT  is the number of probes for history/ time series   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1240 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IHISVR(',I10   ,',1) MUST BE AN INTEGER EGAL to -1      ',/,&
'@      or Zero,    or positive but less than NCAPT =',I10     ,/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@    The calculation could NOT run.                          ',/,&
'@                                                            ',/,&
'@  IHISVR(I,1) is the number of probes for variable  I       ',/,&
'@      (-1 means all probes are used)                        ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1250 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IHISVR(',I10   ,',',I10   ,') MUST BE AN INTEGER        ',/,&
'@    LARGER OR EQUAL TO 1 AND                                ',/,&
'@      LESS OR EQUAL TO NCAPT = ',I10                         ,/,&
'@   IT HAS VALUE =',I10                                       ,/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@  IHISVR(I,j+1) gives the number of the j-ieth probe        ',/,&
'@    to be used with variable number I to post-process       ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1260 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ILISVR(',I10   ,') MUST BE AN INTEGER EQUAL  0 OR  1    ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@  ILISVR(I) tells if variable (I)should be included         ',/,&
'@    in the printed listing                                  ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1270 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ITRSVR(',I10   ,') MUST BE AN INTEGER                   ',/,&
'@    BETWEEN  0 and  NVAR=',I10                               ,/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 1300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@      PARAMETER ',A6,' MUST BE AN INTEGER MULTIPLE OF       ',/,&
'@      THE FOLLOWING INTEGER :                               ',/,&
'@      ',A6,' = ',I10                                         ,/,&
'@      ',A6,' = ',I10                                         ,/,&
'@      ',A6,' = ',I10                                         ,/,&
'@      ',A6,' = ',I10                                         ,/,&
'@   IT HAS VALUE ',A6,' = ', I10                              ,/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@  This parameter tells which extra variables should be      ',/,&
'@    included for post-processing                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER                               ',/,&
'@      STRICTLY POSITIVE AND LESS THAN or EGAL TO',I10        ,/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2005 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    TEMPORAL EXTRAPOLATION OF DENSITY RHO REQUESTED,        ',/,&
'@      BUT IROEXT(IPHAS) = ',I10                              ,/,&
'@    THIS IS INCOMPATIBLE WITH RHO = CONSTANT                ',/,&
'@      IROVAR(IPHAS) = ',I10                                  ,/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@     PARAMETERS OF NUMERICAL SCHEME FOR VARIABLE     ',A8    ,/,&
'@                                                            ',/,&
'@ Parameter               ISTAT     ICONV                    ',/,&
'@ Values   accepted     0 or  1   0 or  1                    ',/,&
'@ Values entered here',I10     ,I10                           ,/,&
'@                                                            ',/,&
'@ PARAMETER                         IDIFF     IDIFFT         ',/,&
'@ Values   accepted               0 or  1   0 or  1          ',/,&
'@ Values entered here',10X     ,I10      ,I10                 ,/,&
'@                                                            ',/,&
'@ PARAMETER               THETAV   BLENCV    ISCHCV    ISSTPC',/,&
'@ Values   accepted     [0.; 1.] [0.; 1.]   0 or  1   0 or  1',/,&
'@ Values entered here      ',E14.5     ,E14.5,I10  ,I10       ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2111 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@   INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME           ',/,&
'@                                                            ',/,&
'@   TIME DISCRETISATION SCHEME for velocity                  ',/,&
'@      THETA DOES NOT HAVE SAME VALUE FOR ALL 3 COMPONENTS   ',/,&
'@                                                            ',/,&
'@ Parameter THETAV              U          V          W      ',/,&
'@ Values  entered here',E10.2,E10.2,E10.2                     ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2112 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME            ',/,&
'@                                                            ',/,&
'@  PARAMETER THETAV FOR PRESSURE MUST BE EQUAL TO    1       ',/,&
'@                                                            ',/,&
'@  It has value ',E14.5                                       ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2113 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME            ',/,&
'@                                                            ',/,&
'@  IN L.E.S.                                                 ',/,&
'@    PARAMETER THETFL FOR TIME DISCRETISATION FOR MASS FLUX  ',/,&
'@     IS DIFFERENT FROM THAT OF PHASE 1 ',E14.5               ,/,&
'@     THETFL WAS SET AT         ',E14.5                       ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2114 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME            ',/,&
'@                                                            ',/,&
'@   PARAMETER ISTMPF IN SCHEME EN TEMPS POUR LE FLUX DE      ',/,&
'@     MASSE  IS DIFFERENT FROM THAT  OF  PHASE 1',I10        ,/, &
'@     ISTMPF IS NOW IMPOSED AS  ',I10                         ,/,&
'@                                                            ',/,&
'@  Computation will NOT proceed                              ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2121 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING  :      WHEN READING INPUT DATA                ',/&
'@    =========                                               ',/,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME                    ',/,&
'@                                                            ',/,&
'@   IN L.E.S.                                                ',/,&
'@   RECOMMENDED VALUE FOR PARAMETER THETAV IN TIME-SCHEME    ',/,&
'@   FOR VARIABLE               ',A8  ,' IS 0.5               ',/,&
'@     THETAV IS NOW IMPOSED AT  ',E14.5                       ,/,&
'@                                                            ',/,&
'@  computation will go on                                    ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@  user interface or usini1.                                 ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2122 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME                    ',/,&
'@                                                            ',/,&
'@   WITH L.E.S.                                              ',/,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER BLENCV  ',/,   &
'@   FOR CONVECTION OF VARIABLE ',A8  ,' EST 1.0              ',/,&
'@     BLENCV IS NOW IMPOSED AS ',E14.5                        ,/,&
'@                                                            ',/,&
'@  computation will go on                                    ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2123 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME                    ',/,&
'@                                                            ',/,&
'@   EN L.E.S.                                                ',/,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER ISSTPC IN SCHEME ',/,&
'@     CONVECTION OF VARIABLE ',A8  ,' IS    1                ',/,&
'@     ISSTPC IS NOW IMPOSED AS  ',I10                         ,/,&
'@                                                            ',/,&
'@  Computation will go on                                    ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2124 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME                    ',/,&
'@                                                            ',/,&
'@   EN L.E.S.                                                ',/,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER ISSTPC IN SCHEME ',/,&
'@     CONVECTION OF VARIABLE   ',A8  ,' IS  0                ',/,&
'@     ISSTPC IS NOW IMPOSED AS  ',I10                         ,/,&
'@                                                            ',/,&
'@  Computation will go on                                    ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2125 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME                    ',/,&
'@                                                            ',/,&
'@   ORDRE 2 EN TEMPS or the                                  ',/,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER NSWRSM FOR       ',/,&
'@        VARIABLE ',A8  ,' IS    ',I10                        ,/,&
'@     NSWRSM IS NOW IMPOSED AS  ',I10                         ,/,&
'@                                                            ',/,&
'@  computation will go on                                    ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2127 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME                    ',/,&
'@                                                            ',/,&
'@   EN L.E.S.                                                ',/,&
'@   THE VALUE RECOMMANDED FOR THE PARAMETER BLENCV IN        ',/,&
'@     CONVECTION OF VARIABLE   ',A8  ,' IS  1.0              ',/,&
'@     BLENCV IS NOW IMPOSED AS  ',E14.5                       ,/,&
'@                                                            ',/,&
'@  Computation will NOT proceed                              ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2131 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@   CHOICE OF TIME-SCHEME                                    ',/,&
'@                                                            ',/,&
'@     TIME-SCHEME FOR VELOCITY IS FIRST ORDER                ',/,&
'@       (THETAV = ',E10.2 ,')                                ',/,&
'@     CERTAIN TERMES ARE HOWEVER SECOND ORDER IN TIME WITH   ',/,&
'@       THE FOLLOWING SETTINGS:                              ',/,&
'@                                                            ',/,&
'@ parameters       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT ',/,&
'@ Values  entered ',6I7                                       ,/,&
'@                                                            ',/,&
'@  computation will go on.                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2132 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@   CHOICE OF TIME-SCHEME                                    ',/,&
'@                                                            ',/,&
'@     TIME-SCHEME FOR VELOCITY IS SECOND ORDER               ',/,&
'@       (THETAV = ',E10.2 ,')                                ',/,&
'@     CERTAIN TERMES ARE HOWEVER FIRST ORDER IN TIME  WITH   ',/,&
'@       THE FOLLOWING SETTINGS:                              ',/,&
'@                                                            ',/,&
'@ parameters       ISTMPF ISNO2T ISTO2T IROEXT IVIEXT ICPEXT ',/,&
'@ Values  entered ',6I7                                       ,/,&
'@                                                            ',/,&
'@  computation will go on.                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2133 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME                    ',/,&
'@                                                            ',/,&
'@   SCALAR   ',I10,' ISSO2T = ',I10                           ,/,&
'@       IS DIFFERENT FROM ISNO2T                             ',/,&
'@     ISNO2T(IPHAS) = ',I10                                   ,/,&
'@                                                            ',/,&
'@  computation will go on                                    ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2134 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@   NON-STANDARD CHOICE WITH  TIME-SCHEME                    ',/,&
'@                                                            ',/,&
'@   SCALAIRE ',I10,' IVSEXT = ',I10                           ,/,&
'@       IS DIFFERENT FROM IVIEXT                             ',/,&
'@     IVIEXT(IPHAS) = ',I10                                   ,/,&
'@                                                            ',/,&
'@  computation will go on                                    ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2135 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME            ',/,&
'@                                                            ',/,&
'@     Specific heat is extrapolated in time with             ',/,&
'@       ICPEXT(IPHAS) = ',I10                                 ,/,&
'@    in which case it should be variable, or                 ',/,&
'@       ICP(IPHAS)    = ',I10                                 ,/,&
'@                                                            ',/,&
'@  Computation will NOT go on                                ',/,&
'@                                                            ',/,&
'@  Verify   the parameters given  via l''interface or usini1 ',/,&
'@    - deactivate xtrapolation of Cp in time                 ',/,&
'@      or                                                    ',/,&
'@    - define Cp as variable                                 ',/,&
'@         (give its variation law in   usphyv)               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2136 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME            ',/,&
'@                                                            ',/,&
'@   Scalar   ISCAL = ',I10                                    ,/,&
'@    Diffusivity   is extrapolated in time wih               ',/,&
'@       IVSEXT(ISCAL) = ',I10                                 ,/,&
'@     it should thus be a variable, or                       ',/,&
'@       IVISLS(ISCAL) = ',I10                                 ,/,&
'@                                                            ',/,&
'@ Computation will  NOT  proceed                             ',/,&
'@                                                            ',/,&
'@  Verify the parameters given via user interface or usini1  ',/,&
'@    - deactivate intepolation in time                       ',/,&
'@                                     for diffusivity        ',/,&
'@      or                                                    ',/,&
'@    - impose diffusivite variable                           ',/,&
'@         (and describe variation law in usphyv)             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2137 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHEN READING INPUT DATA                ',/,&
'@    =========                                               ',/,&
'@   CHOICE INCOMPATIBLE FOR ERROR ESTIMATES                  ',/,&
'@                                                            ',/,&
'@  One or several error estimates are activated for          ',/,&
'@    Navier-Stokes in simulation with frozen velocity field. ',/,&
'@    estimates will not be computed                          ',/,&
'@                                                            ',/,&
'@ Computation will  NOT  proceed                             ',/,&
'@                                                            ',/,&
'@  Verify   the parameters given  via l''interface or        ',/,&
'@    usini1 :                                                ',/,&
'@      desactivate  ERROR ESTIMATES  (ICCVFG)                ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2140 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME   ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  A second ordre time-scheme was requested:                 ',/,&
'@      U,V,W : THETA = ',3E12.4                               ,/,&
'@     Source terme in Navier-Stokes: ISNO2T = ',I10           ,/,&
'@                                    THETSN = ',E12.4         ,/,&
'@      Density                     : IROEXT = ',I10           ,/,&
'@                                    THETRO = ',E12.4         ,/,&
'@      Viscosity                    : IVIEXT = ',I10          ,/,&
'@                                    THETVI = ',E12.4         ,/,&
'@  Current version does not allow this in combination with   ',/,&
'@  one of the following option (which has been activated ):  ',/,&
'@    - Error estimation (IESCAL)                             ',/,&
'@    - reinforced U-P coupling (IPUCOU)                      ',/,&
'@    - specific treatment of hydrostatic pressure            ',/,&
'@      contribution  (IPHYDR et ICALHY)                      ',/,&
'@    - time-step variable with space or iteration or         ',/,&
'@      steady-state   algorithm(IDTVAR)                      ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2141 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME   ',/,&
'@                                                            ',/,&
'@  Pressure-Velocity coupling by fixed-point method          ',/,&
'@    was selected  by setting   NTERUP = ',I10                ,/,&
'@  Current version does not allow this in combination with   ',/,&
'@  one of the following options (which has been activated ): ',/,&
'@    - Error estimation (IESCAL)                             ',/,&
'@    - reinforced U-P coupling (IPUCOU)                      ',/,&
'@    - specific treatment of hydrostatic pressure            ',/,&
'@      contribution  (IPHYDR et ICALHY)                      ',/,&
'@    - time-step variable with space or iteration or         ',/,&
'@      steady-state   algorithm(IDTVAR=-1)                   ',/,&
'@    - compressible module (IPPMOD(ICOMPF)>=0)               ',/,&
'@    - frozen velocity field (ICCVFG=1)                      ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2142 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME   ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  With the k-epsilon turbulence                             ',/,&
'@  model (ITURB = ',I10 ,') solved coupled(IKECOU = ',I10,'):',/,&
'@    the current version does not allow second order         ',/,&
'@    in time resolution of k-epsilon equations in a coupled  ',/,&
'@    manner.                                                 ',/,&
'@                                                            ',/,&
'@   Thus one or more of the values below are not permited    ',/,&
'@                                                            ',/,&
'@       THETST    ISTO2T     THETA K   THETA EPS             ',/,&
'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2143 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME   ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  With the   V2F     turbulence                             ',/,&
'@  model (ITURB = ',I10 ,') solved coupled(IKECOU = ',I10,'):',/,&
'@    the current version does not allow second order         ',/,&
'@    in time resolution of V2F equations in a coupled        ',/,&
'@    manner.                                                 ',/,&
'@                                                            ',/,&
'@   Thus one or more of the values below are not permited    ',/,&
'@                                                            ',/,&
'@       THETST    ISTO2T     THETA K   THETA EPS             ',/,&
'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,&
'@     THETA PHI    THETA FB                                  ',/,&
'@ ',      E12.4,      E12.4                                   ,/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2144 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME   ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  With the k-omega   turbulence                             ',/,&
'@  model (ITURB = ',I10 ,') solved coupled(IKECOU = ',I10,'):',/,&
'@    the current version does not allow second order         ',/,&
'@    in time resolution of k-omega equations in a coupled    ',/,&
'@    manner.                                                 ',/,&
'@                                                            ',/,&
'@   Thus one or more of the values below are not permited    ',/,&
'@                                                            ',/,&
'@       THETST    ISTO2T     THETA K   THETA OMEGA           ',/,&
'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2145 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME   ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  With the Spalart-Allmaras                                 ',/,&
'@  turbulence model (ITURB = ',I10 ,')                       ',/,&
'@    the current version does not allow second order         ',/,&
'@    in time resolution.                                     ',/,&
'@                                                            ',/,&
'@   Thus one or more of the values below are not permited    ',/,&
'@                                                            ',/,&
'@       THETST    ISTO2T     THETA K   THETA OMEGA           ',/,&
'@ ',     E12.4,      I10,      E12.4,      E12.4              ,/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2146 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME   ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  The current version does not allow changing the time      ',/,&
'@  discretisation scheme when a specific physics is active   ',/,&
'@    (combustion, coal, electrical, ...).                    ',/,&
'@                                                            ',/,&
'@  Verify   the parameters given  via interface, usini1      ',/,&
'@  and  usppmo.                                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2147 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME   ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  The source terms stemming from the                        ',/,&
'@    Lagrangien module will not be computed as second order  ',/,&
'@    in this version, despite user settings chosen below     ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@     For Navier-Stokes       For   turbulence               ',/,&
'@       THETSN    ISNO2T      THETST    ISTO2T               ',/,&
'@ ',     E12.4,      I10,      E12.4,      I10                ,/,&
'@                                                            ',/,&
'@  (other termes sources could be second order in time       ',/,&
'@             )                                              ',/,&
'@                                                            ',/,&
'@  Verify   the parameters given  via interface, usini1      ',/,&
'@  and uslag1.                                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2148 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    OPTIONS NOT COMPATIBLE WITH TIME DICRETISATION SCHEME   ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  source terms coming from module ',A11                      ,/,&
'@    will not be computed as second order                    ',/,&
'@    in this version, despite user settings chosen below     ',/,&
'@                                                            ',/,&
'@       Pour le scalaire ',I10                                ,/,&
'@       THETSS    ISSO2T                                     ',/,&
'@ ',     E12.4,      I10                                      ,/,&
'@                                                            ',/,&
'@  (Les autres termes sources pourraient etre traites a      ',/,&
'@   l''ordre 2)                                              ',/,&
'@                                                            ',/,&
'@  Verify   the parameters given  via l''interface, usini1   ',/,&
'@  et ',A6                                                    ,/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2149 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@     ALGORITHME STATIONNAIRE                                ',/,&
'@     COEFFICIENT DE RELAXATION POUR LA VARIABLE ',A8         ,/,&
'@                                                            ',/,&
'@ RELAXV MUST BE A  REAL comprised between    0 et 1         ',/,&
'@ it has here a value of ',E14.5                              ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2150  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@   INCOMPATIBILITE ALGORITHME STATIONNAIRE                  ',/,&
'@                                                            ',/,&
'@   relaxation Coefficient of velocity  phase ',I10           ,/,&
'@      RELAXV should have the same value for all components  ',/,&
'@                                                            ',/,&
'@ PARAMETER RELAXV              U          V          W      ',/,&
'@ Values entered here ',E10.2,E10.2,E10.2                     ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2151  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@   L''ALGORITHME STATIONNAIRE N''EST PAS COMPATIBLE AVEC    ',/,&
'@    LE MODULE LAGRANGIEN QUI EST UNE APPROCHE INSTATIONNAIRE',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@   Integer parameter IILAGR was set to    ',I10              ,/,&
'@    in uslag1 (lagrangian module active for IILAGR>0).      ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2152  format(                                                          &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@   STEADY STATE SOLVER SCHEME IS NOT COMPATIBLE WITH        ',/,&
'@   L.E.S. TRUBULENCE MODELLING WHICH IS TIME-DEPENDENT      ',/,&
'@    BY NATURE                                               ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  Integer parameter ITURB was set to ',I10                   ,/,&
'@      in usini1.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2200 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL  0  OR 1                ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2201 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL  0  OR 1                ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2202 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    HYDROSTATIQUE PRESSURE ONLY VALID FOR SINGLE PHASE FLOW ',/,&
'@    IPHYDR HAS VLAUE ',I10                                   ,/,&
'@    AND NPHAS        ',I10                                   ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
! 2203 format(
!     &'@                                                            ',/,
!     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!     &'@                                                            ',/,
!     &'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,
!     &'@    =========                                               ',/,
!     &'@  GRAVITY IS SET TO ZERO    ',/,
!     &'@  THIS IS NOT COMPATIBLE WITH HYDROSTATIQUE PRESSURE OPTION ',/,
!     &'@                                                            ',/,
!     &'@                                                            ',/,
!     &'@  ICALHY IS EQUAL ',I10                                      ,/,
!     &'@         IS EQUAL ',3E14.5                                   ,/,
!     &'@                                                            ',/,
!     &'@  Computation CAN NOT run                                   ',/,
!     &'@                                                            ',/,
!     &'@ Check the input data given via User Interface or in usini1.',/,
!     &'@                                                            ',/,
!     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!     &'@                                                            ',/)
 2204 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@  INCOMPATIBILITY FOR TIME DISCRETISATION SCHEME            ',/,&
'@                                                            ',/,&
'@  ON DEMANDE LA PRISE UN SCHEMA EN TEMPS FOR    VELOCITY    ',/,&
'@    D''ORDRE 2 AVEC UN PAS DE TEMPS NON CONSTANT, UN        ',/,&
'@    ALGORITHME STATIONNAIRE or the MATRICES POIDS           ',/,&
'@                                                            ',/,&
'@  THETAV IS EQUAL ',E14.5,' FOR    VELOCITY                 ',/,&
'@  ALORS QUE IDTVAR VAUT ',I10                                ,/,&
'@  ET IPUCOU VAUT        ',I10                                ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2205 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EQUAL TO 0, 1, 2, 3 or 4      ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2206 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ANOMAX DOIT ETRE A REAL POSITIVE or NUL AND             ',/,&
'@                             LESS THAN or EQUAL TO PI/2     ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  On demande la reconstruction des gradients par moindres   ',/,&
'@    carres sur voisinage etendu reduit (IMRGRA = ',I10  ,').',/,&
'@    Le critere est base sur l''angle de non orthogonalite   ',/,&
'@    des faces ANOMAX qui doit etre fourni en radians et     ',/,&
'@    compris dans the bornes indiquees ci-dessus.            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2207 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    L''UTILISATION DE LA METHODE DE CALCUL DE GRADIENT PAR  ',/,&
'@      MOINDRES CARRES EST IMPOSSIBLE AVEC                   ',/,&
'@        NPHAS SUPERIEUR A 1   or                            ',/,&
'@        EXTRAG(IPR(1)) DIFFERENT DE 0 ET 1                  ',/,&
'@                                                            ',/,&
'@    ON A ICI                                                ',/,&
'@        IMRGRA         = ',I10                               ,/,&
'@        NPHAS          = ',I10                               ,/,&
'@        EXTRAG(IPR(1)) = ',E14.5                             ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2208 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    EN K-EPS PROD LIN (ITURB=21) ET EN V2F (ITURB=50)       ',/,&
'@    IKECOU DOIT ETRE EGAL A 0                               ',/,&
'@    ITURB  IS EQUAL ',I10                                    ,/,&
'@    IKECOU IS EQUAL ',I10                                    ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2209 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    2u*scales version of WALL FUNCTION (IDEUCH=1 or 2)      ',/,&
'@    EST INCOMPATIBLE AVEC UN CALCUL EN LAMINAIRE, EN        ',/,&
'@    LONGUEUR DE MELANGE or EN L.E.S.                        ',/,&
'@    ON A ICI ITURB=',I10                                     ,/,&
'@         ET IDEUCH=',I10                                     ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2210 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SOLVE STEADY-STATE EQN. OPTION IS NOT COMPATIBLE WITH   ',/,&
'@    COUPLING OF SOURCES TERMES IN K-EPS, V2F or K-OMEGA    ',/, &
'@                                                            ',/,&
'@    WE HAVE  IKECOU=',I10                                    ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2211 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EGAL A 0, 1 or 2             ',/, &
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IMLIGR(',I10   ,') MUST BE AN INTEGER                  ',/, &
'@      LESS THAN or EGAL A 1                                 ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  IMLIGR(I) indique le mode de limitation des gradients     ',/,&
'@    pour la variable I                                      ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2310 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IRCFLU(',I10   ,') MUST BE AN INTEGER EQUAL  0  OR 1    ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  IRCFLU(I) flags if fluxes are reconstructed for            ',/&
'@             variable I                                      ',/&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2311 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IRCFLU(',I10   ,') = ',I10   ,' IS  INCOMPATIBLE WITH   ',/,&
'@    ISCHCV(',I10   ,') = ',I10                               ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  IRCFLU(I) = 0 fluxes are not reconstructed for variable I ',/,&
'@                            ',/,                          &
'@  ISCHCV(I) = 0 (schema SOLU) requests a  reconstruction    ',/,&
'@   for convective fluxes                                    ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2320 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    CLIMGR(',I10   ,') DOIT ETRE A  REAL                    ',/,&
'@      SUPERIEUR or EGAL A 1                                 ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  CLIMGR(I) is the coefficient limiting the gradients       ',/,&
'@    for variable I                                          ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2330 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    EXTRAG(',I10   ,')  MUST BE REAL and equal 0 or 1       ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2331 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    EXTRAG(',I10   ,') MUST BE ZERO                         ',/,&
'@      (non zero values are only possible for pressure       ',/,&
'@                                                            ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2400 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IRESOL(',I10   ,') MUST BE AN INTEGER EQUAL             ',/,&
'@                               to -1 or to     IPOL*1000+ J ',/,&
'@ where IPOL is the order of the preconditionning polynomial ',/,&
'@      and  J    = 0 for conjugate grandient                 ',/,&
'@                = 1 for  JACOBI   (IPOL = 0 DANS CE CAS)    ',/,&
'@                = 2 for  BI-CGSTAB                          ',/,&
'@                = 3 for  GMRES                              ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  IRESOL(I) is the linear system reso methode to use        ',/,&
'@    for  variable I                                         ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2401 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    IDIRCL(',I10   ,') MUST BE AN INTEGER EQUAL  0  OR 1    ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  IDIRCL(I) tells if the diagonal of the matrix for variable',/,&
'@  I should be shifted in the absence of Dirichlet condition ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2410 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    RESOLUTION OF LINEAR SYSTEM COULD FAIL                  ',/,&
'@    IRESOL(',I10   ,') = ',I10                               ,/,&
'@      AND THE VARIABLE IS ADVECTED        (ICONV = ',I10,') ',/,&
'@                                                            ',/,&
'@  The calculation will be launched nevertheless             ',/,&
'@                                                            ',/,&
'@  The chosen linear solver could fail to converge           ',/,&
'@    because of the nature of the problem                    ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2420 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    RISQUE DE PERTE D''INFORMATION EN CALCUL SUITE          ',/,&
'@                                                            ',/,&
'@  The calculation will run.                                  ',/&
'@                                                            ',/,&
'@ A turbulence model was activated by ITURB(IPHAS) = ',I10    ,/,&
'@    but writing to auxiliary restart file was de-activated  ',/,&
'@                                                            ',/,&
'@    ILEAUX = ',I10   ,'    IECAUX = ',I10                    ,/,&
'@  Although this file is not necessary to restart            ',/,&
'@  a computation, it does contain information that avoid     ',/,&
'@   numerical perturbations when restarting a computation    ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2500 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER equal to -1, 0, 1 or 2        ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2510 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' must be a positive real number                   ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2511 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' DOIT ETRE A POSITIVE REAL                        ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2520 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@       DTMIN must be > or =  0.     and                     ',/,&
'@              DTMIN must be  < or = DTMAX                   ',/,&
'@      Here DTMIN = ',E14.5     ,' and DTMAX = ',E14.5        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  When the time-step is non uniforme and constant           ',/,&
'@              DTMIN <= timestep <= DTMAX                    ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2530 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    CDTVAR(',I10   ,') DOIT ETRE A  STRICTLY POSITIVE REAL  ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  CDTVAR(I) multiplyer coefficient applied to the           ',/,&
'@  timestep for the  resolution of variable I.               ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2540 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@  ON DEMANDE UNE LIMITATION DU PAS DE TEMPS LIEE AUX EFFETS ',/,&
'@    DE DENSITE (IPTLRO = ',I10   ,') AVEC UNE OPTION DE     ',/,&
'@    PAS DE TEMPS FIXE (IDTVAR = ',I10   ,')                 ',/,&
'@                                                            ',/,&
'@  Le calcul sera engage, mais le pas de temps ne sera pas   ',/,&
'@    clippe. Le code indiquera juste the eventuels           ',/,&
'@    depassements de critere local.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2541 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@  A time-step reduction in liaison with variable density    ',/,&
'@    effects (IPTLRO = ',I10   ,') was requested while       ',/,&
'@   steady-state algorythm is selected (IDTVAR = ',I10   ,') ',/,&
'@                                                            ',/,&
'@  Computation will run, but option IPTLRO is ignored        ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2600 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER EGAL A 0, 10, 20, 21, 30, 31,',/, &
'@    40, 41, 42, 50, or 60'                                   ,/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2601 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    THE K-OMEGA SST TURBULENCE MODEL IS NOT COMPATIBLE WITH ',/,&
'@    TWO-WAY COUPLING IN LAGRANGIAN MODELLING                ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2605 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    Synthetic Vortex method for LES inlet is not compatible ',/,&
'@    with mutiphase option                                   ',/,&
'@    we have here                                            ',/,&
'@    NPHAS =',I10                                             ,/,&
'@    IVRETX=',I10                                             ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2606 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    Synthetic Vortex method is only for for LES             ',/,&
'@    here we have                                            ',/,&
'@    ITURB(1)=',I10                                           ,/,&
'@    IVRETX  =',I10                                           ,/,&
'@                                                            ',/,&
'@  computation will go on while ignoring keyword  IVRTEX.    ',/,&
'@    (it is reset to       0)                                ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2607 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@   A reduction of the extened neighbourhood was selected    ',/,&
'@   for the caculation of the gradients by least squares.    ',/,&
'@   However this will also be applied to the averaging in    ',/,&
'@    the LES Dynamic model (also selected)                   ',/,&
'@      ITURB =',I10                                           ,/,&
'@      IMRGRA=',I10                                           ,/,&
'@                                                            ',/,&
'@  Computation will run, but                                 ',/,&
'@                                                            ',/,&
'@  averaging of the Smagorinsky constant can be              ',/,&
'@  degraded, as it uses the same reduced neighbourhood.      ',/,&
'@                                                            ',/,&
'@  Recommendation                                            ',/,&
'@    - use extended neighbourhood                            ',/,&
'@                                                            ',/,&
'@      (IMRGRA = 2)                                          ',/,&
'@    - user defines (yourself) the averaging of the dynamic'  ,/,&
'@   Smagorinsky constant via  subroutine USSMAG.             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2610 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    DONNEES INCOHERENTES                                    ',/,&
'@    ',A31,I10                                                ,/,&
'@    ',A31,I10                                                ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
! 2620 FORMAT
!     & (' ATTENTION : ON DEMANDE ',A6,' = ',I10,/,
!     &  '                  AVEC GRAVITE = ',3E14.5)
 2621 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@  Gravity is taken into account                             ',/,&
'@             ',3E14.5                                        ,/,&
'@    without solving for temperature or energy               ',/,&
'@    (ISCALT = ',I10   ,')                                   ',/,&
'@                                                            ',/,&
'@  The calculation will run.                                 ',/,&
'@                                                            ',/,&
'@  The above options are not incompatible                    ',/,&
'@    but gravity is more often activated when density is     ',/,&
'@    variable with temperature (natural convection)          ',/,&
'@   this could be an error                                   ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2622 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@     Gravity is taken into account                          ',/,&
'@   in the turbulence source terms  (',A6,' = ',I10   ,')    ',/,&
'@    without solving for temperature or energy               ',/,&
'@                                                            ',/,&
'@  The calculation will run.                                  ',/&
'@                                                            ',/,&
'@  gravity usualy afects turbulence only via density effects ',/,&
'@    Check that density is variable.                         ',/,&
'@     other than via temperature                             ',/,&
'@  this could be by user defined density                     ',/,&
'@    in subroutine   usphyv.                                 ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2623 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@ Coefficient RELAXV for turbulence variables was modified   ',/,&
'@ although IKECOU is not = 0.      It is in fact to          ',/,&
'@ - for k                   : ',E12.4                         ,/,&
'@ - for  epsilon (or omega) : ',E12.4                         ,/,&
'@                                                            ',/,&
'@ The  modification will be ignored (RELAXV is only useful   ',/,&
'@ if IKECOU=0)                                               ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2624 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@ LCoefficient RELAXV for turbulence variables must          ',/,&
'@ be a REAL comprised between 0 and 1. IT IS EQUAL :         ',/,&
'@ - for  k                  : ',E12.4                         ,/,&
'@ - for  epsilon (ou omega) : ',E12.4                         ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2625 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@ Coefficient RELAXV for the pressure must be a REAL         ',/,&
'@ between  0 et 1. It  IS EQUAL : ',E12.4                     ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2630 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@   Van Driest near wall damping was selected                ',/,&
'@    Van Driest (IDRIES(IPHAS) = ', I10,')'                   ,/,&
'@    with a LES model not compatible (ITURB(IPHAS) = ',I10,')',/,&
'@    (dynamic model or WALE model)                           ',/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2640 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' must be a REAL number in the range  [0;1]        ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2650 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAR   ',A8                                            ,/,&
'@    ',A6,'(',I10   ,') MUST BE AN INTEGER EQUAL  0  OR 1    ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  ICPSYR(I) is a flag for coupling scalar I with            ',/,&
'@    SYRTHES  (solid walls modeller)                         ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2660 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  No coupling with  SYRTHES has been defined                ',/,&
'@  The number of coupled scalars is however     ',I10         ,/,&
'@    (the total number of scalars is',I10   ,')              ',/,&
'@  Verify  the couplage SYRTHES-Noyau.                       ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2661 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE POUR LE COUPLAGE SYRTHES                    ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  A coupling with   SYRTHES was defined                     ',/,&
'@  This requires a coupled scalar (and only one).            ',/,&
'@  The total number of scalars   is ',I10                     ,/,&
'@  The number of coupled scalars is ',I10                     ,/,&
'@  Verify   le couplage SYRTHES-Noyau.                       ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2662 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE for coupling with SYRTHES                   ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  The coupled scalar must be the temperature                ',/,&
'@                                                            ',/,&
'@  Here it is scalar number : ',I10                           ,/,&
'@    which is not temperature because                        ',/,&
'@                    ISCSTH(',I10   ,') = ',I10               ,/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2663 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    INCOHERENCE IN COUPLING WITH SYRTHES                    ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@  For compressible and SYRTHES coupling                     ',/,&
'@    the coupled scalar must be energy.                      ',/,&
'@  here it is scalar                      ',I10               ,/,&
'@    which is not the energy  here since                     ',/,&
'@                    ISCSTH(',I10   ,') = ',I10               ,/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2664 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    CHOIX DU CALCUL DES ESTIMATEURS D''ERREUR               ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@   Flag IESCAL related to                                   ',/,&
'@   error estimate number IEST = ',I10   ,' for              ',/,&
'@    Navier-Stokes MUST BE AN INTEGER egal to 0, 1 or 2.     ',/,&
'@  It IS EQUAL : IESCAL(',I10  ,') = ',I10                    ,/,&
'@                                                            ',/,&
'@  Rq. : the possible values of IEST are    :                ',/,&
'@        IESPRE = ',I10                                       ,/,&
'@        IESDER = ',I10                                       ,/,&
'@        IESCOR = ',I10                                       ,/,&
'@        IESTOT = ',I10                                       ,/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2700 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    Choice of method for computing distance to the wall     ',/,&
'@                                                            ',/,&
'@  ICDPAR MUST BE AN INTEGER  EQUAL TO -2, -1, 1 or 2        ',/,&
'@  IL IS EQUAL ', I10                                         ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2710 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE   A  REAL INCLUD IN SEGMENT       [0.;1.]',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2720 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE   A  REAL SUPERIEUR or EGAL A 1.         ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2730 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE   A  REAL EGAL A 0. or A 1.              ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2740 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE   A  STRICTLY POSITIVE REAL.             ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2750 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER  LESS THAN or EGAL A 1        ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ',A6,'(',I10   ,') MUST BE AN INTEGER EQUAL  0  OR 1    ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@  The multigrid algorithm for the linear system resolution  ',/,&
'@  is not compatible with convected variables.               ',/,&
'@                                                            ',/,&
'@  The calculation could NOT run.                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    VARIABLE ',A8                                            ,/,&
'@    ',A6,'(',I10   ,') MUST BE AN INTEGER                   ',/,&
'@      STRICTLY  POSITIVE                                    ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ',A6,' MUST BE AN INTEGER  STRICTLY  POSITIVE           ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ENTREE DES DONNEES                          ',/,&
'@    =========                                               ',/,&
'@    REFERENCE VELOCITY  UREF WAS NOT DEFINED                ',/,&
'@    or is ill defined  (NEGATIVE value ? ).                 ',/,&
'@   It IS EQUAL ',E14.5                                       ,/,&
'@                                                            ',/,&
'@  calculation can only run if turbulence is defined         ',/,&
'@  via a restart file or                                     ',/,&
'@  suroutine usiniv.                                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 4300 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAR ',A8                                              ,/,&
'@    ISCSTH(',I10   ,') MUST BE AN INTEGER  EGAL A -1,0,1,2,3',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4320 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    ISCAVR(',I10   ,') MUST BE AN INTEGER                   ',/,&
'@      POSITIVE or NUL AND                                   ',/,&
'@      LESS THAN or EGAL A NSCAL = ',I10                      ,/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  If ISCAVR(I) =0 , le scalare  I is not a  variance        ',/,&
'@  If ISCAVR(I) is POSITIVE, scalare I is a variance :       ',/,&
'@    it is the variance of fluctuations of scalaire J        ',/,&
'@    who''s number is   ISCAVR(I)                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4321 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    THE SCALAR ',A8  ,'IS DEFINED AS        FLUCTUATION     ',/,&
'@    OF SCALAR ',A8                                           ,/,&
'@    (ISCAVR(',I10   ,') = ',I10   ,'),                      ',/,&
'@    WHICH ITSELF IS DEFINED AS A FLUCTUATION                ',/,&
'@    (ISCAVR(',I10   ,') = ',I10   ,').                      ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  If ISCAVR(I) is POSITIVE, scalar  I is a  variance :      ',/,&
'@    variance of fluctuations of scalar J                    ',/,&
'@    who''s number is ISCAVR(I) , so we must have            ',/,&
'@    ISCAVR(J) = 0                                           ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4330 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    ICLVFL(',I10   ,') MUST BE AN INTEGER  EGAL A 0, 1 or 2 ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  ICLVFL(I) defines the type of clipping of scalar I        ',/,&
'@    when it is a variance of  fluctuations.                 ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4331 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    ICLVFL(',I10   ,') is only used for  VARIANCES          ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@    BUT THE SCALAR IS NOT A VARIANCE                        ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  ICLVFL(I) flags the type of clipping for scalar I         ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4340 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    VISLS0(',I10   ,') MUST BE   A  REAL POSITIVE           ',/,&
'@      WHILE IVISLS(',I10   ,') = ',I10                       ,/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  VISLS0(I) is the molecular diffusion coefficient  of the  ',/,&
'@    scalar and  MUST BE POSITIVE when ivisls is  different  ',/,&
'@    from 1 (it must then be defined in USPHYV )             ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4350 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    SIGMAS(',I10   ,') MUST BE   A  REAL POSITIVE           ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  SIFMAS(I) is the turbulent Prandtl turbulent              ',/,&
'@    associated to scalar I.                                 ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4360 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    SCAMIN(',I10   ,') IS EQUAL ',E14.5                      ,/,&
'@      AVEC ICLVFL(',I10   ,') = ',I10                        ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  SCAMIN(I) is the  minimale acceptable value for           ',/,&
'@    scalaire I. When this scalar is a variance              ',/,&
'@    (ISCAVR(I) > 0) value of SCAMIN is only used if         ',/,&
'@                  ICLVFL(I) = 2                             ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4361 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    SCAMAX(',I10   ,') IS EQUAL ',E14.5                      ,/,&
'@      AVEC ICLVFL(',I10   ,') = ',I10                        ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  SCAMAX(I) is the maximum acceptable value for             ',/,&
'@    scalar  I. When this is a variance                      ',/,&
'@    (ISCAVR(I) > 0) the value of SCAMAX is only used        ',/,&
'@               if ICLVFL(I) = 2                             ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4370 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    SCAMAX(',I10   ,') MUST BE   A  REAL STRICTLY  POSITIVE ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  SCAMAX(I) is the maximum acceptable value for             ',/,&
'@    scalar   I, which is a variance                         ',/,&
'@  with ICLVFL(I) = 2, value SCAMAX must be                  ',/,&
'@   strictly  positive.                                      ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 4380 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    SCALAIRE ',A8                                            ,/,&
'@    RVARFL(',I10   ,') MUST BE   A  REAL STRICTLY  POSITIVE ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  RVARFL(I) is the coefficient R for the scalar I (which is ',/,&
'@ a variance) related to the dissipation equation sourceterme',/,&
'@    - (1/R) rho scalaire epsilon/k                          ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5003 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    ANGULAR PERIODICITE IS NOT COMPATIBLE WITH THE          ',/,&
'@    ENHANCED PRESSSURE-VELOCITY COUPLING  or ALE METHOD     ',/,&
'@      IN THE CURRENT VERSION                                ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  variable COMMANDE_PERIO was defined in the run-case script',/,&
'@    (periodicity was activated which results in :           ',/,&
'@        IPERIO = ',I10,   ')                                ',/,&
'@    and some peridic boundaries involve rotation.           ',/,&
'@  The flag IPUCOU is defined as  ',I10                       ,/,&
'@   in the interface or usini1 (enhanced coupling for        ',/,&
'@    IPUCOU=1).                                              ',/,&
'@  The ALE fag  IALE is defined as  ',I10                     ,/,&
'@    in the interface or usini1 (method activated if IALE=1) ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    PERIODICITY IS INCOMPATIBLE WITH THE MULTIGRID SOLVER   ',/,&
'@       IN THE CURRENT VERSION                               ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  variable COMMANDE_PERIO was defined in the run-case script',/,&
'@    (periodicity was activated which results in :           ',/,&
'@        IPERIO = ',I10,   ')                                ',/,&
'@  Flag IMGR is defined as 1  for one variable at least      ',/,&
'@    in the interface or usini1.                             ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5005 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    PERIODICITY IS INCOMPATIBLE WITH THIS METHOD FOR        ',/,&
'@    COMPUTING DISTANCE TO WALL IN THE CURRENT VERSION       ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  variable COMMANDE_PERIO was defined in the run-case script',/,&
'@    (periodicity was activated which results in :           ',/,&
'@        IPERIO = ',I10,   ')                                ',/,&
'@  the parameters specified need the calculation of the      ',/,&
'@  distance to the wall (Rij-eps LRR with wall echo term   , ',/,&
'@     van Driest damping   or k-omega SST).                  ',/,&
'@  The method for computing wall distance  :                 ',/,&
'@    ICDPAR = ',I10,   ' is not compatible with              ',/,&
'@     periodicity                                            ',/,&
'@                                                            ',/,&
'@  Recommendation: use ICDPAR = 1 or -1.                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
! 5007 format(
!     &'@                                                            ',/,
!     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!     &'@                                                            ',/,
!     &'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,
!     &'@    =========                                               ',/,
!     &'@    PERIODICITY IS NOT COMPATIBLE WITH RADIATIVE HEAT       ',/,
!     &'@      TRANSFER IN SEMI TRANSPARENT MEDIA                    ',/,
!     &'@      (in the current version)                              ',/,
!     &'@                                                            ',/,
!     &'@   The calculation could NOT run.                           ',/,
!     &'@                                                            ',/,
!     &'@  Flag  IPERIO is equl to ',I10                              ,/,
!     &'@    in usini1 (periodicity actived if IPERIO=1).            ',/,
!     &'@  Flag IIRAYO(',I10,') is equal to ',I10                     ,/,
!     &'@    in usray1.                                              ',/,
!     &'@                                                            ',/,
!     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
!     &'@                                                            ',/)
 5008 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@ ROTATION PERIODICITY IS NOT COMPATIBLE WITH RADIATIVE HEAT ',/,&
'@      TRANSFER IN SEMI TRANSPARENT MEDIA                    ',/,&
'@      (in the current version)                              ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@ variable COMMANDE_PERIO was activated in the runcase script',/,&
'@  Flag  IPERIO is equal to ',I10                             ,/,&
'@              (periodicity actived if IPERIO=1).            ',/,&
'@  Flag IIRAYO is equal to ',I10                              ,/,&
'@    in usray1.                                              ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5009 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    DEFECTS CAN APPEAR WHEN USING A COMBINATION OF',/,    &
'@     ANGULAR PERIODICITE (ROTATION) AND RSTM RIJ-EPSILON.   ',/,&
'@                                                            ',/,&
'@   COMMANDE_PERIO was defined in the run-case script        ',/,&
'@                                                            ',/,&
'@    and IPERIO = ',I10,   ')                                ',/,&
'@    and some periodic boundaries involve rotation           ',/,&
'@      Flag for turb ITURB(IPHAS) is = ',I10                  ,/,&
'@                                                            ',/,&
'@    Job can run.                                            ',/,&
'@                                                            ',/,&
'@    The defects are related to the turbulent transport terms',/,&
'@    in the Re stress equations (equivalent to an anisotropic',/,&
'@    diffusion tensor), but these terms are generaly small   ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 5010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@   WARNING :      WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@   DEFECTS CAN APPEAR WHEN USING A COMBINATION OF           ',/,&
'@     ANGULAR PERIODICITE (ROTATION) AND SECOND ORDER SCHEME ',/,&
'@     IN TIME FOR VELOCITY                                   ',/,&
'@                                                            ',/,&
'@   COMMANDE_PERIO was defined in the run-case script        ',/,&
'@                                                            ',/,&
'@    and IPERIO = ',I10,   ')                                ',/,&
'@  Flags  THETAV for 3 velocity components      Ux, Uy, Uz   ',/,&
'@    of velocity                                             ',/,&
'@    are selected  (in usini1 or by default                  ',/,&
'@    with the following values :'                             ,/,&
'@    THETAV(IU(IPHAS))  THETAV(IV(IPHAS))  THETAV(IW(IPHAS)) ',/,&
'@        ',E14.5     ,'     ',E14.5     ,'     ',E14.5        ,/,&
'@                                                            ',/,&
'@  Job can run.                                              ',/,&
'@     defects are due to the computation of the velocity     ',/,&
'@    gradient tensor at rotation perodicity boundaries       ',/,&
'@    This uses an explicit approach which reduces the order  ',/,&
'@    of the time discretization scheme                       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 6002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    PARALLEL COMPUTING AND LAGRANGIAN PARTICULE TRANSPORT   ',/,&
'@      OR NOT COMPATIBLE IN THE CURRENT VERSION              ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  The present CPU has rank ',I10                             ,/,&
'@  Flag  IILAGR had value',I10                                ,/,&
'@    in uslag1 ( lagrangien module is active for  IILAGR>0). ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6004 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@     PARALLEL COMPUTING AND MUTIGRID SOLVER                 ',/,&
'@                OR NOT COMPATIBLE IN THE CURRENT VERSION    ',/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@   The present CPU has rank ',I10                            ,/,&
'@  Flag IMGR  is set to       1                              ',/,&
'@    for at least one variable in interface or usini1.       ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 6005 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    THIS METHOD FOR COMPUTING WALL DISTANCES                ',/,&
'@    IS NOT COMPATIBLE WITH PARALLEL COMPUTING               ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  The present CPU has rank ',I10                             ,/,&
'@                                                            ',/,&
'@    Wall distance is necessary for RSTM wall echo terms, or,',/,&
'@     van Driest damping or k-omega SST).                    ',/,&
'@  Method for comuting wall distance is defined by           ',/,&
'@    ICDPAR = ',I10,   '                                     ',/,&
'@    l IS NOT COMPATIBLE WITH PARALLEL COMPUTING             ',/,&
'@                                                            ',/,&
'@  Use ICDPAR = 1 or -1.                                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    FLAG FOR ALE  METHOD                                    ',/,&
'@                                                            ',/,&
'@  IALE should be = 0 or 1                                   ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  Verify   the parameters given  in   interface or usalin.  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    NOMBRE D''ITERATIONS D''INITIALISATION DU FLUIDE EN ALE ',/,&
'@                                                            ',/,&
'@  NALINF MUST BE A POSITIVE INTEGER                         ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  Verify   the parameters given  in     interface or usalin.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7020 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    NON VALID COEFFICIENTS IN  NEWMARK METHOD               ',/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@  ALPNMK MUST BE   between  0 and  1                        ',/,&
'@  BETNMK MUST BE   between  0 and 1/2                       ',/,&
'@  GAMNMK MUST BE   between  0 and 1                         ',/,&
'@  We have here:                                             ',/,&
'@                                                            ',/,&
'@       ALPNMK      BETNMK      GAMNMK                       ',/,&
'@ ',     E12.4,      E12.4,      E12.4                        ,/,&
'@                                                            ',/,&
'@  Verify   the parameters given  in interface, usini1       ',/,&
'@    or usalin.                                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7030 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    MAX number of iterations for implicit ALE method        ',/,&
'@                                                            ',/,&
'@  NALIMX MUST BE A POSITIVE INTEGER                         ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  Verify   the parameters given  in interface or usalin.    ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7040 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@    PRECISION DU COUPLAGE IMPLICITE EN ALE                  ',/,&
'@                                                            ',/,&
'@  EPALIM MUST BE A REAL NUMBER,  STRICTLY  POSITIVE         ',/,&
'@   IT HAS VALUE ',E14.5                                      ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  Verify   the parameters given  in interface or usalin.'   ,/, &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 7050 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@   INITIALISATION  ITERATION  FOR ALE                       ',/,&
'@                                                            ',/,&
'@  ITALIN must be =   0 or 1                                 ',/,&
'@   IT HAS VALUE ',I10                                        ,/,&
'@                                                            ',/,&
'@   The calculation could NOT run.                           ',/,&
'@                                                            ',/,&
'@  Verify   the parameters  given  in usalin.                ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 8000 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@  WARNING:   STOP WHILE READING INPUT DATA               ',/,&
'@    =========                                               ',/,&
'@               COMPRESSIBLE FLOW MODULE                     ',/,&
'@    T0 AND  P0 MUST BE STRICTLY POSITIVE REAL NUMBERS       ',/,&
'@    Here they have values:                                  ',/,&
'@                   T0 = ',E14.5                              ,/,&
'@                   P0 = ',E14.5                              ,/,&
'@                                                            ',/,&
'@  Computation CAN NOT run                                   ',/,&
'@                                                            ',/,&
'@ Check the input data given via User Interface or in usini1.',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!===============================================================================
! 8. SORTIE
!===============================================================================

return
end subroutine
