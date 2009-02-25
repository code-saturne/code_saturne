!-------------------------------------------------------------------------------

!     This file is part of the Code_Saturne Kernel, element of the
!     Code_Saturne CFD tool.

!     Copyright (C) 1998-2008 EDF S.A., France

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

                  subroutine mtini1
!================


!===============================================================================
!  FONCTION  :
!  ---------

!         INITIALISATION DES COMMON MATISSE
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
!    nom           !type!mode !                   role                         !
!__________________!____!_____!________________________________________________!
!__________________!____!_____!________________________________________________!

!     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
!            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
!     MODE : <-- donnee, --> resultat, <-> Donnee modifiee
!            --- tableau de travail
!===============================================================================

implicit none

!===============================================================================
!     DONNEES EN COMMON
!===============================================================================

include "paramx.h"
include "numvar.h"
include "optcal.h"
include "ihmpre.h"
include "matiss.h"
include "cstphy.h"
include "cstnum.h"
include "entsor.h"
include "parall.h"

!===============================================================================

! Arguments

! VARIABLES LOCALES

integer          ii    , jj    , kk
integer          iphas
integer          iarret
double precision tsor0
double precision xmapmn, xmapmx
double precision xmapvl
double precision cxsd  , cysd  , por0
double precision ustmat, drfmat
double precision prfmat, erfmat, srfmat

!===============================================================================

!===============================================================================
! 1. PAR DEFAUT MATISSE N'EST PAS ACTIVE
!===============================================================================

!     Ceci doit etre present en dehors du if defined xml :
!       on ne peut activer Matisse que avec IHM et XML

!     Ainsi, dans les autres sous-programmes Matisse, il est inutile de
!       tester si l'interface ou Xml sont disponibles (en effet, si
!       l'interface ou Xml ne sont pas disponibles, mtini1 est le seul
!       sous-programme Matisse dans lequel on entre).

imatis = 0


!===============================================================================
! 2. LECTURE DE LA BALISE MATISSE EN XML (si xml present)
!===============================================================================

! --- Lecture de IMATIS
!     (si IHM presente ; sinon, Matisse ne sera pas actif)

if (iihmpr.eq.1) then
  call csmtpr(imatis)
endif

! --- Si Matisse n'est pas actif, on sort a la fin du if
!     Sinon, on lit les donnees, on les traite, ...

if (imatis.eq.1) then


!===============================================================================
! 3. REMPLISSAGE DES COMMONS A PARTIR DU XML
!===============================================================================

! On remplit les COMMON  de données géométriques et physiques
!  a partir du contenu du fichier xml

! --- /IMTGEO/

  call csgein(nptran, nplgrs, nelgrs, nchest, netran, itypen)

! --- /RMTGEO/

  call csgedb(epregi, epchem, hconve, rconve, hchali, hcheva,     &
              hfttoi, ptrres, frdtra, plgres, epchel, dmcont)

! --- /IMTPHY/

  call csphat(imdcnt, icofor, iconlg, ialveo)

! --- /RMTPHY/

  call csphdb(dtdtmx, puicon, tinit,  tcrit,  emicon, emimur,     &
              hepcnt, dhpcnt, debmas, pdccha, pdcfch, dhchea,     &
              sdchea, pdcche, pdccch, dhches, sdches, pdcalg,     &
              pdcatv, argamt, pdcslg, pdcstv, argavl, amppdc,     &
              dhalve, hreso,  hplen,  dpvent)

! --- Cartes 2D et 3D

  do ii = 1, ncarte
    do jj = 1, nmtdir
      call csnbmp(jj,ii, nzocar(jj,ii))
      if (nzocar(jj,ii) .gt. nzonmx) then
        write(nfecra,9011) nzocar(jj,ii), nzonmx
        call csexit (1)
      endif
    enddo
  enddo

  do ii = 1, ncarte
    do jj = 1, nmtdir
      do kk = 1 , nzocar(jj,ii)
        call csdfmp (kk ,jj, ii, xmapmn, xmapmx, xmapvl)
        vizcar(1, kk, jj, ii) = xmapmn
        vizcar(2, kk, jj, ii) = xmapmx
        if (ii.eq.icpuis) then
          vcarth(kk, jj) = xmapvl
        endif
      enddo
    enddo
  enddo


!===============================================================================
! 4. CALCULS COMPLEMENTAIRES AVANT VERIFICATIONS
!===============================================================================

!     Variables en common

! --- Hauteur d'erosion reduite a un nombre entier de mailles
  hercnt = epchel*int(hepcnt/epchel)


!===============================================================================
! 5. VERIFICATIONS
!===============================================================================

! --- Par defaut, tout est suppose correct (on ne s'arretera donc pas)
  iarret = 0


! --- Matisse n'est pas parallele

  if (irangp.ge.0) then
    iarret = 1
    write(nfecra,9001)
  endif


! --- NSCAUS = 3
!     Matisse a besoin de 3 scalaires et 3 seulement
!     Des developpements supplementaires sont necessaires pour
!       que d'autres scalaires utilisateurs puissent etre pris en compte
!       (exemple : acces a ustssc)

  if(nscaus.ne.3) then
    iarret = 1
    write(nfecra,9021) nscaus
  endif


! --- TCRIT > TINIT (cf calcul de VITREF)
!     TCRIT est la "Temperature d'air de sortie critique en degres C"
!       c'est une estimation de la temperature de l'air en sortie de
!       l'entrepot et donc naturellement strictement superieure a TINIT

  if(tcrit.le.tinit) then
    iarret = 1
    write(nfecra,9031) tcrit, tinit
  endif


! --- RCONVE = 1 si 2D
!     RCONVE est le rapport du convergent represente sur le maillage.
!       Attention, "respresente sur le maillage" est important :
!         il ne s'agit pas du convergent reel.
!       En 2D, il ne peut pas y avoir de convergent (le convergent
!         reduit la dimension X, transverse, des cheminees) et RCONVE
!         est donc necessairement unite.
!       Pour savoir si l'on fait une simulation 2D, on regarde si
!         FRDTRA est non unite (necessaire mais non suffisant : FRDTRA
!         non unite indique que l'on represente une configuration
!         3D par un calcul 2D, i.e. a une seule maille dans la
!         direction X).

  if ( (abs(frdtra-1.d0).ge.epzero).and.                          &
       (abs(rconve-1.d0).ge.epzero)     ) then
    iarret = 1
    write(nfecra,9041) frdtra, rconve
  endif


! --- Si on modelise les panaches, la hauteur doit etre correcte
!     HERCNT est la hauteur reduite a un numbre entier de mailles
!       dans la direction verticale : c'est la hauteur que l'on
!       va reellement prendre en compte dans le calcul. Elle doit
!       etre positive et ne pas depasser la hauteur de l'espace
!       libre au dessus des colis :

!       HERCNT > 0 et HERCNT <= NCHEST*EPCHEL-HRESO

  if(imdcnt.eq.1) then
    if ( (hercnt.le.0.d0).or.                                     &
         (hercnt.gt.nchest*epchel-hreso) ) then
      iarret = 1
      write(nfecra,9051) imdcnt,                                  &
           hepcnt, hercnt,                                        &
           nchest*epchel, nchest, epchel, hreso,                  &
           hercnt, nchest*epchel-hreso
    endif
  endif


! --- ARRET EVENTUEL

  if(iarret.ne.0) then
    call csexit (1)
    !==========
  endif


!===============================================================================
! 6. VALEURS PAR DEFAUT
!===============================================================================

! === On ne traite qu'une seule et unique phase
! ===========================================================

!     On ne le dit qu'une fois, au debut, pour eviter d'en oublier

  iphas = 1


! === Options numeriques
! ===========================================================


! --- IPHYDR = 0 à cause des pertes de charges dans la cheminee
!     Plus exactement, a cause des pertes de charges en
!       sortie, avec presence de gravite.
!     A noter, cependant, qu'il semble que meme avec ICALHY = 0
!       l'option IPHYDR cause des problemes (et que c'était déjà le cas
!       dans les premieres versions de Matisse, puisque l'option etait
!       certes proposee par défaut =1, mais mise systematiquement =0
!       dans tous les cas de calcul).
!     Bien qu'il s'agisse de la valeur par defaut de Code_Saturne 1.2,
!       on laisse l'option ici dans la mesure ou des tests
!       supplementaires seraient necessaires pour comprendre pourquoi
!       on rencontre des difficultes avec les configurations etudiees
!       dans Matisse.

  iphydr = 0


! --- Par contre, l'extrapolation de la pression fournit des resultats
!       satisfaisants. Elle est utile et ne fait pas apparaitre de
!       probleme sur les cas testes.

  extrag(ipr(iphas)) = 1.d0


! --- IPUCOU prend sa valeur par defaut de Code_Saturne 1.2
!     L'option est cependant conservee ici, dans la mesure ou elle
!       etait prise egale a 1 auparavant et ou cela fait une difference
!       (toute petite, certes, mais on garde ici la ligne pour memoire).

  ipucou = 0


! --- Pas de convection pour les temperatures solides
!     Pas de phenomene diffusif pour les parois (le rayonnement est
!       traite par terme source)

  iconv(isca(itpcmt)) = 0
  iconv(isca(itppmt)) = 0
  idiff(isca(itppmt)) = 0


! --- Schema convectif : robuste

  blencv(iu(iphas)) = 0.0d0
  blencv(iv(iphas)) = 0.0d0
  blencv(iw(iphas)) = 0.0d0
  if(nscaus.ge.1) then
    do ii = 1, nscaus
      blencv(isca(ii)) = 0.0d0
    enddo
  endif


! === Proprietes physiques et pas de temps
! ===========================================================


! --- Connectivite du rayonnement a calculer
!     (au premier passage dans mttssc)

  icnrok = 0


! --- Masse volumique, chaleur massique

!     On prend comme reference :
!     - la masse volumique RRFMAT a la temperature TRFMAT en degres C
!     - on utilise la loi des gaz parfaits
!     - on fixe CP0

  ro0(iphas) = (trfmat + tkelvi)*rrfmat /(tinit + tkelvi)
  cp0(iphas) = crfmat

!     On fixe P0 a 1 atm et PRED0 a 0
!     La reference de pression est prise au niveau de l'alimentation
  p0(iphas) = 1.013d5
  pred0(iphas) = 0.d0
  xyzp0(1,iphas) = 0.d0
  xyzp0(2,iphas) = 0.d0
  xyzp0(3,iphas) = hchali
  ixyzp0(iphas) = 1


! --- Vitesse de reference
!     pour le calcul de la viscosite turbulente et
!       pour les pertes de charge, le Richardson, ...
!     la vitesse de reference n'est pas utilisee pour le calcul du pas
!       de temps : le pas de temps est pris ici uniforme et constant ;
!       s'il etait envisage de le prendre variable en temps, pour
!       l'adapter automatiquement a la configuration finale du
!       calcul une fois l'etat stationnaire etabli, il conviendrait
!       alors d'utiliser la vitesse de reference pour calculer une
!       limite basee par exemple sur un nombre de Courant.


!   - En convection naturelle
!       l'echelle de vitesse est calculee par analyse dimensionnelle
!         comme VITREF = P / (E S), a partir :
!       . de l'energie (E, en Joule/m3) que l'air initialement a TINIT
!          peut extraire s'il atteint la temperature critique de
!          sortie TCRIT : E = RO0(IPHAS)*CP0(IPHAS)*(TSOR0-TINIT)
!       . de la puissance (P, en Joule/s) d'une ligne (x cst) de
!           NPLGRS conteneurs, P = PUICON * NPLGRS
!       . de la surface verticale (S en m2) de la coupe transverse
!           d'une ligne de conteneurs, S = PTRRES*EPCHEL*NCHEST

  if(icofor.eq.0)then

    tsor0  = tcrit
    prfmat = puicon * nplgrs
    erfmat = ro0(iphas)*cp0(iphas)*(tsor0-tinit)
    srfmat = ptrres*epchel*nchest

    vitref = prfmat/(erfmat*srfmat)


!   - En convection forcee
!       La vitesse de reference est calculee comme le rapport du debit
!         reel total (DEBMAS/RO0(IPHAS)) a la surface de la zone de
!         stockage (NPTRAN*PTRRES*NCHEST*EPCHEL).
!       Elle est corrigee par le facteur de reduction transverse
!         FRDTRA du maillage par rapport a la realite afin d'obtenir
!         une vitesse representative de la vitesse reelle.

  else

    vitref =                                                      &
         debmas/(ro0(iphas)*frdtra*nptran*ptrres*nchest*epchel)

! TSOR0 inutile ici ; on conserve la formule pour memoire
!     TSOR0= TINIT + NPTRAN*NPLGRS*PUICON/(DEBMAS/FRDTRA)/CP0(IPHAS)

  endif


! --- Viscosite dynamique totale (modelisee, constante et uniforme)

!     Le tableau VISCL0 contient la viscosite dynamique totale
!       (moleculaire + turbulente). On evalue d'abord la viscosite
!       turbulente puis on ajoute la viscosite moleculaire XMUMAT.

!     La viscosite dynamique turbulente est modelisee sous la forme
!       mu_t = rho * u*_ref * D_ref, formule dans laquelle :
!     . u*_ref est la vitesse de frottement deduite de la vitesse de
!              reference VITREF calculee ci-dessus et d'une intensite
!              turbulente imposee RTURB0
!     . D_ref  est une distance de reference basee sur l'encombrement
!              du milieu

!     La valeur de u*_ref est calculee a partir de l'energie cinetique
!       turbulente k par u*_ref = Cmu**(1/4) k**(1/2). La valeur de k
!       se deduit de l'intensite turbulente I, supposee connue, par
!       k = 3/2 (I V_ref)**2 avec V_ref la vitesse de reference
!       VITREF calculee precedemment.

!     La valeur de D_ref est calculee par D_ref = 0.2(Pt - d) ou Pt est
!       le pas transverse du reseau et d le diametre des conteneurs.

!     Ainsi, mu_t = rho * u*_ref * L_ref
!            . u*_ref = Cmu**(1/4) * (3/2)**(1/2) * I * V_ref
!            . L_ref  = 0.2 * (Pt - d)

!     On ajoute la viscosite moleculaire XMUMAT a mu_t pour obtenir
!       la viscosite dynamique totale VISCL0

!     VISCL0 est recalcule dans mtphyv (VITREF est mis a jour au cours
!       du calcul)


!   - Calcul de la viscosite turbulente

  ustmat = cmu**0.25d0 * sqrt(1.5d0) * (rturb0/100.d0) * vitref
  drfmat = 0.2d0 * (ptrres-dmcont)
  viscl0(iphas) = ro0(iphas) * ustmat * drfmat

!   - Ajout de la viscosite moleculaire

  viscl0(iphas) = viscl0(iphas) + xmumat


! --- Pas de temps

!   - Pas de temps uniforme en espace et constant en temps
!       Initialement, on souhaitait adopter dans Matisse un pas de
!         temps variable en temps, mais l'option n'a jamais ete
!         effectivement utilisee.
!       Quels que soient les choix par defaut de Code_Saturne, on
!         conserve IDTVAR ici pour memoire, au cas ou l'on reprendrait
!         les tests plus tard.

  idtvar = 0

!   - Calcul de la porosite verticale liee a la presence d'un
!       conteneur :
!       POR0 = section horizontale passante /section horizontale totale
!       (on suppose le conteneur cylindrique à base circulaire,
!        d'axe vertical)
!       Utilise pour le calcul du pas de temps ci-dessous

  cxsd = ptrres/dmcont
  cysd = plgres/dmcont
  por0 = 1.d0 - pi/(cxsd*cysd*4.d0)

!   - On calcule le pas de temps par analyse dimensionnelle pour que
!       le volume d'air avoisinant un conteneur ne voie pas sa
!       temperature croitre de plus de DTDTMX degres en un pas de
!       temps.
!     On a donc DTREF = E * V / P avec :
!       . E = RO0(IPHAS)*CP0(IPHAS)*DTDTMX, en Joule/m3, l'energie
!           maximale que l'air peut recevoir par pas de temps,
!       . V = (PLGRES*PTRRES*EPCHEL*NCHEST)*POR0, en m3, le volume
!           d'air entourant un conteneur.
!       . P = PUICON, en Joule/s, la puissance d'un conteneur

  dtref = ro0(iphas)*cp0(iphas)*dtdtmx                            &
       * (plgres*ptrres*epchel*nchest)*por0 / puicon


! --- Diffusivite des scalaires
!     Une diffusivite de reference est fournie ici pour tous les
!       scalaires, mais la diffusivite est variable (voir mtphyv).
!     Il s'agit plutot d'une securite.
!     Pour la temperature de l'air ambiant (scalaire ITAAMT),
!       l'initialisation est ici de type "Prandtl=1" (on considere
!       un Prandtl "efficace" incluant la turbulence)
!     Pour la temperature de peau des parois (scalaires ITPPMT),
!       il n'y a pas de phenomene de diffusion pris en compte :
!       peu importe donc la valeur adoptee.

  if(nscaus.gt.0) then

!       On boucle sur les scalaires utilisateurs :
    do ii = 1, nscaus
!         Pour les scalaires qui ne sont pas des variances
      if(iscavr(ii).le.0) then
!           On definit la diffusivite
        visls0(ii) = viscl0(iphsca(ii))
      endif
    enddo

  endif

!     Fin du test IF (IMATIS.EQ.1)
endif


!--------
! FORMATS
!--------

 9001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : MATISSE                                     ',/,&
'@    =========                                               ',/,&
'@      ARRET SUR LE PARALLELISME DANS MTINI1                 ',/,&
'@                                                            ',/,&
'@    Le parallelisme n''est pas prevu avec Matisse           ',/,&
'@      or le calcul semble avoir ete lance en parallele.     ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,&
'@    =========                                               ',/,&
'@      ARRET SUR LE NOMBRE DE ZONES DES CARTES 2D ET 3D.     ',/,&
'@                                                            ',/,&
'@    Le nombre de zones demande est                 ',I10   ,  /,&
'@    Le nombre maximal de zones autorise est NZONMX ',I10   ,  /,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9021 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,&
'@    =========                                               ',/,&
'@      NOMBRE DE SCALAIRES INCORRECT.                        ',/,&
'@                                                            ',/,&
'@    NSCAUS represente                                       ',/,&
'@      le nombre de scalaires requis pour Matisse.           ',/,&
'@                                                            ',/,&
'@    NSCAUS doit etre exactement egal a 3.                   ',/,&
'@                                                            ',/,&
'@    Il vaut ici                                             ',/,&
'@      NSCAUS = ',I10   ,'                                   ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Contacter l''equipe de developpement.                   ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9031 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,&
'@    =========                                               ',/,&
'@      ECART TCRIT-TINIT NEGATIF OU NUL.                     ',/,&
'@                                                            ',/,&
'@    TCRIT represente                                        ',/,&
'@      la temperature d''air de sortie critique en degres C. ',/,&
'@    TINIT represente                                        ',/,&
'@      la temperature d''air en entree          en degres C. ',/,&
'@                                                            ',/,&
'@    TCRIT est une estimation de la temperature de l''air    ',/,&
'@      en sortie de l''entrepot.                             ',/,&
'@    TCRIT doit donc naturellement etre strictement          ',/,&
'@      superieure a TINIT.                                   ',/,&
'@                                                            ',/,&
'@    Or, les donnees saisies sont telles que                 ',/,&
'@      TCRIT = ',E12.5   ,'                                  ',/,&
'@      TINIT = ',E12.5   ,'                                  ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier les donnes saisies pour que TCRIT > TINIT.     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9041 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,&
'@    =========                                               ',/,&
'@      RAPPORT DE CONVERGENT NON UNITE EN 2D.                ',/,&
'@                                                            ',/,&
'@    RCONVE represente                                       ',/,&
'@      le rapport du convergent represente sur le maillage.  ',/,&
'@                                                            ',/,&
'@    RCONVE doit etre unite pour un calcul bidimensionnel.   ',/,&
'@                                                            ',/,&
'@    Ici, le calcul est apparemment bidimensionnel puisque   ',/,&
'@      le facteur de reduction transverse saisi n''est pas   ',/,&
'@      unite : FRDTRA = ',E12.5   ,'                         ',/,&
'@    Or, la valeur saisie pour le rapport du convergent      ',/,&
'@      n''est pas unite : RCONVE = ',E12.5   ,'              ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Modifier la valeur du facteur de reduction transverse   ',/,&
'@      ou du rapport de convergent.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9051 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET MATISSE (MTINI1)                      ',/,&
'@    =========                                               ',/,&
'@      HAUTEUR D''EROSION DES PANACHES INADAPTEE.            ',/,&
'@                                                            ',/,&
'@    La modelisation des panaches est activee                ',/,&
'@      avec IMDCNT = ',I10   ,'                              ',/,&
'@    La hauteur d''erosion des panaches reduite a un nombre  ',/,&
'@      entier de mailles doit etre strictement positive et   ',/,&
'@      inferieure ou egale a la distance separant le plafond ',/,&
'@      du sommet du reseau de colis.                         ',/,&
'@                                                            ',/,&
'@    La hauteur d''erosion saisie est   HEPCNT = ',E12.5      ,/,&
'@                                                            ',/,&
'@    La hauteur d''erosion reduite a un nombre entier de     ',/,&
'@      mailles est                      HERCNT = ',E12.5      ,/,&
'@    Le plafond est a la hauteur NCHEST*EPCHEL = ',E12.5      ,/,&
'@      avec NCHEST = ',I10   ,' et EPCHEL = ',E12.5   ,'     ',/,&
'@    La hauteur du reseau de colis est  HRESO  = ',E12.5      ,/,&
'@                                                            ',/,&
'@    L''inegalite suivante n''est pas verifiee :             ',/,&
'@             0 <    HERCNT    <= (NCHEST*EPCHEL - HRESO)    ',/,&
'@                 ',E12.5   ,'           ',E12.5   ,'        ',/,&
'@                                                            ',/,&
'@    Le calcul ne sera pas execute.                          ',/,&
'@                                                            ',/,&
'@    Descativer la modelisation des panaches ou modifier la  ',/,&
'@      hauteur d''erosion prescrite.                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


return
end
