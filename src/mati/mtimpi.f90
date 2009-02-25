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

                  subroutine mtimpi
!================


!===============================================================================
!  FONCTION  :
!  ---------

!         AFFICHAGE DES DONNEES DE CALCUL MATISSE
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
include "ihmpre.h"
include "matiss.h"
include "cstphy.h"
include "entsor.h"
include "optcal.h"

!===============================================================================

! Arguments

! VARIABLES LOCALES

character*15     name
integer          kechrg, kergrs, keclgr, kergch, keciel
integer          ii    , iphas
double precision djechr, djercl, djeclr, djerch
double precision hbdtoi
double precision tsor0 , dbm   , cxsd
double precision ureel , betmat, hrfmat, richar

!===============================================================================
!===============================================================================


!===============================================================================
! RECAPITULATIF DANS LE FICHIER resuMatisse.txt
!===============================================================================


! --- Ouverture du fichier
!       (il sera ferme au dernier pas de temps par mtproj)

NAME='resuMatisse'
open(unit=impmat,file=name,                                       &
         FORM='FORMATTED', STATUS='UNKNOWN', ERR=900)
goto 950

!   - En cas d'erreur : message et stop

 900  write (0, 9998) name
write (nfecra, 9999) name
call csexit (1)
!==========
 950  continue


! --- Ecriture des donnees geometriques

!     On recupere les donnees non stockees dans les COMMON
!       et utiles seulement pour l'affichage

call csmhdb(djechr, djercl, djeclr, djerch,                       &
            kechrg, kergrs, keclgr, kergch,                       &
            hbdtoi, keciel)

!     On imprime

write(impmat,1001)
write(impmat,1002)
write(impmat,1001)

if (itypen.eq.1) then
   WRITE(IMPMAT,2096) ' Emm  '
else
   WRITE(IMPMAT,2096) ' Vault'
endif

if (ialveo.eq.0) then
   write(impmat,2097) hreso
else
   write(impmat,2098) hreso
   write(impmat,2099) hplen
endif

write(impmat,2100) epregi
write(impmat,2101) epchem
write(impmat,2102) djechr
write(impmat,2103) kechrg
write(impmat,2104) djercl
write(impmat,2105) kergrs
if (itypen.eq.0) then
  write(impmat,2106) djeclr
  write(impmat,2107) keclgr
  write(impmat,2108) djerch
  write(impmat,2109) kergch
endif
if (itypen.eq.1) then
  write(impmat,2124) hbdtoi
  write(impmat,2125) hfttoi
  write(impmat,2126) keciel
endif

write(impmat,2110) hconve
write(impmat,2111) rconve
write(impmat,2112) hchali
write(impmat,2113) hcheva
write(impmat,2114) ptrres
write(impmat,2115) nptran
write(impmat,2116) netran
write(impmat,2117) frdtra
write(impmat,2118) plgres
write(impmat,2119) nplgrs
write(impmat,2120) nelgrs
write(impmat,2121) epchel
write(impmat,2122) nchest
write(impmat,2123) dmcont


! --- Ecriture des parametres de calcul

write(impmat,*)
write(impmat,1001)
write(impmat,1003)
write(impmat,1001)

write(impmat,3099) ntmabs
write(impmat,3100) iphydr
write(impmat,1151)
write(impmat,3101) dtdtmx


! --- Chargement thermique

write(impmat,*)
write(impmat,1001)
write(impmat,1004)
write(impmat,1001)

write(impmat,3102) imdcnt
write(impmat,1151)
write(impmat,3103) puicon
write(impmat,3104) tinit
write(impmat,3105) tcrit
write(impmat,3106) emicon
write(impmat,3107) emimur


! --- Chargement hydraulique

write(impmat,*)
write(impmat,1001)
write(impmat,1005)
write(impmat,1001)

write(impmat,4107) icofor
if (iconlg.eq.1) then
  WRITE(IMPMAT,4108) ' En ligne        '
else
  WRITE(IMPMAT,4108) ' Pas triangulaire'
endif
write(impmat,4109) ialveo
write(impmat,4110) debmas
write(impmat,4111) pdccha
write(impmat,4112) pdcfch
write(impmat,4113) dhchea
write(impmat,4114) sdchea
write(impmat,4115) pdcche
write(impmat,4116) pdccch
write(impmat,4117) dhches
write(impmat,4118) sdches
write(impmat,4119) pdcalg
write(impmat,4120) pdcatv
write(impmat,4121) argamt
write(impmat,4122) pdcslg
write(impmat,4123) pdcstv
write(impmat,4124) argavl
write(impmat,4125) amppdc
write(impmat,4126) dpvent
if (ialveo.eq.1) then
  write(impmat,4127) dhalve
endif


! --- Pertes de charge : cartes

write(impmat,*)
write(impmat,1001)
write(impmat,1006)
write(impmat,1001)

write(impmat,5100)
write(impmat,5000)
do ii = 1, nzocar(irange,icpdce)
  write(impmat,5003) vizcar(1, ii, irange, icpdce),               &
                     vizcar(2, ii, irange, icpdce)
enddo
write(impmat,5001)
do ii = 1, nzocar(iligne,icpdce)
  write(impmat,5003) vizcar(1, ii, iligne, icpdce),               &
                     vizcar(2, ii, iligne, icpdce)
enddo
write(impmat,5002)
do ii = 1, nzocar(ialtit,icpdce)
  write(impmat,5003) vizcar(1, ii, ialtit, icpdce),               &
                     vizcar(2, ii, ialtit, icpdce)
enddo

write(impmat,5101)
write(impmat,5000)
do ii = 1, nzocar(irange,icpdcs)
  write(impmat,5003) vizcar(1, ii, irange, icpdcs),               &
                     vizcar(2, ii, irange, icpdcs)
enddo
write(impmat,5001)
do ii = 1, nzocar(iligne,icpdcs)
  write(impmat,5003) vizcar(1, ii, iligne, icpdcs),               &
                     vizcar(2, ii, iligne, icpdcs)
enddo
write(impmat,5002)
do ii = 1, nzocar(ialtit,icpdcs)
  write(impmat,5003) vizcar(1, ii, ialtit, icpdcs),               &
                     vizcar(2, ii, ialtit, icpdcs)
enddo

write(impmat,5102)
write(impmat,5000)
do ii = 1, nzocar(irange,icpdcr)
  write(impmat,5003) vizcar(1, ii, irange, icpdcr),               &
                     vizcar(2, ii, irange, icpdcr)
enddo
write(impmat,5001)
do ii = 1, nzocar(iligne,icpdcr)
  write(impmat,5003) vizcar(1, ii, iligne, icpdcr),               &
                     vizcar(2, ii, iligne, icpdcr)
enddo


! --- Sources de chaleur homogeneisees : cartes

write(impmat,*)
write(impmat,1001)
write(impmat,1007)
write(impmat,1001)

write(impmat,5000)
do ii = 1, nzocar(irange,icpuis)
  write(impmat,5004) vizcar(1, ii, irange, icpuis),               &
                     vizcar(2, ii, irange, icpuis),               &
                     vcarth(ii, irange)

enddo
write(impmat,5001)
do ii = 1, nzocar(iligne,icpuis)
  write(impmat,5004) vizcar(1, ii, iligne, icpuis),               &
                     vizcar(2, ii, iligne, icpuis),               &
                     vcarth(ii, iligne)

enddo
write(impmat,5002)
do ii = 1, nzocar(ialtit,icpuis)
  write(impmat,5004) vizcar(1, ii, ialtit, icpuis),               &
                     vizcar(2, ii, ialtit, icpuis),               &
                     vcarth(ii, ialtit)

enddo


! --- Affichage banniere de resultats

write(impmat,*)
write(impmat,1001)
write(impmat,1008)
write(impmat,1001)


! --- Affichage Richardson en convection forcee
!       selon la formule Ri = g beta DeltaT H_ref / U**2 avec :
!       . g = -GZ gravite verticale, valeur positive
!       . beta = 1/T en gaz parfait, ici 1/((TINIT+TSOR0)*0.5 + TKELVI)
!       . DeltaT = TSOR0-TINIT
!       . H_ref = EPCHEL*NCHEST hauteur de la zone de stockage
!       . U = UREEL vitesse reelle horizontale dans la zone de stockage

!       en outre :
!       . TSOR0 = temperature de sortie evaluee a partir de la
!           puissance des conteneurs, ie telle que on ait equilibre entre
!           debit * Cp (TSOR0-TINIT) en Joule/s et
!           NPTRAN*NPLGRS*PUICON en Joule/s
!       . le debit est calcule a partir du debit reel en corrigeant
!           par le facteur d'echelle transverse du modele par rapport au
!           cas reel FRDTRA.

!     En convection naturelle, le Richardson est affiche par mttsns

if (icofor.eq.1) then

  iphas  = 1
  dbm    = debmas/frdtra
  tsor0  = tinit + nptran*nplgrs*puicon/(dbm*cp0(iphas))

  cxsd   = ptrres/dmcont
  ureel  = vitref/(cxsd-1.d0)*cxsd

  betmat = 1.d0/((tinit+tsor0)*0.5d0 + tkelvi)

  hrfmat = epchel*nchest

  richar = -gz*betmat*(tsor0-tinit)*hrfmat/(ureel**2)
  write(impmat,1011) richar

endif


! --- Des resultats seront affiches par mttsns et mtproj
!       le fichier n'est donc pas ferme

!===============================================================================
! FORMATS
!===============================================================================

! --- Texte

 1001 format(74('-'))
 1002 format(7(' '),'Geometrie :')
 1003 format(7(' '),'Parametres de calcul :')
 1004 format(7(' '),'Chargement thermique :')
 1005 format(7(' '),'Chargement hydraulique :')
 1006 format(7(' '),'Zones de pertes de charges :')
 1007 format(7(' '),'Repartition des sources de chaleurs',              &
  ' homogeneisees:')
 1008 format(7(' '),'Resultats :')
 1011 format(' Nombre de Richardson                                   ',&
 ':', E12.5)

 1151 format(' (1 oui, 0 non)')

! --- Geometrie

 2096 format(' Concept d''entrepot                                   ', &
 '  :',A6)
 2097 format(' Hauteur du réseau de conteneurs                       ', &
 ' :',E12.5,' m')
 2098 format(' Hauteur max des alvéoles                              ', &
 ' :',E12.5,' m')
 2099 format(' Hauteur min des alvéoles                              ', &
 ' :',E12.5,' m')
 2100 format(' Epaisseur des registres/cloisons amont et aval        ', &
 ' :',E12.5,' m')
 2101 format(' Epaisseur des cheminees                               ', &
 ' :',E12.5,' m')
 2102 format(' Jeu amont entre cheminee/registre                     ', &
 ' :',E12.5,' m')
 2103 format(' Nb d''elements entre cheminee/registre amont          ', &
 '  :',I12)
 2104 format(' Jeu entre registre amont/colis                        ', &
 ' :',E12.5,' m')
 2105 format(' Nb d''elements entre registre amont/reseau de colis   ', &
 '  :',I12)
 2106 format(' Jeu entre colis/registre aval                         ', &
 ' :',E12.5,' m')
 2107 format(' Nb d''elements entre colis/registre aval              ', &
 '  :',I12)
 2108 format(' Jeu aval entre registre/cheminee                      ', &
 ' :',E12.5,' m')
 2109 format(' Nb d''elements entre registre/cheminee aval           ', &
 '  :',I12)
 2110 format(' Hauteur du convergent                                 ', &
 ' :',E12.5,' m')
 2111 format(' Rapport du convergent                                 ', &
 ' :',E12.5,' m')
 2112 format(' Hauteur de la cheminee d''alimentation                ', &
 '  :',E12.5,' m')
 2113 format(' Hauteur de la cheminee d''evacuation                  ', &
 '  :',E12.5,' m')
 2114 format(' Pas transversal du reseau de conteneur                ', &
 ' :',E12.5,' m')
 2115 format(' Nombre de pas d''espace transversal                   ', &
 '  :',I12)
 2116 format(' Nombre d''elements par pas transversal                ', &
 '  :',I12)
 2117 format(' Facteur de reduction transversal du modele/reel       ', &
 ' :',E12.5)
 2118 format(' Pas longitudinal du reseau de conteneur               ', &
 ' :',E12.5,' m')
 2119 format(' Nombre de pas d''espace longitudinal                  ', &
 '  :',I12)
 2120 format(' Nombre d''elements par pas longitudinal               ', &
 '  :',I12)
 2121 format(' Epaisseur d''une couche d''element (zone stockage)    ', &
 '   :',E12.5,' m')
 2122 format(' Nombre de couche d''element dans la zone stockage     ', &
 '  :',I12)
 2123 format(' Diametre des conteneurs                               ', &
 ' :',E12.5,' m')
 2124 format(' Hauteur du bord de toit                               ', &
 ' :',E12.5,' m')
 2125 format(' Hauteur du faite de toit                              ', &
 ' :',E12.5,' m')
 2126 format(' Nombre de couches d''elements du ciel d''entrepot     ', &
 '   :',I12)

! --- Parametres de calcul

 3099 format(' Nombre de pas de temps                                ', &
 ' :',I12)
 3100 format(' Prise en compte de la pression hydrostatique          ', &
 ' :',I12)
 3101 format(' delta temperature max/pas de temps                    ', &
 ' :',E12.5)
 3102 format(' Modelisation des panaches de convection naturelle     ', &
 ' :',I12)
 3103 format(' Puissance d''un conteneur                             ', &
 '  :',E12.5,' W')
 3104 format(' Temperature d''air en entree                          ', &
 '  :',E12.5,' °C')
 3105 format(' Temperature d''air de sortie critique                 ', &
 '  :',E12.5,' °C')
 3106 format(' Emissivite des conteneurs                             ', &
 ' :',E12.5)
 3107 format(' Emissivite des murs                                   ', &
 ' :',E12.5)

! --- Chargement hydraulique

 4107 format(' Regime hydraulique de circulation forcee (1 oui, 0 non', &
 '):',I12)
 4108 format(' Reseau de conteneur                                   ', &
 ' :',A17)
 4109 format(' Entreposage en alveole (1 oui, 0 non)                 ', &
 ' :',I12)
 4110 format(' Debit de circulation forcee                           ', &
 ' :',E12.5, ' kg/s')
 4111 format(' Perte de charge du diffuseur de cheminee d''alimentati', &
 'on:',E12.5)
 4112 format(' Perte de charge du filtre de cheminee d''alimentation ', &
 '  :',E12.5)
 4113 format(' Diametre hydraulique de la cheminee d''alimentation   ', &
 '  :',E12.5, ' m')
 4114 format(' Surface debitante de la cheminee d''alimentation      ', &
 '  :',E12.5, ' m²')
 4115 format(' Perte de charge du diffuseur de cheminee d''evacuation', &
 '  :',E12.5)
 4116 format(' Perte de charge du clapet de cheminee d''evacuation   ', &
 '  :',E12.5)
 4117 format(' Diametre hydraulique de la cheminee d''evacuation     ', &
 '  :',E12.5, ' m')
 4118 format(' Surface debitante de la cheminee d''evacuation        ', &
 '  :',E12.5, ' m²')
 4119 format(' Perte de charge porte d''entree AMONT longitudinale   ', &
 '  :',E12.5)
 4120 format(' Perte de charge porte d''entree AMONT transversale    ', &
 '  :',E12.5)
 4121 format(' Angle d''inclinaison du registre AMONT (degre)        ', &
 '  :',E12.5, ' °')
 4122 format(' Perte de charge porte de sortie AVAL longitudinale    ', &
 ' :',E12.5)
 4123 format(' Perte de charge porte de sortie AVAL transversale     ', &
 ' :',E12.5)
 4124 format(' Angle d''inclinaison du registre AMONT (degre)        ', &
 '  :',E12.5, ' °')
 4125 format(' Amplificateur des pertes de charge de reseau          ', &
 ' :',E12.5)
 4126 format(' Differentiel de pression entree/sortie                ', &
 ' :',E12.5, ' Pa')
 4127 format(' Diametre hydraulique de l''alveole                    ', &
 '  :',E12.5, ' m')


! --- Cartes

 5000 format('   - Direction longitudinale',                            &
       ' (distance en nombre de rangees)')
 5001 format('   - Direction transversale',                             &
       ' (distance en nombre de lignes)')
 5002 format('   - Direction verticale',                                &
       ' (hauteur en nombre de mailles)')

 5003 format('       * min :', E12.5,' max :', E12.5)
 5004 format('       * min :', E12.5,' max :', E12.5,' val :', E12.5)

 5100 format(' Entree')
 5101 format(' Sortie')
 5102 format(' Réseau')

! --- Erreurs

 9998 format(/,                                                   &
'Code_Saturne : Erreur d''initialisation :',/,              &
'Impossible d''ouvrir le fichier : ',A,/)
 9999 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET MATISSE (MTIMPI)                      ',/,&
'@    =========                                               ',/,&
'@      ARRET A L''OUVERTURE DU FICHIER RESUME MATISSE        ',/,&
'@                                                            ',/,&
'@  Le fichier ',A15,' ne peut etre ouvert.                   ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute.                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)


return
end
