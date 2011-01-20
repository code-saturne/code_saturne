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

subroutine letgeo &
!================

 ( ndim   , ncelet , ncel   , nfac   , nfabor , nfml   , nprfml , &
   nnod   , lndfac , lndfbr ,                                     &
   ntetra , npyram , nprism , nhexae , inodal ,                   &
   ifacel , ifabor , ifmfbr , ifmcel , iprfml ,                   &
   icotet , icopyr , icopri , icohex ,                            &
   ipnfac , nodfac , ipnfbr , nodfbr ,                            &
   xyzcen , surfac , surfbo , cdgfac , cdgfbo , xyznod )

!===============================================================================
!  FONCTION :
!  --------

! LECTURE DES TABLEAUX ENTITES GEOMETRIQUES
!ARGU
!-------------------------------------------------------------------------------
! Arguments
!__________________.____._____.________________________________________________.
! name             !type!mode ! role                                           !
!__________________!____!_____!________________________________________________!
! ndim             ! i  ! <-- ! spatial dimension                              !
! ncelet           ! i  ! <-- ! number of extended (real + ghost) cells        !
! ncel             ! i  ! <-- ! number of cells                                !
! nfac             ! i  ! <-- ! number of interior faces                       !
! nfabor           ! i  ! <-- ! number of boundary faces                       !
! nfml             ! i  ! <-- ! number of families (group classes)             !
! nprfml           ! i  ! <-- ! number of properties per family (group class)  !
! nnod             ! i  ! <-- ! number of vertices                             !
! lndfac           ! i  ! <-- ! size of nodfac indexed array                   !
! lndfbr           ! i  ! <-- ! size of nodfbr indexed array                   !
! ntetra           ! e  ! <-- ! nombre de tetraedres du maillage               !
! npyram           ! e  ! <-- ! nombre de pyramides  du maillage               !
! nprism           ! e  ! <-- ! nombre de prismes    du maillage               !
! nhexae           ! e  ! <-- ! nombre de hexaedres  du maillage               !
! inodal           ! e  ! --> ! indique si l'on doit lire la                   !
!                  !    !     !  connectivite nodale pour le                   !
!                  !    !     !  post traitement                               !
! ifacel           ! te ! --> ! elements voisins d'une face interne            !
! ifabor           ! te ! --> ! element  voisin  d'une face de bord            !
! (nfabor)         !    !     !                                                !
! ifmfbr           ! te ! --> ! numero de famille d'une face de bord           !
! (nfabor)         !    !     !                                                !
! ifmcel           ! te ! --> ! numero de famille d'une cellule                !
! (ncelet)         !    !     !                                                !
! iprfml           ! te ! --> ! proprietes d'une famille                       !
! nfml  ,nprfml    !    !     !                                                !
! icotet           ! te ! --> ! connectivite tetraedres-noeuds                 !
! icopyr           ! te ! --> ! connectivite pyramides-noeuds                  !
! icopri           ! te ! --> ! connectivite prismes-noeuds                    !
! icohex           ! te ! --> ! connectivite hexaedres-noeuds                  !
! ipnfac           ! te ! --- ! position du premier noeud de chaque            !
!   (lndfac)       !    !     !  face interne dans nodfac (optionnel)          !
! nodfac           ! te ! --- ! connectivite faces internes/noeuds             !
!   (nfac+1)       !    !     !  (optionnel)                                   !
! ipnfbr           ! te ! --- ! position du premier noeud de chaque            !
!   (lndfbr)       !    !     !  face de bord dans nodfbr (optionnel)          !
! nodfbr           ! te ! --- ! connectivite faces de bord/noeuds              !
!   (nfabor+1)     !    !     !  (optionnel)                                   !
! xyzcen           ! tr ! --> ! point associes aux volumes de control          !
! (ndim,ncelet     !    !     !                                                !
! surfac           ! tr ! --> ! vecteur surface des faces internes             !
! (ndim,nfac)      !    !     !                                                !
! surfbo           ! tr ! --> ! vecteur surface des faces de bord              !
! (ndim,nfabor)    !    !     !                                                !
! cdgfac           ! tr ! --> ! centre de gravite des faces internes           !
! (ndim,nfac)      !    !     !                                                !
! cdgfbo           ! tr ! --> ! centre de gravite des faces de bord            !
! (ndim,nfabor)    !    !     !                                                !
! xyznod           ! tr ! --- ! coordonnes des noeuds                          !
! (ndim,nnod)      !    !     !                                                !
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
use entsor

!===============================================================================

implicit none

! Arguments

integer          ndim   , ncelet , ncel   , nfac   , nfabor
integer          nfml   , nprfml
integer          nnod   , lndfac , lndfbr
integer          ntetra , npyram , nprism , nhexae , inodal

integer          ifacel(2,nfac) , ifabor(nfabor)
integer          ifmfbr(nfabor) , ifmcel(ncelet)
integer          iprfml(nfml,nprfml)
integer          icotet(4,ntetra), icopyr(5,npyram)
integer          icopri(6,nprism), icohex(8,nhexae)
integer          ipnfac(nfac+1), nodfac(lndfac)
integer          ipnfbr(nfabor+1), nodfbr(lndfbr)

double precision xyzcen(ndim,ncelet)
double precision surfac(ndim,nfac)
double precision surfbo(ndim,nfabor)
double precision cdgfac(ndim,nfac)
double precision cdgfbo(ndim,nfabor)
double precision xyznod(ndim,nnod)

! Local variables

integer          kel, iel, idim, kface, iface, kfafbr, kprffb
integer          iok, iok1, iok2, ind
integer          jdim,jcel ,jfac  ,jfabor
integer          jprffb, jfafbr
integer          jpoint, jtetra, jpyram, jprism, jhexae
integer          ifac, nn, n1, n2, ip, ii1, ii2
integer          nbfac1, nbfac2

!===============================================================================

iok = 0

inodal = 0

!===============================================================================
! 1.  OUVERTURE
!===============================================================================

open (file=ficgeo, unit=impgeo, form='formatted')


!===============================================================================
! 2.  LECTURE DES DIMENSIONS ET VERIFICATIONS
!===============================================================================

read (impgeo,   *)
read (impgeo,   *)
jdim = 3
read (impgeo,1100)  jcel, jfac, jfabor, jpoint

if (jdim.ne.ndim.or.jcel .ne.ncel ) then
  write(nfecra,2010) '  jdim',jdim,'  jcel',jcel ,                &
                     '  ndim',ndim,'  ncel',ncel
  iok = iok + 1
endif
if (jfac  .ne.nfac  .or.jfabor.ne.nfabor) then
  write(nfecra,2010) '  jfac',jfac  ,'jfabor',jfabor,             &
                     '  nfac',nfac  ,'nfabor',nfabor
  iok = iok + 1
endif
if (jpoint.ne.nnod) then
  write(nfecra,2011) 'jpoint',jpoint,'nnod  ',nnod
  iok = iok + 1
endif

read (impgeo,   *)
read (impgeo,   *)
read (impgeo,1100) jtetra,jpyram,jprism,jhexae

if (jtetra.ne.ntetra.or.jpyram.ne.npyram) then
  write(nfecra,2010) 'jtetra',jtetra,'jpyram',jpyram,             &
                     'ntetra',ntetra,'npyram',npyram
  iok = iok + 1
endif
if (jprism.ne.nprism.or.jhexae.ne.nhexae) then
  write(nfecra,2010) 'jprism',jprism,'jhexae',jhexae,             &
                     'nprism',nprism,'nhexae',nhexae
  iok = iok + 1
endif

read (impgeo,   *)
read (impgeo,   *)
read (impgeo,1100)  jprffb, jfafbr

if (jprffb.ne.nprfml.or.jfafbr.ne.nfml  ) then
  write(nfecra,2010) 'jprffb',jprffb,'jfafbr',jfafbr,             &
                     'nprfml',nprfml,'nfml  ',nfml
  iok = iok + 1
endif

if (iok.ne.0) then
  write(nfecra,9999)
  call csexit (1)
endif


!===============================================================================
! 3.  LECTURE DES TABLEAUX PRINCIPAUX
!===============================================================================

read (impgeo,   *)
read (impgeo,   *)
do kface = 1, nfac
  read (impgeo,1100) iface,(ifacel(iel,kface),iel=1,2)
enddo

read (impgeo,   *)
read (impgeo,   *)
do kface = 1, nfabor
  read (impgeo,1100) iface,ifabor(kface)
enddo

read (impgeo,   *)
read (impgeo,   *)
do kel = 1, ncel
  read (impgeo,1120) iel,xyzcen(1,kel),xyzcen(2,kel),xyzcen(3,kel)
enddo

read (impgeo,   *)
read (impgeo,   *)
do kface = 1, nfac
  read (impgeo,1120) iface,(surfac(idim,kface),idim=1,ndim)
enddo

read (impgeo,   *)
read (impgeo,   *)
do kface = 1, nfabor
  read (impgeo,1120) iface,(surfbo(idim,kface),idim=1,ndim)
enddo

read (impgeo,   *)
read (impgeo,   *)
do kface = 1, nfac
  read (impgeo,1120) iface,(cdgfac(idim,kface),idim=1,ndim)
enddo

read (impgeo,   *)
read (impgeo,   *)
do kface = 1, nfabor
  read (impgeo,1120) iface,(cdgfbo(idim,kface),idim=1,ndim)
enddo

read (impgeo,   *)
read (impgeo,   *)
do kface = 1, nfabor
  read (impgeo,1100) iface,ifmfbr(kface)
enddo

read (impgeo,   *)
read (impgeo,   *)
do kfafbr = 1, nfml
  read (impgeo,1100)                                              &
       jfafbr,(iprfml(jfafbr,kprffb),kprffb=1,nprfml)
enddo

!     IFMCEL NON LU ICI, ON INITIALISE A LA PLACE

do kel = 1, ncel
  ifmcel(kel) = 0
enddo

!===============================================================================
! 4.  LECTURE DES POINTS SUPPORT ET DE LA CONNECTIVITE NODALE
!===============================================================================

!     Selon si l'on dispose des connectivites faces->sommets ou non,
!     on pourra reconstruire une connectivite nodale pour le post
!     traitement. Si celui-ci est demande et que cette connectivite
!     n'est pas disponible, on utilise la connectivite nodale
!     fournie dans le fichier solcom.

if (lndfac.eq.0 .and. lndfbr.eq.0 .and. ichrvl.gt.0) then
  inodal = 1
else
  inodal = 0
endif

!     Coordonnees des sommets
!     -----------------------

!     Selon si l'on dispose des connectivites faces->sommets ou non,
!     le tableau XYZNOD passe en argument appartient a la connectivite
!     complete ou au maillage nodal destine au post traitement.

read(impgeo,*)
read(impgeo,*)

if ((lndfac.gt.0 .or. lndfbr.gt.0) .or. inodal.eq.1) then

  do n1 = 1,nnod
    read(impgeo,1120) nn,(xyznod(ip,nn),ip=1,3)
  enddo

else

  do n1 = 1,nnod
    read(impgeo,*)
  enddo

endif

!     Connectivites cellules - sommets
!     --------------------------------

!     On n'utilise ces connectivites que lorsque l'on ne dispose pas
!     des connectivites faces -> sommets, et que l'on ne pourra donc
!     pas reconstruire une connectivite nodale pour le post traitement.

read (impgeo,*)
read (impgeo,*)

if (inodal .eq. 1) then

  do n1 = 1, ntetra
    read(impgeo,1100) nn,(icotet(ip,n1),ip=1,4)
  enddo
  do n1 = 1, npyram
    read(impgeo,1100) nn,(icopyr(ip,n1),ip=1,5)
  enddo
  do n1 = 1, nprism
    read(impgeo,1100) nn,(icopri(ip,n1),ip=1,6)
  enddo
  do n1 = 1, nhexae
    read(impgeo,1100) nn,(icohex(ip,n1),ip=1,8)
  enddo

else

  do n1 = 1, ntetra
    read(impgeo,*)
  enddo
  do n1 = 1, npyram
    read(impgeo,*)
  enddo
  do n1 = 1, nprism
    read(impgeo,*)
  enddo
  do n1 = 1, nhexae
    read(impgeo,*)
  enddo

endif

!===============================================================================
! 5.  LECTURE DES CONNECTIVITES FACE-SOMMET
!===============================================================================

!     On peut avoir besoin de la connectivite faces -> sommets pour
!       les modules rayonnement et Lagrangien. Les tableaux
!       correspondants se trouvent en fin de fichier, après
!       les coordonnees des sommets.

!     Si l'on n'a pas detecte la presence de ces tableaux
!       supplementaires, on sort immediatement.

if (lndfac.eq.0 .and. lndfbr.eq.0) return

!     Remarque : ces rubriques ont deja ete lues dans LEDGEO, pour
!                determiner les dimensions. On ne prend donc pas de
!                precaution quant a leur existence.

! Connectivite faces de bord - sommets
ii2 = 0
read(impgeo,*,err=9000)
read(impgeo,*,err=9000)
do n1 = 1,nfabor
  read(impgeo,'(20i10)',err=9000)                                 &
       nn, ipnfbr(nn), (nodfbr(ii1),ii1=ii2+1,ii2+ipnfbr(nn))
  ii2 = ii2 + ipnfbr(nn)
enddo

! Connectivite faces internes - sommets
ii2 = 0
read(impgeo,*,err=9000)
read(impgeo,*,err=9000)
do n1 = 1,nfac
  read(impgeo,'(20i10)',err=9000)                                 &
       nn, ipnfac(nn), (nodfac(ii1),ii1=ii2+1,ii2+ipnfac(nn))
  ii2 = ii2 + ipnfac(nn)
enddo

!===============================================================================
! 6. CREATION DES TABLEAUX IPNFAC et IPNFBR
!===============================================================================

! --> On calcule les pointeurs des tableaux NODFAC et NODFBR

nbfac2 = ipnfac(1)
ipnfac(1) = 1
do ifac = 2,nfac
  nbfac1 = nbfac2
  nbfac2 = ipnfac(ifac)
  ipnfac(ifac) = ipnfac(ifac-1) + nbfac1
enddo
ipnfac(nfac+1) = ipnfac(nfac) + nbfac2

nbfac2 = ipnfbr(1)
ipnfbr(1) = 1
do ifac = 2, nfabor
  nbfac1 = nbfac2
  nbfac2 = ipnfbr(ifac)
  ipnfbr(ifac) = ipnfbr(ifac-1) + nbfac1
enddo
ipnfbr(nfabor+1) = ipnfbr(nfabor) + nbfac2

!===============================================================================
! 7.  CONTROLES SIMPLES
!===============================================================================

! --> Verifications de IPNFAC, IPNFBR

iok = 0
if (ipnfac(nfac+1).ne.lndfac+1) then
  write(nfecra,2001)nfac  ,lndfac+1,ipnfac(nfac+1)
  iok = iok + 1
endif
if (ipnfbr(nfabor+1).ne.lndfbr+1) then
  write(nfecra,2002)nfabor,lndfbr+1,ipnfbr(nfabor+1)
  iok = iok + 1
endif
if(iok.ne.0) then
  call csexit (1)
endif


! --> connectivite faces internes - sommets
iok1=0
do ifac = 1,nfac
  ind = 0
  do n1 = ipnfac(ifac),ipnfac(ifac+1)-1
    do n2 = ipnfac(ifac),ipnfac(ifac+1)-1
      if (n1.ne.n2 .and.                                          &
         ((xyznod(1,nodfac(n1)).eq.xyznod(1,nodfac(n2))) .and.    &
          (xyznod(2,nodfac(n1)).eq.xyznod(2,nodfac(n2))) .and.    &
          (xyznod(3,nodfac(n1)).eq.xyznod(3,nodfac(n2))) )) then
        ind=1
      endif
    enddo
  enddo
  if (ind.ne.0) then
    do n1 = ipnfac(ifac),ipnfac(ifac+1)-1
      write(nfecra,3011)nodfac(n1),                               &
                        n1,xyznod(1,nodfac(n1)),                  &
                        n1,xyznod(2,nodfac(n1)),                  &
                        n1,xyznod(3,nodfac(n1))
    enddo
    iok1=iok1+1
  endif
enddo

! --> connectivite faces de bord - sommets
iok2=0
do ifac = 1,nfabor
  ind = 0
  do n1 = ipnfbr(ifac),ipnfbr(ifac+1)-1
    do n2 = ipnfbr(ifac),ipnfbr(ifac+1)-1
      if (n1.ne.n2 .and.                                          &
         ((xyznod(1,nodfbr(n1)).eq.xyznod(1,nodfbr(n2))) .and.    &
          (xyznod(2,nodfbr(n1)).eq.xyznod(2,nodfbr(n2))) .and.    &
          (xyznod(3,nodfbr(n1)).eq.xyznod(3,nodfbr(n2))) )) then
        ind=1
      endif
    enddo
  enddo
  if (ind.ne.0) then
    do n1 = ipnfbr(ifac),ipnfbr(ifac+1)-1
      write(nfecra,3012)nodfbr(n1),                               &
                        n1,xyznod(1,nodfbr(n1)),                  &
                        n1,xyznod(2,nodfbr(n1)),                  &
                        n1,xyznod(3,nodfbr(n1))
    enddo
    iok2=iok2+1
  endif
enddo

! --> stop si erreur connectivites

if(iok1.ne.0) then
  write(nfecra,3001)
endif
if(iok2.ne.0) then
  write(nfecra,3002)
endif
if(iok1.ne.0.or.iok2.ne.0) then
  call csexit (1)
endif

!===============================================================================
! 8.  FERMETURE
!===============================================================================

close(impgeo)

return

!===============================================================================
! 9. SORTIES SUR ERREUR DE LECTURE
!===============================================================================

! ---> Stop si le maillage lu est incomplet ou erreur a la lecture.

 9000 write(nfecra,9100)
call csexit (1)

!---> FORMATS

#if defined(_CS_LANG_FR)

 1100 format(20i10)
 1120 format(i10, 5e23.15)

 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA CREATION DES CONNECTIVITES       ',/,&
'@    =========   (LETGEO).                                   ',/,&
'@                                                            ',/,&
'@    Le nombre de faces internes est NFAC   = ',I10           ,/,&
'@    IPNFAC(NFAC  +1)devrait valoir ',I10                     ,/,&
'@                           il vaut ',I10                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@  Verifier LETGEO.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA CREATION DES CONNECTIVITES       ',/,&
'@    =========   (LETGEO).                                   ',/,&
'@                                                            ',/,&
'@    Le nombre de faces de bord  est NFABOR = ',I10           ,/,&
'@    IPNFBR(NFABOR+1)devrait valoir ',I10                     ,/,&
'@                           il vaut ',I10                     ,/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@  Verifier LETGEO.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,&
'@    =========                                               ',/,&
'@      INCOHERENCES RENCONTREES                              ',/,&
'@                                                            ',/,&
'@      Les dimensions lues dans letgeo ne sont pas en accord ',/,&
'@        avec les dimensions lues dans ledgeo.               ',/,&
'@      Lecture de  ',A6,' = ',I10   ,', ',A6,' = ',I10        ,/,&
'@        alors que ',A6,' = ',I10   ,', ',A6,' = ',I10        ,/,&
'@                                                            ',/,&
'@      Le calcul ne peut etre execute.                       ',/,&
'@                                                            ',/,&
'@      Contacter l''assistance.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,&
'@    =========                                               ',/,&
'@      INCOHERENCES RENCONTREES                              ',/,&
'@                                                            ',/,&
'@      Les dimensions lues dans letgeo ne sont pas en accord ',/,&
'@        avec les dimensions lues dans ledgeo.               ',/,&
'@      Lecture de  ',A6,' = ',I10                             ,/,&
'@        alors que ',A6,' = ',I10                             ,/,&
'@                                                            ',/,&
'@      Le calcul ne peut etre execute.                       ',/,&
'@                                                            ',/,&
'@      Contacter l''assistance.                              ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA CREATION DES CONNECTIVITES       ',/,&
'@    =========   (LETGEO)                                    ',/,&
'@                                                            ',/,&
'@    Le test a echoue sur la connectivite                    ',/,&
'@         faces internes -> points supports                  ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@  Verifier LETGEO.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA CREATION DES CONNECTIVITES       ',/,&
'@    =========   (LETGEO)                                    ',/,&
'@                                                            ',/,&
'@    Le test a echoue sur la connectivite                    ',/,&
'@         faces de bord  -> points supports                  ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier le maillage.                                     ',/,&
'@  Verifier LETGEO.                                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3011 format(                                                           &
'@                                                            ',/,&
'@ @@ ATTENTION : LETGEO - CONNECTIVITE : NODFAC = ',I10       ,/,&
'@                   X(',I10,') = ',E23.15                     ,/,&
'@                   Y(',I10,') = ',E23.15                     ,/,&
'@                   Z(',I10,') = ',E23.15                       )
 3012 format(                                                           &
'@                                                            ',/,&
'@ @@ ATTENTION : LETGEO - CONNECTIVITE : NODFBR = ',I10       ,/,&
'@                   X(',I10,') = ',E23.15                     ,/,&
'@                   Y(',I10,') = ',E23.15                     ,/,&
'@                   Z(',I10,') = ',E23.15                       )

 9100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU MAILLAGE              ',/,&
'@    =========                                               ',/,&
'@                                                            ',/,&
'@    LE MAILLAGE LU (.slc ou .tlc) EST INCOMPLET             ',/,&
'@                                                            ',/,&
'@  Le fichier de maillage ne contient pas les informations   ',/,&
'@    relatives a la connectivite faces-noeuds, ou contient   ',/,&
'@    des informations incompletes.                           ',/,&
'@  Ces informations sont necessaires pour effectuer un calcul',/,&
'@    en utilisant le module Lagrangien ou Rayonnement        ',/,&
'@                                                            ',/,&
'@  Le calcul ne peut etre execute.                           ',/,&
'@                                                            ',/,&
'@  Verifier le fichier de maillage.                          ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ ATTENTION : ARRET A LA LECTURE DU FICHIER GEOMETRIE     ',/,&
'@    =========                                               ',/,&
'@      INCOHERENCES RENCONTREES                              ',/,&
'@                                                            ',/,&
'@  Le calcul ne sera pas execute (',I10,' erreurs).          ',/,&
'@                                                            ',/,&
'@  Se reporter aux impressions precedentes pour plus de      ',/,&
'@    renseignements.                                         ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#else

 1100 format(20i10)
 1120 format(i10, 5e23.15)

 2001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE CONNECTIVITY CREATION (LETGEO)    ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    The number of interior faces is NFAC   = ',I10           ,/,&
'@    IPNFAC(NFAC  +1) should be ',I10                         ,/,&
'@                         it is ',I10                         ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh.                                          ',/,&
'@  Verify LETGEO.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 2002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE CONNECTIVITY CREATION (LETGEO)    ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    The number of interior faces is NFABOR = ',I10           ,/,&
'@    IPNFAC(NFABOR+1) should be ',I10                         ,/,&
'@                         it is ',I10                         ,/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh.                                          ',/,&
'@  Verify LETGEO.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2010 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE READING THE GEOMETRY FILE          ',/,&
'@    ========                                                ',/,&
'@      THERE ARE SOME INCOHERENCIES                          ',/,&
'@                                                            ',/,&
'@      The dimensions read in letgeo do not comply with      ',/,&
'@        the dimensions read in ledgeo.                      ',/,&
'@      Reading   ',A6,' = ',I10   ,', ',A6,' = ',I10          ,/,&
'@        whereas ',A6,' = ',I10   ,', ',A6,' = ',I10          ,/,&
'@                                                            ',/,&
'@      The calculation will not be run.                      ',/,&
'@                                                            ',/,&
'@      Contact the support.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 2011 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE READING THE GEOMETRY FILE          ',/,&
'@    ========                                                ',/,&
'@      THERE ARE SOME INCOHERENCIES                          ',/,&
'@                                                            ',/,&
'@      The dimensions read in letgeo do not comply with      ',/,&
'@        the dimensions read in ledgeo.                      ',/,&
'@      Reading   ',A6,' = ',I10                               ,/,&
'@        whereas ',A6,' = ',I10                               ,/,&
'@                                                            ',/,&
'@      The calculation will not be run.                      ',/,&
'@                                                            ',/,&
'@      Contact the support.                                  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3001 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE CONNECTIVITY CREATION (LETGEO)    ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    The test failed on the interior faces -> vertices       ',/,&
'@      connectivity.                                         ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh.                                          ',/,&
'@  Verify LETGEO.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 3002 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT IN THE CONNECTIVITY CREATION (LETGEO)    ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    The test failed on the boundary faces -> vertices       ',/,&
'@      connectivity.                                         ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh.                                          ',/,&
'@  Verify LETGEO.                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

 3011 format(                                                           &
'@                                                            ',/,&
'@ @@ WARNING: LETGEO - CONNECTIVITY : NODFAC = ',I10          ,/,&
'@                   X(',I10,') = ',E23.15                     ,/,&
'@                   Y(',I10,') = ',E23.15                     ,/,&
'@                   Z(',I10,') = ',E23.15                       )
 3012 format(                                                           &
'@                                                            ',/,&
'@ @@ WARNING: LETGEO - CONNECTIVITY : NODFBR = ',I10          ,/,&
'@                   X(',I10,') = ',E23.15                     ,/,&
'@                   Y(',I10,') = ',E23.15                     ,/,&
'@                   Z(',I10,') = ',E23.15                       )

 9100 format(                                                           &
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE READING THE MESH                   ',/,&
'@    ========                                                ',/,&
'@                                                            ',/,&
'@    THE MESH (.slc or .tlc) IS INCOMPLETE                   ',/,&
'@                                                            ',/,&
'@  The mesh file does not contain the information relative   ',/,&
'@    to the face-vertices connectivity, or contains          ',/,&
'@    incomplete informations.                                ',/,&
'@  These informations are mandatory to run a simulation      ',/,&
'@    using the Lagrangian or radiative transfer module.      ',/,&
'@                                                            ',/,&
'@  The calculation will not be run.                          ',/,&
'@                                                            ',/,&
'@  Verify the mesh file.                                     ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)
 9999 format(                                                           &
'@                                                            ',/,&
'@                                                            ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/,&
'@ @@ WARNING: ABORT WHILE READING THE GEOMETRY FILE          ',/,&
'@    ========                                                ',/,&
'@      THERE ARE SOME INCOHERENCIES                          ',/,&
'@                                                            ',/,&
'@  The calculation will not be run (',I10,' errors).         ',/,&
'@                                                            ',/,&
'@  Refer to the previous warnings for further informations.  ',/,&
'@                                                            ',/,&
'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,&
'@                                                            ',/)

#endif

!----
! FIN
!----

end subroutine
