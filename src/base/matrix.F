c@a
c@versb
C-----------------------------------------------------------------------
C
C     This file is part of the Code_Saturne Kernel, element of the
C     Code_Saturne CFD tool.
C
C     Copyright (C) 1998-2008 EDF S.A., France
C
C     contact: saturne-support@edf.fr
C
C     The Code_Saturne Kernel is free software; you can redistribute it
C     and/or modify it under the terms of the GNU General Public License
C     as published by the Free Software Foundation; either version 2 of
C     the License, or (at your option) any later version.
C
C     The Code_Saturne Kernel is distributed in the hope that it will be
C     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
C     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C     GNU General Public License for more details.
C
C     You should have received a copy of the GNU General Public License
C     along with the Code_Saturne Kernel; if not, write to the
C     Free Software Foundation, Inc.,
C     51 Franklin St, Fifth Floor,
C     Boston, MA  02110-1301  USA
C
C-----------------------------------------------------------------------
c@verse
                        SUBROUTINE MATRIX
C                       *****************
C     -------------------------------------------------------------
     & ( NCELET , NCEL   , NFAC   , NFABOR ,
     &   ICONVP , IDIFFP , NDIRCP , ISYM   , NFECRA ,
     &   THETAP ,
     &   IFACEL , IFABOR ,
     &   COEFBP , ROVSDT , FLUMAS , FLUMAB , VISCF  , VISCB  ,
     &   DA     , XA     )
C     -------------------------------------------------------------
C***********************************************************************
C FONCTION :
C ----------
c@foncb
CFONC
CFONC CONSTRUCTION DE LA MATRICE DE CONVECTION UPWIND/DIFFUSION/TS
CFONC
CFONC      IL EST INTERDIT DE MODIFIER ROVSDT ICI
CFONC
c@fonce
C
C-----------------------------------------------------------------------
c@argub
CARGU                             ARGUMENTS
CARGU .______________.____._____.______________________________________.
CARGU !    NOM       !TYPE!MODE !                   ROLE               !
CARGU !______________!____!_____!______________________________________!
CARGU ! NCELET       ! E  !  -> ! NOMBRE D'ELEMENTS HALO COMPRIS       !
CARGU ! NCEL         ! E  !  -> ! NOMBRE D'ELEMENTS ACTIFS             !
CARGU ! NFAC         ! E  !  -> ! NOMBRE DE FACES INTERNES             !
CARGU ! NFABOR       ! E  !  -> ! NOMBRE DE FACES DE BORD              !
CARGU ! ICONVP       ! E  !  -> ! INDICATEUR = 1 CONVECTION, 0 SINON   !
CARGU ! IDIFFP       ! E  !  -> ! INDICATEUR = 1 DIFFUSION , 0 SINON   !
CARGU ! NDIRCP       ! E  !  -> ! INDICATEUR = 0 SI DECALAGE DIAGONALE !
CARGU ! ISYM         ! E  !  -> ! INDICATEUR = 1 MATRICE SYMETRIQUE    !
CARGU !              !    !     !              2 MATRICE NON SYMETRIQUE!
CARGU ! THETAP       ! R  !  -> ! COEFFICIENT DE PONDERATION POUR LE   !
CARGU !              !    !     ! THETA-SCHEMA (ON NE L'UTILISE POUR LE!
CARGU !              !    !     ! MOMENT QUE POUR U,V,W ET LES SCALAIRE!
CARGU !              !    !     ! - THETAP = 0.5 CORRESPOND A UN SCHEMA!
CARGU !              !    !     !   TOTALEMENT CENTRE EN TEMPS (MIXAGE !
CARGU !              !    !     !   ENTRE CRANK-NICOLSON ET ADAMS-     !
CARGU !              !    !     !   BASHFORTH)                         !
CARGU ! IFACEL(2,NFAC! TE !  -> ! No DES ELTS VOISINS D'UNE FACE INTERN!
CARGU ! IFABOR(NFABOR! TE !  -> ! No DE L'ELT VOISIN D'UNE FACE DE BORD!
CARGU ! COEFBP(NFABOR! TR !  -> ! TAB B DES CL POUR LA VAR CONSIDEREE  !
CARGU ! FLUMAS(NFAC) ! TR !  -> ! FLUX DE MASSE AUX FACES INTERNES     !
CARGU ! FLUMAB(NFABOR! TR !  -> ! FLUX DE MASSE AUX FACES DE BORD      !
CARGU ! VISCF(NFAC)  ! TR !  -> ! VISC*SURFACE/DIST AUX FACES INTERNES !
CARGU ! VISCB(NFABOR ! TR !  -> ! VISC*SURFACE/DIST AUX FACES DE BORD  !
CARGU ! DA (NCELET   ! TR ! <-  ! PARTIE DIAGONALE DE LA MATRICE       !
CARGU ! XA (NFAC,*)  ! TR ! <-  ! EXTRA  DIAGONALE DE LA MATRICE       !
CARGU !______________!____!_____!______________________________________!
c@argue
C
c@commb
CCOMM                             COMMONS
CCOMM .______________.____._____.______________________________________.
CCOMM !    NOM       !TYPE!MODE !                   ROLE               !
CCOMM !______________!____!_____!______________________________________!
CCOMM !______________!____!_____!______________________________________!
c@comme
C
C     TYPE : E (ENTIER), R (REEL), A (ALPHANUMERIQUE), T (TABLEAU)
C            L (LOGIQUE)   .. ET TYPES COMPOSES (EX : TR TABLEAU REEL)
C     MODE : -> DONNEE, <- RESULTAT, <-> DONNEE MODIFIEE,
C            - TABLEAU DE TRAVAIL
C***********************************************************************
C
      IMPLICIT NONE
C
C***********************************************************************
C     DONNEES EN COMMON
C***********************************************************************
C
      INCLUDE "vector.h"
C
C***********************************************************************
C
C
C ARGUMENTS
C
      INTEGER          NCELET , NCEL   , NFAC   , NFABOR
      INTEGER          ICONVP , IDIFFP , NDIRCP , ISYM
      INTEGER          NFECRA
      DOUBLE PRECISION THETAP
C
      INTEGER          IFACEL(2,NFAC), IFABOR(NFABOR)
      DOUBLE PRECISION COEFBP(NFABOR), ROVSDT(NCELET)
      DOUBLE PRECISION FLUMAS(NFAC), FLUMAB(NFABOR)
      DOUBLE PRECISION VISCF(NFAC), VISCB(NFABOR)
      DOUBLE PRECISION DA(NCELET ),XA(NFAC ,ISYM)
C
C VARIABLES LOCALES
C
      INTEGER          IFAC,II,JJ,IEL
      DOUBLE PRECISION FLUI,FLUJ,EPSI
C
C***********************************************************************
C
C=======================================================================
C 1. INITIALISATION
C=======================================================================
C
      IF(ISYM.NE.1.AND.ISYM.NE.2) THEN
        WRITE(NFECRA,1000) ISYM
        CALL CSEXIT (1)
      ENDIF
C
      EPSI = 1.D-7
C
      DO IEL = 1, NCEL
        DA(IEL) = ROVSDT(IEL)
      ENDDO
      IF(NCELET.GT.NCEL) THEN
        DO IEL = NCEL+1, NCELET
          DA(IEL) = 0.D0
        ENDDO
      ENDIF
C
      IF(ISYM.EQ.2) THEN
        DO IFAC = 1, NFAC
          XA(IFAC,1) = 0.D0
          XA(IFAC,2) = 0.D0
        ENDDO
      ELSE
        DO IFAC = 1, NFAC
          XA(IFAC,1) = 0.D0
        ENDDO
      ENDIF
C
C=======================================================================
C 2.    CALCUL DES TERMES EXTRADIAGONAUX
C=======================================================================
C
      IF(ISYM.EQ.2) THEN
C
        DO IFAC = 1,NFAC
          FLUI = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
          FLUJ =-0.5D0*( FLUMAS(IFAC) +ABS(FLUMAS(IFAC)) )
          XA(IFAC,1) = THETAP*(ICONVP*FLUI -IDIFFP*VISCF(IFAC))
          XA(IFAC,2) = THETAP*(ICONVP*FLUJ -IDIFFP*VISCF(IFAC))
        ENDDO
C
      ELSE
C
        DO IFAC = 1,NFAC
          FLUI = 0.5D0*( FLUMAS(IFAC) -ABS(FLUMAS(IFAC)) )
          XA(IFAC,1) = THETAP*(ICONVP*FLUI -IDIFFP*VISCF(IFAC))
        ENDDO
C
      ENDIF
C
C=======================================================================
C 3.     CONTRIBUTION DES TERMES X-TRADIAGONAUX A LA DIAGONALE
C=======================================================================
C
      IF(ISYM.EQ.2) THEN
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            DA(II) = DA(II) -XA(IFAC,2)
            DA(JJ) = DA(JJ) -XA(IFAC,1)
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            DA(II) = DA(II) -XA(IFAC,2)
            DA(JJ) = DA(JJ) -XA(IFAC,1)
          ENDDO
C
        ENDIF
C
      ELSE
C
        IF (IVECTI.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            DA(II) = DA(II) -XA(IFAC,1)
            DA(JJ) = DA(JJ) -XA(IFAC,1)
          ENDDO
C
        ELSE
C
C VECTORISATION NON FORCEE
          DO IFAC = 1,NFAC
            II = IFACEL(1,IFAC)
            JJ = IFACEL(2,IFAC)
            DA(II) = DA(II) -XA(IFAC,1)
            DA(JJ) = DA(JJ) -XA(IFAC,1)
          ENDDO
C
        ENDIF
C
      ENDIF
C
C=======================================================================
C 4.     CONTRIBUTION DES FACETTES DE BORDS A LA DIAGONALE
C=======================================================================
C
      IF (IVECTB.EQ.1) THEN
C
!OCL NOVREC,VRL(16)
        DO IFAC=1,NFABOR
          II = IFABOR(IFAC)
          FLUI = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
          FLUJ =-0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
          DA(II) = DA(II) + THETAP*(
     &                     ICONVP*(-FLUJ + FLUI*COEFBP(IFAC) )
     &                    +IDIFFP*VISCB(IFAC)*(1.D0-COEFBP(IFAC))
     &                             )
        ENDDO
C
      ELSE
C
C VECTORISATION NON FORCEE
        DO IFAC=1,NFABOR
          II = IFABOR(IFAC)
          FLUI = 0.5D0*( FLUMAB(IFAC) -ABS(FLUMAB(IFAC)) )
          FLUJ =-0.5D0*( FLUMAB(IFAC) +ABS(FLUMAB(IFAC)) )
          DA(II) = DA(II) + THETAP*(
     &                     ICONVP*(-FLUJ + FLUI*COEFBP(IFAC) )
     &                    +IDIFFP*VISCB(IFAC)*(1.D0-COEFBP(IFAC))
     &                             )
        ENDDO
C
      ENDIF
C
C=======================================================================
C 5.  NON PRESENCE DE PTS DIRICHLET --> LEGER RENFORCEMENT DE LA
C     DIAGONALE POUR DECALER LE SPECTRE DES VALEURS PROPRES
C=======================================================================
C     (si IDIRCL=0, on a force NDIRCP a valoir au moins 1 pour ne pas
C      decaler la diagonale)
C
      IF ( NDIRCP.LE.0 ) THEN
        DO IEL=1,NCEL
          DA(IEL) = (1.D0+EPSI)*DA(IEL)
        ENDDO
      ENDIF
C
C--------
C FORMATS
C--------
C
 1000 FORMAT(
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/,
     &'@ @@ ATTENTION : ARRET DANS matrix                           ',/,
     &'@    *********                                               ',/,
     &'@     APPEL DE matrix              AVEC ISYM   = ',I10        ,/,
     &'@                                                            ',/,
     &'@  Le calcul ne peut pas etre execute.                       ',/,
     &'@                                                            ',/,
     &'@  Contacter l''assistance.                                  ',/,
     &'@                                                            ',/,
     &'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',/,
     &'@                                                            ',/)
C
C----
C FIN
C----
C
      RETURN
C
      END
c@z
