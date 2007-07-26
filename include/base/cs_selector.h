#ifndef __CS_SELECTOR_H__
#define __CS_SELECTOR_H__

/*============================================================================
 * Structure principale associée à un maillage
 *
 *  Bibliothèque : Code_Saturne 1.3                    Copyright EDF 1999-2006
 *============================================================================*/

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfbr, CSGFBR)
(
 const char          *const fstr,          /* <-- Fortran string */
 int                 *const len,           /* <-- String Length  */
 int                 *const faces_number,   /* --> faces number */
 int                 *const faces         /* --> faces  */
 CS_ARGF_SUPP_CHAINE
 );

/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfac, CSGFAC)
(
 const char          *const fstr,          /* <-- Fortran string */
 int                 *const len,           /* <-- String Length  */
 int                 *const faces_number,   /* --> faces number */
 int                 *const faces         /* --> faces  */
 CS_ARGF_SUPP_CHAINE
 );

/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgcel, CSGCEL)
(
 const char          *const fstr,          /* <-- Fortran string */
 int                 *const len,           /* <-- String Length  */
 int                 *const cells_number,   /* --> cells number */
 int                 *const cells         /* --> cells  */
 CS_ARGF_SUPP_CHAINE
 );


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* __CS_SELECTOR_H__ */
