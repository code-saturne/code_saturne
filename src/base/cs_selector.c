#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include <bft_mem_usage.h>
#include <bft_mem.h>
#include <bft_error.h>
#include <bft_printf.h>
#include <bft_sys_info.h>
#include <bft_timer.h>

/* pour les API Fortran : a deplacer */
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"
#include "cs_selector.h"

#include "fvm_selector.h"

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */

/* API Fortran : a deplacer */

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
)
{
  char * cpyfstr;
  int i;
  int lenwithoutblank;

  /* Initialization */
  *faces_number = 0;

  /* Copy fstr without last blanks  */
  BFT_MALLOC(cpyfstr, *len + 1, char);
  lenwithoutblank = *len - 1;

  while(fstr[lenwithoutblank--] == ' ' &&  lenwithoutblank >= 0);

  if (lenwithoutblank < -1) {
    BFT_FREE(cpyfstr);
    return;
  }
  else
    lenwithoutblank += 2;

  for(i = 0 ; i < lenwithoutblank; i++)
    cpyfstr[i] = fstr[i];
  cpyfstr[lenwithoutblank] = '\0';

  BFT_REALLOC(cpyfstr, lenwithoutblank + 1, char);

  /* Get faces with C string */

  fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                        cpyfstr,
                        faces_number,
                        faces);
  BFT_FREE(cpyfstr);
}

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
)
{
  char * cpyfstr;
  int i;
  int lenwithoutblank;

  /* Initialization */
  *faces_number = 0;

  /* Copy fstr without last blanks  */
  BFT_MALLOC(cpyfstr, *len + 1, char);
  lenwithoutblank = *len - 1;

  while(fstr[lenwithoutblank--] == ' ' &&  lenwithoutblank >= 0);

  if (lenwithoutblank < -1) {
    BFT_FREE(cpyfstr);
    return;
  }
  else
    lenwithoutblank += 2;

  for(i = 0 ; i < lenwithoutblank; i++)
    cpyfstr[i] = fstr[i];
  cpyfstr[lenwithoutblank] = '\0';

  BFT_REALLOC(cpyfstr, lenwithoutblank + 1, char);

  /* Get faces with C string */

  fvm_selector_get_list(cs_glob_mesh->select_i_faces,
                        cpyfstr,
                        faces_number,
                        faces);
  BFT_FREE(cpyfstr);
}

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
)
{
  char * cpyfstr;
  int i;
  int lenwithoutblank;

  /* Initialization */
  *cells_number = 0;

  /* Copy fstr without last blanks  */
  BFT_MALLOC(cpyfstr, *len + 1, char);
  lenwithoutblank = *len - 1;

  while(cpyfstr[lenwithoutblank--] == ' ' &&  lenwithoutblank >= 0);

  if (lenwithoutblank < -1) {
    BFT_FREE(cpyfstr);
    return;
  }
  else
    lenwithoutblank += 2;

  for(i = 0 ; i < lenwithoutblank; i++)
    cpyfstr[i] = fstr[i];
  cpyfstr[lenwithoutblank] = '\0';

  BFT_REALLOC(cpyfstr, lenwithoutblank + 1, char);

  /* Get cells with C string */
  fvm_selector_get_list(cs_glob_mesh->select_cells,
                        cpyfstr,
                        cells_number,
                        cells);

  BFT_FREE(cpyfstr);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
