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
 const char          *const fstr,     /* <-- Fortran string */
 int                 *const len,      /* <-- String Length  */
 int                 *const n_faces,  /* --> number of faces */
 int                 *const faces     /* --> faces */
 CS_ARGF_SUPP_CHAINE
)
{
  char * cpyfstr;
  int i, c_id;
  int lenwithoutblank;

  /* Initialization */
  *n_faces = 0;

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

  c_id = fvm_selector_get_list(cs_glob_mesh->select_b_faces,
                               cpyfstr,
                               n_faces,
                               faces);

  if (fvm_selector_n_missing(cs_glob_mesh->select_b_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_b_faces, c_id, 0);
      bft_error(__FILE__, __LINE__, 0,
                _("Le groupe ou attribut \"%s\" figurant dans le\n"
                  "critère de sélection:\n"
                  "\"%s\"\ne correspond à aucune face de bord."),
                missing, cpyfstr);
  }

  BFT_FREE(cpyfstr);
}

/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgfac, CSGFAC)
(
 const char          *const fstr,      /* <-- Fortran string */
 int                 *const len,       /* <-- String Length  */
 int                 *const n_faces,   /* --> number of faces */
 int                 *const faces      /* --> faces  */
 CS_ARGF_SUPP_CHAINE
)
{
  char * cpyfstr;
  int i, c_id;
  int lenwithoutblank;

  /* Initialization */
  *n_faces = 0;

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

  c_id = fvm_selector_get_list(cs_glob_mesh->select_i_faces,
                               cpyfstr,
                               n_faces,
                               faces);

  if (fvm_selector_n_missing(cs_glob_mesh->select_i_faces, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_i_faces, c_id, 0);
      bft_error(__FILE__, __LINE__, 0,
                _("Le groupe ou attribut \"%s\" figurant dans le\n"
                  "critère de sélection:\n"
                  "\"%s\"\ne correspond à aucune face de interne."),
                missing, cpyfstr);
  }

  BFT_FREE(cpyfstr);
}

/*----------------------------------------------------------------------------
 * Give the list of the faces checking the constraint defined by
 * the Fortran string "fstr"
 *----------------------------------------------------------------------------*/

void CS_PROCF(csgcel, CSGCEL)
(
 const char          *const fstr,     /* <-- Fortran string */
 int                 *const len,      /* <-- String Length  */
 int                 *const n_cells,  /* --> number of cells */
 int                 *const cells     /* --> cells  */
 CS_ARGF_SUPP_CHAINE
)
{
  char * cpyfstr;
  int i, c_id;
  int lenwithoutblank;

  /* Initialization */
  *n_cells = 0;

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
  c_id = fvm_selector_get_list(cs_glob_mesh->select_cells,
                               cpyfstr,
                               n_cells,
                               cells);

  if (fvm_selector_n_missing(cs_glob_mesh->select_cells, c_id) > 0) {
    const char *missing
      = fvm_selector_get_missing(cs_glob_mesh->select_cells, c_id, 0);
      bft_error(__FILE__, __LINE__, 0,
                _("Le groupe ou attribut \"%s\" figurant dans le\n"
                  "critère de sélection:\n"
                  "\"%s\"\ne correspond à aucune cellule."),
                missing, cpyfstr);
  }

  BFT_FREE(cpyfstr);
}


#ifdef __cplusplus
}
#endif /* __cplusplus */
