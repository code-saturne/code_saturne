#ifndef CS_PHYSICAL_UNITS_H
#define CS_PHYSICAL_UNITS_H

/*============================================================================
 * Physical units management and catalogue
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stddef.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_defs.h"

/*============================================================================
 * Public function prototypes
 *============================================================================*/

#if defined(__cplusplus)
extern "C" {
#endif

/*----------------------------------------------------------------------------
 * Define or register a unit symbol with an optional description.
 *
 * If the unit symbol already exists, its description is updated if non-null.
 * If symbol is NULL or empty, NULL is returned.
 *
 * parameters:
 *   symbol      <-- unit symbol (e.g. "m/s", "Pa", "kg/m3")
 *   description <-- optional unit description (e.g. "meter per second")
 *
 * returns:
 *   pointer to persistent interned unit symbol string, or NULL
 *----------------------------------------------------------------------------*/

const char *
cs_physical_unit_define(const char  *symbol,
                        const char  *description);

/*----------------------------------------------------------------------------
 * Return the description associated with a unit symbol.
 *
 * parameters:
 *   symbol <-- unit symbol
 *
 * returns:
 *   pointer to unit description, or NULL if none registered
 *----------------------------------------------------------------------------*/

const char *
cs_physical_unit_get_description(const char  *symbol);

/*----------------------------------------------------------------------------
 * Return the default physical unit symbol matching a canonical name.
 *
 * The canonical name can be a standard field name (e.g. "velocity",
 * "pressure", "density", "temperature", "k", "epsilon", ...) or a fluid
 * property / reference constant name (e.g. "ro0", "viscl0", "cp0", "p0",
 * "t0", "lambda0", ...).
 *
 * parameters:
 *   name <-- canonical quantity name
 *
 * returns:
 *   pointer to persistent interned unit symbol string, or NULL if unknown
 *----------------------------------------------------------------------------*/

const char *
cs_physical_unit_get_default(const char  *name);

/*----------------------------------------------------------------------------
 * Initialize physical units map and populate with standard catalogue.
 *
 * Called automatically upon first access if not already initialized.
 *----------------------------------------------------------------------------*/

void
cs_physical_units_init(void);

/*----------------------------------------------------------------------------
 * Finalize physical units dictionary and free allocated resources.
 *----------------------------------------------------------------------------*/

void
cs_physical_units_finalize(void);

#if defined(__cplusplus)
}
#endif

#endif /* CS_PHYSICAL_UNITS_H */
