/*============================================================================
 *
 *     This file is part of the Code_Saturne Kernel, element of the
 *     Code_Saturne CFD tool.
 *
 *     Copyright (C) 1998-2008 EDF S.A., France
 *
 *     contact: saturne-support@edf.fr
 *
 *     The Code_Saturne Kernel is free software; you can redistribute it
 *     and/or modify it under the terms of the GNU General Public License
 *     as published by the Free Software Foundation; either version 2 of
 *     the License, or (at your option) any later version.
 *
 *     The Code_Saturne Kernel is distributed in the hope that it will be
 *     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
 *     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with the Code_Saturne Kernel; if not, write to the
 *     Free Software Foundation, Inc.,
 *     51 Franklin St, Fifth Floor,
 *     Boston, MA  02110-1301  USA
 *
 *============================================================================*/

#ifndef __CS_CALCIUM_H__
#define __CS_CALCIUM_H__

/*============================================================================
 * Basic CALCIUM-mappable functions for code coupling using SALOME's YACS
 *============================================================================*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_base.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Function name mappings */

#define cp_cd    cs_calcium_connect
#define cp_fin   cs_calcium_disconnect

/* Instance continuation directive */

#define CS_CALCIUM_CONTINUE  20
#define CS_CALCIUM_STOP      22

/* Maximum string Lengths (based on CALCIUM's limits) */

#define CS_CALCIUM_INSTANCE_LEN 72
#define CS_CALCIUM_VARIABLE_LEN 144

/*=============================================================================
 * Type Definitions
 *============================================================================*/

/* CALCIUM Variable type dependency */

typedef enum {

  CALCIUM_time,         /* Physical time */
  CALCIUM_iteration     /* Iteration number */

} cs_calcium_timedep_t;

/* CALCIUM Variable type dependency */

typedef enum {

  CALCIUM_continue,     /* Use last values after disconnect */
  CALCIUM_stop          /* Stop after disconnect */

} cs_calcium_continuation_t;

/*----------------------------------------------------------------------------
 * Function pointer types
 *----------------------------------------------------------------------------*/

typedef int
(cs_calcium_connect_t) (void  *component,
                        char  *s);

typedef int
(cs_calcium_disconnect_t) (void  *component,
                           int    cont);

typedef int
(cs_calcium_read_int_t)(void    *component,
                        int      time_dep,
                        float   *min_time,
                        float   *max_time,
                        int     *iteration,
                        char    *var_name,
                        int      n_val_max,
                        int     *n_val_read,
                        int      val[]);

typedef int
(cs_calcium_read_float_t)(void    *component,
                          int      time_dep,
                          float   *min_time,
                          float   *max_time,
                          int     *iteration,
                          char    *var_name,
                          int      n_val_max,
                          int     *n_val_read,
                          float    val[]);

typedef int
(cs_calcium_read_double_t)(void    *component,
                           int      time_dep,
                           double  *min_time,
                           double  *max_time,
                           int     *iteration,
                           char    *var_name,
                           int      n_val_max,
                           int     *n_val_read,
                           double   val[]);

typedef int
(cs_calcium_write_int_t)(void    *component,
                         int      time_dep,
                         float    cur_time,
                         int      iteration,
                         char    *var_name,
                         int      n_val,
                         int      val[]);

typedef int
(cs_calcium_write_float_t)(void    *component,
                           int      time_dep,
                           float    cur_time,
                           int      iteration,
                           char    *var_name,
                           int      n_val,
                           float    val[]);

typedef int
(cs_calcium_write_double_t)(void    *component,
                            int      time_dep,
                            double   cur_time,
                            int      iteration,
                            char    *var_name,
                            int      n_val,
                            double   val[]);

/*=============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Connection
 *
 * parameters:
 *   comp_id <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   s       --> name of calling instance
 *               (CS_CALCIUM_INSTANCE_LEN chars max)
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_connect(int   comp_id,
                   char *s);

/*----------------------------------------------------------------------------
 * Disconnection
 *
 * parameters:
 *   comp_id <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   cont    --> continuation directive (continue with last values or stop)
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_disconnect(int                       comp_id,
                      cs_calcium_continuation_t cont);

/*----------------------------------------------------------------------------
 * Read values, blocking until they are available.
 *
 * parameters:
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   min_time   <-> lower bound of read interval
 *   max_time   <-- upper bound of read interval
 *   iteration  <-> iteration number of read
 *   var_name   <-- name of the variable to read
 *   n_val_max  <-- maximum number of values to read
 *   n_val_read <-- maximum number of values to read
 *   val        --> values read
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_read_int(int                    comp_id,
                    cs_calcium_timedep_t   time_dep,
                    double                *min_time,
                    double                *max_time,
                    int                   *iteration,
                    const char            *var_name,
                    int                    n_val_max,
                    int                   *n_val_read,
                    int                    val[]);

int
cs_calcium_read_float(int                    comp_id,
                      cs_calcium_timedep_t   time_dep,
                      double                *min_time,
                      double                *max_time,
                      int                   *iteration,
                      const char            *var_name,
                      int                    n_val_max,
                      int                   *n_val_read,
                      float                  val[]);

int
cs_calcium_read_double(int                    comp_id,
                       cs_calcium_timedep_t   time_dep,
                       double                *min_time,
                       double                *max_time,
                       int                   *iteration,
                       const char            *var_name,
                       int                    n_val_max,
                       int                   *n_val_read,
                       double                 val[]);

/*----------------------------------------------------------------------------
 * Write values.
 *
 * parameters:
 *   comp_id    <-- id of component to connect (0 to n-1, Code_Saturne local)
 *   time_dep   <-- type of time dependency (time or iteration)
 *   cur_time   <-- current time
 *   iteration  <-- iteration number
 *   var_name   <-- name of the variable to read
 *   n_val      <-- number of values to read
 *   val        <-- values written
 *
 * returns:
 *   0 in case of success, error code otherwise
 *----------------------------------------------------------------------------*/

int
cs_calcium_write_int(int                    comp_id,
                     cs_calcium_timedep_t   time_dep,
                     double                 cur_time,
                     int                    iteration,
                     const char            *var_name,
                     int                    n_val,
                     const int              val[]);

int
cs_calcium_write_float(int                    comp_id,
                       cs_calcium_timedep_t   time_dep,
                       double                 cur_time,
                       int                    iteration,
                       const char            *var_name,
                       int                    n_val,
                       const float            val[]);

int
cs_calcium_write_double(int                    comp_id,
                        cs_calcium_timedep_t   time_dep,
                        double                 cur_time,
                        int                    iteration,
                        const char            *var_name,
                        int                    n_val,
                        const double           val[]);

/*----------------------------------------------------------------------------
 * Assign a component and its id
 *
 * parameters:
 *   comp_id <-- id of component (0 to n-1, Code_Saturne local)
 *   comp    <-- pointer to component
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_component(int    comp_id,
                         void  *comp);

/*----------------------------------------------------------------------------
 * Set proxy IPC communication mode
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_comm_proxy(void);

/*----------------------------------------------------------------------------
 * Set the CALCIUM-mappable function's verbosity
 *
 * parameters:
 *   n_echo <-- verbosity (none if -1, headers if 0,
 *              headers + n first and last elements if > 0.
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_verbosity(int  n_echo);

/*----------------------------------------------------------------------------
 * Set connect and disconnect functions
 *
 * parameters:
 *   cp_cd_func  <-- pointer to cp_cd function or equivalent
 *   cp_fin_func <-- pointer to cp_fin function or equivalent
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_connection_funcs(cs_calcium_connect_t     *cp_cd_func,
                                cs_calcium_disconnect_t  *cp_fin_func);

/*----------------------------------------------------------------------------
 * Set read, non-blocking read, and write functions for integers
 *
 * parameters:
 *   cp_len_func  <-- pointer to cp_len function or equivalent
 *   cp_een_func  <-- pointer to cp_een function or equivalent
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_int_rw_funcs(cs_calcium_read_int_t     *cp_len_func,
                            cs_calcium_write_int_t    *cp_een_func);

/*----------------------------------------------------------------------------
 * Set read, non-blocking read, and write functions for floats
 *
 * parameters:
 *   cp_lre_func  <-- pointer to cp_lre function or equivalent
 *   cp_ere_func  <-- pointer to cp_ere function or equivalent
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_float_rw_funcs(cs_calcium_read_float_t     *cp_lre_func,
                              cs_calcium_write_float_t    *cp_ere_func);

/*----------------------------------------------------------------------------
 * Set read, non-blocking read, and write functions for doubles
 *
 * parameters:
 *   cp_ldb_func  <-- pointer to cp_ldb function or equivalent
 *   cp_edb_func  <-- pointer to cp_edb function or equivalent
 *----------------------------------------------------------------------------*/

void
cs_calcium_set_double_rw_funcs(cs_calcium_read_double_t     *cp_ldb_func,
                               cs_calcium_write_double_t    *cp_edb_func);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_CALCIUM_H__ */

