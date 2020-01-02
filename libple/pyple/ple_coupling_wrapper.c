/*============================================================================
 * Wrapper file used to generate the Python bindings for ple_coupling functions
 *============================================================================*/

/*
  This file is part of the "Parallel Location and Exchange" library,
  intended to provide mesh or particle-based code coupling services.

  Copyright (C) 2005-2020  EDF S.A.

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with this library; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include <Python.h>
#include "mpi4py/mpi4py.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "ple_defs.h"
#include "ple_config_defs.h"
#include "ple_coupling.h"

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#if 0
} /* Fake brace to force Emacs auto-indentation back to column 0 */
#endif
#endif /* __cplusplus */


/*============================================================================
 * Static local variables
 *============================================================================*/

static int                      _n_sets = 0;
static ple_coupling_mpi_set_t **_mpi_sets;

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Python bindings for the C functions which will be called by Python
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Build a group id within a communicator based on its name
 *
 * \param  [in]  Comm       mpi4py communicator.
 * \param  [in]  group_name name associated with current group.
 *
 * \return id associated with local name.
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_name_to_id(PyObject *self, PyObject *args)
{
  PyObject *py_comm = NULL;
  MPI_Comm *comm_p  = NULL;

  const char *group_name = "";
  int         group_name_size = 0;

  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "Os#:ple_coupling_mpi_name_to_id",
                        &py_comm, &group_name, &group_name_size))
    return NULL;

  /* Retrieve the C MPI_Comm from the mpi4py.Comm */
  comm_p = PyMPIComm_Get(py_comm);
  if (comm_p == NULL)
    return NULL;

  /* Call the C function */
  int new_id = ple_coupling_mpi_name_to_id(*comm_p, group_name);

  /* Return a Python value */
  return Py_BuildValue("i", new_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Discover other applications in a set with a common communicator.
 *
 * \param[in] sync_flag 1 if application is to be synchronized at each
 *                      time step, 0 if independent from others.
 * \param[in] app_type  name of current application type (software name).
 * \param[in] app_name  name of current application (data/case name).
 * \param[in] base_comm communicator associated with all applications.
 * \param[in] app_comm  communicator associated with local application.
 *
 * \return PLE coupling MPI set info id within the static local array
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_create(PyObject *self, PyObject *args)
{
  int sync_flags = 0;

  const char *app_type = "";
  const char *app_name = "";

  PyObject *py_base_comm = NULL;
  PyObject *py_app_comm  = NULL;

  MPI_Comm *base_comm_p = NULL;
  MPI_Comm *app_comm_p  = NULL;

  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "issOO:ple_coupling_mpi_set_create",
                        &sync_flags,
                        &app_type,
                        &app_name,
                        &py_base_comm,
                        &py_app_comm))
    return NULL;

  /* Retrieve the C MPI_Comm from the mpi4py.Comm */
  base_comm_p = PyMPIComm_Get(py_base_comm);
  app_comm_p  = PyMPIComm_Get(py_app_comm);

  if (base_comm_p == NULL || app_comm_p == NULL)
    return NULL;

  /* Call the C function */
  if (_n_sets == 0)
    PLE_MALLOC(_mpi_sets, 1, ple_coupling_mpi_set_t *);
  else
    PLE_REALLOC(_mpi_sets, _n_sets + 1, ple_coupling_mpi_set_t *);

  _mpi_sets[_n_sets] = ple_coupling_mpi_set_create(sync_flags,
                                                   app_type,
                                                   app_name,
                                                   *base_comm_p,
                                                   *app_comm_p);

  _n_sets++;

  int ret_index = _n_sets - 1;

  return Py_BuildValue("i", ret_index);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free an PLE coupling MPI set info structure.
 *
 * \param[in] index id of the set to free
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_destroy(PyObject *self, PyObject *args)
{
  int set_index = 0;
  if (!PyArg_ParseTuple(args, "i:ple_coupling_mpi_set_destroy", &set_index))
    return NULL;

  ple_coupling_mpi_set_destroy(&_mpi_sets[set_index]);

  Py_INCREF(Py_None);
  return Py_None;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return application information in set's common communicator.
 *
 * \param[in] index   index of the PLE coupling MPI set info structure.
 * \param[in] app_id  application id
 *
 * \return application information structure as a dictionnary.
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_get_info(PyObject *self, PyObject *args)
{
  int set_index = 0;
  int app_id    = 0;

  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "ii:ple_coupling_mpi_set_get_info",
                        &set_index, &app_id))
    return NULL;

  /* Call the C function */
  ple_coupling_mpi_set_info_t ai =
    ple_coupling_mpi_set_get_info(_mpi_sets[set_index], app_id);

  /* Return a Python dictionary */
  return Py_BuildValue("{sisisissss}",
                       "status", ai.status,
                       "root_rank", ai.root_rank,
                       "n_ranks", ai.n_ranks,
                       "app_type", ai.app_type,
                       "app_name", ai.app_name);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the number of applications in a coupled set.
 *
 * \param[in] index   index of the PLE coupling MPI set info structure.
 *
 * \return number of application in set's common communicator.
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_n_apps(PyObject *self, PyObject *args)
{
  int set_index = 0;

  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "i:ple_coupling_mpi_set_n_apps", &set_index))
    return NULL;

  /* Call the C function */
  int n_apps = ple_coupling_mpi_set_n_apps(_mpi_sets[set_index]);

  /* Return a Python value */
  return Py_BuildValue("i", n_apps);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Return the id of the local application in a coupled set.
 *
 * \param[in] index   index of the PLE coupling MPI set info structure.
 *
 * \return id of the local application in set's common communicator.
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_get_app_id(PyObject *self, PyObject *args)
{
  int set_index = 0;

  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "i:ple_coupling_mpi_set_get_app_id", &set_index))
    return NULL;

  /* Call the C function */
  int app_id = ple_coupling_mpi_set_get_app_id(_mpi_sets[set_index]);

  /* Return a Python value */
  return Py_BuildValue("i", app_id);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Synchronize applications in a set.
 *
 * \param[in] index       index of the PLE coupling MPI set info structure.
 * \param[in] sync_flags  app sync flags
 * \param[in] time_step   app time step
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_synchronize(PyObject *self, PyObject *args)
{
  int    set_index = 0;
  int    sync_flag = 0;
  double time_step = 0.;

  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "iid:ple_coupling_mpi_set_synchronize",
                        &set_index, &sync_flag, &time_step))
    return NULL;

  /* Call the C function */
  ple_coupling_mpi_set_synchronize(_mpi_sets[set_index], sync_flag, time_step);

  /* Return a Python value */
  Py_INCREF(Py_None);
  return Py_None;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get status of applications in a set.
 *
 * \param[in] index       index of the PLE coupling MPI set info structure.
 *
 * \return a python list of status flags
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_get_status(PyObject *self, PyObject *args)
{
  int set_index = 0;
  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "i:ple_coupling_mpi_set_get_status", &set_index))
    return NULL;

  /* Call the C function */
  const int *set_status = ple_coupling_mpi_set_get_status(_mpi_sets[set_index]);
  const int n_apps = ple_coupling_mpi_set_n_apps(_mpi_sets[set_index]);

  PyObject *status_list = PyList_New(n_apps);
  for (int i = 0; i < n_apps; i++)
    PyList_SetItem(status_list, i, set_status[i]);

  /* Return a Python list */
  return status_list;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Get time steps in a set.
 *
 * \param[in] index       index of the PLE coupling MPI set info structure.
 *
 * \return a python list of the set's time steps
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_get_timestep(PyObject *self, PyObject *args)
{
  int set_index = 0;
  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "i:ple_coupling_mpi_set_get_timestep", &set_index))
    return NULL;

  /* Call the C function */
  const double *set_dt = ple_coupling_mpi_set_get_timestep(_mpi_sets[set_index]);
  const int n_apps = ple_coupling_mpi_set_n_apps(_mpi_sets[set_index]);

  PyObject *dt_list = PyList_New(n_apps);
  for (int i = 0; i < n_apps; i++)
    PyList_SetItem(dt_list, i, PyFloat_FromDouble(set_dt[i]));

  /* Return a Python list */
  return dt_list;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Dump printout of an PLE coupling MPI set info structure.
 *
 * \param[in] index       index of the PLE coupling MPI set info structure.
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_dump(PyObject *self, PyObject *args)
{

  int set_index = 0;
  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "i:ple_coupling_mpi_set_dump", &set_index))
    return NULL;

  /* Call the C function */
  ple_coupling_mpi_set_dump(_mpi_sets[set_index]);

  /* Return a Python value */
  Py_INCREF(Py_None);
  return Py_None;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Create an intracommunicator from a local and distant communicator
 *        within a base communicator.
 *
 * \param[in] base_comm     communicator associated with both applications
 * \param[in] app_comm      communicator associated with local application
 * \param[in] distant_root  rank of distant group leader in base_comm
 * \param[in] new_comm      pointer to new communicator
 * \param[in] local_range   first and past-the last ranks of local application
 *                          in new communicator
 * \param[in] distant_range first and past-the last ranks of distant
 *                          application in new communicator
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_intracomm_create(PyObject *self, PyObject *args)
{

  PyObject *py_base_comm = NULL;
  PyObject *py_app_comm  = NULL;
  PyObject *py_new_comm  = NULL;

  int distant_root = 0;

  PyObject *llist = NULL;
  PyObject *dlist = NULL;

  /* We verify that the correct arguments are provided.
   * format is of the form "var1_typevar2_type..varn_type:function_name"
   * function_name is used so that in case of a traceback error Python will
   * know it failed here.
   * If a CPython function returns NULL, it means it failed.
   */
  if (!PyArg_ParseTuple(args, "OOiOO!O!:ple_coupling_mpi_intracomm_create",
                        &py_base_comm, &py_app_comm,
                        &distant_root, &py_new_comm,
                        &PyList_Type, &llist,
                        &PyList_Type, &dlist))
    return NULL;

  MPI_Comm *base_comm_p = NULL;
  MPI_Comm *app_comm_p  = NULL;
  MPI_Comm *new_comm_p  = NULL;

  base_comm_p = PyMPIComm_Get(py_base_comm);
  app_comm_p  = PyMPIComm_Get(py_app_comm);
  new_comm_p  = PyMPIComm_Get(py_new_comm);

  int lranks[2] = {-1, -1};
  int dranks[2] = {-1, -1};

  /* Call the C function */
  ple_coupling_mpi_intracomm_create(*base_comm_p,
                                    *app_comm_p,
                                    distant_root,
                                    new_comm_p,
                                    lranks,
                                    dranks);

  for (int i = 0; i < 2; i++) {
    PyList_SetItem(llist, i, Py_BuildValue("i",lranks[i]));
    PyList_SetItem(dlist, i, Py_BuildValue("i",dranks[i]));
  }

  /* Return a Python value */
  Py_INCREF(Py_None);
  return Py_None;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Make the PLE_COUPLING_* flags available in python
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_get_masks(PyObject *self, PyObject *args)
{
  /* Return a Python dictionary */
  return Py_BuildValue("{sisisisisisisisisisisisisisisi}",
                       "INIT", PLE_COUPLING_INIT,
                       "NO_SYNC", PLE_COUPLING_NO_SYNC,
                       "STOP", PLE_COUPLING_STOP,
                       "LAST", PLE_COUPLING_LAST,
                       "NEW_ITERATION", PLE_COUPLING_NEW_ITERATION,
                       "REDO_ITERATION", PLE_COUPLING_REDO_ITERATION,
                       "TS_MIN", PLE_COUPLING_TS_MIN,
                       "TS_LEADER", PLE_COUPLING_TS_LEADER,
                       "UNSTEADY", PLE_COUPLING_UNSTEADY,
                       "STEADY", PLE_COUPLING_STEADY,
                       "CONVERGED", PLE_COUPLING_CONVERGED,
                       "USER_1", PLE_COUPLING_USER_1,
                       "USER_2", PLE_COUPLING_USER_2,
                       "USER_3", PLE_COUPLING_USER_3,
                       "USER_4", PLE_COUPLING_USER_4);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Free all the ple_coupling_mpi_set_t pointers created.
 */
/*----------------------------------------------------------------------------*/

static PyObject *
pyple_coupling_mpi_set_destroy_all(PyObject *self, PyObject *args)
{

  for (int i = 0; i < _n_sets; i++)
    ple_coupling_mpi_set_destroy(&_mpi_sets[i]);

  PLE_FREE(_mpi_sets);

  /* Return a Python value */
  Py_INCREF(Py_None);
  return Py_None;

}

/*============================================================================
 * Define the list of methods available for Python.
 * Two first arguments are the function's Python name and a pointer to
 * the CPython function which will be called.
 *============================================================================*/

static struct PyMethodDef
pyple_coupling_methods[] = {
  {"coupling_mpi_name_to_id",
   (PyCFunction)pyple_coupling_mpi_name_to_id,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_create",
   (PyCFunction)pyple_coupling_mpi_set_create,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_destroy",
   (PyCFunction)pyple_coupling_mpi_set_destroy,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_get_info",
   (PyCFunction)pyple_coupling_mpi_set_get_info,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_n_apps",
   (PyCFunction)pyple_coupling_mpi_set_n_apps,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_get_app_id",
   (PyCFunction)pyple_coupling_mpi_set_get_app_id,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_synchronize",
   (PyCFunction)pyple_coupling_mpi_set_synchronize,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_get_status",
   (PyCFunction)pyple_coupling_mpi_set_get_status,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_get_timestep",
   (PyCFunction)pyple_coupling_mpi_set_get_timestep,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_dump",
   (PyCFunction)pyple_coupling_mpi_set_dump,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_intracomm_create",
   (PyCFunction)pyple_coupling_mpi_intracomm_create,
   METH_VARARGS,
   NULL},
  {"coupling_get_masks",
   (PyCFunction)pyple_coupling_get_masks,
   METH_VARARGS,
   NULL},
  {"coupling_mpi_set_destroy_all",
   (PyCFunction)pyple_coupling_mpi_set_destroy_all,
   METH_VARARGS,
   NULL},
  {NULL, NULL, 0, NULL}
};

/*============================================================================
 * Define the python module specific functions needed for it to be imported
 *============================================================================*/

#if PY_MAJOR_VERSION < 3
/* --- Python 2 --- */

PyMODINIT_FUNC initlibpyplecoupling(void)
{
  PyObject *m = NULL;

  /* Initialize mpi4py C-API */
  if (import_mpi4py() < 0) goto bad;

  /* Module initialization  */
  m = Py_InitModule("libpyplecoupling", pyple_coupling_methods);
  if (m == NULL) goto bad;

  return;

 bad:
  return;
}

#else
/* --- Python 3 --- */

static struct PyModuleDef libpyplecoupling_module = {
  PyModuleDef_HEAD_INIT,
  "libpyplecoupling",       /* m_name */
  NULL,                     /* m_doc */
  -1,                       /* m_size */
  pyple_coupling_methods    /* m_methods */,
  NULL,                     /* m_reload */
  NULL,                     /* m_traverse */
  NULL,                     /* m_clear */
  NULL                      /* m_free */
};

PyMODINIT_FUNC
PyInit_libpyplecoupling(void)
{
  PyObject *m = NULL;

  /* Initialize mpi4py's C-API */
  if (import_mpi4py() < 0) goto bad;

  /* Module initialization  */
  m = PyModule_Create(&libpyplecoupling_module);
  if (m == NULL) goto bad;

  return m;

 bad:
  return NULL;
}

#endif

/*----------------------------------------------------------------------------*/

#ifdef __cplusplus
}
#endif /* __cplusplus */
