#ifndef __CS_PARAMEDMEM_PRIV_HXX__
#define __CS_PARAMEDMEM_PRIV_HXX__

/*============================================================================
 * MEDCoupling ParaMESH/ParaFIELD wrapper functions.
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2026 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "base/cs_medcoupling_mesh.hxx"
#include "base/cs_paramedmem_coupling.h"

#if defined(HAVE_PARAMEDMEM)

#include <MEDCouplingField.hxx>
#include <MEDCouplingFieldDouble.hxx>

#include <InterpKernelDEC.hxx>
#include <ParaFIELD.hxx>
#include <ParaMESH.hxx>

using namespace MEDCoupling;
#endif

#define USE_PARAFIELD 1
/* 1 to use <ParaFIELD> fields,                                            \
 0 to directly use          MEDCouplingFieldDouble */

/*----------------------------------------------------------------------------*/

/*============================================================================
 * Structure definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * MEDCoupling writer/reader structure
 *----------------------------------------------------------------------------*/

struct _cs_paramedmem_coupling_t {
  /* Coupling Name */
  std::string _name; /* Coupling name */

  /* Current app name */
  ple_coupling_mpi_set_info_t apps[2];

#if defined(HAVE_PARAMEDMEM)

  cs_medcoupling_mesh_t *mesh;

  ParaMESH *para_mesh; /* Associated ParaMESH structure. */

  InterpKernelDEC *dec; /* Data Exchange Channel */

#if USE_PARAFIELD == 1
  std::vector<ParaFIELD *> fields;
#else
  std::vector<MEDCouplingFieldDouble *> fields;
#endif

#else

  void *para_mesh;
  void *dec;
  void *fields;

#endif

#if defined(HAVE_MPI)
  MPI_Comm comm;
#else
  int comm;
#endif

  bool dec_synced;

#ifdef __cplusplus

private:
  void
  _error_without_paramedmem() const
  {
#if !defined(HAVE_PARAMEDMEM)

    bft_error(__FILE__,
              __LINE__,
              0,
              _("Error: %s cannot be called without "
                "MEDCoupling MPI support."),
              __func__);
#endif
  };

#if defined(HAVE_PARAMEDMEM)

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Get a field by its name
   *
   * \param[in] name          name of field
   *
   * \return medcoupling field
   */
  /*----------------------------------------------------------------------------*/

  MEDCouplingFieldDouble *
  _get_field(const char *name) const;

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Generate mesh structure from user's defintion
   *
   * \param[in] select_criteria   selection criteria (string)
   * \param[in] elt_dim           mesh dimension (2 or 3)
   */
  /*----------------------------------------------------------------------------*/

  void
  _generate_coupling_mesh(const char *select_criteria, int elt_dim);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Generate mesh structure from user's defintion
   *
   * \param[in] n_elts   local number of elements
   * \param[in] elt_ids  list of local elements
   * \param[in] elt_dim  dimension of elements (2: faces, 3: cells)
   */
  /*----------------------------------------------------------------------------*/

  void
  _generate_coupling_mesh_from_ids(cs_lnum_t       n_elts,
                                   const cs_lnum_t elt_ids[],
                                   int             elt_dim);

#endif /* HAVE_PARAMEDMEM */

public:
  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Define coupled mesh based on a selection criteria
   *
   * \param[in] sel_crit  geometrical selection criteria (string)
   * \param[in] elt_dim   dimension of coupled elements
   */
  /*----------------------------------------------------------------------------*/

  void
  add_mesh_from_criteria(const char *sel_crit, int elt_dim);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Define coupled mesh based on a cs_zone_t pointer
   *
   * \param[in] zone  pointer to cs_zone_t struct
   */
  /*----------------------------------------------------------------------------*/

  void
  add_mesh_from_zone(const cs_zone_t *zone);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Define coupled mesh based on a cs_zone_t pointer
   *
   * \param[in] n_elts   local number of elements
   * \param[in] elt_ids  list of local elements
   * \param[in] elt_dim  dimension of elements (2: faces, 3: cells)
   */
  /*----------------------------------------------------------------------------*/

  void
  add_mesh_from_ids(cs_lnum_t n_elts, const cs_lnum_t elt_ids[], int elt_dim);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Get name of the coupling
   *
   * \return name of the coupling
   */
  /*----------------------------------------------------------------------------*/

  std::string
  getName() const
  {
    return _name;
  }

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Get number of elements of coupled mesh
   *
   * \return number of elements in mesh associated to coupling
   */
  /*----------------------------------------------------------------------------*/

  cs_lnum_t
  get_n_elts() const
  {
    cs_lnum_t retval = 0;

#if !defined(HAVE_PARAMEDMEM)

    this->_error_without_paramedmem();

#else

    retval = this->mesh->n_elts;

#endif

    return retval;
  }

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Get indirection list for elements in coupled mesh
   *
   * \return cs_lnum_t pointer to indirection list
   */
  /*----------------------------------------------------------------------------*/

  const cs_lnum_t *
  get_elt_list() const
  {
    const cs_lnum_t *retval = nullptr;

#if !defined(HAVE_PARAMEDMEM)

    this->_error_without_paramedmem();

#else

    retval = this->mesh->elt_list;

#endif

    return retval;
  }

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Get number of vertices of coupled mesh
   *
   * \return number of elements in mesh associated to coupling
   */
  /*----------------------------------------------------------------------------*/

  cs_lnum_t
  get_n_vertices() const
  {
    cs_lnum_t retval = 0;

#if !defined(HAVE_PARAMEDMEM)

    this->_error_without_paramedmem();

#else

    retval = cs_medcoupling_mesh_get_n_vertices(this->mesh);

#endif

    return retval;
  }

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Get indirection list for vertices in coupled mesh.
   *
   * \return pointer to indirection list; null if locally contiguous or empty
   */
  /*----------------------------------------------------------------------------*/

  const cs_lnum_t *
  get_vertex_list() const
  {
    const cs_lnum_t *retval = nullptr;

#if !defined(HAVE_PARAMEDMEM)

    this->_error_without_paramedmem();

#else

    retval = this->mesh->vtx_list;

#endif

    return retval;
  }

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Add a coupled field
   *
   * \param[in] name          name of field
   * \param[in] dim           field dimension
   * \param[in] field_nature  field nature flag
   * \param[in] space_discr   field space discretisation (nodes or cells)
   * \param[in] time_discr    field coupling time discretisation
   *
   * \return index of field within the storing vector
   */
  /*----------------------------------------------------------------------------*/

  int
  add_field(const char              *name,
            int                      dim,
            cs_medcpl_field_nature_t field_nature,
            cs_medcpl_space_discr_t  space_discr,
            cs_medcpl_time_discr_t   time_discr);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Add a coupled field based on a cs_field_t pointer
   *
   * \param[in] f           pointer to cs_field_t struct
   * \param[in] fn          field nature flag
   * \param[in] time_discr  field coupling time discretisation
   *
   * \return index of field within the storing vector
   */
  /*----------------------------------------------------------------------------*/

  int
  add_field(const cs_field_t        *f,
            cs_medcpl_field_nature_t fn,
            cs_medcpl_time_discr_t   td);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Assign values based on parent mesh location to associated
   *        ParaFIELD objects.
   *
   * \param[in]  name    name of field
   * \param[in]  values  array of values to write
   *                     (defined on parent mesh location)
   */
  /*----------------------------------------------------------------------------*/

  void
  set_values(const char *name, const double values[]);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Assign values based on mesh location corresponding to coupled
   *        elements (and associated ParaMESH) to associated ParaFIELD objects.
   *
   * If the whole mesh is coupled, the behavior is the sames as that of
   * \ref set_values.
   *
   * \param[in]  name    name of field
   * \param[in]  values  array of values to write
   *                     (defined on selected mesh subset)
   */
  /*----------------------------------------------------------------------------*/

  void
  set_values_l(const char *name, const double values[]);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Copy values from associated ParaFIELD object to array defined
   *        parent mesh location.
   *
   * \param[in]  c       pointer to cs_paramedmem_coupling_t structure
   * \param[in]  name    name of field
   * \param[in]  values  array in which values will be stored
   */
  /*----------------------------------------------------------------------------*/

  void
  get_values(const char *name, double values[]) const;

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Copy values from associated ParaFIELD structure to array defined
   *        on mesh location corresponding to coupled elements
   *        (and associated ParaMESH).
   *
   * If the whole mesh is coupled, the behavior is the sames as that of
   * \ref get_values.
   *
   * \param[in]  name    name of field
   * \param[in]  values  array in which values will be stored
   */
  /*----------------------------------------------------------------------------*/

  void
  get_values_l(const char *name, double values[]) const;

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Send values of field attached to DEC
   *
   */
  /*----------------------------------------------------------------------------*/

  void
  send_data() const;

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Send values of a field. If vals pointer is non-null,
   * values are updated before send
   *
   * \param[in] name  name of field
   * \param[in] vals  array of values to write
   */
  /*----------------------------------------------------------------------------*/

  void
  send_data(const char *name, const double *vals);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Send values of a field. If vals pointer is non-null,
   * values are updated before send
   *
   * \param[in] name  name of field
   * \param[in] vals  array of values to write
   */
  /*----------------------------------------------------------------------------*/

  void
  send_data_l(const char *name, const double *vals);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Recieve values of field attached to DEC
   *
   */
  /*----------------------------------------------------------------------------*/

  void
  recv_data();

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Recieve values of a field.
   *
   * \param[in] name  name of field
   * \param[in] vals  array of values to read
   */
  /*----------------------------------------------------------------------------*/

  void
  recv_data(const char *name, double *vals);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Receive values of a field.
   *
   * \param[in] name  name of field
   * \param[in] vals  array of values to read
   */
  /*----------------------------------------------------------------------------*/

  void
  recv_data_l(const char *name, double *vals);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Sync the coupling's the DEC
   *
   */
  /*----------------------------------------------------------------------------*/

  void
  sync_dec();

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Attach a field to the DEC for send operation using its index
   *
   * \param[in] field_id  index of field in storing vector
   */
  /*----------------------------------------------------------------------------*/

  void
  attach_field_by_id(int field_id);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Attach a field to the DEC for send operation using its name
   *
   * \param[in] name  name of field (string)
   */
  /*----------------------------------------------------------------------------*/

  void
  attach_field_by_name(const char *name);

  /*----------------------------------------------------------------------------*/
  /*!
   * \brief Log ParaMEDMEM coupling setup information
   *
   */
  /*----------------------------------------------------------------------------*/

  void
  log() const;

#endif // cplusplus
};

#endif /* __CS_PARAMEDMEM_PRIV_HXX__ */
