/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to Catalyst objects
 *============================================================================*/

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2024 EDF S.A.

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

/*----------------------------------------------------------------------------*/

/* On glibc-based systems, define _GNU_SOURCE so as to enable
   modification of floating-point error exceptions handling;
   _GNU_SOURCE must be defined before including any headers, to ensure
   the correct feature macros are defined first. */

#if defined(__linux__) || defined(__linux) || defined(linux)
#define CS_FPE_TRAP
#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif
#endif

#include "cs_defs.h"

#if defined(HAVE_CATALYST)

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Catalyst and VTK library headers
 *----------------------------------------------------------------------------*/

#include <vtkCellType.h>

#include <vtkCellData.h>
#include <vtkCharArray.h>
#include <vtkFloatArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkIdList.h>
#include <vtkIdTypeArray.h>
#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkMPI.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkUnstructuredGrid.h>
#if defined(HAVE_VTKPVVERSION_H)
#include <vtkPVVersion.h>  // Preferred from version 5.10 on.
#else
#include <vtkPVConfig.h>
#endif

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#if defined(HAVE_VTKCPPYTHONSCRIPTV2PIPELINE_H)
#include <vtkCPPythonScriptV2Pipeline.h>
#endif
#include <vtkSmartPointer.h>

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkDoubleArray.h>

#include <vtkFileOutputWindow.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"
#include "bft_printf.h"

#include "fvm_defs.h"
#include "fvm_convert_array.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_writer_priv.h"

#include "cs_block_dist.h"
#include "cs_file.h"
#include "cs_parall.h"
#include "cs_part_to_block.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_catalyst.h"

/*=============================================================================
 * Local Macro Definitions
 *============================================================================*/

/* Macro for version test (tests are done on ParaView and not VTK version,
   as ParaView usually includes its own VTK, and ParaView 5.0 seems to
   indicate VTK 7.1 just like 5.1, but that version did not contain
   SetTypedTuple or have SetTupleValue deprecated for vtkTypedDataArray) */

#define CS_PV_VERSION  (PARAVIEW_VERSION_MAJOR*10 + PARAVIEW_VERSION_MINOR)

/*----------------------------------------------------------------------------
 * Catalyst field structure
 *----------------------------------------------------------------------------*/

typedef struct {

  int                   mesh_id;       /* Associated mesh structure id */
  int                   dim;           /* Field dimension */
  vtkUnstructuredGrid  *f;             /* Pointer to VTK writer fields */

} fvm_catalyst_field_t;

/*----------------------------------------------------------------------------
 * Catalyst writer/reader structure
 *----------------------------------------------------------------------------*/

typedef struct {

  char                       *name;            /* Writer name */

  int                         rank;            /* Rank of current process
                                                  in communicator */
  int                         n_ranks;         /* Size of communicator */

  vtkMultiBlockDataSet       *mb;              /* Associated dataset */

  fvm_writer_time_dep_t       time_dependency; /* Mesh time dependency */

  int                         n_fields;        /* Number of fields */
  fvm_catalyst_field_t      **fields;          /* Array of field helper
                                                  structures */

  int                         time_step;       /* Latest time step */
  double                      time_value;      /* Latest time value */

  char                       *input_name;      /* input name, or NULL for
                                                  default */
  bool                        private_comm;    /* Use private communicator */
  bool                        ensight_names;   /* Use EnSight rules for
                                                  field names */

  bool                        modified;        /* Has output been added since
                                                  last coprocessing ? */

} fvm_to_catalyst_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

int _n_writers = 0;

int _n_scripts = 0;                         /* Number of scripts */
char **_scripts = NULL;                     /* List of scripts */

vtkCPProcessor  *_processor = NULL;         /* Co processor */

#if defined(HAVE_MPI)
MPI_Comm  _comm = MPI_COMM_NULL;            /* Associated communicator */
MPI_Comm  _reference_comm = MPI_COMM_NULL;  /* Reference communicator */
#endif

/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Check if a script is a Catalyst script
 *
 * The script only does cursory checks, so may return false positives
 * in some cases, but using "incorrect" scripts for Catalyst only
 * leads to extra warnings, which should be ok pending a more general
 * file checker.
 *
 * parameters:
 *   path <-- scripts path
 *
 * returns:
 *   1 if script is a Catalyst V1, 2 if it seems to be a V2 script,
 *   script, 0 otherwise
 *----------------------------------------------------------------------------*/

static int
_check_script_is_catalyst(const char  *path)
{
  assert(path != NULL);
  int retval = 0;

  FILE *fp = fopen(path, "r");

  if (fp == NULL)
    return retval;

  int checks[] = {false, false, false};
  int n_checks = 0;
  bool import_catalyst = false;

  const char *check_strings[]
    = {"CreateCoProcessor(",
       "RequestDataDescription(",
       "DoCoProcessing("};

  /* Note: we could simply check for "from paraview import coprocessing"
     in most cases but in cases of non-autogenerated scripts, this might
     appear in a different form, so we limit ourselves to the main
     required methods */

  while (true) {
    char buffer[1024];
    char *e = buffer+1024;

    char *s = fgets(buffer, 1024, fp);
    if (s == NULL) break;

    while (s < e && *s == ' ' && *s == '\t')  /* skip initial whitespace */
      s++;

    if (strncmp(s, "def", 3) == 0) {
      s += 3;
      if (*s == ' ' || *s == '\t') {
        while (*s == ' ' && *s == '\t' && *s != '\0')
          s++;
        /* Remove whitespace */
        size_t l = strlen(s);
        size_t i, j;
        for (i = 0, j = 0 ; i < l ; i++) {
          if (s[i] != ' ' && s[i] != '\t')
            s[j++] = s[i];
        }

        for (i = 0; i < 3; i++) {
          if (strncmp(s, check_strings[i], strlen(check_strings[i])) == 0) {
            if (checks[i] == false) {
              checks[i] = true;
              n_checks += 1;
            }
          }
        }

      }
    }
    else if (strncmp(s, "from", 4) == 0) {

      /* cleanup whitespace */
      int n_space = 0;
      int i = 0;
      for (int j = 0; j < 1024 && s+j < e && s[j] != '\0'; j++) {
        if (s[j] == ' ' || s[j] == '\t') {
          if (n_space < 1)
            s[i++] = ' ';
          n_space += 1;
        }
        else {
          s[i++] = s[j];
          n_space = 0;
        }
      }

      if (strncmp(s, "from paraview import catalyst", 29) == 0) {
        import_catalyst = true;
      }

    }

    if (n_checks == 3) {
      retval = 1;
      break;
    }
    else if (import_catalyst) {
      retval = 2;
      break;
    }

    /* Check for end of line; if not present, continue reading from buffer */

    while (s < e && *s != '\0' && *s != '\n')
      s++;
    while (s >= e) {
      s = fgets(buffer, 1024, fp);
      if (s == NULL) break;
      while (s < e && *s != '\0' && *s != '\n')
        s++;
    }
  }

  fclose(fp);

  return retval;
}

/*----------------------------------------------------------------------------
 * Add a Catalyst V2 pipeline if not already present.
 *
 * parameters:
 *   path <-- V2 pipeline file path
 *
 * returns:
 *   id of pipeline script in list, or -1 if not valid
 *----------------------------------------------------------------------------*/

static int
_add_v2_pipeline(const char  *path)
{
  assert(path != NULL);

  /* Check that we did not already add this file */

  for (int i = 0; i < _n_scripts; i++) {
    if (strcmp(_scripts[i], path) == 0)
      return i;
  }

#if defined(HAVE_VTKCPPYTHONSCRIPTV2PIPELINE_H)

  /* Create Catalyst pipeline and check the file is valid */

  vtkNew<vtkCPPythonScriptV2Pipeline> pipeline;
  if (pipeline->Initialize(path) == false) {
    bft_printf(_("\nFile \"%s\"\n"
                 "does not seem to contain a Catalyst Python pipeline.\n"),
               path);
    return -1;
  }

  int id = _n_scripts;
  _n_scripts += 1;

  BFT_REALLOC(_scripts, _n_scripts, char *);

  size_t l = strlen(path);
  BFT_MALLOC(_scripts[id], l+1, char);
  strncpy(_scripts[id], path, l);
  _scripts[id][l] = '\0';

  /* pipeline->SetGhostLevel(1); */
  _processor->AddPipeline(pipeline);

  return id;

#else

  bft_printf
    (_("\nFile \"%s\"\n"
       "seems to contain a Catalyst V2 Python pipeline but the version.\n"
       "of Catalyst linked with only supports version 1 pipelines.\n\n"
       "You may need to use a more recent Catalyst version or downgrade\n"
       "to earlier V1 pipelines.\n"),
     path);

  return -1;

#endif /* defined(HAVE_VTKCPPYTHONSCRIPTV2PIPELINE_H) */

}

/*----------------------------------------------------------------------------
 * Add a Catalyst script if not already present.
 *
 * parameters:
 *   path <-- scripts path
 *
 * returns:
 *   id of script in list, or -1 if not valid
 *----------------------------------------------------------------------------*/

static int
_add_script(const char         *path)
{
  assert(path != NULL);

  int is_catalyst = 0;
  int rank = 0, n_ranks = 1;

#if defined(HAVE_MPI)
  if (_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(_comm, &rank);
    MPI_Comm_size(_comm, &n_ranks);
  }
#endif

  if (rank < 1)
    is_catalyst = _check_script_is_catalyst(path);

#if defined(HAVE_MPI)
  if (n_ranks > 1)
    MPI_Bcast(&is_catalyst, 1, MPI_INT, 0, _comm);
#endif

  if (is_catalyst < 1)
    return -1;

  if (is_catalyst == 2) {
    _add_v2_pipeline(path);
    return 2;
  }

  for (int i = 0; i < _n_scripts; i++) {
    if (strcmp(_scripts[i], path) == 0)
      return i;
  }

  int id = _n_scripts;
  _n_scripts += 1;

  BFT_REALLOC(_scripts, _n_scripts, char *);

  size_t l = strlen(path);
  BFT_MALLOC(_scripts[id], l+1, char);
  strncpy(_scripts[id], path, l);
  _scripts[id][l] = '\0';

  /* Add Catalyst pipeline */

  vtkNew<vtkCPPythonScriptPipeline> pipeline;
  if (pipeline->Initialize(path) == 0)
    bft_error(__FILE__, __LINE__, 0,
              _("Error initializing pipeline from \"%s\""), path);

  /* pipeline->SetGhostLevel(1); */
  _processor->AddPipeline(pipeline);

  return id;
}

/*----------------------------------------------------------------------------
 * Add Catalyst scripts from directoty if not already present.
 *
 * Currently assumes all Python files in the given directory
 * are Catalyst scripts.
 *
 * parameters:
 *   dir_path <-- directory path
 *----------------------------------------------------------------------------*/

static void
_add_dir_scripts(const char  *dir_path)
{
  char **dir_files = cs_file_listdir(dir_path);

  for (int i = 0; dir_files[i] != NULL; i++) {

    const char *file_name = dir_files[i];
    const char *ext = NULL;
    int l_ext = 0;

    /* Find extension */
    for (int j = strlen(file_name) - 1; j > -1; j--) {
      l_ext++;
      if (file_name[j] == '.') {
        ext = file_name + j;
        break;
      }
    }
    if (ext == NULL) {
      BFT_FREE(dir_files[i]);
      continue;
    }

    /* Filter: Python files only */
    if (l_ext == 3 && strncmp(ext, ".py", 3) == 0) {
      char *tmp_name = NULL;
      BFT_MALLOC(tmp_name,
                 strlen(dir_path) + 1 + strlen(file_name) + 1,
                 char);
      sprintf(tmp_name, "%s/%s", dir_path, file_name);
      _add_script(tmp_name);
      BFT_FREE(tmp_name);
    }

#if defined(HAVE_VTKCPPYTHONSCRIPTV2PIPELINE_H)

    /* Filter: Catalyst V2 pipeline might be in ".zip" form */
    else if (l_ext == 4 && strncmp(ext, ".zip", 4) == 0) {
      char *tmp_name = NULL;
      BFT_MALLOC(tmp_name,
                 strlen(dir_path) + 1 + strlen(file_name) + 1,
                 char);
      sprintf(tmp_name, "%s/%s", dir_path, file_name);
      _add_v2_pipeline(tmp_name);
      BFT_FREE(tmp_name);
    }
#endif /* defined(HAVE_VTKCPPYTHONSCRIPTV2PIPELINE_H) */

    BFT_FREE(dir_files[i]);
  }

  BFT_FREE(dir_files);
}

/*----------------------------------------------------------------------------
 * Initialize coprocessor.
 *
 * parameters:
 *   private_comm <-- if true, use dedicated communicator
 *   comm         <-- associated MPI communicator.
 *----------------------------------------------------------------------------*/

static void
#if defined(HAVE_MPI)
_init_coprocessor(bool      private_comm,
                  MPI_Comm  comm)
#else
_init_coprocessor(void)
#endif
{
  int mpi_flag = 0;
  int mpi_rank = -1;

#if defined(HAVE_MPI)

  MPI_Initialized(&mpi_flag);

  if (mpi_flag) {
    if (comm != _reference_comm) {
      if (comm != MPI_COMM_NULL && _reference_comm != MPI_COMM_NULL)
        bft_error(__FILE__, __LINE__, 0,
                  _("All Catalyst writers must use the same MPI communicator"));
    }
  }

#endif

  if (_processor == NULL) {

    _processor = vtkCPProcessor::New();

#if defined(HAVE_MPI)

    if (mpi_flag) {

      _reference_comm = comm;
      if (private_comm && comm != MPI_COMM_NULL)
        MPI_Comm_dup(comm, &(_comm));
      else
        _comm = comm;

      if (comm != MPI_COMM_NULL)
        MPI_Comm_rank(_comm, &mpi_rank);

      vtkMPICommunicatorOpaqueComm vtk_comm
        = vtkMPICommunicatorOpaqueComm(&_comm);
      _processor->Initialize(vtk_comm);

    }

#endif

    if (!mpi_flag)
      _processor->Initialize();

    vtkFileOutputWindow *log_output = vtkFileOutputWindow::New();
    if (mpi_rank < 1)
      log_output->SetFileName("./catalyst.log");
    else
      log_output->SetFileName("/dev/null");
    vtkFileOutputWindow::SetInstance(log_output);

  }
}

/*----------------------------------------------------------------------------
 * Finalize coprocessor.
 *
 * parameters:
 *   private_comm <-- if true, use dedicated communicator
 *   comm         <-- associated MPI communicator.
 *----------------------------------------------------------------------------*/

static void
_free_coprocessor(void)
{
  if (_processor != NULL && _n_writers < 2) {

    _processor->Finalize();

    /* Workaround for segmentation fault observed on a Red Hat 8.3 system
       with gcc 8.3.1, at coprocessor destruction, not observed on other
       machines. */

    bool cp_delete = true;

    if (_n_scripts > 0) {
      const char *s = getenv("CS_PV_CP_DELETE_CRASH_WORKAROUND");
      if (s != NULL) {
        if (atoi(s) > 0)
          cp_delete = false;
      }
    }
    if (cp_delete)
      _processor->Delete();

    _processor = NULL;

    for (int i = 0; i < _n_scripts; i++)
      BFT_FREE(_scripts[i]);
    BFT_FREE(_scripts);

#if defined(HAVE_MPI)
    {
      if (_comm != _reference_comm && _comm != MPI_COMM_NULL)
        MPI_Comm_free(&_comm);
    }
#endif

  }
}

/*----------------------------------------------------------------------------
 * Return the Catalyst mesh id associated with a given mesh name,
 * or -1 if no association found.
 *
 * parameters:
 *   writer    <-- writer structure
 *   mesh_name <-- mesh name
 *
 * returns:
 *    Catalyst mesh id, or -1 if Catalyst mesh name is not associated
 *    with this writer structure in FVM
 *----------------------------------------------------------------------------*/

static int
_get_catalyst_mesh_id(fvm_to_catalyst_t  *writer,
                      const char         *mesh_name)
{
  int i;
  int retval = -1;

  assert(writer != NULL);

  int nb = writer->mb->GetNumberOfBlocks();
  for (i = 0; i < nb; i++) {
    if (strcmp
         (mesh_name,
          writer->mb->GetMetaData(i)->Get(vtkCompositeDataSet::NAME())) == 0)
      break;
  }

  if (i < nb)
    retval = i;

  return retval;
}

/*----------------------------------------------------------------------------
 * Create a Catalyst mesh structure.
 *
 * parameters:
 *   writer    <-- Catalyst structure.
 *   mesh      <-- FVM mesh  structure.
 *
 * returns:
 *   Catalyst mesh id.
 *----------------------------------------------------------------------------*/

static int
_add_catalyst_mesh(fvm_to_catalyst_t  *writer,
                   const fvm_nodal_t  *mesh)
{
  assert(writer != NULL);

  /* Add a new Catalyst mesh structure */

  int id = writer->mb->GetNumberOfBlocks();
  vtkSmartPointer<vtkUnstructuredGrid> ugrid
    = vtkSmartPointer<vtkUnstructuredGrid>::New();

  writer->mb->SetBlock(id, ugrid);
  writer->mb->GetMetaData(id)->Set(vtkCompositeDataSet::NAME(), mesh->name);

  return id;
}

/*----------------------------------------------------------------------------
 * Define VTK geometrical element type according to FVM element type
 *
 * parameters:
 *   fvm_elt_type <-- pointer to fvm element type.
 *
 * return:
 *   med geometrical element type.
 *----------------------------------------------------------------------------*/

static VTKCellType
_get_norm_elt_type(const fvm_element_t fvm_elt_type)
{
  VTKCellType  norm_elt_type;

  switch (fvm_elt_type) {

  case FVM_EDGE:
    norm_elt_type = VTK_LINE;
    break;

  case FVM_FACE_TRIA:
    norm_elt_type = VTK_TRIANGLE;
    break;

  case FVM_FACE_QUAD:
    norm_elt_type = VTK_QUAD;
    break;

  case FVM_FACE_POLY:
    norm_elt_type = VTK_POLYGON;
    break;

  case FVM_CELL_TETRA:
    norm_elt_type = VTK_TETRA;
    break;

  case FVM_CELL_PYRAM:
    norm_elt_type = VTK_PYRAMID;
    break;

  case FVM_CELL_PRISM:
    norm_elt_type = VTK_WEDGE;
    break;

  case FVM_CELL_HEXA:
    norm_elt_type = VTK_HEXAHEDRON;
    break;

  case FVM_CELL_POLY:
    norm_elt_type = VTK_POLYHEDRON;
    break;

  default:
    norm_elt_type = VTK_EMPTY_CELL;
    bft_error(__FILE__, __LINE__, 0,
              "_get_norm_elt_type(): "
              "No association with VTK element type has been found\n"
              "FVM element type: \"%i\"\n",
              (int)fvm_elt_type);

  } /* End of switch on element type */

  return norm_elt_type;
}

/*----------------------------------------------------------------------------
 * Get vertex order to describe Catalyst element type.
 *
 * parameters:
 *   norm_elt_type  <-- Catalyst element type.
 *   vertex_order  --> Pointer to vertex order array (0 to n-1).
 *----------------------------------------------------------------------------*/

static void
_get_vertex_order(VTKCellType   norm_elt_type,
                  int          *vertex_order)
{
  switch(norm_elt_type) {

  case VTK_LINE:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    break;

  case VTK_TRIANGLE:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    break;

  case VTK_QUAD:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    break;

  case VTK_TETRA:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    break;

  case VTK_PYRAMID:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    vertex_order[4] = 4;
    break;

  case VTK_WEDGE:
    vertex_order[0] = 0;
    vertex_order[1] = 2;
    vertex_order[2] = 1;
    vertex_order[3] = 3;
    vertex_order[4] = 5;
    vertex_order[5] = 4;
    break;

  case VTK_HEXAHEDRON:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    vertex_order[4] = 4;
    vertex_order[5] = 5;
    vertex_order[6] = 6;
    vertex_order[7] = 7;
    break;

  case VTK_POLYGON:
    vertex_order[0] = -1;
    break;

  case VTK_POLYHEDRON:
    vertex_order[0] = -1;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "_get_vertex_order(): No associated FVM element type known\n"
              "VTK element type: \"%i\"\n",
              (int)norm_elt_type);
  }

  return;
}

/*----------------------------------------------------------------------------
 * Return the Catalyst field id associated with given mesh and field names,
 * or -1 if no association found.
 *
 * parameters:
 *   writer    <-- Catalyst writer structure.
 *   fieldname <-- input fieldname.
 *   mesh_id   <-- id of associated mesh in structure.
 *   dim       <-- number of field components.
 *   location  <-- mesh location (cells or vertices).
 *   type      <-- field type (cells, nodes)
 *   td        <-- time discretization type
 *
 * returns
 *   field id in writer structure, or -1 if not found
 *----------------------------------------------------------------------------*/

static int
_get_catalyst_field_id(fvm_to_catalyst_t         *writer,
                       const char                *fieldname,
                       int                        mesh_id,
                       int                        dim,
                       cs_datatype_t              datatype,
                       fvm_writer_var_loc_t       location)
{
  CS_UNUSED(dim);
  CS_UNUSED(datatype);

  int i;

  for (i = 0; i < writer->n_fields; ++i){

    vtkFieldData *fData_ptr ;

    if (writer->fields[i]->mesh_id != mesh_id)
      continue;

    if (location == FVM_WRITER_PER_NODE)
      fData_ptr = (writer->fields[i])->f->GetPointData();
    else
      fData_ptr = writer->fields[i]->f->GetCellData();

    if (fData_ptr->HasArray(fieldname)){
      break;
    }
  }

  if (i < writer->n_fields)
    return i;

  return -1;
}

/*----------------------------------------------------------------------------
 * Create a Catalyst field structure.
 *
 * parameters:
 *   writer     <-- Catalyst writer structure.
 *   fieldname  <-- input fieldname.
 *   mesh_id    <-- id of associated mesh in structure.
 *   dim        <-- number of field components.
 *   location   <-- mesh location (cells or vertices).
 *
 * returns
 *   field id in writer structure
 *----------------------------------------------------------------------------*/

static int
_add_catalyst_field(fvm_to_catalyst_t         *writer,
                    const char                *fieldname,
                    int                        mesh_id,
                    int                        dim,
                    cs_datatype_t              datatype,
                    fvm_writer_var_loc_t       location)
{
  CS_UNUSED(datatype);

  int f_id = writer->n_fields;

  BFT_REALLOC(writer->fields,
              writer->n_fields+ 1,
              fvm_catalyst_field_t *);

  vtkUnstructuredGrid *f = NULL;

  const int dest_dim = (dim == 6) ? 9 : dim;

  if (writer->mb->GetMetaData(mesh_id) != NULL) {

    f = vtkUnstructuredGrid::SafeDownCast(vtkDataSet::SafeDownCast
          (writer->mb->GetBlock(mesh_id)));

    vtkNew<vtkDoubleArray> tmp;
    tmp->SetName(fieldname);

    tmp->SetNumberOfComponents(dest_dim);

    if (location == FVM_WRITER_PER_NODE) {
      tmp->SetNumberOfTuples(f->GetNumberOfPoints());
      // f->GetPointData()->AllocateArrays(f->GetNumberOfPoints());
      vtkDataSetAttributes::SafeDownCast
        (f->GetPointData())->AddArray(vtkAbstractArray::SafeDownCast(tmp));
    }
    else {
      tmp->SetNumberOfTuples(f->GetNumberOfCells());
      // f->GetCellData()->AllocateArrays(f->GetNumberOfCells());
      /* force to AddArray, due to protection of SetArray */
      vtkDataSetAttributes::SafeDownCast
        (f->GetCellData())->AddArray(vtkAbstractArray::SafeDownCast(tmp));
    }

  }

  BFT_REALLOC(writer->fields,
              writer->n_fields+ 1,
              fvm_catalyst_field_t *);

  BFT_MALLOC(writer->fields[f_id], 1, fvm_catalyst_field_t);

  writer->fields[f_id]->mesh_id = mesh_id;

  writer->fields[f_id]->dim = dim;
  writer->fields[f_id]->f = f;

  writer->n_fields++;

  return f_id;
}

/*----------------------------------------------------------------------------
 * Write vertex coordinates to VTK.
 *
 * parameters:
 *   mesh        <-- pointer to nodal mesh structure
 *   vtk_mesh    <-- pointer to VTK Mesh object
 *----------------------------------------------------------------------------*/

static void
_export_vertex_coords(const fvm_nodal_t        *mesh,
                      vtkUnstructuredGrid      *vtk_mesh)
{
  cs_lnum_t   i, j;
  size_t stride;

  double point[3];

  const double  *vertex_coords = mesh->vertex_coords;
  const cs_lnum_t  n_vertices = mesh->n_vertices;

  /* Vertex coordinates */
  /*--------------------*/

  stride = (size_t)(mesh->dim);

  vtkNew<vtkPoints> points;

  points->Allocate(mesh->n_vertices);

  if (mesh->parent_vertex_id != NULL) {
    const cs_lnum_t  *parent_vertex_id = mesh->parent_vertex_id;
    for (i = 0; i < n_vertices; i++) {
      for (j = 0; j < mesh->dim; j++)
        point[j] = vertex_coords[parent_vertex_id[i]*stride + j];
      for (; j < 3; j++)
        point[j] = 0.;
      points->InsertNextPoint(point[0], point[1], point[2]);
    }
  }
  else {
    for (i = 0; i < n_vertices; i++) {
      for (j = 0; j < mesh->dim; j++)
        point[j] = vertex_coords[i*stride + j];
      for (; j < 3; j++)
        point[j] = 0.;
      points->InsertNextPoint(point[0], point[1], point[2]);
    }
  }

  if (mesh->global_vertex_num != NULL) {

    const cs_gnum_t *g_vtx_num
      = fvm_io_num_get_global_num(mesh->global_vertex_num);

    vtkNew<vtkIdTypeArray> g_vtx_id;

    g_vtx_id->SetNumberOfComponents(1);
    g_vtx_id->SetName("GlobalNodeIds");
    g_vtx_id->SetNumberOfTuples(n_vertices);

    for (i = 0; i < n_vertices; i++) {
      vtkIdType ii = g_vtx_num[i]-1;
#if CS_PV_VERSION < 51
      g_vtx_id->SetTupleValue(i, &ii);
#else
      g_vtx_id->SetTypedTuple(i, &ii);
#endif
    }

    vtk_mesh->GetPointData()->SetGlobalIds(g_vtx_id);

  }

  vtk_mesh->SetPoints(points);
}

/*----------------------------------------------------------------------------
 * Write strided connectivity block to VTK.
 *
 * The connectivity on input may use 1 to n numbering, so it is shifted
 * by -1 here.
 *
 * parameters:
 *   type     <-- FVM element type
 *   n_elts   <-- number of elements in block
 *   connect  <-- connectivity array
 *   vtk_mesh <-> pointer to VTK mesh object
 *----------------------------------------------------------------------------*/

static void
_write_connect_block(fvm_element_t         type,
                     cs_lnum_t             n_elts,
                     const cs_lnum_t       connect[],
                     vtkUnstructuredGrid  *vtk_mesh)
{
  int vertex_order[8];
  cs_lnum_t  i;
  int  j;

  vtkIdType *vtx_ids = new vtkIdType[8];

  const int  stride = fvm_nodal_n_vertices_element[type];
  VTKCellType vtk_type = _get_norm_elt_type(type);

  _get_vertex_order(vtk_type, vertex_order);

  assert(vtk_mesh != NULL);

  for (i = 0; i < n_elts; i++) {
    for (j = 0; j < stride; j++)
      vtx_ids[j] = connect[i*stride + vertex_order[j]] - 1;
    vtk_mesh->InsertNextCell(vtk_type, stride, vtx_ids);
  }

  delete [] vtx_ids;
}

/*----------------------------------------------------------------------------
 * Write polygons from a nodal mesh to VTK.
 *
 * parameters:
 *   export_section <-- pointer to Catalyst section helper structure
 *   vtk_mesh       <-> pointer to VTK mesh object
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polygons(const fvm_nodal_section_t  *section,
                       vtkUnstructuredGrid        *vtk_mesh)
{
  cs_lnum_t   i, j;

  int vtx_ids_size = 8;
  vtkIdType *vtx_ids = NULL;

  BFT_MALLOC(vtx_ids, vtx_ids_size, vtkIdType);

  /* Loop on all polygonal faces */
  /*-----------------------------*/

  for (i = 0; i < section->n_elements; i++) {

    int k = 0;

    int face_size = section->vertex_index[i+1] - section->vertex_index[i];
    while (vtx_ids_size < face_size) {
      vtx_ids_size *= 2;
      BFT_REALLOC(vtx_ids, vtx_ids_size, vtkIdType);
    }

    for (j = section->vertex_index[i];
         j < section->vertex_index[i+1];
         j++)
      vtx_ids[k++] = section->vertex_num[j] - 1;

    vtk_mesh->InsertNextCell(VTK_POLYGON, k, vtx_ids);

  } /* End of loop on polygonal faces */

  BFT_FREE(vtx_ids);
}

/*----------------------------------------------------------------------------
 * Write polyhedra from a nodal mesh to VTK.
 *
 * parameters:
 *   n_vertices <-- number of vertices in FVM mesh
 *   section    <-- pointer to Catalyst section helper structure
 *   vtk_mesh   <-> VTK mesh object
 *----------------------------------------------------------------------------*/

static void
_export_nodal_polyhedra(cs_lnum_t                   n_vertices,
                        const fvm_nodal_section_t  *section,
                        vtkUnstructuredGrid        *vtk_mesh)
{
  int  face_sgn;
  int  buf_size = 8;
  cs_lnum_t  i, j, k, l;

  cs_lnum_t  face_length, face_id;

  int *vtx_marker = NULL;
  vtkIdType *face_array = NULL;
  vtkIdType *vtx_ids = NULL;

  BFT_MALLOC(vtx_marker, n_vertices, int);
  BFT_MALLOC(face_array, buf_size, vtkIdType);
  BFT_MALLOC(vtx_ids, buf_size, vtkIdType);

  for (i = 0; i < n_vertices; i++)
    vtx_marker[i] = -1;

  /* Write cell/vertex connectivity */
  /*--------------------------------*/

  for (i = 0; i < section->n_elements; i++) {

    int  m = 0;
    int  cell_vtx_count = 0;
    int  n_elt_faces = section->face_index[i+1] - section->face_index[i];

    for (j = section->face_index[i];
         j < section->face_index[i+1];
         j++) {

      if (section->face_num[j] > 0) {
        face_id = section->face_num[j] - 1;
        face_sgn = 1;
      }
      else {
        face_id = -section->face_num[j] - 1;
        face_sgn = -1;
      }

      face_length = (  section->vertex_index[face_id+1]
                     - section->vertex_index[face_id]);

      while (m + face_length + 1 > buf_size) {
        buf_size *= 2;
        BFT_REALLOC(face_array, buf_size, vtkIdType);
        BFT_REALLOC(vtx_ids, buf_size, vtkIdType);
      }

      face_array[m++] = face_length;

      for (k = 0; k < face_length; k++) {
        int vtx_id;
        l =    section->vertex_index[face_id]
            + (face_length + (k*face_sgn))%face_length;
        vtx_id = section->vertex_num[l] - 1;
        face_array[m++] = vtx_id;
        if (vtx_marker[vtx_id] < i) {
          vtx_ids[cell_vtx_count] = vtx_id;
          vtx_marker[vtx_id] = i;
          cell_vtx_count += 1;
        }
      }

    } /* End of loop on cell faces */

    vtk_mesh->InsertNextCell(VTK_POLYHEDRON, cell_vtx_count,
                             vtx_ids, n_elt_faces, face_array);

  } /* End of loop on polyhedral cells */

  BFT_FREE(vtx_ids);
  BFT_FREE(face_array);
  BFT_FREE(vtx_marker);
}

/*----------------------------------------------------------------------------
 * Write field values associated with nodal values of a nodal mesh to VTK.
 *
 * Output fields ar either scalar or 3d vectors or scalars, and are
 * non interlaced. Input arrays may be less than 2d, in which case the z
 * values are set to 0, and may be interlaced or not.
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   dim              <-- field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- input data type (output is real)
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated vtkUnstructuredGrid object
 *----------------------------------------------------------------------------*/

static void
_export_field_values_n(const fvm_nodal_t    *mesh,
                       const char           *fieldname,
                       int                   dim,
                       cs_interlace_t        interlace,
                       int                   n_parent_lists,
                       const cs_lnum_t       parent_num_shift[],
                       cs_datatype_t         datatype,
                       const void           *const field_values[],
                       vtkUnstructuredGrid  *f)
{
  assert(f != NULL);

  double *values = NULL;

  const int dest_dim = (dim == 6) ? 9 : dim;

  values = vtkDoubleArray::SafeDownCast
    (f->GetPointData()->GetArray(fieldname))
    ->WritePointer(0, dest_dim*f->GetNumberOfPoints());

  fvm_convert_array(dim,
                    0,
                    dest_dim,
                    0,
                    mesh->n_vertices,
                    interlace,
                    datatype,
                    CS_DOUBLE,
                    n_parent_lists,
                    parent_num_shift,
                    mesh->parent_vertex_id,
                    field_values,
                    values);

  /* Special case for symmetric tensors */

  if (dim == 6) {
    for (cs_lnum_t i = 0; i < mesh->n_vertices; i++) {
      values[9*i + 8] = values[9*i + 2];
      values[9*i + 7] = values[9*i + 4];
      values[9*i + 6] = values[9*i + 5];
      values[9*i + 4] = values[9*i + 1];
      values[9*i + 2] = values[9*i + 5];
      values[9*i + 1] = values[9*i + 3];
      values[9*i + 5] = values[9*i + 7];
    }
  }

}

/*----------------------------------------------------------------------------
 * Write field values associated with element values of a nodal mesh to VTK.
 *
 * Output fields are non interlaced. Input arrays may be interlaced or not.
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   dim              <-- field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated object handle
 *----------------------------------------------------------------------------*/

static void
_export_field_values_e(const fvm_nodal_t         *mesh,
                       const char                *fieldname,
                       int                        dim,
                       cs_interlace_t             interlace,
                       int                        n_parent_lists,
                       const cs_lnum_t            parent_num_shift[],
                       cs_datatype_t              datatype,
                       const void          *const field_values[],
                       vtkUnstructuredGrid       *f)
{
  assert(f != NULL);

  int  section_id;
  double  *values = NULL;

  const int dest_dim = (dim == 6) ? 9 : dim;

  values = vtkDoubleArray::SafeDownCast
    (f->GetCellData()->GetArray(fieldname))
    ->WritePointer(0, dest_dim*f->GetNumberOfCells());

  /* Distribute partition to block values */

  cs_lnum_t start_id = 0;
  cs_lnum_t src_shift = 0;

  /* loop on sections which should be appended */

  const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

    assert(values != NULL || section->n_elements == 0);

    fvm_convert_array(dim,
                      0,
                      dest_dim,
                      src_shift,
                      section->n_elements + src_shift,
                      interlace,
                      datatype,
                      CS_DOUBLE,
                      n_parent_lists,
                      parent_num_shift,
                      section->parent_element_id,
                      field_values,
                      values + start_id);

    start_id += section->n_elements*dest_dim;
    if (n_parent_lists == 0)
      src_shift += section->n_elements;

  }

  /* Special case for symmetric tensors */

  if (dim == 6) {

    cs_lnum_t n_elts = f->GetNumberOfCells();

    assert(values != NULL || n_elts == 0);

    for (cs_lnum_t i = 0; i < n_elts; i++) {
      values[9*i + 8] = values[9*i + 2];
      values[9*i + 7] = values[9*i + 4];
      values[9*i + 6] = values[9*i + 5];
      values[9*i + 4] = values[9*i + 1];
      values[9*i + 2] = values[9*i + 5];
      values[9*i + 1] = values[9*i + 3];
      values[9*i + 5] = values[9*i + 7];
    }
  }
}

/*============================================================================
 * Public function definitions
 *============================================================================*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------
 * Initialize FVM to Catalyst object writer.
 *
 * Options are:
 *   private_comm        use private MPI communicator (default: false)
 *   names=<fmt>         use same naming rules as <fmt> format
 *                       (default: ensight)
 *   input_name=<name>   define input name (default: writer name)
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separaed, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque Catalyst writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_catalyst_init_writer(const char             *name,
                            const char             *path,
                            const char             *options,
                            fvm_writer_time_dep_t   time_dependency,
                            MPI_Comm                comm)
#else
void *
fvm_to_catalyst_init_writer(const char             *name,
                            const char             *path,
                            const char             *options,
                            fvm_writer_time_dep_t   time_dependency)
#endif
{
  CS_UNUSED(path);

  fvm_to_catalyst_t  *w = NULL;

  bool private_comm = false;

  /* Initialize writer */

  BFT_MALLOC(w, 1, fvm_to_catalyst_t);

  w->rank = 0;
  w->n_ranks = 1;

  w->time_dependency = time_dependency;

  w->mb = vtkMultiBlockDataSet::New();

  w->n_fields  = 0;
  w->fields = NULL;

  w->time_step  = -1;
  w->time_value = 0.0;

  w->ensight_names = true;
  w->input_name = NULL;

  /* Writer name */

  if (name != NULL) {
    BFT_MALLOC(w->name, strlen(name) + 1, char);
    strcpy(w->name, name);
  }
  else {
    const char _name[] = "catalyst";
    BFT_MALLOC(w->name, strlen(_name) + 1, char);
    strcpy(w->name, _name);
  }

  /* Parse options */

  if (options != NULL) {

    int i1, i2, l_opt;
    int l_tot = strlen(options);

    i1 = 0; i2 = 0;
    while (i1 < l_tot) {

      for (i2 = i1; i2 < l_tot && options[i2] != ' '; i2++);
      l_opt = i2 - i1;

      if ((l_opt == 12) && (strncmp(options + i1, "private_comm", l_opt) == 0))
        private_comm = true;
      else if ((l_opt > 6) && (strncmp(options + i1, "names=", 6) == 0)) {
        if ((l_opt == 6+7) && (strncmp(options + i1 + 6, "ensight", 7) == 0))
          w->ensight_names = true;
        else
          w->ensight_names = false;
      }
      else if ((l_opt > 11) && (strncmp(options + i1, "input_name=", 11) == 0)) {
        int l = strlen(options + i1 + 11);
        BFT_MALLOC(w->input_name, l+1, char);
        strncpy(w->input_name, options + i1 + 11, l);
        w->input_name[l] = '\0';
      }

      for (i1 = i2 + 1; i1 < l_tot && options[i1] == ' '; i1++);

    }

  }

#if CS_PV_VERSION < 55
  if (w->input_name == NULL) {
    char n[] = "input";
    BFT_MALLOC(w->input_name, strlen(n)+1, char);
    strcpy(w->input_name, n);
  }
#endif

  /* Ensure coprocessor is built */

  _n_writers += 1;

#if defined(HAVE_MPI)
  _init_coprocessor(private_comm, comm);
#else
  _init_coprocessor();
#endif

  /* Parallel parameters */

#if defined(HAVE_MPI)
  if (_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(_comm, &(w->rank));
    MPI_Comm_size(_comm, &(w->n_ranks));
  }
#endif

  if (_n_scripts < 1)
    _add_dir_scripts(".");

  w->modified = true;

  return w;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to Catalyst object writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_catalyst_finalize_writer(void  *this_writer_p)
{
  int i;

  fvm_to_catalyst_t  *w = (fvm_to_catalyst_t *)this_writer_p;

  assert(w != NULL);

  /* Write output if not done already */

  fvm_to_catalyst_flush(this_writer_p);

  /* Free structures */

  BFT_FREE(w->name);

  /* Free vtkUnstructuredGrid and field structures
     (reference counters should go to 0) */

  _free_coprocessor();

  _n_writers -= 1;

  w->mb->Delete();

  for (i = 0; i < w->n_fields; i++) {
    w->fields[i]->f = NULL; // delete w->fields[i]->f;
    BFT_FREE(w->fields[i]);
  }

  BFT_FREE(w->fields);

  /* Free fvm_to_catalyst_t structure */

  BFT_FREE(w);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with an Catalyst geometry.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst_set_mesh_time(void    *this_writer_p,
                              int      time_step,
                              double   time_value)
{
  fvm_to_catalyst_t  *w = (fvm_to_catalyst_t *)this_writer_p;

  int _time_step = (time_step > -1) ? time_step : 0;
  double _time_value = (time_value > 0.0) ? time_value : 0.0;

  if (_time_step > w->time_step) {
    if (   w->time_dependency == FVM_WRITER_TRANSIENT_CONNECT
        && _time_step > w->time_step) {
      for (int i = 0; i < w->n_fields; i++) {
        w->fields[i]->f = NULL; // delete w->fields[i]->f;
        BFT_FREE(w->fields[i]);
      }
      BFT_FREE(w->fields);
      w->n_fields = 0;
      w->mb->Delete();
      w->mb = vtkMultiBlockDataSet::New();
    }
    w->time_step = _time_step;
    assert(time_value >= w->time_value);
    w->time_value = _time_value;
  }
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a Catalyst object
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst_export_nodal(void               *this_writer_p,
                             const fvm_nodal_t  *mesh)
{
  int  mesh_id, section_id;

  fvm_to_catalyst_t  *w = (fvm_to_catalyst_t *)this_writer_p;

  const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);

  /* Initialization */
  /*----------------*/

  /* Get matching mesh */

  mesh_id = _get_catalyst_mesh_id(w, mesh->name);

  if (mesh_id < 0)
    mesh_id = _add_catalyst_mesh(w, mesh);

  vtkUnstructuredGrid *ugrid
    = vtkUnstructuredGrid::SafeDownCast
        (vtkDataSet::SafeDownCast(w->mb->GetBlock(mesh_id)));

  /* Vertex coordinates */
  /*--------------------*/

  _export_vertex_coords(mesh, ugrid);

  /* Element connectivity size */
  /*---------------------------*/

  cs_lnum_t  n_elts = 0;

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

    n_elts += section->n_elements;

  } /* End of loop on sections */

  if (n_elts > 0)
    ugrid->Allocate(n_elts);

  /* Element connectivity */
  /*----------------------*/

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

    if (section->stride > 0)
      _write_connect_block(section->type,
                           section->n_elements,
                           section->vertex_num,
                           ugrid);

    else if (section->type == FVM_FACE_POLY)
      _export_nodal_polygons(section, ugrid);

    else if (section->type == FVM_CELL_POLY)
      _export_nodal_polyhedra(mesh->n_vertices, section, ugrid);

  } /* End of loop on sections */

  w->modified = true;
}

/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a Catalyst object.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- variable definition location (nodes or elements)
 *   dimension        <-- variable dimension (0: constant, 1: scalar,
 *                        3: vector, 6: sym. tensor, 9: asym. tensor)
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *   field_values     <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst_export_field(void                  *this_writer_p,
                             const fvm_nodal_t     *mesh,
                             const char            *name,
                             fvm_writer_var_loc_t   location,
                             int                    dimension,
                             cs_interlace_t         interlace,
                             int                    n_parent_lists,
                             const cs_lnum_t        parent_num_shift[],
                             cs_datatype_t          datatype,
                             int                    time_step,
                             double                 time_value,
                             const void      *const field_values[])
{
  int  mesh_id, field_id;
  char _name[128];

  fvm_to_catalyst_t *w = (fvm_to_catalyst_t *)this_writer_p;

  /* Initialization */
  /*----------------*/

  mesh_id = _get_catalyst_mesh_id(w, mesh->name);

  strncpy(_name, name, 127);
  _name[127] = '\0';
  if (w->ensight_names) {
    for (int i = 0; i < 127 && _name[i] != '\0'; i++) {
      switch (_name[i]) {
      case '(':
      case ')':
      case ']':
      case '[':
      case '+':
      case '-':
      case '@':
      case ' ':
      case '\t':
      case '!':
      case '#':
      case '*':
      case '^':
      case '$':
      case '/':
        _name[i] = '_';
        break;
      default:
        break;
      }
      if (_name[i] == ' ')
        _name[i] = '_';
    }
  }

  if (mesh_id < 0) {
    mesh_id = _add_catalyst_mesh(w, mesh);
    fvm_to_catalyst_export_nodal(w, mesh);
  }

  int _time_step = (time_step > -1) ? time_step : 0;
  double _time_value = (time_value > 0.0) ? time_value : 0.0;
  if (_time_step > w->time_step) {
    w->time_step = _time_step;
    assert(time_value >= w->time_value);
    w->time_value = _time_value;
  }

  /* Get field id */

  field_id = _get_catalyst_field_id(w,
                                    _name,
                                    mesh_id,
                                    dimension,
                                    datatype,
                                    location);

  if (field_id < 0)
    field_id = _add_catalyst_field(w,
                                   _name,
                                   mesh_id,
                                   dimension,
                                   datatype,
                                   location);

  vtkUnstructuredGrid  *f = w->fields[field_id]->f;

  /* Per node variable */
  /*-------------------*/

  if (location == FVM_WRITER_PER_NODE)
    _export_field_values_n(mesh,
                           _name,
                           dimension,
                           interlace,
                           n_parent_lists,
                           parent_num_shift,
                           datatype,
                           field_values,
                           f);


  /* Per element variable */
  /*----------------------*/

  else if (location == FVM_WRITER_PER_ELEMENT)
    _export_field_values_e(mesh,
                           _name,
                           dimension,
                           interlace,
                           n_parent_lists,
                           parent_num_shift,
                           datatype,
                           field_values,
                           f);

  /* Update field status */
  /*---------------------*/

  fvm_to_catalyst_set_mesh_time(w, time_step, time_value);

  w->modified = true;
}

/*----------------------------------------------------------------------------
 * Flush files associated with a given writer.
 *
 * In this case, the effective call to coprocessing is done.
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *----------------------------------------------------------------------------*/

void
fvm_to_catalyst_flush(void  *this_writer_p)
{
  fvm_to_catalyst_t *w = (fvm_to_catalyst_t *)this_writer_p;

  vtkNew<vtkCPDataDescription> dataDescription;
  if (w->input_name != NULL)
    dataDescription->AddInput(w->input_name);
  else
    dataDescription->AddInput(w->name);
  dataDescription->SetTimeData(w->time_value, w->time_step);

  if (_processor->RequestDataDescription(dataDescription) != 0 && w->modified) {
    int n = dataDescription->GetNumberOfInputDescriptions();
    if (n == 1)
      dataDescription->GetInputDescription(0)->SetGrid(w->mb);
    else {
      if (w->input_name != NULL)
        dataDescription->GetInputDescriptionByName(w->input_name)->SetGrid(w->mb);
      else
        dataDescription->GetInputDescriptionByName(w->name)->SetGrid(w->mb);
    }

    _processor->CoProcess(dataDescription);
    w->modified = false;
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* HAVE_CATALYST */
