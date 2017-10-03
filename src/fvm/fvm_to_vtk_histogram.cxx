/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to png histogram files
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Catalyst and VTK library headers
 *----------------------------------------------------------------------------*/

#if defined(HAVE_CATALYST)

#include <vtkAxis.h>
#include <vtkContextView.h>
#include <vtkContextScene.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkMath.h>
#include <vtkNew.h>
#include <vtkPNGWriter.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTable.h>
#include <vtkTextProperty.h>
#include <vtkVector.h>
#include <vtkWindowToImageFilter.h>
#include <vtkChartXY.h>
#include <vtkPlot.h>
#include <vtkIntArray.h>
#include <vtkPlot.h>
#include <vtkStringArray.h>

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_to_histogram.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_vtk_histogram.h"

/*----------------------------------------------------------------------------*/

using namespace std;

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * Display histograms in png format.
 *
 * parameters:
 *  var_min  <--  minimum variable value
 *  var_max  <--  maximum variable value
 *  count    <--  count for each histogram slice
 *  w        <--> histogram writer
 *  var_name <--  name of the variable
 */
/*----------------------------------------------------------------------------*/

void
fvm_to_vtk_display_histogram_png(cs_real_t                   var_min,
                                 cs_real_t                   var_max,
                                 cs_gnum_t                   count[],
                                 fvm_to_histogram_writer_t  *w,
                                 char                       *var_name)
{
  double var_step = CS_ABS(var_max - var_min) / w->n_sub;

  /* If non-zero histogram */
  if (var_step > 0) {

    VTK_CREATE(vtkContextView, view);
    VTK_CREATE(vtkChartXY, chart);
    view->GetScene()->AddItem(chart);

    /* Create a table containing the count array */
    VTK_CREATE(vtkTable, table);

    VTK_CREATE(vtkDoubleArray, arrSteps);
    arrSteps->SetName(var_name);
    table->AddColumn(arrSteps);

    VTK_CREATE(vtkIntArray, arrCount);
    arrCount->SetName("Number of elements");
    table->AddColumn(arrCount);

    table->SetNumberOfRows(w->n_sub);

    for (int i = 0; i < w->n_sub; i++) {
      table->SetValue(i, 0, var_min + (i+0.5)*var_step);
      table->SetValue(i, 1, count[i]);
    }

    /* Add plot, setting the colors etc */
    vtkPlot *line = 0;
    line = chart->AddPlot(vtkChart::BAR);
    line->SetInputData(table, 0, 1);
    line->SetColor(0, 0, 255, 255);
    chart->SetBarWidthFraction(0.5);

    VTK_CREATE(vtkStringArray, labels);
    labels->SetNumberOfValues(w->n_sub);
    char buffer[128];
    for (int i = 0; i < w->n_sub; i++) {
      sprintf(buffer, "%.4e",var_min + (i+0.5)*var_step);
      labels->SetValue(i, buffer);
    }
    line->SetIndexedLabels(labels);

    /* Set X-axis properties */
    vtkAxis *XAxis = line->GetXAxis();
    XAxis->SetTitle(var_name);
    XAxis->SetCustomTickPositions(arrSteps);
    XAxis->SetNotation(vtkAxis::PRINTF_NOTATION);
    XAxis->SetLabelFormat("%.3e");
    XAxis->GetTitleProperties()->BoldOff();
    XAxis->GetTitleProperties()->SetFontSize(14);
    line->SetXAxis(XAxis);

    /* Set Y-axis properties */
    vtkAxis *YAxis = line->GetYAxis();
    YAxis->SetTitle("Number of Elements");
    YAxis->GetTitleProperties()->BoldOff();
    YAxis->GetTitleProperties()->SetFontSize(14);
    YAxis->SetTickLabelAlgorithm(vtkAxis::TICK_WILKINSON_EXTENDED);
    line->SetYAxis(YAxis);

    /* Render window */
    vtkNew<vtkRenderWindow> renderWindow;
    renderWindow->AddRenderer(view->GetRenderer());
    renderWindow->SetSize(650, 650);
    renderWindow->OffScreenRenderingOn();
    renderWindow->Render();

    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    windowToImageFilter->SetInput(renderWindow.Get());

    /* Write the png file */
    vtkNew<vtkPNGWriter> writer;
    writer->SetFileName(w->file_name);
    writer->SetInputConnection(windowToImageFilter->GetOutputPort());
    writer->Write();
  }
}

/*----------------------------------------------------------------------------*/

#endif /* HAVE_CATALYST */

END_C_DECLS
