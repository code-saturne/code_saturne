# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2022 EDF S.A.
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation; either version 2 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
# Street, Fifth Floor, Boston, MA 02110-1301, USA.

#-------------------------------------------------------------------------------

"""
This module defines the values of reference.

This module contains the following classes and function:
- MobileMeshView
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore    import *
from code_saturne.gui.base.QtGui     import *
from code_saturne.gui.base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.base.Toolbox   import GuiParam
from code_saturne.gui.case.MobileMeshForm  import Ui_MobileMeshForm
from code_saturne.gui.base.QtPage    import IntValidator,  ComboModel, from_qvariant
from code_saturne.model.MobileMeshModel import MobileMeshModel

from code_saturne.gui.case.QMegEditorView import QMegEditorView
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MobileMeshView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class MobileMeshView(QWidget, Ui_MobileMeshForm):
    """
    Class to open Page.
    """
    viscosity_iso = """# Viscosity of the mesh allows to control the deformation
# of the mesh. Viscosity must be greater than zero.
# It could be isotropic (the same for all directions) or
# orthotropic.
#
# In the following example, a hight value of viscosity
# is imposed around a mobile cylinder.
# The viscosity is specfied for all cells
# on the initial mesh before any deformation.
#
xr2 = 1.5^2;
xcen = 5.0;
ycen = 0.;
zcen = 6.0;
xray2 = (x-xcen)^2 + (y-ycen)^2 + (z-zcen)^2;
mesh_viscosity = 1;
if (xray2 < xr2) mesh_viscosity = 1e10;
"""

    viscosity_ortho = """# Viscosity of the mesh allows to control the deformation
# of the mesh. Viscosity must be greater than zero.
# It could be isotropic (the same for all directions) or
# orthotropic.
#
# In the following example, a hight value of viscosity
# is imposed around a mobile cylinder.
# The viscosity is specfied for all cells
# on the initial mesh before any deformation.
#
xr2 = 1.5^2;
xcen = 5.0;
ycen = 0.;
zcen = 6.0;
xray2 = (x-xcen)^2 + (y-ycen)^2 + (z-zcen)^2;
mesh_viscosity[X] = 1;
mesh_viscosity[Y] = 1;
mesh_viscosity[Z] = 1;
if (xray2 < xr2) {
    mesh_viscosity[X] = 1e10;
    mesh_viscosity[Y] = 1e10;
    mesh_viscosity[Z] = 1e10;
}
"""

    def __init__(self, parent, case, browser):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_MobileMeshForm.__init__(self)
        self.setupUi(self)

        self.case = case
        self.case.undoStopGlobal()
        self.mdl = MobileMeshModel(self.case)
        self.browser = browser
        self.notebook = NotebookModel(self.case)

        # Combo model VISCOSITY
        self.modelVISCOSITY = ComboModel(self.comboBoxVISCOSITY,2,1)

        self.modelVISCOSITY.addItem(self.tr("isotropic"), 'isotrop')
        self.modelVISCOSITY.addItem(self.tr("orthotropic"), 'orthotrop')

        # Connections
        self.lineEditNALINF.textChanged[str].connect(self.slotNalinf)
        self.comboBoxVISCOSITY.activated[str].connect(self.slotViscosityType)
        self.pushButtonFormula.clicked.connect(self.slotFormula)

        # Validators
        validatorNALINF = IntValidator(self.lineEditNALINF, min=0)
        self.lineEditNALINF.setValidator(validatorNALINF)

        # Settings
        nalinf = self.mdl.getSubIterations()
        self.lineEditNALINF.setText(str(nalinf))
        value = self.mdl.getViscosity()
        self.modelVISCOSITY.setItem(str_model=value)
        exp = self.mdl.getFormula()
        if exp:
            self.pushButtonFormula.setStyleSheet("background-color: green")
            self.pushButtonFormula.setToolTip(exp)
        else:
            self.pushButtonFormula.setStyleSheet("background-color: red")

        self.case.undoStartGlobal()


    @pyqtSlot(str)
    def slotNalinf(self, text):
        """
        Input viscosity type of mesh : isotrop or orthotrop.
        """
        if self.lineEditNALINF.validator().state == QValidator.Acceptable:
            nalinf = from_qvariant(text, int)
            self.mdl.setSubIterations(nalinf)


    @pyqtSlot(str)
    def slotViscosityType(self, text):
        """
        Input viscosity type of mesh : isotrop or orthotrop.
        """
        self.viscosity_type = self.modelVISCOSITY.dicoV2M[str(text)]
        visco = self.viscosity_type
        self.mdl.setViscosity(visco)
        exp = self.mdl.getFormula()
        if exp:
            self.pushButtonFormula.setStyleSheet("background-color: green")
            self.pushButtonFormula.setToolTip(exp)
        else:
            self.pushButtonFormula.setStyleSheet("background-color: red")
        return visco


    @pyqtSlot()
    def slotFormula(self):
        """
        Run formula editor.
        """
        exp = self.mdl.getFormula()

        exp, req, sca, symb = self.mdl.getFormulaViscComponents()

        if self.mdl.getViscosity() == 'isotrop':
            exa = MobileMeshView.viscosity_iso
        else:
            exa = MobileMeshView.viscosity_ortho

        dialog = QMegEditorView(parent = self,
                                function_type = 'vol',
                                zone_name     = "all_cells",
                                variable_name = "mesh_viscosity",
                                expression    = exp,
                                required      = req,
                                symbols       = symb,
                                known_fields  = [],
                                examples      = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMobileMeshView -> %s" % str(result))
            self.mdl.setFormula(str(result))
            self.pushButtonFormula.setStyleSheet("background-color: green")
            self.pushButtonFormula.setToolTip(result)


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
