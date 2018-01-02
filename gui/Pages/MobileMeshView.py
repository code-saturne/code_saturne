# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2018 EDF S.A.
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

from code_saturne.Base.QtCore    import *
from code_saturne.Base.QtGui     import *
from code_saturne.Base.QtWidgets import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.Base.Toolbox   import GuiParam
from code_saturne.Pages.MobileMeshForm  import Ui_MobileMeshForm
from code_saturne.Base.QtPage    import IntValidator,  ComboModel, from_qvariant
from code_saturne.Pages.MobileMeshModel import MobileMeshModel

from code_saturne.Pages.QMeiEditorView import QMeiEditorView
from code_saturne.Pages.NotebookModel import NotebookModel

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
mesh_viscosity_1 = 1;
if (xray2 < xr2) mesh_viscosity_1 = 1e10;
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
mesh_viscosity_1 = 1;
mesh_viscosity_2 = 1;
mesh_viscosity_3 = 1;
if (xray2 < xr2) {
    mesh_viscosity_1 = 1e10;
    mesh_viscosity_2 = 1e10;
    mesh_viscosity_3 = 1e10;
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
        self.groupBoxALE.clicked.connect(self.slotMethod)
        self.lineEditNALINF.textChanged[str].connect(self.slotNalinf)
        self.comboBoxVISCOSITY.activated[str].connect(self.slotViscosityType)
        self.pushButtonFormula.clicked.connect(self.slotFormula)

        # Validators
        validatorNALINF = IntValidator(self.lineEditNALINF, min=0)
        self.lineEditNALINF.setValidator(validatorNALINF)

        if self.mdl.getMethod() == 'on':
            self.groupBoxALE.setChecked(True)
            checked = True
        else:
            self.groupBoxALE.setChecked(False)
            checked = False

        self.slotMethod(checked)

        # Enable / disable formula state
        self.pushButtonFormula.setStyleSheet("background-color: None")

        self.case.undoStartGlobal()


    @pyqtSlot(bool)
    def slotMethod(self, checked):
        """
        Private slot.

        Activates ALE method.

        @type checked: C{True} or C{False}
        @param checked: if C{True}, shows the QGroupBox ALE parameters
        """
        self.groupBoxALE.setFlat(not checked)
        if checked:
            self.frame.show()
            self.mdl.setMethod ("on")
            nalinf = self.mdl.getSubIterations()
            self.lineEditNALINF.setText(str(nalinf))
            value = self.mdl.getViscosity()
            self.modelVISCOSITY.setItem(str_model=value)
        else:
            self.frame.hide()
            self.mdl.setMethod("off")
        exp = self.mdl.getFormula()
        if exp:
            self.pushButtonFormula.setStyleSheet("background-color: green")
            self.pushButtonFormula.setToolTip(exp)
        else:
            self.pushButtonFormula.setStyleSheet("background-color: red")


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

        if self.mdl.getViscosity() == 'isotrop':
            if not exp:
                exp = "mesh_viscosity_1 = 1;"
            req = [('mesh_viscosity_1', 'mesh viscosity')]
            exa = MobileMeshView.viscosity_iso
        else:
            if not exp:
                exp = "mesh_viscosity_1 = 1;\nmesh_viscosity_2 = 1;\nmesh_viscosity_3 = 1;"
            req = [('mesh_viscosity_1', 'mesh viscosity X'),
                   ('mesh_viscosity_2', 'mesh viscosity Y'),
                   ('mesh_viscosity_3', 'mesh viscosity Z')]
            exa = MobileMeshView.viscosity_ortho

        symb = [('x', "X cell's gravity center"),
                ('y', "Y cell's gravity center"),
                ('z', "Z cell's gravity center"),
                ('dt', 'time step'),
                ('t', 'current time'),
                ('iter', 'number of iteration')]

        for (nme, val) in self.notebook.getNotebookList():
            symb.append((nme, 'value (notebook) = ' + str(val)))

        dialog = QMeiEditorView(self,
                                check_syntax = self.case['package'].get_check_syntax(),
                                expression = exp,
                                required   = req,
                                symbols    = symb,
                                examples   = exa)
        if dialog.exec_():
            result = dialog.get_result()
            log.debug("slotFormulaMobileMeshView -> %s" % str(result))
            self.mdl.setFormula(str(result))
            self.pushButtonFormula.setStyleSheet("background-color: green")
            self.pushButtonFormula.setToolTip(result)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
