# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2014 EDF S.A.
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
This module contains the following class:
- SolutionVerifView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, logging
import string, shutil

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from Base.QtPage import ComboModel, Py2
from Base.XMLengine import Case
from Base.XMLinitialize import XMLinit
from Pages.SolutionVerifForm import Ui_SolutionVerifForm
from Pages.MeshQualityCriteriaLogDialogForm import Ui_MeshQualityCriteriaLogDialogForm
from Pages.SolutionDomainModel import SolutionDomainModel
from Pages.OutputControlModel import OutputControlModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("SolutionVerifView")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

class MeshQualityCriteriaLogDialogView(QDialog, Ui_MeshQualityCriteriaLogDialogForm):
    """
    Advanced dialog
    """
    def __init__(self, parent, case, case2):
        """
        Constructor
        """
        QDialog.__init__(self, parent)

        Ui_MeshQualityCriteriaLogDialogForm.__init__(self)
        self.setupUi(self)

        self.setWindowTitle(self.tr("Run mesh quality criteria"))
        self.pushButton.setEnabled(False)

        self.case = case
        self.case2 = case2
        self.case.undoStopGlobal()
        self.mdl = SolutionDomainModel(self.case)
        self.out2 = OutputControlModel(self.case2)

        self.proc = QProcess()
        self.connect(self.proc, SIGNAL('readyReadStandardOutput()'), self.__readFromStdout)
        self.connect(self.proc, SIGNAL('readyReadStandardError()'), self.__readFromStderr)
        self.procErrorFlag = False

        self.cwd = os.getcwd()

        self.exec_dir = os.path.join(self.case['resu_path'], 'check_mesh')
        if os.path.isdir(self.exec_dir):
            shutil.rmtree(self.exec_dir)

        os.mkdir(self.exec_dir)
        os.chdir(self.exec_dir)

        self.fmt = OutputControlModel(self.case).getWriterFormat("-1").lower()
        self.out2.setWriterLabel("-1", "quality")
        self.out2.setWriterFormat("-1", self.fmt)

        # Prepare preprocessing

        mesh_input = self.mdl.getMeshInput()

        if mesh_input:

            if not os.path.isabs(mesh_input):
                mesh_input = os.path.join(self.case['case_path'], mesh_input)
            try:
                os.symlink(mesh_input, 'mesh_input')
            except AttributeError:
                shutil.copy2(mesh_input, 'mesh_input')

            self.__csProcess()

        else:

            self.preprocess_cmd = []
            nodeList = self.mdl.node_meshes.xmlGetNodeList('mesh', 'name')

            if len(nodeList) > 1:
                os.mkdir('mesh_input')

            for meshNode in nodeList:

                cmd = []

                mesh   = meshNode['name']
                format = meshNode['format']
                path   = meshNode['path']
                if path != None:
                    mesh = os.path.join(path, mesh)
                if not os.path.isabs(mesh) and self.case['mesh_path'] != None:
                    mesh = os.path.join(self.case['mesh_path'], mesh)
                if meshNode['num']:
                    cmd = cmd + ['--num', meshNode['num']]
                if meshNode['reorient'] == 'on':
                    cmd.append('--reorient')
                if meshNode['grp_fac']:
                    cmd = cmd + ['--grp-fac', meshNode['grp_fac']]
                if meshNode['grp_cel']:
                    cmd = cmd + ['--grp-cel', meshNode['grp_cel']]

                # Define postprocessing output for errors and warnings.

                if self.fmt in ('med', 'cgns'):
                    cmd = cmd + ['--post-error', self.fmt]
                else:
                    cmd = cmd + ['--post-error', 'ensight']

                str_add = '_%02d' % (len(self.preprocess_cmd)+1)

                cmd = cmd + ['--case', 'preprocess'+str_add]
                if len(nodeList) > 1:
                    cmd = cmd + ['--out',
                                 os.path.join('mesh_input', 'mesh'+str_add)]

                cmd.append(mesh)

                log.debug("ecs_cmd = %s" % str(cmd))

                self.preprocess_cmd.append(cmd)

            self.__preProcess()

        self.case.undoStartGlobal()


    def __preProcess(self):


        # Modify the PATH for relocatable installation

        import sys

        if self.case['package'].config.features['relocatable'] == "yes":
            if sys.platform.startswith("win"):
                env = QProcessEnvironment.systemEnvironment()
                saved_path = env.value('PATH')
                env.insert("PATH",
                           self.case['package'].get_dir('bindir') + ";" + \
                           env.value("PATH"))
                self.proc.setProcessEnvironment(env)

        # Prepare command

        args = self.preprocess_cmd.pop(0)
        nodelist = self.mdl.node_meshes.xmlGetNodeList('mesh', 'name')

        if len(self.preprocess_cmd) < len(nodelist):
            self.disconnect(self.proc,
                            SIGNAL('finished(int, QProcess::ExitStatus)'),
                            self.__preProcess)

        self.proc.start(str(self.case['package'].get_preprocessor()),
                        [str(s) for s in args])

        # Run Preprocessor

        next_task = None

        if len(self.preprocess_cmd) > 0:
            self.connect(self.proc,
                         SIGNAL('finished(int, QProcess::ExitStatus)'),
                         self.__preProcess)
        else:
            self.connect(self.proc,
                         SIGNAL('finished(int, QProcess::ExitStatus)'),
                         self.__ecsPostTreatment)

        # Reset the PATH to its previous value

        if self.case['package'].config.features['relocatable'] == "yes":
            if sys.platform.startswith("win"):
                env.insert("PATH", saved_path)
                self.proc.setProcessEnvironment(env)


    def __ecsPostTreatment(self):
        if self.proc.exitStatus() == QProcess.NormalExit and not self.procErrorFlag:
            mesh_input = os.path.join(self.exec_dir, 'mesh_input')
            if os.path.isfile(mesh_input) or os.path.isdir(mesh_input):
                self.__csProcess()
            else:
                self.__finished()

        else:
            self.__finished()


    def __csProcess(self):

        # Modify the PATH for relocatable installation

        import sys

        if self.case['package'].config.features['relocatable'] == "yes":
            if sys.platform.startswith("win"):
                env = QProcessEnvironment.systemEnvironment()
                saved_path = env.value('PATH')
                env.insert("PATH",
                           self.case['package'].get_dir('bindir') + ";" + \
                           env.value("PATH"))
                self.proc.setProcessEnvironment(env)

        # Run Kernel
        self.disconnect(self.proc,
                        SIGNAL('finished(int, QProcess::ExitStatus)'),
                        self.__ecsPostTreatment)

        self.case2.xmlSaveDocument()
        args = ['--quality', '--log', '0', '--param', self.case2['xmlfile']]

        self.proc.start(str(self.case['package'].get_solver()),
                        [str(s) for s in args])

        self.connect(self.proc,
                     SIGNAL('finished(int, QProcess::ExitStatus)'),
                     self.__csPostTreatment)

        # Reset the PATH to its previous value

        if self.case['package'].config.features['relocatable'] == "yes":
            if sys.platform.startswith("win"):
                env.insert("PATH", saved_path)
                self.proc.setProcessEnvironment(env)


    def __csPostTreatment(self):

        # Cleanup
        mesh_input = os.path.join(self.exec_dir, 'mesh_input')
        if os.path.isdir(mesh_input) and not os.path.islink(mesh_input):
            shutil.rmtree(mesh_input)
        else:
            os.remove(mesh_input)
        if os.path.isfile('mesh_output'):
            os.remove('mesh_output')
        os.remove('cs_cmd')

        self.__saveLog()
        self.__finished()


    def __readFromStdout(self):
        """
        Private slot to handle the readyReadStandardOutput signal of the process.
        """
        if self.proc is None:
            return
        self.proc.setReadChannel(QProcess.StandardOutput)

        while self.proc and self.proc.canReadLine():
            ba = self.proc.readLine()
            if ba.isNull(): return
            s = (ba.data()).decode("utf-8")[:-1]
            self.logText.append(s)


    def __readFromStderr(self):
        """
        Private slot to handle the readyReadStandardError signal of the process.
        """
        if self.proc is None:
            return
        self.proc.setReadChannel(QProcess.StandardError)

        while self.proc and self.proc.canReadLine():
            ba = self.proc.readLine()
            if ba.isNull(): return
            s = (ba.data()).decode("utf-8")[:-1]
            self.logText.append('<font color="red">' + s + '</font>')
            self.procErrorFlag = True


    def __saveLog(self):
        if Py2:
            logFile = open(os.path.join(self.exec_dir, 'check_mesh.log'), 'w')
        else:
            logFile = open(os.path.join(self.exec_dir, 'check_mesh.log'), 'w', encoding='utf-8')
        logFile.write(str(self.logText.toPlainText().toAscii()))
        logFile.close()


    def __finished(self):
        os.chdir(self.cwd)
        self.pushButton.setEnabled(True)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class SolutionVerifView(QWidget, Ui_SolutionVerifForm):
    """
    """
    def __init__(self, parent, case):
        """
        Constructor
        """
        QWidget.__init__(self, parent)

        Ui_SolutionVerifForm.__init__(self)
        self.setupUi(self)

        self.parent = parent
        self.case = case
        self.case.undoStopGlobal()
        self.mdl = SolutionDomainModel(self.case)
        self.out = OutputControlModel(self.case)

        self.case2 = Case(package = self.case['package'], file_name = None)
        XMLinit(self.case2).initialize()
        self.case2['xmlfile'] = 'cs_cmd'
        self.case2['salome'] = self.case['salome']

        self.node_domain  = self.case.xmlGetNode('solution_domain')
        faces_cutting = self.node_domain.xmlGetNode('faces_cutting')
        joining = self.node_domain.xmlGetNode('joining')
        periodicity = self.node_domain.xmlGetNode('periodicity')

        sd_node = self.case2.xmlGetNode('solution_domain')
        if faces_cutting != None:
            if (faces_cutting)['status'] == 'on':
                sd_node.xmlInitNode('faces_cutting',
                                    status='on').xmlChildsCopy(faces_cutting)
        if joining != None:
            sd_node.xmlInitNode('joining').xmlChildsCopy(joining)
        if periodicity != None:
            sd_node.xmlInitNode('periodicity').xmlChildsCopy(periodicity)

        self.out2 = OutputControlModel(self.case2)

        # combo models
        self.modelFMTCHR         = ComboModel(self.comboBoxFMTCHR, 3, 1)
        self.modelFormat         = ComboModel(self.comboBoxFormat, 2, 1)
        self.modelPolygon        = ComboModel(self.comboBoxPolygon, 3, 1)
        self.modelPolyhedra      = ComboModel(self.comboBoxPolyhedra, 3, 1)

        self.modelFMTCHR.addItem(self.tr("EnSight Gold"), 'ensight')
        self.modelFMTCHR.addItem(self.tr("MED"), 'med')
        self.modelFMTCHR.addItem(self.tr("CGNS"), 'cgns')
        self.modelFMTCHR.addItem(self.tr("Catalyst"), 'catalyst')
        self.modelFMTCHR.addItem(self.tr("CCM-IO"), 'ccm')

        import cs_config
        cfg = cs_config.config()
        if cfg.libs['med'].have == "no":
            self.comboBoxFMTCHR.setItemData(1, QColor(Qt.red), Qt.TextColorRole);
        if cfg.libs['cgns'].have == "no":
            self.comboBoxFMTCHR.setItemData(2, QColor(Qt.red), Qt.TextColorRole);
        if cfg.libs['catalyst'].have == "no":
            self.comboBoxFMTCHR.setItemData(3, QColor(Qt.red), Qt.TextColorRole);
        if cfg.libs['ccm'].have == "no":
            self.comboBoxFMTCHR.setItemData(4, QColor(Qt.red), Qt.TextColorRole);

        self.modelFormat.addItem(self.tr("binary"), 'binary')
        self.modelFormat.addItem(self.tr("text"), 'text')

        self.modelPolygon.addItem(self.tr("display"), 'display')
        self.modelPolygon.addItem(self.tr("discard"), 'discard_polygons')
        self.modelPolygon.addItem(self.tr("subdivide"), 'divide_polygons')

        self.modelPolyhedra.addItem(self.tr("display"), 'display')
        self.modelPolyhedra.addItem(self.tr("discard"), 'discard_polyhedra')
        self.modelPolyhedra.addItem(self.tr("subdivide"), 'divide_polyhedra')

        # connections

        self.connect(self.comboBoxFMTCHR, SIGNAL("activated(const QString&)"), self.slotOutputFormat)
        self.connect(self.comboBoxFormat, SIGNAL("activated(const QString&)"), self.slotOutputOptions)
        self.connect(self.comboBoxPolygon, SIGNAL("activated(const QString&)"), self.slotOutputOptions)
        self.connect(self.comboBoxPolyhedra, SIGNAL("activated(const QString&)"), self.slotOutputOptions)
        self.connect(self.checkBoxBigEndian, SIGNAL("clicked()"), self.slotOutputOptions)
        self.connect(self.toolButtonBatch, SIGNAL("clicked()"), self.slotMeshChecking)

        # INITIALISATIONS

        # 1 - Values of post processing's format

        fmt = self.out.getWriterFormat("-1")
        self.modelFMTCHR.setItem(str_model=fmt)
        line = self.out.getWriterOptions("-1")
        self.__updateOptionsFormat(line)

        if not (self.mdl.getMeshList() or self.mdl.getMeshInput()):
            self.toolButtonBatch.setEnabled(False)

        self.case.undoStartGlobal()


    @pyqtSignature("const QString &")
    def slotOutputFormat(self, text):
        """
        Input format of post-processing
        """
        format = self.modelFMTCHR.dicoV2M[str(text)]

        if self.out.getWriterFormat("-1") != format:
            self.out.setWriterFormat("-1",format)
            l = self.out.defaultWriterValues()['options']
            self.out.setWriterOptions("-1",l)

        if self.out2.getWriterFormat("-1") != format:
            self.out2.setWriterFormat("-1",format)
            l = self.out2.defaultWriterValues()['options']
            self.out2.setWriterOptions("-1",l)
            self.__updateOptionsFormat(l)


    @pyqtSignature("")
    def slotOutputOptions(self):
        """
        Create format's command line options
        """
        line = []
        opt_format = self.modelFormat.dicoV2M[str(self.comboBoxFormat.currentText())]
        line.append(opt_format)

        if self.checkBoxBigEndian.isChecked():
            line.append('big_endian')

        opt_polygon = self.modelPolygon.dicoV2M[str(self.comboBoxPolygon.currentText())]
        opt_polyhed = self.modelPolyhedra.dicoV2M[str(self.comboBoxPolyhedra.currentText())]
        if opt_polygon != 'display': line.append(opt_polygon)
        if opt_polyhed != 'display': line.append(opt_polyhed)

        l = string.join(line, ',')
        log.debug("slotOutputOptions-> OPTCHR = %s" % l)
        self.out.setWriterOptions("-1",l)
        self.out2.setWriterOptions("-1",l)


    def __updateOptionsFormat(self, line):
        """
        Update command-line options at each modification of
        post processing format
        """
        lst = line.split(',')
        format = self.modelFMTCHR.dicoV2M[str(self.comboBoxFMTCHR.currentText())]
        log.debug("__updateOptionsFormat-> FMTCHR = %s" % format)
        log.debug("__updateOptionsFormat-> OPTCHR = %s" % line)

        # update widgets from the options list

        for opt in lst:

            if opt == 'binary' or opt == 'text' :
                self.modelFormat.setItem(str_model=opt)

            if opt == 'discard_polygons' or opt == 'divide_polygons':
                self.modelPolygon.setItem(str_model=opt)

            if opt == 'discard_polyhedra' or opt == 'divide_polyhedra':
                self.modelPolyhedra.setItem(str_model=opt)

            if format == 'ensight':
                if opt == 'big_endian':
                    self.checkBoxBigEndian.setChecked(True)

        if 'discard_polygons' not in lst and 'divide_polygons' not in lst:
            self.modelPolygon.setItem(str_model="display")
        if 'discard_polyhedra' not in lst and 'divide_polyhedra' not in lst:
            self.modelPolyhedra.setItem(str_model="display")
        if 'big_endian' not in lst:
            self.checkBoxBigEndian.setChecked(False)

        # enable and disable options related to the format

        self.modelPolygon.enableItem(str_model='discard_polygons')
        self.modelPolygon.enableItem(str_model='divide_polygons')
        self.modelPolyhedra.enableItem(str_model='discard_polyhedra')
        self.modelPolyhedra.enableItem(str_model='divide_polyhedra')
        self.comboBoxPolygon.setEnabled(True)
        self.comboBoxPolyhedra.setEnabled(True)

        if format != "ensight":
            if format == "cgns":
                self.modelPolyhedra.setItem(str_model='divide_polyhedra')
                self.modelPolyhedra.disableItem(str_model='display')
            elif format in ["catalyst", "ccm"]:
                self.modelPolyhedra.setItem(str_model='display')
                self.modelPolygon.setItem(str_model='display')
                self.comboBoxPolygon.setEnabled(False)
                self.comboBoxPolyhedra.setEnabled(False)
            self.modelFormat.setItem(str_model="binary")
            self.modelFormat.disableItem(str_model='text')
            self.labelBigEndian.setEnabled(False)
            self.checkBoxBigEndian.setEnabled(False)
        else:
            self.modelFormat.enableItem(str_model='text')
            self.comboBoxFormat.setEnabled(True)
            self.labelBigEndian.setEnabled(True)
            self.checkBoxBigEndian.setEnabled(True)


    def __setButtonEnabled(self):
        """
        Block the QButton during the display of the dialog.
        """
        try:
            self.toolButtonBatch.setEnabled(not self.toolButtonBatch.isEnabled())
        except:
            pass


    def slotMeshChecking(self):
        """
        """
        self.__setButtonEnabled()
        dialog = MeshQualityCriteriaLogDialogView(self.parent, self.case, self.case2)
        dialog.show()
        self.connect(dialog, SIGNAL("accepted()"), self.__setButtonEnabled)
        self.connect(dialog, SIGNAL("rejected()"), self.__setButtonEnabled)


    def tr(self, text):
        """
        Translation
        """
        return text

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
