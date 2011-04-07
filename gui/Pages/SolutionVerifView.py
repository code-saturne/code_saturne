# -*- coding: utf-8 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2010 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     The Code_Saturne User Interface is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     The Code_Saturne User Interface is distributed in the hope that it will be
#     useful, but WITHOUT ANY WARRANTY; without even the implied warranty
#     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with the Code_Saturne Kernel; if not, write to the
#     Free Software Foundation, Inc.,
#     51 Franklin St, Fifth Floor,
#     Boston, MA  02110-1301  USA
#
#-------------------------------------------------------------------------------

"""
This module contains the following class:
- SolutionVerifView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import os, logging, subprocess
import string, shutil, cStringIO

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from Base.Toolbox import GuiParam
from Base.QtPage import ComboModel
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
        self.cs = self.case['package'].get_solver()
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

        self.fmt = string.split(self.__getPostCommand())[0]

        # Prepare preprocessing

        mesh_input = None
        node = self.mdl.node_meshes.xmlGetNode('mesh_input', 'path')
        if node:
            mesh_input = node['path']
            if mesh_input:
                if not os.path.isabs(mesh_input):
                    mesh_input = os.path.join(self.case['case_path'], mesh_input)

        if mesh_input:
            try:
                os.symlink(mesh_input, 'mesh_input')
            except AttributeError:
                shutil.copy2(mesh_input, 'mesh_input')

        else:

            self.preprocess_cmd = []
            nodeList = self.mdl.node_meshes.xmlGetNodeList('mesh', 'name')

            if len(nodeList) > 1:
                os.mkdir('mesh_input')

            for meshNode in nodeList:

                cmd = self.case['package'].get_preprocessor()

                name   = meshNode['name']
                format = meshNode['format']
                mesh = self.case['mesh_path'] + '/' + name
                cmd += ' --mesh ' + mesh
                if meshNode['num']:
                    cmd += ' --num ' + meshNode['num']
                if meshNode['reorient'] == 'on':
                    cmd += ' --reorient'
                if meshNode['grp_fac']:
                    cmd += ' --grp-fac ' + meshNode['grp_fac']
                if meshNode['grp_cel']:
                    cmd += ' --grp-cel ' + meshNode['grp_cel']

                cmd += ' ' + self.__getPostCommand()

                # Limit postprocessing output to errors and info.

                cmd += ' --info'

                cmd += ' --case preprocess'
                if len(nodeList) > 1:
                    str_add = '_%02d' % (len(self.preprocess_cmd) + 1)
                    cmd += str_add
                    cmd += ' --out ' + os.path.join('mesh_input', 'mesh' + str_add)
                else:
                    cmd += ' --out mesh_input'
                log.debug("ecs_cmd = %s" % cmd)

                self.preprocess_cmd.append(cmd)

            self.__preProcess()


    def __preProcess(self):

        cmd = self.preprocess_cmd.pop(0)
        nodelist = self.mdl.node_meshes.xmlGetNodeList('mesh', 'name')

        if len(self.preprocess_cmd) < len(nodelist):
            self.disconnect(self.proc,
                            SIGNAL('finished(int, QProcess::ExitStatus)'),
                            self.__preProcess)

        # Run Preprocessor

        self.proc.start(cmd)

        next_task = None

        if len(self.preprocess_cmd) > 0:
            self.connect(self.proc,
                         SIGNAL('finished(int, QProcess::ExitStatus)'),
                         self.__preProcess)
        else:
            self.connect(self.proc,
                         SIGNAL('finished(int, QProcess::ExitStatus)'),
                         self.__ecsPostTreatment)


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
        # Run Kernel
        self.disconnect(self.proc,
                        SIGNAL('finished(int, QProcess::ExitStatus)'),
                        self.__ecsPostTreatment)

        self.case2.xmlSaveDocument()
        args = ' --quality --log 0 --param ' + self.case2['xmlfile']
        self.proc.start(self.cs + args)

        self.connect(self.proc,
                     SIGNAL('finished(int, QProcess::ExitStatus)'),
                     self.__csPostTreatment)


    def __csPostTreatment(self):
        if self.proc.exitStatus() == QProcess.NormalExit and not self.procErrorFlag:

            try:

                if self.fmt == "--ensight":

                    os.rename(os.path.join(self.exec_dir, 'chr.ensight'),
                              os.path.join(self.exec_dir, 'quality.ensight'))

                    os.chdir(os.path.join(self.exec_dir, 'quality.ensight'))

                    for src in os.listdir(os.getcwd()):
                        if src[:4] == "chr.":
                            dst = src.replace("chr.", "quality.")
                            os.rename(src, dst)

                    os.rename('CHR.case', 'QUALITY.case')

                    out = cStringIO.StringIO()
                    f = open('QUALITY.case')
                    for line in f:
                        out.write(line.replace('chr', 'quality'))
                    f.close()
                    out2 = open('QUALITY.case', 'w')
                    out2.write(out.getvalue())
                    out2.close()

                    os.chdir(self.exec_dir)

                elif self.fmt == "--med":
                    os.rename('chr.med', 'QUALITY.med')

                elif self.fmt == "--cgns":
                    os.rename('chr.cgns', 'QUALITY.cgns')

            except OSError: # file to rename might not exist
                pass

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


    def __getPostCommand(self):
        """
        Return the preprocessor argument for postprocessing.
        """
        format = self.out2.getPostProFormat()
        if format == "EnSight":
            l = " --ensight"
        elif format == "MED":
            l = " --med"
        elif format == "CGNS":
            l = " --cgns"
        return l


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
            str = QString()
            s = QString(str.fromUtf8(ba.data()))[:-1]
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
            str = QString()
            s = QString(str.fromUtf8(ba.data()))[:-1]
            self.logText.append(s.prepend('<font color="red">').append('</font>'))
            self.procErrorFlag = True


    def __saveLog(self):
        logFile = open(os.path.join(self.exec_dir, 'check_mesh.log'), 'w')
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
        self.mdl = SolutionDomainModel(self.case)
        self.out = OutputControlModel(self.case)

        self.case2 = Case(package = self.case['package'], file_name = None)
        XMLinit(self.case2)
        self.case2['xmlfile'] = 'cs_cmd'
        self.case2['salome'] = self.case['salome']

        faces_cutting = self.case.xmlGetNode('faces_cutting')
        joining = self.case.xmlGetNode('joining')
        periodicity = self.case.xmlGetNode('periodicity')

        sd_node = self.case2.xmlGetNode('solution_domain')
        if faces_cutting != None:
            if (faces_cutting)['status'] == 'on':
                sd_node.xmlInitNode('faces_cutting',
                                    status='on').xmlChildsCopy(faces_cutting)
        if joining != None:
            sd_node.xmlInitNode('joining').xmlChildsCopy(joining)
        if periodicity != None:
            sd_node.xmlInitNode('solution_domain').xmlChildsCopy(periodicity)

        self.out2 = OutputControlModel(self.case2)

        # combo models
        self.modelFMTCHR         = ComboModel(self.comboBoxFMTCHR, 3, 1)
        self.modelFormat         = ComboModel(self.comboBoxFormat, 2, 1)
        self.modelPolygon        = ComboModel(self.comboBoxPolygon, 3, 1)
        self.modelPolyhedra      = ComboModel(self.comboBoxPolyhedra, 3, 1)

        self.modelFMTCHR.addItem(self.tr("EnSight Gold"), 'EnSight')
        self.modelFMTCHR.addItem(self.tr("MED"), 'MED')
        self.modelFMTCHR.addItem(self.tr("CGNS"), 'CGNS')

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

        fmt = self.out.getPostProFormat()
        self.modelFMTCHR.setItem(str_model=fmt)
        line = self.out.getPostProOptionsFormat()
        self.__updateOptionsFormat(line)

        if not self.mdl.getMeshList():
            self.toolButtonBatch.setEnabled(False)


    @pyqtSignature("const QString &")
    def slotOutputFormat(self, text):
        """
        Input format of post-processing
        """
        format = self.modelFMTCHR.dicoV2M[str(text)]

        if self.out.getPostProFormat() != format:
            self.out.setPostProFormat(format)
            l = self.out.defaultInitialValues()['postprocessing_options']
            self.out.setPostProOptionsFormat(l)

        if self.out2.getPostProFormat() != format:
            self.out2.setPostProFormat(format)
            l = self.out2.defaultInitialValues()['postprocessing_options']
            self.out2.setPostProOptionsFormat(l)
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
        self.out.setPostProOptionsFormat(l)
        self.out2.setPostProOptionsFormat(l)


    def __updateOptionsFormat(self, line):
        """
        Update command-line options at each modification of
        post processing format
        """
        list = string.split(line, ',')
        format = self.modelFMTCHR.dicoV2M[str(self.comboBoxFMTCHR.currentText())]
        log.debug("__updateOptionsFormat-> FMTCHR = %s" % format)
        log.debug("__updateOptionsFormat-> OPTCHR = %s" % line)

        # update widgets from the options list

        for opt in list:

            if opt == 'binary' or opt == 'text' :
                self.modelFormat.setItem(str_model=opt)

            if opt == 'discard_polygons' or opt == 'divide_polygons':
                self.modelPolygon.setItem(str_model=opt)

            if opt == 'discard_polyhedra' or opt == 'divide_polyhedra':
                self.modelPolyhedra.setItem(str_model=opt)

            if format == 'EnSight':
                if opt == 'big_endian':
                    self.checkBoxBigEndian.setChecked(True)

        if 'discard_polygons' not in list and 'divide_polygons' not in list:
            self.modelPolygon.setItem(str_model="display")
        if 'discard_polyhedra' not in list and 'divide_polyhedra' not in list:
            self.modelPolyhedra.setItem(str_model="display")
        if 'big_endian' not in list:
            self.checkBoxBigEndian.setChecked(False)

        # enable and disable options related to the format

        if format != "EnSight":

            if format == "CGNS":
                self.modelPolyhedra.setItem(str_model='divide_polyhedra')
                self.modelPolyhedra.disableItem(str_model='display')

            self.modelFormat.setItem(str_model="binary")
            self.modelFormat.disableItem(str_model='text')
            self.labelBigEndian.setEnabled(False)
            self.checkBoxBigEndian.setEnabled(False)
        else:
            self.modelFormat.enableItem(str_model='text')
            self.comboBoxFormat.setEnabled(True)
            self.labelBigEndian.setEnabled(True)
            self.checkBoxBigEndian.setEnabled(True)
            self.modelPolyhedra.enableItem(str_model='display')
            self.comboBoxPolyhedra.setEnabled(True)


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
