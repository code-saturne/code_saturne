# -*- coding: iso-8859-1 -*-
#
#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne User Interface, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2009 EDF S.A., France
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
This module contains the following classes and function:
- MatisseCustomView
"""

#-------------------------------------------------------------------------------
# Standard modules
#-------------------------------------------------------------------------------

import string
import logging

#-------------------------------------------------------------------------------
# Third-party modules
#-------------------------------------------------------------------------------

from PyQt4.QtCore import *
from PyQt4.QtGui  import *

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from MatisseCustomForm import Ui_MatisseCustomForm
import Base.QtPage as QtPage

import Pages.MatisseTypeModel as MatisseType
import Pages.MatisseGeomModel as MatisseGeom

from Pages.MatisseRangeDescriptionModel import MatisseRangeDescriptionModel
from Pages.MatisseNetworkModel import MatisseNetworkModel
from Pages.MatisseThermicModel import MatisseThermicModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("MatisseCustomView")
log.setLevel(logging.DEBUG)

#-------------------------------------------------------------------------------
# StandarItemModel class
#-------------------------------------------------------------------------------

_EPSILON = 1e-6

class StandardItemModelMatisseCustom(QStandardItemModel):
    """
    QStandardItemModel associated with the QtableView for common Matisse operations.

    Behavior depends on the value of attribute tagName (self.tagName).
    """

    def __init__(self, parent, case, tagName, dico):
        """
        """
        QStandardItemModel.__init__(self)

        self.parent = parent
        self.case = case
        self.tagName = tagName
        self.dico = dico
        self.dataMatisse = []
        self.setColumnCount(len(self.dico['headers']))
        self.__initXMLModel()
        self.__initData()


    def __initXMLModel(self):
        """
        """
        if self.tagName in ["inlet_range_line", "inlet_range_height",
                            "outlet_range_line", "outlet_range_height",
                            "network_line", "network_row"]:

            self.default_row = ["default", 0., 0., 0., 0.]
            self.columns_editable = [0, 1, 2]
            self.column_step_min = 3
            self.column_step_max = 4

            # XML Model
            if self.tagName in ["network_line", "network_row"]:

                self.model = MatisseNetworkModel(self.case)

                model_geom = MatisseGeom.MatisseGeomModel(self.case)
                self.lineStep = model_geom.getMatisseGeomDoubleVar('ptrres')
                self.rowStep = model_geom.getMatisseGeomDoubleVar('plgres')
                self.heightStep = model_geom.getMatisseGeomDoubleVar('epchel')
                self.lineMax = model_geom.getMatisseGeomDoubleVar('nptran')
                self.rowMax = model_geom.getMatisseGeomDoubleVar('nplgrs')
                self.heightMax = model_geom.getMatisseGeomDoubleVar('nchest')

                model_mat_type = MatisseType.MatisseTypeModel(self.case)
                self.alveoStat = model_mat_type.node_alveo['status']

            elif self.tagName in ["inlet_range_line", "inlet_range_height"]:

                self.model = MatisseRangeDescriptionModel(self.case, 'inlet_range')

                model_geom = MatisseGeom.MatisseGeomModel(self.case)
                self.lineStep = model_geom.getMatisseGeomDoubleVar('ptrres')
                self.heightStep = model_geom.getMatisseGeomDoubleVar('epchel')
                self.lineMax = model_geom.getMatisseGeomDoubleVar('nptran')
                self.heightMax = model_geom.getMatisseGeomDoubleVar('nchest')

            elif self.tagName in ["outlet_range_line", "outlet_range_height"]:

                self.model = MatisseRangeDescriptionModel(self.case, 'outlet_range')

                model_geom = MatisseGeom.MatisseGeomModel(self.case)
                self.lineStep = model_geom.getMatisseGeomDoubleVar('ptrres')
                self.heightStep = model_geom.getMatisseGeomDoubleVar('epchel')
                self.lineMax = model_geom.getMatisseGeomDoubleVar('nptran')
                self.heightMax = model_geom.getMatisseGeomDoubleVar('nchest')

        else:
        #elif self.tagName in ["thermal_line", "thermal_row", "thermal_height"]:

            # Thermic load
            self.areatype = self.dico['areatype']

            self.default_row = ["default", 0., 0., 0., 0., 0.]
            self.columns_editable = [0, 1, 2, 3]
            self.column_step_min = 4
            self.column_step_max = 5

            # XML Model
            self.model = MatisseThermicModel(self.case)

            self.model_mat_type = MatisseType.MatisseTypeModel(self.case)
            model_geom = MatisseGeom.MatisseGeomModel(self.case)
            self.lineStep = model_geom.getMatisseGeomDoubleVar('ptrres')
            self.rowStep = model_geom.getMatisseGeomDoubleVar('plgres')
            self.heightStep = model_geom.getMatisseGeomDoubleVar('epchel')

            self.lineMax = model_geom.getMatisseGeomDoubleVar('nptran')
            self.rowMax = model_geom.getMatisseGeomDoubleVar('nplgrs')
            self.heightMax = model_geom.getMatisseGeomDoubleVar('nchest')

        #
        self.areatype = self.dico['areatype']
        if self.areatype == 'line' :
            self.step = self.lineStep
        elif self.areatype == 'height' :
            self.step = self.heightStep
        elif self.areatype == 'row' :
            self.step = self.rowStep


    def __initData(self):
        """
        Load previous values
        """
        if self.tagName in ["inlet_range_line", "inlet_range_height",
                            "outlet_range_line", "outlet_range_height",
                            "network_line", "network_row"]:

            llabel, lbmin, lbmax  = self.model.GetAreas(self.areatype)

            n = len(llabel)
            if ((n != len(lbmin)) or (n != len(lbmax))):
                print "XML format error : bad definition of <area> "
                sys.exit(1)

            for area in range(0,len(llabel)):
                row = self.default_row
                row[0] = llabel[area]
                row[1] = float(lbmin[area])
                row[2] = float(lbmax[area])
                row[self.column_step_min] = self.step * row[1]
                row[self.column_step_max] = self.step * row[2]
                self.dataMatisse.append(row)
                nrows = self.rowCount()
                self.setRowCount(nrows+1)

        elif self.tagName in ["thermal_line", "thermal_row", "thermal_height"]:

            llabel, lbmin, lbmax, lval  = self.model.GetAreas(self.areatype)

            n = len(llabel)
            if ((n != len(lbmin)) or (n != len(lbmax)) or (n != len(lval))):
                print "XML format error : bad definition of <area> "
                sys.exit(1)

            # Default values
            if len(llabel) == 0 :
                dlabel, dmin, dmax, dval = self.model.DefaultArea(self.areatype)
                self.model.NewArea(self.areatype, dlabel, dmin, dmax, dval)
                llabel.append(dlabel)
                lbmin.append(dmin)
                lbmax.append(dmax)
                lval.append(dval)

            for area in range(0,len(llabel)):
                row = self.default_row
                row[0] = llabel[area]
                row[1] = float(lbmin[area])
                row[2] = float(lbmax[area])
                row[3] = float(lval[area])
                row[self.column_step_min] = self.step * row[1]
                row[self.column_step_max] = self.step * row[2]
                self.dataMatisse.append(row)
                nrows = self.rowCount()
                self.setRowCount(nrows+1)


    def data(self, index, role):
        if not index.isValid():
            return QVariant()
        if role == Qt.DisplayRole:
            row = index.row()
            col = index.column()
            return QVariant(self.dataMatisse[row][col])
        return QVariant()


    def flags(self, index):
        if not index.isValid():
            return Qt.ItemIsEnabled
        else:
            return Qt.ItemIsEnabled | Qt.ItemIsSelectable


    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return QVariant(self.dico['headers'][section])
        return QVariant()


    def setData(self, index, value, role):
        self.emit(SIGNAL("dataChanged(const QModelIndex &, const QModelIndex &)"), index, index)
        return True


    def addRow(self, label, sbmin, sbmax, svalue):
        """
        Add a row in the model.
        """
        log.debug("addRow")
        bool_create = False
        if self.tagName in ["inlet_range_line", "inlet_range_height",
                            "outlet_range_line", "outlet_range_height",
                            "network_line", "network_row"]:

            bool_create = self.__addRowRangeDescriptionNetwork(label, sbmin, sbmax, self.areatype, self.step)

        if self.tagName in ["thermal_line", "thermal_row", "thermal_height"]:

            bool_create = self.__addRowThermal(label, sbmin, sbmax, svalue, self.areatype, self.step)


        if bool_create:
            row = self.default_row
            row[0] = label
            row[1] = sbmin
            row[2] = sbmax
            row[self.column_step_min] = self.step * row[1]
            row[self.column_step_max] = self.step * row[2]
            self.dataMatisse.append(row)
            nrows = self.rowCount()
            self.setRowCount(nrows+1)


    def editRow(self, row, new_label, new_bmin, new_bmax, new_value):
        """
        Edit a row in the model
        """
        if self.tagName in ["inlet_range_line", "inlet_range_height",
                            "outlet_range_line", "outlet_range_height",
                            "network_line", "network_row"]:

            [label, sbmin, sbmax, sbmin2, sbmax2] = self.dataMatisse[row]
            bool_modify = self.__editRangeDescription(new_label, new_bmin, new_bmax,
                                                      label, sbmin, sbmax,
                                                      self.areatype, self.step)

        if self.tagName in ["thermal_line", "thermal_row", "thermal_height"]:
            [label, sbmin, sbmax, svalue, sbmin2, sbmax2] = self.dataMatisse[row]
            bool_modify = self.__editRangeThermal(new_label, new_bmin, new_bmax, new_value,
                                                  label, sbmin, sbmax, svalue,
                                                  self.areatype, self.step)


        if bool_modify:
            pass # TODO

    def deleteRow(self, row):
        """
        Delete a row in the model.
        """
        log.debug("deleteRow")

        if self.tagName in ["inlet_range_line", "inlet_range_height",
                            "outlet_range_line", "outlet_range_height",
                            "network_line", "network_row"]:

            [label, sbmin, sbmax, sbmin2, sbmax2] = self.dataMatisse[row]
            bool_cancel = self.__deleteRowRangeDescription(label, sbmin, sbmax, self.areatype)

        if self.tagName in ["thermal_line", "thermal_row", "thermal_height"]:
            [label, sbmin, sbmax, svalue, sbmin2, sbmax2] = self.dataMatisse[row]
            bool_cancel = self.__deleteRowThermal(label, sbmin, sbmax, svalue, self.areatype)



        if bool_cancel:
            del self.dataMatisse[row]
            nrows = self.rowCount()
            self.setRowCount(nrows-1)


    # RangeDescription & Network
    # --------------------------

    def __addRowRangeDescriptionNetwork(self, label, sbmin, sbmax, areatype, step):
        """
        Add a row : test for RangeDescription & Network
        """
        log.debug("__addRowRangeDescriptionNetwork")

        if not label:
            label = "Default"

        try :
            bmin = float(sbmin)
            bmax = float(sbmax)
        except:
            title = self.tr("Warning")
            msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
            QMessageBox.warning(self.parent, title, msg)
            return

        llabel, lbmin, lbmax  = self.model.GetAreas(areatype)
        create = True
        #
        # Bounds check 1/2
        if (bmin > bmax) or (bmin < 0):
            create = False
        if (areatype == 'line'):
            if (bmax > self.lineMax) :
                create = False
        elif (areatype == 'height'):
            if (bmax > self.heightMax) :
                create = False
        elif (areatype == 'row'):
            if (bmax > self.rowMax) :
                create = False

        for area in range(0,len(llabel)) :
            if ((label == llabel[area]) and
                (abs(bmin - float(lbmin[area])) < _EPSILON) and
                (abs(bmax - float(lbmax[area])) < _EPSILON)):

                create = False

            #
            # Bounds check 2/2
            if (((bmin > float(lbmin[area])) and (bmin < float(lbmax[area]))) or
                ((bmax > float(lbmin[area])) and (bmax < float(lbmax[area])))) :
                create = False

        if create :
            #name = self.insertHlist(h, label, sbmin, sbmax, step)
            self.model.NewArea(areatype, label, bmin, bmax)
        else :
            title = self.tr("Warning")
            msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
            QMessageBox.warning(self.parent, title, msg)

        return create


    def __editRowRangeDescriptionNetwork(self, new_label, new_bmin, new_bmax, label, sbmin, sbmax, areatype):
        """
        Edit row : test for RangeDescription & Network
        """

        if not new_label:
            new_label = "Default"

        if areatype == 'line' :
            step = self.lineStep
        elif areatype == 'height' :
            step = self.heightStep
        elif areatype == 'row' :
            step = self.rowStep

        #if len(self.currentEntry) == 1 : # TODO ONE SELECTION
        if 1 :

##             entry = self.currentEntry
##             if h == self.h :
##                 label, sbmin, sbmax = self.select.areaInfo(entry)
##             elif h == self.h3 :
##                 label, sbmin, sbmax = self.select3.areaInfo(entry)

            bmin = float(sbmin)
            bmax = float(sbmax)

            llabel, lbmin, lbmax  = self.model.GetAreas(areatype)

            #
            # Bounds check
            modify = True

            #
            # Type check
            try :
                fnew_bmin=float(new_bmin)
                fnew_bmax=float(new_bmax)
            except :

                title = self.tr("Warning")
                msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
                QMessageBox.warning(self.parent, title, msg)
                return

            if (fnew_bmin > fnew_bmax) or (fnew_bmin < 0):
                modify = False
            if (areatype == 'line'):
                if fnew_bmax > self.lineMax :
                    modify = False
            elif (areatype == 'height'):
                if (fnew_bmax > self.heightMax) :
                    modify = False

            if not modify :
                title = self.tr("Warning")
                msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
                QMessageBox.warning(self.parent, title, msg)

            else:
                for area in range(0,len(llabel)) :
                    if ((label == llabel[area]) and
                        (abs(bmin - float(lbmin[area])) <= _EPSILON) and
                        (abs(bmax - float(lbmax[area])) <= _EPSILON)):


                        if ((label != new_label) or
                            (abs(bmax - fnew_bmax) >= _EPSILON) or
                            (abs(bmin - fnew_bmin) >= _EPSILON)) :

                            for area1 in range(0,len(llabel)) :
                                if area1 != area :
                                    if (((fnew_bmin > float(lbmin[area1])) and (fnew_bmin < float(lbmax[area1]))) or
                                        ((fnew_bmax > float(lbmin[area1])) and (fnew_bmax < float(lbmax[area1])))) :
                                        modify = False

                            if modify :
                                self.model.SetArea(areatype, area, new_label, new_bmin, new_bmax)
                                #self.replaceHlist(h, entry, new_label, new_bmin, new_bmax, step)
                            else:
                                title = self.tr("Warning")
                                msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
                                QMessageBox.warning(self.parent, title, msg)

        else: # # TODO MULTI-SELECTION
            for entry in self.currentEntry:

##                 if h == self.h :
##                     label, sbmin, sbmax = self.select.areaInfo(entry)
##                 elif h == self.h3 :
##                     label, sbmin, sbmax = self.select3.areaInfo(entry)

                bmin = float(sbmin)
                bmax = float(sbmax)

                llabel, lbmin, lbmax  = self.model.GetAreas(areatype)

                for area in range(0,len(llabel)) :
                    if ((label == llabel[area]) or
                        (abs(bmin - float(lbmin[area])) <= _EPSILON) and
                        (abs(bmax - float(lbmax[area])) <= _EPSILON)):

                        if (label != new_label) :

                            if (new_label != _MULTISEL) :
                                self.model.SetArea(areatype, area, new_label, None , None)
                                #self.replaceHlist( h, entry, new_label, bmin, bmax,step)

                            else :
                                pass


    def __deleteRowRangeDescriptionNetwork(self, label, sbmin, sbmax, areatype):
        """
        Delete row : test for RangeDescription & Network
        """
        log.debug("__deleteRowRangeDescriptionNetwork")

        bmin = float(sbmin)
        bmax = float(sbmax)

        llabel, lbmin, lbmax  = self.model.GetAreas(areatype)

        cancel = False
        for area in range(0,len(llabel)) :
            if ((label == llabel[area]) and
                (abs(bmin - float(lbmin[area])) < _EPSILON) and
                (abs(bmax - float(lbmax[area])) < _EPSILON)):

                self.model.EraseArea(areatype, area)
                cancel = True

        return cancel

##             # Delete boundary region from the Hlist
##             #
##             if not cancel :
##                 print "Entry is in Hlist but not in XML"
##                 sys.exit(1)
##             else :
## #
## #    Pb : si on ne supprime pas le dernier elt : lors de la creation d'une nouvelle hlist
## #         il recupere un nom deja donne !!!
## ##                if h == self.h :
## ##                    self.entriesNumberH1 = self.entriesNumberH1 -1
## ##                elif  h == self.h3 :
## ##                    self.entriesNumberH3 = self.entriesNumberH3 -1
##                 h.hlist.delete_entry(entry)



    # Thermal
    # -------

    def __addRowRangeThermal(self, label, sbmin, sbmax, svalue, areatype, step):
        """
        Add a row : test for Thermal load
        """
        log.debug("__addRowThermal")

        # Reading  and add to Hlist
        #
        if not label:
            label = "Default"

        try :
            bmin = float(sbmin)
            bmax = float(sbmax)
            value = float(svalue)
        except:
            title = self.tr("Warning")
            msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
            QMessageBox.warning(self.parent, title, msg)
            return

        llabel, lbmin, lbmax, lval  = self.model.GetAreas(areatype)

        create = True
        #
        # Bounds check
        if (bmin > bmax) or (bmin < 0):
            create = False
        if (areatype == 'line'):
            if (bmax > self.lineMax) :
                create = False
        elif (areatype == 'row'):
            if (bmax > self.rowMax) :
                create = False
        elif (areatype == 'height'):
            if (bmax > self.heightMax) :
                create = False

        for area in range(0,len(llabel)) :
            #
            # No duplicate
            if ((label == llabel[area]) and
                (abs(bmin - float(lbmin[area])) < _EPSILON) and
                (abs(bmax - float(lbmax[area])) < _EPSILON) and
                (abs(value - float(lval[area])) < _EPSILON)) :
                create = False

            if (((bmin > float(lbmin[area])) and (bmin < float(lbmax[area]))) or
                ((bmax > float(lbmin[area])) and (bmax < float(lbmax[area])))) :
                create = False


        if create :
            #name = self.insertHlist(h, label, sbmin, sbmax, svalue ,step)
            self.model.NewArea(areatype, label, bmin, bmax, value)
        else :
            title = self.tr("Warning")
            msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
            QMessageBox.warning(self.parent, title, msg)

        return create


    def __editRowThermal(self, new_label, new_bmin, new_bmax, label, sbmin, sbmax, areatype):
        """
        Edit row : test for Thermal load
        """
        log.debug("__editRowThermal")


        if not new_label:
            new_label = "Default"

        if areatype == 'line' :
            step = self.lineStep
        elif areatype == 'row' :
            step = self.rowStep
        elif areatype == 'height' :
            step = self.heightStep

        if 1:
        #if len(self.currentEntry) == 1 :

            bmin = float(sbmin)
            bmax = float(sbmax)
            value = float(svalue)

            llabel, lbmin, lbmax, lval  = self.model.GetAreas(areatype)

            #
            # Bounds check
            modify = True

            #
            # Type check
            try :
                fnew_bmin=float(new_bmin)
                fnew_bmax=float(new_bmax)
                fnew_value=float(new_value)
            except :
                title = self.tr("Warning")
                msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
                QMessageBox.warning(self.parent, title, msg)
                return

            if (fnew_bmin > fnew_bmax) or (fnew_bmin < 0):
                modify = False
            if (areatype == 'line'):
                if fnew_bmax > self.lineMax :
                    modify = False
            elif (areatype == 'height'):
                if (fnew_bmax > self.heightMax) :
                    modify = False
            elif (areatype == 'row'):
                if (fnew_bmax > self.rowMax) :
                    modify = False

            if not modify :
                title = self.tr("Warning")
                msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
                QMessageBox.warning(self.parent, title, msg)
            else:
                for area in range(0,len(llabel)) :
                    if ((label == llabel[area]) and
                        (abs(bmin - float(lbmin[area])) <= _EPSILON) and
                        (abs(bmax - float(lbmax[area])) <= _EPSILON) and
                        (abs(value - float(lval[area])) <= _EPSILON)) :

                        if ((label != new_label) or
                            (abs(bmax - fnew_bmax) >= _EPSILON) or
                            (abs(bmin - fnew_bmin) >= _EPSILON) or
                            (abs(value - fnew_value) >= _EPSILON)) :

                            for area1 in range(0,len(llabel)) :
                                if area1 != area :
                                    if (((fnew_bmin > float(lbmin[area1])) and (fnew_bmin < float(lbmax[area1]))) or
                                        ((fnew_bmax > float(lbmin[area1])) and (fnew_bmax < float(lbmax[area1])))) :
                                        modify = False

                            if modify :
                                self.model.SetArea(areatype, area, new_label, new_bmin, new_bmax, new_value)
                                #self.replaceHlist(h, entry, new_label, new_bmin, new_bmax, new_value, step)
                            else:
                                title = self.tr("Warning")
                                msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
                                QMessageBox.warning(self.parent, title, msg)


        else :
            for entry in self.currentEntry:

                bmin = float(sbmin)
                bmax = float(sbmax)
                value = float(svalue)

                llabel, lbmin, lbmax, lval  = self.model.GetAreas(areatype)

                #
                # Type check
                try :
                    fnew_value=float(new_value)
                except :
                    if new_value != _MULTISEL :
                        title = self.tr("Warning")
                        msg   = self.tr("'%1' bad definition: see the bounds definition").arg(label)
                        QMessageBox.warning(self.parent, title, msg)

                        return

                for area in range(0,len(llabel)) :
                    if ((label == llabel[area]) or
                        (abs(bmin - float(lbmin[area])) <= _EPSILON) and
                        (abs(bmax - float(lbmax[area])) <= _EPSILON) and
                        (abs(value - float(lval[area])) <= _EPSILON)) :

                        if ((label != new_label) or new_value == _MULTISEL or
                            (abs(value - fnew_value) >= _EPSILON)) :

                            if (new_label == _MULTISEL) and (new_value != _MULTISEL):
                                self.model.SetArea(areatype, area, None, None , None , new_value)
                                #self.replaceHlist( h, entry, label, bmin, bmax, new_value, step)

                            elif (new_value == _MULTISEL) and (new_label != _MULTISEL) :
                                self.model.SetArea(areatype, area, new_label , None , None , None)
                                #self.replaceHlist( h, entry, new_label, bmin, bmax, value, step)

                            elif (new_value != _MULTISEL) and (new_label != _MULTISEL) :
                                self.model.SetArea(areatype, area, new_label, None , None , new_value)
                                #self.replaceHlist( h, entry, new_label, bmin, bmax, new_value, step)

                            else :
                                pass


    def __deleteRowThermal(self, label, sbmin, sbmax, svalue, areatype):
        """
        Delete row : test for Thermal
        """
        log.debug("__deleteRowThermal")

        llabel, lbmin, lbmax, lval  = self.model.GetAreas(areatype)
        cancel = False
        for area in range(0,len(llabel)) :
            if ((label == llabel[area]) and
                (abs(bmin - float(lbmin[area])) < _EPSILON) and
                (abs(bmax - float(lbmax[area])) < _EPSILON) and
                (abs(value - float(lval[area])) < _EPSILON)) :
                self.model.EraseArea(areatype, area)
                cancel = True

        return cancel

#-------------------------------------------------------------------------------
# Main class
#-------------------------------------------------------------------------------

class MatisseCustomView(QWidget, Ui_MatisseCustomForm):
    """
    """

    def __init__(self, *args):
        """
        Constructor
        """
        QWidget.__init__(self, *args)
        Ui_MatisseCustomForm.__init__(self)
        self.setupUi(self)


    def initWidget(self, case, tagName):
        """
        Method to initialize the widget.
        Must be explicitly called.
        """
        self.case = case
        self.tagName = tagName
        self._createWidgets()


    def getLabel(self):
        """
        Return the name of the boundary condition. It 's not allowed to have
        blank or not ASCII caracters in this name.
        """
        label = str(self.lineEditLabel.text())
        return  string.join(string.split(label), '_')


    def getBmin(self):
        """
        """
        bmin, ok = self.lineEditMin.text().toFloat()
        return bmin


    def getBmax(self):
        """
        """
        bmax, ok = self.lineEditMax.text().toFloat()
        return bmax


    def getValue(self):
        """
        """
        value, ok = self.lineEditValue.text().toFloat()
        return value


    def createItem(self):
        label = self.getLabel()
        sbmin = self.getBmin()
        sbmax = self.getBmax()
        svalue = self.getValue()
        self.modelMatisseCustom.addRow(label, sbmin, sbmax, svalue)


    def modifyItem(self):
        new_label = self.getLabel()
        new_bmin = self.getBmin()
        new_bmax = self.getBmax()
        new_value = self.getValue()
        nrows = 0
        tab_rows = []
        selectionModel = self.tableView.selectionModel()
        for index in selectionModel.selectedRows():
            tab_rows.append(index.row())
            nrows = nrows + 1
        log.debug("modifyItem nrows = %i tab_rows = %s "%(nrows, str(tab_rows)))
        for row in tab_rows:
            self.modelMatisseCustom.editRow(row, new_label, new_bmin, new_bmax, new_value)


    def cancelItem(self):
        nrows = 0
        tab_rows = []
        selectionModel = self.tableView.selectionModel()
        for index in selectionModel.selectedRows():
            tab_rows.append(index.row())
            nrows = nrows + 1
        log.debug("cancelItem nrows = %i tab_rows = %s "%(nrows, str(tab_rows)))
        for row in tab_rows:
            self.modelMatisseCustom.deleteRow(row)


    def _createWidgets(self):
        """
        Create the Page layout.
        """
        if self.tagName == "":
            log.debug("_createWidgets self.tagName is not set!")
            sys.exit(1)

        self.__initDico()
        dico = self.dico_custom[self.tagName]

        # Update title and labels in the view
        self.groupBox.setTitle(dico['title'])
        self.labelLabel.setText(QString(dico['label']))
        self.labelMin.setText(QString(dico['min']))
        self.labelMax.setText(QString(dico['max']))
        self.labelValue.hide()
        self.lineEditValue.hide()

        if 'value' in dico.keys():
            self.labelValue.setText(QString(dico['value']))
            self.labelValue.show()
            self.lineEditValue.show()

        # Set model for tableView
        self.modelMatisseCustom = StandardItemModelMatisseCustom(self, self.case, self.tagName, dico)
        self.tableView.setModel(self.modelMatisseCustom)
        self.tableView.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.tableView.setSelectionMode(QAbstractItemView.ExtendedSelection)

        # Connections
        self.connect(self.pushButtonNew,    SIGNAL("clicked()"), self.createItem)
        self.connect(self.pushButtonModify, SIGNAL("clicked()"), self.modifyItem)
        self.connect(self.pushButtonDelete, SIGNAL("clicked()"), self.cancelItem)

        # validators
        regExp = QRegExp("[_A-Za-z0-9]*") # QRegExp("^all[ ]*$|^[0-9\ ]*$")
        validatorLabel = QRegExpValidator(regExp, self.lineEditLabel)
        self.lineEditLabel.setValidator(validatorLabel)

        validatorFloat = QtPage.DoubleValidator(self.lineEditMin)
        self.lineEditMin.setValidator(validatorFloat)
        self.lineEditMax.setValidator(validatorFloat)
        self.lineEditValue.setValidator(validatorFloat)


    def tr(self, text):
        """
        Translation
        """
        return text


    def __initDico(self):
        """
        Set a dictionnary which contains field names.
        """

        line   = { 'areatype' : "line",
                   'title' : self.tr("headloss in line") ,
                   'label' : self.tr("Label"),
                   'min' : self.tr("Distance min \n  (in lines)"),
                   'max' : self.tr("Distance max \n  (in lines)"),
                   'headers' : [self.tr("Label"),
                                self.tr("Distance min \n  (in lines)"),
                                self.tr("Distance max \n  (in lines)"),
                                self.tr("Distance min \n      (m)"),
                                self.tr("Distance max \n      (m)")]
                   }

        height = { 'areatype' : "height",
                   'title' : self.tr("headloss in height") ,
                   'label' : self.tr("Label"),
                   'min' : self.tr("Height min \n  (in cells)"),
                   'max' : self.tr("Height max \n  (in cells)"),
                   'headers' : [self.tr("Label"),
                                self.tr("Height min \n  (in cells)"),
                                self.tr("Height max \n  (in cells)"),
                                self.tr("Height min \n     (m)"),
                                self.tr("Height max \n     (m)")]
                   }

        row    = { 'areatype' : "row",
                   'title' : self.tr("Headloss in row") ,
                   'label' : self.tr("Label"),
                   'min' : self.tr("Distance min \n (in rows)"),
                   'max' : self.tr("Distance max \n (in rows)"),
                   'headers' : [self.tr("Label"),
                                self.tr("Distance min \n (in rows)"),
                                self.tr("Distance max \n (in rows)"),
                                self.tr("Distance min \n      (m)"),
                                self.tr("Distance max \n      (m)")]
                   }

        dico_custom = {}
        dico_custom['inlet_range_line']    = line
        dico_custom['inlet_range_height']  = height
        dico_custom['outlet_range_line']   = line
        dico_custom['outlet_range_height'] = height
        dico_custom['network_line']        = line
        dico_custom['network_row']         = row
        dico_custom['thermal_line']        = line
        dico_custom['thermal_height']      = height
        dico_custom['thermal_row']         = row

        self.dico_custom = dico_custom

#-------------------------------------------------------------------------------
# Testing part
#-------------------------------------------------------------------------------

if __name__ == "__main__":
    pass

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------