# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of code_saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2024 EDF S.A.
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
This module defines the Immersed boundaries model data management.

This module contains the following classes and function:
- ImmersedBoundariesModel
"""

#-------------------------------------------------------------------------------
# Library modules
#-------------------------------------------------------------------------------

import sys, unittest

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Variables, Model
from code_saturne.model.XMLmodel import ModelTest
from code_saturne.model.NotebookModel import NotebookModel

#-------------------------------------------------------------------------------
# Cathare coupling model class
#-------------------------------------------------------------------------------

class ImmersedBoundariesModel(Variables, Model):
    """
    Manage the input/output markups in the xml doc
    """

    # ----------------------------------
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case

        self.__node_models = self.case.xmlGetNode('thermophysical_models')
        self.__node_ibm    = self.__node_models.xmlInitNode('immersed_boundaries')
        self.notebook = NotebookModel(self.case)
    # ----------------------------------

    # ----------------------------------
    def defaultValues(self):
        """
        Return in a dictionnary which contains default values.
        """
        default = {}
        default['OnOff']                = 'off'
        default['ibm_dimension']        = "3D computation"
        default['object_name']          = "object"
        default['object_method']        = "defined by function"
        default['object_init']      = "0.0"
        default['file_name']        = ""

        #Volumic zone
        default['object_is_fsi']         = "off"
        default['object_moving']         = "fixed"
        default['object_conjugate_heat_transfer']  = "off"
        default['object_initialization']  = "off"
        default['object_physical_properties']  = "off"
        default['object_thermal_source_term']  = "off"

        #FSI
        default['block_rotationX']   = "off"
        default['block_rotationY']   = "off"
        default['block_rotationZ']   = "off"
        default['block_displacementX']  = "off"
        default['block_displacementY']  = "off"
        default['block_displacementZ']  = "off"
        default['object_nb_cycle_min'] = 3
        default['object_nb_cycle_max'] = 12
        default['object_conv_criteria'] = 0.01

        #physical properties
        default['object_modeling_mode'] = "density" #or "mass"
        default['object_density'] = 0.0
        default['object_mass']    = 0.0

        default['object_stiffness'] = 0.0
        default['object_damping']   = 0.0
        default['object_specific_heat'] = 4000.0
        default['object_thermal_conductivity'] = 1e-5

        default['object_rho_property_mode'] = "constant" #user_law
        default['object_mass_property_mode'] = "constant"
        default['object_dam_property_mode'] = "constant"
        default['object_sti_property_mode'] = "constant"
        default['object_cp_property_mode']  = "constant"
        default['object_al_property_mode']  = "constant"

        #Boundary condition
        default['object_boundary_condition_nature'] = 'Wall'
         #off (no thermal), temperature, temperature_formula, flux, flux_formula
        default['object_boundary_energy_mode'] = 'off'
        default['object_boundary_velocity_mode'] = 'wall_law'
        default['object_boundary_heat_flux'] = 0.0
        default['object_boundary_temperature'] = 0.0

        return default
    # ----------------------------------

    # ----------------------------------
    def getNumberOfObjects(self):

        return len(self.__node_ibm.xmlGetNodeList('ibm_object'))
    # ----------------------------------

    # ----------------------------------
    def getNumberOfSTLObjectSeedPoints(self, num_obj):
        self.isLowerOrEqual(num_obj, self.getNumberOfObjects())
        node_obj = self.__node_ibm.xmlGetNodeList('ibm_object')[num_obj-1]

        return len(node_obj.xmlGetNodeList('STL_exterior_point'))
    # ----------------------------------

    # ----------------------------------
    def getObjectsNodeList(self):

        return self.__node_ibm.xmlGetNodeList('ibm_object')
    # ----------------------------------

    # ----------------------------------
    def getObjectsNameList(self):
        """
        return list of objects' names
        """
        objectsList = []
        num = 0
        for node in self.__node_ibm.xmlGetNodeList('ibm_object'):
            objectsList.append(self.getObjectName(num+1))
            num = num +1

        return objectsList
    # ----------------------------------

    # ----------------------------------

    def getObjectNumFromName(self, object_name):
        """
        return id of object from his label
        """
        num = 1
        for node in self.__node_ibm.xmlGetNodeList('ibm_object'):
            name = self.getObjectName(num)
            if (name == object_name):
                return num
            else:
                num = num + 1

        raise ValueError("Unknown immersed object_name = {}".format(object_name))


    # ----------------------------------

    # ----------------------------------

    def setOnOff(self, state):

        self.__node_ibm.xmlSetData('ibm_state', state)


    def getOnOff(self):

        state = self.__node_ibm.xmlGetString('ibm_state')
        if state is None or state == '':
            state = self.defaultValues()['OnOff']

        return state
    # ----------------------------------

    def setIBMDim(self, dimension):
        self.isStr(dimension)
        self.__node_ibm.xmlSetData('ibm_dimension', dimension)


    def getIBMDim(self):

        dimension = self.__node_ibm.xmlGetString('ibm_dimension')
        if dimension is None or dimension == '':
            dimension = self.defaultValues()['ibm_dimension']

        return dimension

    # ----------------------------------

    #------------------------------------------------------------------
    # Helper function
    #------------------------------------------------------------------

    def __getStringData(self, index, name, setFunction, tag=None):
        """
        Get string value from xml file.
        """
        self.isLowerOrEqual(index+1, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[index]
        value = node.xmlGetString(name)
        return self.__getDefaultDataIfNone(index, value, name, setFunction, tag)

    def __getIntData(self, index, name, setFunction):
        """
        Get int value from xml file.
        """
        self.isLowerOrEqual(index+1, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[index]
        value = node.xmlGetInt(name)
        return self.__getDefaultDataIfNone(index, value, name, setFunction)


    def __getDefaultDataIfNone(self, index, value, name, setFunction, tag=None):
        """
        Get default value if value is none.
        """
        if value is None or value == "":
            value = self.defaultValues()[name]
            if tag is not None:
                setFunction(index+1, value, tag)
            else:
                setFunction(index+1, value)
        return value

    # ----------------------------------
    @Variables.undoGlobal
    def addObject(self, name, method, is_fsi, moving_type):

        num = self.getNumberOfObjects()

        node_new = self.__node_ibm.xmlAddChild('ibm_object')

        num += 1

        self.setObjectName(num, name)
        self.setObjectMethod(num, method)
        self.setObjectFSI(num, is_fsi)
        self.setObjectMoving(num, moving_type)

        return num


    @Variables.undoLocal
    def deleteObject(self, num):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node_list = self.__node_ibm.xmlGetNodeList('ibm_object')
        node = node_list[num-1]
        node.xmlRemoveNode()
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectName(self, num, name):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(name)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_name', name)

    @Variables.noUndo
    def getObjectName(self, num):

        return self.__getStringData(num-1, 'object_name',
                                    self.setObjectName)
    # ----------------------------------

    # ----------------------------------

    @Variables.undoLocal
    def setObjectMethod(self, num, method):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(method)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_method', method)

    @Variables.noUndo
    def getObjectMethod(self, num):

        return self.__getStringData(num-1, 'object_method',
                                    self.setObjectMethod)

    # ----------------------------------

    # ----------------------------------

    @Variables.undoLocal
    def setObjectFSI(self, num, is_fsi):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(is_fsi)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_is_fsi', is_fsi)

    @Variables.noUndo
    def getObjectFSI(self, num):

        return self.__getStringData(num-1, 'object_is_fsi',
                                    self.setObjectFSI)
    # ----------------------------------

    # ----------------------------------

    @Variables.undoLocal
    def setObjectMoving(self, num, moving):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(moving)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_moving', moving)


    @Variables.noUndo
    def getObjectMoving(self, num):

        return self.__getStringData(num-1, 'object_moving',
                                    self.setObjectMoving)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectDensity(self, num, density):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_density', density)


    @Variables.noUndo
    def getObjectDensity(self, num):

        return self.__getStringData(num-1, 'object_density',
                                    self.setObjectDensity)
    # ----------------------------------

    # ----------------------------------

    @Variables.undoLocal
    def setObjectMass(self, num, mass):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_mass', mass)


    @Variables.noUndo
    def getObjectMass(self, num):

        return self.__getStringData(num-1, 'object_mass',
                                    self.setObjectMass)
    # ----------------------------------

    # ----------------------------------

    @Variables.undoLocal
    def setObjectStiffness(self, num, stiffness):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_stiffness', stiffness)


    @Variables.noUndo
    def getObjectStiffness(self, num):

        return self.__getStringData(num-1, 'object_stiffness',
                                    self.setObjectStiffness)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectDamping(self, num, damping):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_damping', damping)


    @Variables.noUndo
    def getObjectDamping(self, num):

        return self.__getStringData(num-1, 'object_damping',
                                    self.setObjectDamping)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectSpecificHeat(self, num, cp):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_specific_heat', cp)


    @Variables.noUndo
    def getObjectSpecificHeat(self, num):

        return self.__getStringData(num-1, 'object_specific_heat',
                                    self.setObjectSpecificHeat)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectThermalConductivity(self, num, al):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_thermal_conductivity', al)


    @Variables.noUndo
    def getObjectThermalConductivity(self, num):

        return self.__getStringData(num-1, 'object_thermal_conductivity',
                                    self.setObjectThermalConductivity)
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectPropertyMode(self, num, mode, tag):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]

        if (tag == 'density'):
            node.xmlSetData('object_rho_property_mode', mode)
        elif (tag == 'stiffness'):
            node.xmlSetData('object_sti_property_mode', mode)
        elif (tag == 'damping'):
            node.xmlSetData('object_dam_property_mode', mode)
        elif (tag == 'specific_heat'):
            node.xmlSetData('object_cp_property_mode', mode)
        elif (tag == 'thermal_conductivity'):
            node.xmlSetData('object_al_property_mode', mode)
        elif (tag == 'mass'):
            node.xmlSetData('object_mass_property_mode', mode)
        else:
            raise ValueError("Unknown immersed tag = {} for solid object property.".format(tag))


    @Variables.noUndo
    def getObjectPropertyMode(self, num, tag):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]

        if (tag == 'density'):
            return self.__getStringData(num-1, 'object_rho_property_mode',
                                        self.setObjectPropertyMode, tag)
        elif (tag == 'stiffness'):
            return self.__getStringData(num-1, 'object_sti_property_mode',
                                        self.setObjectPropertyMode, tag)
        elif (tag == 'damping'):
            return self.__getStringData(num-1, 'object_dam_property_mode',
                                        self.setObjectPropertyMode, tag)
        elif (tag == 'specific_heat'):
            return self.__getStringData(num-1, 'object_cp_property_mode',
                                        self.setObjectPropertyMode, tag)
        elif (tag == 'thermal_conductivity'):
            return self.__getStringData(num-1, 'object_al_property_mode',
                                        self.setObjectPropertyMode, tag)
        elif (tag == 'mass'):
            return self.__getStringData(num-1, 'object_mass_property_mode',
                                        self.setObjectPropertyMode, tag)
        else:
            raise ValueError("Unknown immersed tag = {} for solid object property.".format(tag))

    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectEqPosition(self, num, xeq=None, yeq=None, zeq=None):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_eq = node.xmlInitNode('equilibrium_position')
        if xeq:
            node_eq.xmlSetData('xeq',xeq)
        if yeq:
            node_eq.xmlSetData('yeq',yeq)
        if zeq:
            node_eq.xmlSetData('zeq',zeq)


    @Variables.noUndo
    def getObjectEqPosition(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_eq = node.xmlGetChildNode('equilibrium_position')

        if node_eq is None:
            self.setObjectEqPosition(num,
                                     self.defaultValues()['object_init'],
                                     self.defaultValues()['object_init'],
                                     self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']

        else:
            return node_eq.xmlGetString('xeq'), \
                   node_eq.xmlGetString('yeq'), \
                   node_eq.xmlGetString('zeq')
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectInitPosition(self, num, xini=None, yini=None, zini=None):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_ini = node.xmlInitNode('initial_position')
        if xini:
            node_ini.xmlSetData('xini', xini)
        if yini:
            node_ini.xmlSetData('yini', yini)
        if zini:
            node_ini.xmlSetData('zini', zini)


    @Variables.noUndo
    def getObjectInitPosition(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_ini = node.xmlGetChildNode('initial_position')

        if node_ini is None:
            self.setObjectInitPosition(num,
                                       self.defaultValues()['object_init'],
                                       self.defaultValues()['object_init'],
                                       self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']
        else:
            return node_ini.xmlGetString('xini'), \
                   node_ini.xmlGetString('yini'), \
                   node_ini.xmlGetString('zini')
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectInitVel(self, num, vx=None, vy=None, vz=None):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_vel = node.xmlInitNode('initial_velocity')
        if vx:
            node_vel.xmlSetData('vx', vx)
        if vy:
            node_vel.xmlSetData('vy', vy)
        if vz:
            node_vel.xmlSetData('vz', vz)


    @Variables.noUndo
    def getObjectInitVel(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_vel = node.xmlGetChildNode('initial_velocity')

        if node_vel is None:
            self.setObjectInitVel(num,
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']
        else:
            return node_vel.xmlGetString('vx'), \
                   node_vel.xmlGetString('vy'), \
                   node_vel.xmlGetString('vz')
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectAngularVel(self, num, wx=None, wy=None, wz=None):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_ang_vel = node.xmlInitNode('angular_velocity')
        if wx:
            node_ang_vel.xmlSetData('wx', wx)
        if wy:
            node_ang_vel.xmlSetData('wy', wy)
        if wz:
            node_ang_vel.xmlSetData('wz', wz)


    @Variables.noUndo
    def getObjectAngularVel(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_ang_vel = node.xmlGetChildNode('angular_velocity')

        if node_ang_vel is None:
            self.setObjectAngularVel(num,
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']
        else:
            return node_ang_vel.xmlGetString('wx'), \
                   node_ang_vel.xmlGetString('wy'), \
                   node_ang_vel.xmlGetString('wz')
    # ----------------------------------

    # ----------------------------------
    @Variables.undoLocal
    def setObjectInitAcc(self, num, ax=None, ay=None, az=None):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_acc = node.xmlInitNode('initial_acceleration')
        if ax:
            node_acc.xmlSetData('ax', ax)
        if ay:
            node_acc.xmlSetData('ay', ay)
        if az:
            node_acc.xmlSetData('az', az)


    @Variables.noUndo
    def getObjectInitAcc(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_vel = node.xmlGetChildNode('initial_acceleration')

        if node_vel is None:
            self.setObjectInitAcc(num,
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'],
                                  self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']

        else:
            return node_vel.xmlGetString('ax'), \
                   node_vel.xmlGetString('ay'), \
                   node_vel.xmlGetString('az')
    # ----------------------------------

    # ----------------------------------

    @Variables.undoLocal
    def setObjectInitialAngles(self, num, theta_x=None, theta_y=None, theta_z=None):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_init_angles = node.xmlInitNode('initial_angles')
        if theta_x:
            node_init_angles.xmlSetData('theta_x', theta_x)
        if theta_y:
            node_init_angles.xmlSetData('theta_y', theta_y)
        if theta_z:
            node_init_angles.xmlSetData('theta_z', theta_z)


    @Variables.noUndo
    def getObjectInitialAngles(self, num):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_init_angles = node.xmlGetChildNode('initial_angles')

        if node_init_angles is None:
            self.setObjectInitialAngles(num,
                                        self.defaultValues()['object_init'],
                                        self.defaultValues()['object_init'],
                                        self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']
        else:
            return node_init_angles.xmlGetString('theta_x'), \
                   node_init_angles.xmlGetString('theta_y'), \
                   node_init_angles.xmlGetString('theta_z')
    # ----------------------------------

    # ----------------------------------

    def deleteObjectFormula(self, objId):

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        node.xmlRemoveChild('explicit_formula')


    def setObjectFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('explicit_formula')

        n.xmlSetTextNode(formula)


    def getObjectFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('explicit_formula')

        return formula

    def setObjectImposedMovingFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('imposed_moving_formula')

        n.xmlSetTextNode(formula)


    def getObjectImposedMovingFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('imposed_moving_formula')

        return formula

    def setObjectTemperatureFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('init_temperature_formula')

        n.xmlSetTextNode(formula)


    def getObjectTemperatureFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('init_temperature_formula')

        return formula


    def setObjectRhoFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('rho_formula')

        n.xmlSetTextNode(formula)


    def getObjectRhoFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('rho_formula')

        return formula


    def setObjectCpFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('cp_formula')

        n.xmlSetTextNode(formula)


    def getObjectCpFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('cp_formula')

        return formula


    def setObjectAlFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('al_formula')

        n.xmlSetTextNode(formula)


    def getObjectAlFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('al_formula')

        return formula


    def setObjectInertiaFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('inertia_formula')

        n.xmlSetTextNode(formula)


    def getObjectInertiaFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('inertia_formula')

        return formula


    def setObjectBoundaryTemperatureFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('b_temp_formula')

        n.xmlSetTextNode(formula)


    def getObjectBoundaryTemperatureFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('b_temp_formula')

        return formula


    def setObjectBoundaryThermalFluxFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('b_flux_formula')

        n.xmlSetTextNode(formula)


    def getObjectBoundaryThermalFluxFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('b_flux_formula')

        return formula


    def setObjectThermalSourceTermFormula(self, objId, formula):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        n = node.xmlInitChildNode('thermal_st_formula')

        n.xmlSetTextNode(formula)


    def getObjectThermalSourceTermFormula(self, objId):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[objId]
        formula = node.xmlGetString('thermal_st_formula')

        return formula


    def getIBMFormulaComponents(self, objId):

        exp = self.getObjectFormula(objId)
        if not exp:
            exp = """indicator = 1;\n"""

        req = [('indicator', 'Solid object indicator (1 is solid, 0 is fluid)')]
        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('time', 'current time'),
               ('cog_fsi_x', '1st component of the center of gravity (only for FSI objects)'),
               ('cog_fsi_y', '2nt component of the center of gravity (only for FSI objects)'),
               ('cog_fsi_z', '3rd component of the center of gravity (only for FSI objects)'),
               ('rot_fsi', '3x3 rotation matrix with rot[i][j], \
               i,j = 0,1,2 indexing (only for FSI objects)')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym


    def getImposedMovingFormulaComponents(self, objId):

        exp = self.getObjectImposedMovingFormula(objId)
        if not exp:
            exp = (
                "porous_velocity[0] = 0.;\n"
                "porous_velocity[1] = 0.;\n"
                "porous_velocity[2] = 0.;"
            )

        req = [('porous_velocity[0]', '1st component of the porous velocity'),
               ('porous_velocity[1]', '2nd component of the porous velocity'),
               ('porous_velocity[2]', '3rd component of the porous velocity')]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('time', 'current time'),
               ('temperature', 'temperature of the immersed object')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym
    # ----------------------------------


    def getFormulaTemperature(self, objId):

        exp = self.getObjectTemperatureFormula(objId)
        if not exp:
            exp = """temperature = 273.15;\n"""

        req = [('temperature', 'temperature of the immersed object')]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('time', 'current time')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym
    # ----------------------------------


    def getFormulaRho(self, objId):

        exp = self.getObjectRhoFormula(objId)
        if not exp:
            exp = (
                "rho = 1.8;"
            )

        req = [('rho', 'Density of the immersed object')]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('time', 'current time'),
               ('temperature', 'temperature of the immersed object')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym
    # ----------------------------------


    def getFormulaCp(self, objId):

        exp = self.getObjectCpFormula(objId)
        if not exp:
            exp = (
                "cp = 4000.;"
            )

        req = [('cp', "Specific heat of the immersed object")]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('time', 'current time'),
               ('temperature', 'temperature of the immersed object')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym
    # ----------------------------------


    def getFormulaAl(self, objId):

        exp = self.getObjectAlFormula(objId)
        if not exp:
            exp = (
                "lambda = 1.e-5;"
            )

        req = [('lambda', 'Thermal conductivity of the immersed object')]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('time', 'current time'),
               ('temperature', 'temperature of the immersed object')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym
    # ----------------------------------



    def getFormulaThermalSourceTerm(self, objId):

        exp = self.getObjectThermalSourceTermFormula(objId)
        if not exp:
            exp = (
                "#Explicit part (W/m^3)\n"
                "st_exp = 0.;\n\n"
                "#Implicit part (W/m^3/K)\n"
                "st_imp = 0.;"
            )

        req = [('st_exp', 'Explicit thermal source term (W/m^3)'),
               ('st_imp', 'Implicit thermal source term (W/m^3/K)')]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('time', 'current time')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym
    # ----------------------------------



    def getFormulaInertia(self, objId):

        exp = self.getObjectInertiaFormula(objId)
        if not exp:
            exp = (
                "i11 = 0.0;\n"
                "i22 = 0.0;\n"
                "i33 = 0.0;\n"
                "i12 = 0.0;\n"
                "i13 = 0.0;\n"
                "i23 = 0.0;\n"
            )

        req = [('i11', 'inertia matrix of the immersed object (1,1)'),
               ('i22', 'inertia matrix of the immersed object (2,2)'),
               ('i33', 'inertia matrix of the immersed object (3,3)'),
               ('i12', 'inertia matrix of the immersed object (1,2)'),
               ('i13', 'inertia matrix of the immersed object (1,3)'),
               ('i23', 'inertia matrix of the immersed object (2,3)')]

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('time', 'current time')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym
    # ----------------------------------


    def getFormulaBoundaryEnergy(self, objId, energy_mode):

        exp = ''
        if energy_mode == "flux_formula":
            exp = self.getObjectBoundaryThermalFluxFormula(objId)
            req = [('flux', 'Heat flux [W/m2]')]
            if not exp:
                exp = (
                    "flux = 1000.;"
                )
        elif energy_mode == "temperature_formula":
            exp = self.getObjectBoundaryTemperatureFormula(objId)
            req = [('temperature', 'Temperature [K]')]

            if not exp:
                exp = (
                    "temperature = 293.15;"
                )

        sym = [('x', 'cell center coordinate'),
               ('y', 'cell center coordinate'),
               ('z', 'cell center coordinate'),
               ('time', 'current time')]

        for (name, val) in self.notebook.getNotebookList():
            sym.append((name, 'value (notebook) = ' + str(val)))

        return exp, req, sym
    # ----------------------------------


    @Variables.undoLocal
    def setMesh(self, num, file_name):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(file_name)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('file_name', file_name)

    def getMesh(self, num):

        return self.__getStringData(num-1, 'file_name',
                                    self.setMesh)

    def setObjectBlockRX(self, num, val):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(val)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('block_rotationX', val)

    @Variables.noUndo
    def getObjectBlockRX(self, num):

        return self.__getStringData(num-1, 'block_rotationX',
                                    self.setObjectBlockRX)

    def setObjectBlockRY(self, num, val):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(val)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('block_rotationY', val)

    @Variables.noUndo
    def getObjectBlockRY(self, num):

        return self.__getStringData(num-1, 'block_rotationY',
                                    self.setObjectBlockRY)

    def setObjectBlockRZ(self, num, val):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(val)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('block_rotationZ', val)

    @Variables.noUndo
    def getObjectBlockRZ(self, num):

        return self.__getStringData(num-1, 'block_rotationZ',
                                    self.setObjectBlockRZ)

    def setObjectBlockDX(self, num, val):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(val)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('block_displacementX', val)

    @Variables.noUndo
    def getObjectBlockDX(self, num):

        return self.__getStringData(num-1, 'block_displacementX',
                                    self.setObjectBlockDX)

    def setObjectBlockDY(self, num, val):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(val)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('block_displacementY', val)

    @Variables.noUndo
    def getObjectBlockDY(self, num):

        return self.__getStringData(num-1, 'block_displacementY',
                                    self.setObjectBlockDY)

    def setObjectBlockDZ(self, num, val):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(val)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('block_displacementZ', val)

    @Variables.noUndo
    def getObjectBlockDZ(self, num):

        return self.__getStringData(num-1, 'block_displacementZ',
                                    self.setObjectBlockDZ)





    # ----------------------------------

    @Variables.undoLocal
    def setObjectInit(self, num, is_init):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(is_init)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_initialization', is_init)

    @Variables.noUndo
    def getObjectInit(self, num):

        return self.__getStringData(num-1, 'object_initialization',
                                    self.setObjectInit)
    # ----------------------------------


    # ----------------------------------

    @Variables.undoLocal
    def setObjectPhysicalProperties(self, num, physical_pro):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(physical_pro)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_physical_properties', physical_pro)

    @Variables.noUndo
    def getObjectPhysicalProperties(self, num):

        return self.__getStringData(num-1, 'object_physical_properties',
                                    self.setObjectPhysicalProperties)
    # ----------------------------------


    @Variables.undoLocal
    def setObjectThermalSourceTerm(self, num, st):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(st)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_thermal_source_term', st)

    @Variables.noUndo
    def getObjectThermalSourceTerm(self, num):

        return self.__getStringData(num-1, 'object_thermal_source_term',
                                    self.setObjectThermalSourceTerm)
    # ----------------------------------


    @Variables.undoLocal
    def setObjectCHT(self, num, is_CHT):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(is_CHT)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_conjugate_heat_transfer', is_CHT)

    @Variables.noUndo
    def getObjectCHT(self, num):

        return self.__getStringData(num-1, 'object_conjugate_heat_transfer',
                                    self.setObjectCHT)
    # ----------------------------------


    @Variables.undoLocal
    def setObjectBoundaryConditionNature(self, num, bc_nature):
        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(bc_nature)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_boundary_condition_nature', bc_nature)


    @Variables.noUndo
    def getObjectBoundaryConditionNature(self, num):

        return self.__getStringData(num-1, 'object_boundary_condition_nature',
                                    self.setObjectBoundaryConditionNature)


    # ----------------------------------
    @Variables.undoLocal
    def setObjectBoundaryTemperature(self, num, temp):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_boundary_temperature', temp)


    @Variables.noUndo
    def getObjectBoundaryTemperature(self, num):

        return self.__getStringData(num-1, 'object_boundary_temperature',
                                    self.setObjectBoundaryTemperature)

    # ----------------------------------


    # ----------------------------------
    @Variables.undoLocal
    def setObjectBoundaryThermalFlux(self, num, flux):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_boundary_heat_flux', flux)


    @Variables.noUndo
    def getObjectBoundaryThermalFlux(self, num):

        return self.__getStringData(num-1, 'object_boundary_heat_flux',
                                    self.setObjectBoundaryThermalFlux)

    # ----------------------------------


    @Variables.undoLocal
    def setObjectBoundaryEnergyMode(self, num, mode):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(mode)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_boundary_energy_mode', mode)

    @Variables.noUndo
    def getObjectBoundaryEnergyMode(self, num):

        return self.__getStringData(num-1, 'object_boundary_energy_mode',
                                    self.setObjectBoundaryEnergyMode)
    # ----------------------------------

    # ----------------------------------


    @Variables.undoLocal
    def setObjectBoundaryVelocityMode(self, num, mode):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(mode)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_boundary_velocity_mode', mode)

    @Variables.noUndo
    def getObjectBoundaryVelocityMode(self, num):

        return self.__getStringData(num-1, 'object_boundary_velocity_mode',
                                    self.setObjectBoundaryVelocityMode)
    # ----------------------------------

    # ----------------------------------


    @Variables.undoLocal
    def setObjectModelingMode(self, num, mode):

        self.isLowerOrEqual(num, self.getNumberOfObjects())
        self.isStr(mode)

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_modeling_mode', mode)

    @Variables.noUndo
    def getObjectModelingMode(self, num):

        return self.__getStringData(num-1, 'object_modeling_mode',
                                    self.setObjectModelingMode)
    # ----------------------------------

    # ----------------------------------


    @Variables.undoLocal
    def setObjectMinCycle(self, num, min_cycle):
        self.isLowerOrEqual(num, self.getNumberOfObjects())

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_nb_cycle_min', min_cycle)

    @Variables.noUndo
    def getObjectMinCycle(self, num):

        return self.__getStringData(num-1, 'object_nb_cycle_min',
                                    self.setObjectMinCycle)
    # ----------------------------------

    # ----------------------------------


    @Variables.undoLocal
    def setObjectMaxCycle(self, num, max_cycle):
        self.isLowerOrEqual(num, self.getNumberOfObjects())

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_nb_cycle_max', max_cycle)

    @Variables.noUndo
    def getObjectMaxCycle(self, num):

        return self.__getStringData(num-1, 'object_nb_cycle_max',
                                    self.setObjectMaxCycle)
    # ----------------------------------

    # ----------------------------------


    @Variables.undoLocal
    def setObjectCvCriteria(self, num, criteria):
        self.isLowerOrEqual(num, self.getNumberOfObjects())

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node.xmlSetData('object_conv_criteria', criteria)

    @Variables.noUndo
    def getObjectCvCriteria(self, num):

        return self.__getStringData(num-1, 'object_conv_criteria',
                                    self.setObjectCvCriteria)
    # ----------------------------------

    # ----------------------------------


    @Variables.undoLocal
    def addSTLSeedPoint(self, num_obj, num_pt, x, y, z):
        self.isLowerOrEqual(num_obj, self.getNumberOfObjects())
        self.isLowerOrEqual(num_pt, self.getNumberOfSTLObjectSeedPoints(num_obj))

        node_obj = self.__node_ibm.xmlGetNodeList('ibm_object')[num_obj-1]

        num_pt += 1

        node_stl_point = node_obj.xmlInitNode('STL_exterior_point', id=num_pt)

        self.setSTLObjectSeedPoint(num_obj, num_pt, x, y, z)


    @Variables.noUndo
    def deleteSTLSeedPoint(self, num_obj, num_pt):
        self.isLowerOrEqual(num_obj, self.getNumberOfObjects())
        self.isLowerOrEqual(num_pt, self.getNumberOfSTLObjectSeedPoints(num_obj))

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num_obj-1]
        node.xmlRemoveChild('STL_exterior_point', id=num_pt)


    # ----------------------------------

    # ----------------------------------


    @Variables.undoLocal
    def setSTLObjectSeedPoint(self, num_obj, num_pt, x=None, y=None, z=None):
        self.isLowerOrEqual(num_obj, self.getNumberOfObjects())
        self.isLowerOrEqual(num_pt, self.getNumberOfSTLObjectSeedPoints(num_obj))

        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num_obj-1]
        node_stl_point = node.xmlInitNode('STL_exterior_point', id=num_pt)

        if x or x == 0:
            node_stl_point.xmlSetData('x', x)
        if y or y == 0:
            node_stl_point.xmlSetData('y', y)
        if z or z == 0:
            node_stl_point.xmlSetData('z', z)


    # ----------------------------------

    # ----------------------------------


    @Variables.noUndo
    def getSTLObjectSeedPoint(self, num, num_pt):
        node = self.__node_ibm.xmlGetNodeList('ibm_object')[num-1]
        node_stl_point = node.xmlGetChildNode('STL_exterior_point', id=num_pt)

        if node_stl_point is None:
            self.setSTLObjectSeedPoint(num,
                                       num_pt,
                                       self.defaultValues()['object_init'],
                                       self.defaultValues()['object_init'],
                                       self.defaultValues()['object_init'])
            return self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init'], \
                   self.defaultValues()['object_init']
        else:
            return node_stl_point.xmlGetString('x'), \
                   node_stl_point.xmlGetString('y'), \
                   node_stl_point.xmlGetString('z')


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
