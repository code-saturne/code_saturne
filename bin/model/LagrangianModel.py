# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2020 EDF S.A.
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
This module defines the lagrangian two phase flow modelling management.

This module contains the following classes and function:
- LagrangianModel
- LagrangianTestCase
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, unittest, logging

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.model.Common import *
from code_saturne.model.XMLvariables import Model, Variables
from code_saturne.model.XMLmodel import ModelTest
from code_saturne.model.CoalCombustionModel import CoalCombustionModel

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LagrangianModel")
log.setLevel(GuiParam.DEBUG)

#-------------------------------------------------------------------------------
# lagrangian model class
#-------------------------------------------------------------------------------

class LagrangianModel(Model):
    """
    Manage the input/output markups in the xml doc about Lagrangian module.
    """
    def __init__(self, case):
        """
        Constructor.
        """
        self.case = case
        self.node_lagr = self.case.root().xmlInitNode('lagrangian', 'model')
        self.__lagrangianCouplingMode = ('off', 'one_way', 'two_way', 'frozen')
        self.__particlesModels = ("thermal", "coal", "off")


    def defaultParticlesValues(self):
        """
        Return a dictionnary which contains default values.
        """
        default = {}
        default['model']                               = "off"
        default['restart']                             = "off"
        default['carrier_field_stationary']            = "off"
        default['deposition_submodel']                 = "off"
        default['particles_models']                    = "off"
        default['evaporation']                         = "off"
        default['break_up']                            = "off"
        default['coal_fouling']                        = "off"
        default['threshold_temperature']               = 600.
        default['critical_viscosity']                  = 10000.
        default['fouling_coefficient_1']               = 0.316608
        default['fouling_coefficient_2']               = -1.6786
        default['iteration_start']                     = 1
        default['thermal']                             = "off"
        default['dynamic']                             = "off"
        default['mass']                                = "off"
        default['scheme_order']                        = 1
        default['turbulent_dispersion']                = "on"
        default['fluid_particles_turbulent_diffusion'] = "off"
        default['complete_model_iteration']            = 0
        default['complete_model_direction']            = 4

        return default


    def lagrangianStatus(self):
        """
        Return a tuple with the lagrangian status allowed.
        """
        from code_saturne.model.TurbulenceModel import TurbulenceModel
        model = TurbulenceModel(self.case).getTurbulenceModel()
        del TurbulenceModel
        if model not in ('off',
                         'k-epsilon',
                         'k-epsilon-PL',
                         'Rij-epsilon',
                         'Rij-SSG',
                         'Rij-EBRSM',
                         'v2f-BL-v2/k',
                         'k-omega-SST',
                         'Spalart-Allmaras'):
            return ('off',)
        else:
            return self.__lagrangianStatus


    def particlesModels(self):
        """
        Return the list of models associated to the particles.
        """
        return self.__particlesModels


    def lagrangianCouplingMode(self):
        """
        Return all defined lagrangian models in a tuple.
        """
        return self.__lagrangianCouplingMode


    @Variables.undoGlobal
    def setLagrangianModel(self, mdl):
        """
        Update the lagrangian module model markup.
        """
        self.isInList(mdl, self.__lagrangianCouplingMode)

        old_mdl = self.node_lagr['model']
        self.node_lagr['model'] = mdl
        from code_saturne.model.OutputControlModel import OutputControlModel
        if mdl != 'off':
            OutputControlModel(self.case).addDefaultLagrangianWriter()
            OutputControlModel(self.case).addDefaultLagrangianMesh()
            self.node_out_lag = self.node_lagr.xmlInitChildNode('output')
        elif old_mdl and old_mdl != 'off':
            OutputControlModel(self.case).deleteMesh(['-3'])
            OutputControlModel(self.case).deleteWriter(['-3', '-4'])
        del OutputControlModel

        if mdl == 'two_way':
            node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        else:
            node_2way = self.node_lagr.xmlGetChildNode('two_way_coupling')
            if node_2way:
                node_2way.xmlRemoveNode()


    @Variables.noUndo
    def getLagrangianModel(self):
        """
        Return the status for lagrangian module markup from the XML document.
        """
        mdl = self.node_lagr['model']
        if mdl == "":
            mdl = self.defaultParticlesValues()['model']
            self.setLagrangianModel(mdl)
        return mdl


    @Variables.undoLocal
    def setRestart(self, status):
        """
        Update the restart status markup from the XML document.
        """
        self.isOnOff(status)
        node_restart = self.node_lagr.xmlInitChildNode('restart', 'status')
        node_restart['status'] = status


    @Variables.noUndo
    def getRestart(self):
        """
        Return status of restart file.
        """
        node_restart = self.node_lagr.xmlInitChildNode('restart', 'status')
        status = node_restart['status']
        if not status:
            status = self.defaultParticlesValues()['restart']
            self.setRestart(status)
        return status


    @Variables.undoGlobal
    def setCarrierFlowStationary(self, status):
        """
        Update the status for steady flow markup from the XML document.
        """
        self.isOnOff(status)
        node_steady = self.node_lagr.xmlInitChildNode('carrier_field_stationary', 'status')

        if not (self.getLagrangianModel() == "frozen"):
            node_steady['status'] = status
        else:
            node_steady['status'] = "on"


    @Variables.noUndo
    def getCarrierFlowStationary(self):
        """
        Return status of steady (on) or unsteady (off) state
        of the continuous phase flow.
        """
        node_steady = self.node_lagr.xmlInitChildNode('carrier_field_stationary', 'status')
        status = node_steady['status']
        if not status:
            status = self.defaultParticlesValues()['carrier_field_stationary']
            self.setCarrierFlowStationary(status)
        return status


    @Variables.undoLocal
    def setDepositionSubmodel(self, status):
        """
        Update the status for the particle deposition submodel.
        """
        self.isOnOff(status)
        node_deposition = self.node_lagr.xmlInitChildNode('deposition_submodel', 'status')
        node_deposition['status'] = status


    @Variables.noUndo
    def getDepositionSubmodel(self):
        """
        Return status for the particle deposition submodel.
        """
        node_deposition = self.node_lagr.xmlInitChildNode('deposition_submodel', 'status')
        status = node_deposition['status']
        if not status:
            status = self.defaultParticlesValues()['deposition_submodel']
            self.setDepositionSubmodel(status)
        return status


    @Variables.undoGlobal
    def setParticlesModel(self, model):
        """
        Update the particles model markup from the XML document.
        """
        self.isInList(model, self.__particlesModels)
        node_model = self.node_lagr.xmlInitChildNode('particles_models', 'model')
        node_model['model'] = model

        if model == "off":
            node_model.xmlRemoveChild('thermal')
            node_model.xmlRemoveChild('evaporation')
            node_model.xmlRemoveChild('break_up')
            node_model.xmlRemoveChild('coal_fouling')
        elif model == "thermal":
            node_model.xmlRemoveChild('coal_fouling')
        elif model == "coal":
            node_model.xmlRemoveChild('thermal')
            node_model.xmlRemoveChild('evaporation')
            node_model.xmlRemoveChild('break_up')
            self.getCoalFouling()


    def __nodeParticlesModel(self):
        """
        Return xml node for particles model. Create node if it does not exists.
        """
        node_model = self.node_lagr.xmlInitChildNode('particles_models', 'model')
        model = node_model['model']
        if model not in self.__particlesModels:
            model = self.defaultParticlesValues()['particles_models']
            self.setParticlesModel(model)
        return node_model


    @Variables.noUndo
    def getParticlesModel(self):
        """
        Return the current particles model.
        """
        return self.__nodeParticlesModel()['model']


    @Variables.undoLocal
    def setHeating(self, status):
        """
        Update the status for the activation of an evolution equation on the particle temperature.
        """
        self.isOnOff(status)
        node_model = self.__nodeParticlesModel()
        node_thermal = node_model.xmlInitChildNode('thermal', 'status')
        node_thermal['status'] = status


    @Variables.noUndo
    def getHeating(self):
        """
        Return status for the activation of an evolution equation on the particle temperature.
        """
        node_model = self.__nodeParticlesModel()
        node_temp = node_model.xmlInitChildNode('thermal', 'status')
        status = node_temp['status']
        if not status:
            status = self.defaultParticlesValues()['thermal']
            self.setHeating(status)
        return status


    @Variables.undoLocal
    def setEvaporation(self, status):
        """
        Update the status for the activation of an evolution equation on the particle temperature.
        """
        self.isOnOff(status)
        node_model = self.__nodeParticlesModel()
        node_mass = node_model.xmlInitChildNode('evaporation', 'status')
        node_mass['status'] = status


    @Variables.noUndo
    def getEvaporation(self):
        """
        Return status for the activation of an evolution equation on the particle temperature.
        """
        node_model = self.__nodeParticlesModel()
        node_mass = node_model.xmlInitChildNode('evaporation', 'status')
        status = node_mass['status']
        if not status:
            status = self.defaultParticlesValues()['evaporation']
            self.setEvaporation(status)
        return status


    @Variables.undoGlobal
    def setCoalFouling(self, status):
        """
        Update the status for coal particle fouling.
        """
        self.isOnOff(status)
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        node_coal['status'] = status
        if status == "on":
            self.coalCombustionModel = CoalCombustionModel(self.case)
            for icoal in range(self.coalCombustionModel.getCoalNumber()):
                self.getThresholdTemperatureOfFouling(icoal+1)
                self.getCriticalViscosityOfFouling(icoal+1)
                self.getCoef1OfFouling(icoal+1)
                self.getCoef2OfFouling(icoal+1)
        elif status == "off":
            node_coal.xmlRemoveChild('threshold_temperature')
            node_coal.xmlRemoveChild('critical_viscosity')
            node_coal.xmlRemoveChild('fouling_coefficient_1')
            node_coal.xmlRemoveChild('fouling_coefficient_2')


    @Variables.noUndo
    def getCoalFouling(self):
        """
        Return status for coal particle fouling.
        """
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        status = node_coal['status']
        if not status:
            status = self.defaultParticlesValues()['coal_fouling']
            self.setCoalFouling(status)
        return status


    @Variables.undoLocal
    def setThresholdTemperatureOfFouling(self, icoal, value):
        """
        Update the value for the threshold temperature for the coal specified.
        """
        self.isFloat(value)
        self.isGreater(value, 0.)
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        node_temp = node_coal.xmlInitChildNode('threshold_temperature', coal=icoal)
        node_temp.xmlSetTextNode(str(value))


    @Variables.noUndo
    def getThresholdTemperatureOfFouling(self, icoal):
        """
        Return the value for the threshold temperature for the specified coal.
        """
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        value = node_coal.xmlGetDouble('threshold_temperature', coal=icoal)
        if value == None:
            value = self.defaultParticlesValues()['threshold_temperature']
            self.setThresholdTemperatureOfFouling(icoal, value)
        return value


    @Variables.undoLocal
    def setCriticalViscosityOfFouling(self, icoal, value):
        """
        Update the value for the critical viscosity for the coal specified.
        """
        self.isFloat(value)
        self.isGreater(value, 0.)
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        node_visc = node_coal.xmlInitChildNode('critical_viscosity', coal=icoal)
        node_visc.xmlSetTextNode(str(value))


    @Variables.noUndo
    def getCriticalViscosityOfFouling(self, icoal):
        """
        Return the value for the critical viscosity for the coal specified.
        """
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        value = node_coal.xmlGetDouble('critical_viscosity', coal=icoal)
        if value == None:
            value = self.defaultParticlesValues()['critical_viscosity']
            self.setCriticalViscosityOfFouling(icoal, value)
        return value


    @Variables.undoLocal
    def setCoef1OfFouling(self, icoal, value):
        """
        Update the value for the coefficient for the coal specified.
        """
        self.isFloat(value)
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        node_visc = node_coal.xmlInitChildNode('fouling_coefficient_1', coal=icoal)
        node_visc.xmlSetTextNode(str(value))


    @Variables.noUndo
    def getCoef1OfFouling(self, icoal):
        """
        Return the value for the coefficient for the coal specified.
        """
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        value = node_coal.xmlGetDouble('fouling_coefficient_1', coal=icoal)
        if not value:
            value = self.defaultParticlesValues()['fouling_coefficient_1']
            self.setCoef1OfFouling(icoal, value)
        return value


    @Variables.undoLocal
    def setCoef2OfFouling(self, icoal, value):
        """
        Update the value for the coefficient for the coal specified.
        """
        self.isFloat(value)
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        node_visc = node_coal.xmlInitChildNode('fouling_coefficient_2', coal=icoal)
        node_visc.xmlSetTextNode(str(value))


    @Variables.noUndo
    def getCoef2OfFouling(self, icoal):
        """
        Return the value for the coefficient for the coal specified.
        """
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        value = node_coal.xmlGetDouble('fouling_coefficient_2', coal=icoal)
        if not value:
            value = self.defaultParticlesValues()['fouling_coefficient_2']
            self.setCoef2OfFouling(icoal, value)
        return value


    @Variables.undoLocal
    def set2WayCouplingStartIteration(self, value):
        """
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_2way.xmlSetData('iteration_start', value)


    @Variables.noUndo
    def get2WayCouplingStartIteration(self):
        """
        """
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        niter = node_2way.xmlGetInt('iteration_start')
        if niter == None:
            niter = self.defaultParticlesValues()['iteration_start']
            self.set2WayCouplingStartIteration(niter)
        return niter


    @Variables.undoLocal
    def set2WayCouplingDynamic(self, status):
        """
        Update the status for 2 way coupling on continuous phase dynamic markup from the XML document.
        """
        self.isOnOff(status)
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_dyn = node_2way.xmlInitChildNode('dynamic')
        node_dyn['status'] = status


    @Variables.noUndo
    def get2WayCouplingDynamic(self):
        """
        Return status of 2 way coupling on continuous phase dynamic.
        """
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_dyn = node_2way.xmlInitChildNode('dynamic')
        status = node_dyn['status']
        if not status:
            status = self.defaultParticlesValues()['dynamic']
            self.set2WayCouplingDynamic(status)
        return status


    @Variables.undoLocal
    def set2WayCouplingMass(self, status):
        """
        Update the status markup for 2 way coupling on mass from the XML document.
        """
        self.isOnOff(status)
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_mass = node_2way.xmlInitChildNode('mass')
        node_mass['status'] = status


    @Variables.noUndo
    def get2WayCouplingMass(self):
        """
        Return status of 2 way coupling on mass.
        """
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_mass = node_2way.xmlInitChildNode('mass')
        status = node_mass['status']
        if not status:
            status = self.defaultParticlesValues()['mass']
            self.set2WayCouplingMass(status)
        return status


    @Variables.undoLocal
    def set2WayCouplingTemperature(self, status):
        """
        Update the status markup for 2 way coupling on temperature from the XML document.
        """
        self.isOnOff(status)
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_temp = node_2way.xmlInitChildNode('thermal')
        node_temp['status'] = status


    @Variables.noUndo
    def get2WayCouplingTemperature(self):
        """
        Return status of 2 way coupling on temperature.
        """
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_temp = node_2way.xmlInitChildNode('thermal')
        status = node_temp['status']
        if not status:
            status = self.defaultParticlesValues()['thermal']
            self.set2WayCouplingTemperature(status)
        return status


    @Variables.undoLocal
    def setSchemeOrder(self, value):
        """
        Update value for scheme order.
        """
        self.isInt(value)
        self.isInList(value, (1,2))
        node_order = self.node_lagr.xmlInitChildNode('scheme_order', 'choice')
        node_order['choice'] = value


    @Variables.noUndo
    def getSchemeOrder(self):
        """
        Return value for scheme order.
        """
        node_order = self.node_lagr.xmlInitChildNode('scheme_order', 'choice')
        if node_order:
            val = node_order['choice']
            if val == "":
                val = self.defaultParticlesValues()['scheme_order']
                self.setSchemeOrder(val)
        return val


    @Variables.undoLocal
    def setTurbulentDispersion(self, status):
        """
        Update the status markup for turbulent dispersion status from the XML document.
        """
        self.isOnOff(status)
        node_turb = self.node_lagr.xmlInitChildNode('turbulent_dispersion', 'status')
        node_turb['status'] = status


    @Variables.noUndo
    def getTurbulentDispersion(self):
        """
        Return status of turbulent dispersion status.
        """
        node_turb = self.node_lagr.xmlInitChildNode('turbulent_dispersion', 'status')
        status = node_turb['status']
        if not status:
            status = self.defaultParticlesValues()['turbulent_dispersion']
            self.setTurbulentDispersion(status)
        return status


    @Variables.undoLocal
    def setTurbulentDiffusion(self, status):
        """
        Update the status markup for turbulent diffusion status from the XML document.
        """
        self.isOnOff(status)
        node_turb = self.node_lagr.xmlInitChildNode('fluid_particles_turbulent_diffusion', 'status')
        node_turb['status'] = status


    @Variables.noUndo
    def getTurbulentDiffusion(self):
        """
        Return status of turbulent diffusion status.
        """
        node_turb = self.node_lagr.xmlInitChildNode('fluid_particles_turbulent_diffusion', 'status')
        status = node_turb['status']
        if not status:
            status = self.defaultParticlesValues()['fluid_particles_turbulent_diffusion']
            self.setTurbulentDiffusion(status)
        return status


    @Variables.undoLocal
    def setCompleteModelStartIteration(self, iteration):
        """
        Set value for complete model start iteration.
        """
        self.isInt(iteration)
        self.isGreaterOrEqual(iteration, 0)
        self.node_lagr.xmlSetData('complete_model', iteration)


    @Variables.noUndo
    def getCompleteModelStartIteration(self):
        """
        Return value for complete model iteration.
        """
        iteration = self.node_lagr.xmlGetInt('complete_model')
        if iteration == None:
            iteration = self.defaultParticlesValues()['complete_model_iteration']
            self.setCompleteModelStartIteration(iteration)
        return iteration


    @Variables.undoLocal
    def setCompleteModelDirection(self, value):
        """
        Set value for complete model direction.
        """
        self.isInt(value)
        self.isInList(value, (1,2,3,4))
        node_direction = self.node_lagr.xmlInitChildNode('complete_model_direction', 'choice')
        node_direction['choice'] = value


    @Variables.noUndo
    def getCompleteModelDirection(self):
        """
        Return value for complete model direction.
        """
        node_direction = self.node_lagr.xmlInitChildNode('complete_model_direction', 'choice')
        if node_direction:
            val = node_direction['choice']
            if val == "":
                val = self.defaultParticlesValues()['complete_model_direction']
                self.setCompleteModelDirection(val)
        return val

#-------------------------------------------------------------------------------
# Lagrangian test case
#-------------------------------------------------------------------------------

class LagrangianTestCase(ModelTest):
    """
    Unittest.
    """
    def checkLagrangianInstantiation(self):
        """Check whether the LagrangianModel class could be instantiated"""
        model = None
        model = LagrangianModel(self.case)

        assert model != None, 'Could not instantiate LagrangianModel'

        doc = """<lagrangian model="off"/>"""

        assert model.node_lagr == self.xmlNodeFromString(doc), \
               'Could not instantiate LagrangianModel'


    def checkRestart(self):
        """Check whether the restart method could be set and get."""
        mdl = LagrangianModel(self.case)
        status = mdl.getRestart()

        assert status == mdl.defaultParticlesValues()['restart'], \
            'Could not get default values for restart status'

        mdl.setRestart('on')
        doc = """<lagrangian model="off">
                      <restart status="on"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set values for restart status'


    def checkCarrierFlowStationary(self):
        """Check whether the stationary behavior of the carrier flow
        could be set and get."""
        mdl = LagrangianModel(self.case)

        status = mdl.getCarrierFlowStationary()

        assert status == mdl.defaultParticlesValues()['carrier_field_stationary'], \
            'Could not get default values for stationary \
            behavior of the carrier flow status'

        mdl.setCarrierFlowStationary('on')
        doc = """<lagrangian model="off">
                     <carrier_field_stationary status="on"/>
                     <coupling_mode model="one_way"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set default values for stationary \
            behavior of the carrier flow status'


    def checkDepositionSubmodel(self):
        """Check whether the particle deposition model could be set and get."""
        mdl = LagrangianModel(self.case)
        status = mdl.getDepositionSubmodel()

        assert status == mdl.defaultParticlesValues()['deposition_submodel'], \
            'Could not get default values for particle deposition model status'

        mdl.setDepositionSubmodel('on')
        doc = """<lagrangian model="off">
                      <deposition_submodel status="on"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set values for particle deposition model status'


    def checkParticlesModel(self):
        """Check whether the particles models could be set and get."""
        mdl = LagrangianModel(self.case)
        mdl.setParticlesModel('thermal')
        doc = """<lagrangian model="off">
                     <particles_models model="thermal"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc) ,\
            'Could not set values for particles model (thermal)'

        mdl.setParticlesModel('coal')
        doc = """<lagrangian model="off">
                     <particles_models model="coal">
                         <coal_fouling status="off"/>
                     </particles_models>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc) ,\
            'Could not set values for particles model (coal)'

        for name in mdl.particlesModels():
            mdl.setParticlesModel(name)
            name2 = mdl.getParticlesModel()
            assert name == name2 ,\
                   'Could not use the get/setParticlesModel method for model name %s '%name


    def checkBreakUp(self):
        """Check whether the break-up model could be set and get."""
        mdl = LagrangianModel(self.case)
        mdl.setParticlesModel('thermal')
        status = mdl.getBreakUp()

        assert status == mdl.defaultParticlesValues()['break_up'], \
            'Could not get default values for break-up model'

        mdl.setBreakUp('on')
        doc = """<lagrangian model="off">
                     <particles_models model="thermal">
                         <break_up status="on"/>
                     </particles_models>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set value for break-up model'


    def checkHeating(self):
        """Check whether the heating model could be set and get."""
        mdl = LagrangianModel(self.case)
        mdl.setParticlesModel('thermal')
        status = mdl.getHeating()

        assert status == mdl.defaultParticlesValues()['thermal'], \
            'Could not get default values for heating model'

        mdl.setHeating('on')
        doc = """<lagrangian model="off">
                     <particles_models model="thermal">
                         <thermal status="on"/>
                     </particles_models>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set value for heating model'


    def checkEvaporation(self):
        """Check whether the heating model could be set and get."""
        mdl = LagrangianModel(self.case)
        mdl.setParticlesModel('thermal')
        status = mdl.getEvaporation()

        assert status == mdl.defaultParticlesValues()['evaporation'], \
            'Could not get default values for heating model'

        mdl.setEvaporation('on')
        doc = """<lagrangian model="off">
                     <particles_models model="thermal">
                         <evaporation status="on"/>
                     </particles_models>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set value for heating model'


    def checkCoalFouling(self):
        """Check whether the heating model could be set and get."""
        mdl = LagrangianModel(self.case)
        mdl.setParticlesModel('coal')
        status = mdl.getCoalFouling()

        assert status == mdl.defaultParticlesValues()['coal_fouling'], \
            'Could not get default values for coal fouling model'

        mdl.setCoalFouling('on')
        doc = """<lagrangian model="off">
                     <particles_models model="coal">
                         <coal_fouling status="on">
                            <threshold_temperature coal="1">
                                600.0
                            </threshold_temperature>
                            <critical_viscosity coal="1">
                                10000.0
                            </critical_viscosity>
                            <fouling_coefficient_1 coal="1">
                                0.316608
                            </fouling_coefficient_1>
                            <fouling_coefficient_2 coal="1">
                                -1.6786
                            </fouling_coefficient_2>
                        </coal_fouling>
                     </particles_models>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set value for coal fouling model'

        mdl.setThresholdTemperatureOfFouling(1, 800)
        mdl.setCriticalViscosityOfFouling(1, 12300)
        mdl.setCoef1OfFouling(1, -0.001)
        mdl.setCoef2OfFouling(1, 0.002)
        doc = """<lagrangian model="off">
                    <particles_models model="coal">
                        <coal_fouling status="on">
                            <threshold_temperature coal="1">
                                800
                            </threshold_temperature>
                            <critical_viscosity coal="1">
                                12300
                            </critical_viscosity>
                            <fouling_coefficient_1 coal="1">
                                -0.001
                            </fouling_coefficient_1>
                            <fouling_coefficient_2 coal="1">
                                0.002
                            </fouling_coefficient_2>
                        </coal_fouling>
                    </particles_models>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set values for coal fouling model'


    def check2WayCouplingStartIteration(self):
        """Check whether the 2 way coupling start iteration could be set and get."""
        mdl = LagrangianModel(self.case)
        n = mdl.get2WayCouplingStartIteration()

        assert n == mdl.defaultParticlesValues()['iteration_start'], \
            'Could not set default values for 2 way coupling start iteration'

        mdl.set2WayCouplingStartIteration(123456)
        doc = """<lagrangian model="off">
                    <two_way_coupling>
                        <iteration_start>
                            123456
                        </iteration_start>
                    </two_way_coupling>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not get default values for 2 way coupling start iteration'


    def check2WayCouplingDynamic(self):
        """Check whether the 2 way coupling dynamic could be set and get."""
        mdl = LagrangianModel(self.case)
        n = mdl.get2WayCouplingDynamic()

        assert n == mdl.defaultParticlesValues()['dynamic'], \
            'Could not set default values for 2 way coupling dynamic'

        mdl.set2WayCouplingDynamic("on")
        doc = """<lagrangian model="off">
                    <two_way_coupling>
                        <dynamic status="on"/>
                    </two_way_coupling>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not get default values for 2 way coupling dynamic'


    def check2WayCouplingMass(self):
        """Check whether the 2 way coupling mass could be set and get."""
        mdl = LagrangianModel(self.case)
        n = mdl.get2WayCouplingMass()

        assert n == mdl.defaultParticlesValues()['mass'], \
            'Could not set default values for 2 way coupling mass'

        mdl.set2WayCouplingMass("on")
        doc = """<lagrangian model="off">
                    <two_way_coupling>
                        <mass status="on"/>
                    </two_way_coupling>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not get default values for 2 way coupling mass'


    def check2WayCouplingThermal(self):
        """Check whether the 2 way coupling mass could be set and get."""
        mdl = LagrangianModel(self.case)
        n = mdl.get2WayCouplingTemperature()

        assert n == mdl.defaultParticlesValues()['thermal'], \
            'Could not set default values for 2 way coupling thermal'

        mdl.set2WayCouplingTemperature("on")
        doc = """<lagrangian model="off">
                    <two_way_coupling>
                        <thermal status="on"/>
                    </two_way_coupling>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not get default values for 2 way coupling thermal'


    def checkSchemeOrder(self):
        """Check whether the scheme order could be set and get."""
        mdl = LagrangianModel(self.case)
        n = mdl.getSchemeOrder()

        assert n == mdl.defaultParticlesValues()['scheme_order'], \
            'Could not set default values for scheme order.'

        mdl.setSchemeOrder(1)
        doc = """<lagrangian model="off">
                      <scheme_order choice="1"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
             'Could not get default values for scheme order.'


    def checkTurbulentDispersion(self):
        """Check whether the turbulent dispersion could be set and get."""
        mdl = LagrangianModel(self.case)
        status = mdl.getTurbulentDispersion()

        assert status == mdl.defaultParticlesValues()['turbulent_dispersion'], \
            'Could not get default values for turbulent dispersion status.'

        mdl.setTurbulentDispersion('on')
        doc = """<lagrangian model="off">
                     <turbulent_dispersion status="on"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set default values for turbulent dispersion status.'


    def checkTurbulentDiffusion(self):
        """Check whether the turbulent diffusion could be set and get."""
        mdl = LagrangianModel(self.case)
        status = mdl.getTurbulentDiffusion()

        assert status == mdl.defaultParticlesValues()['fluid_particles_turbulent_diffusion'], \
            'Could not get default values for turbulent diffusion status.'

        mdl.setTurbulentDiffusion('on')
        doc = """<lagrangian model="off">
                     <fluid_particles_turbulent_diffusion status="on"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set default values for turbulent diffusion status.'


    def checkCompleteModelStartIteration(self):
        """Check whether the complete model start iteration could be set and get."""
        mdl = LagrangianModel(self.case)
        status = mdl.getCompleteModelStartIteration()

        assert status == mdl.defaultParticlesValues()['complete_model_iteration'], \
            'Could not get default values for complete model start iteration.'

        mdl.setCompleteModelStartIteration(1234)
        doc = """<lagrangian model="off">
                     <complete_model>1234</complete_model>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set default values for complete model start iteration.'


    def checkCompleteModelDirection(self):
        """Check whether the complete model direction could be set and get."""
        mdl = LagrangianModel(self.case)
        status = mdl.getCompleteModelDirection()

        assert status == mdl.defaultParticlesValues()['complete_model_direction'], \
            'Could not get default values for complete model direction.'

        mdl.setCompleteModelDirection(2)
        doc = """<lagrangian model="off">
                    <complete_model_direction choice="2"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set default values for complete model direction.'


def suite():
    testSuite = unittest.makeSuite(LagrangianTestCase, "check")
    return testSuite


def runTest():
    print(os.path.basename(__file__))
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
