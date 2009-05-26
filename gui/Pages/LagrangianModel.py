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

from Base.Common import *
import Base.Toolbox as Tool
from Base.XMLvariables import Model
from Base.XMLmodel import ModelTest
import Pages.CoalThermoChemistry as CoalThermoChemistry

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger("LagrangianModel")
log.setLevel(Tool.GuiParam.DEBUG)

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
        self.getLagrangianStatus()
        self.__lagrangianStatus = ('off', 'on')
        self.__lagrangianCouplingMode = ('one_way', 'two_way', 'frozen')
        self.__particlesModels = ("thermal", "coal", "off")


    def defaultParticlesValues(self):
        """
        Return a dictionnary which contains default values.
        """
        default = {}
        default['model'] = "off"
        default['coupling_mode'] = "one_way"
        default['restart'] = "off"
        default['carrier_field_stationary'] = "off"
        default['particles_max_number'] = 10000
        default['continuous_injection'] = "off"
        default['particles_models'] = "off"
        default['thermal'] = "off"
        #default['particle_temperature'] = 700.
        #default['particle_specific_heat'] = 5200.
        default['evaporation'] = "off"
        default['break_up'] = "off"
        default['coal_fouling'] = "off"
        default['threshold_temperature'] = 600.
        default['critical_viscosity'] = 10000.
        default['fouling_coefficient_1'] = 0.316608
        default['fouling_coefficient_2'] = -1.6786
        default['iteration_start'] = 1
        default['thermal'] = "off"
        default['dynamic'] = "off"
        default['mass'] = "off"
        default['scheme_order'] = 2
        default['turbulent_dispersion'] = "on"
        default['fluid_particles_turbulent_diffusion'] = "off"
        default['complete_model_iteration'] = 0
        default['complete_model_direction'] = 1
        return default


    def lagrangianStatus(self):
        """
        Return a tuple with the lagrangian status allowed.
        """
        from TurbulenceModel import TurbulenceModel
        model = TurbulenceModel(self.case).getTurbulenceModel()
        del TurbulenceModel
        if model not in ('off',
                         'k-epsilon',
                         'k-epsilon-PL',
                         'Rij-epsilon',
                         'Rij-SSG',
                         'v2f-phi',
                         'k-omega-SST'):
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


    def setLagrangianStatus(self, status):
        """
        Update the lagrangian module status markup.
        """
        self.isOnOff(status)
        self.node_lagr['model'] = status

        # WARNING: the 'coal_lagr' model is deprecated.
#        if status == 'off':
#            import CoalCombustionModel
#            coal = CoalCombustionModel.CoalCombustionModel(self.case).getCoalCombustionModel()
#            # WARNING: the 'coal_lagr' model is deprecated.
#            if coal == 'coal_lagr':
#                CoalCombustionModel.CoalCombustionModel(self.case).setCoalCombustion('off')
#            del CoalCombustionModel


    def getLagrangianStatus(self):
        """
        Return the status for lagrangian module markup from the XML document.
        """
        status = self.node_lagr['model']
        if status == "":
            status = self.defaultParticlesValues()['model']
            self.setLagrangianStatus(status)
        return status


    def setCouplingMode(self, model):
        """
        Update the lagrangian model markup from the XML document.
        """
        self.isInList(model, self.__lagrangianCouplingMode)
        node_coupling = self.node_lagr.xmlInitChildNode('coupling_mode', 'model')
        node_coupling['model'] = model

        if model == 'two_way':
            node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')


    def getCouplingMode(self):
        """
        Return the current lagrangian model.
        """
        node_coupling = self.node_lagr.xmlInitChildNode('coupling_mode', 'model')
        model = node_coupling['model']
        if model not in self.__lagrangianCouplingMode:
            model = self.defaultParticlesValues()['coupling_mode']
            self.setCouplingMode(model)
        return model


    def setRestart(self, status):
        """
        Update the restart status markup from the XML document.
        """
        self.isOnOff(status)
        node_restart = self.node_lagr.xmlInitChildNode('restart', 'status')
        node_restart['status'] = status


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


    def setCarrierFlowStationary(self, status):
        """
        Update the status for steady flow markup from the XML document.
        """
        self.isOnOff(status)
        node_steady = self.node_lagr.xmlInitChildNode('carrier_field_stationary', 'status')
        if not (self.getCouplingMode() == "frozen" and status == "off"):
            node_steady['status'] = status


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


    def setMaxNumber(self, value):
        """
        Update value for maximum number of particles allowed
        simultaneously in the calculation domain.
        """
        self.isInt(value)
        self.isGreater(value, 0)
        self.node_lagr.xmlSetData('particles_max_number', value)


    def getMaxNumber(self):
        """
        Return the value for maximum number of particles allowed
        simultaneously in the calculation domain.
        """
        nbpmax = self.node_lagr.xmlGetInt('particles_max_number')
        if nbpmax == None:
            nbpmax = self.defaultParticlesValues()['particles_max_number']
            self.setMaxNumber(nbpmax)
        return nbpmax


    def setContinuousInjection(self, status):
        """
        Update the status for continuous injection of particles.
        """
        self.isOnOff(status)
        node_injection = self.node_lagr.xmlInitChildNode('continuous_injection', 'status')
        node_injection['status'] = status


    def getContinuousInjection(self):
        """
        Return status for continuous injection of particles.
        """
        node_injection = self.node_lagr.xmlInitChildNode('continuous_injection', 'status')
        status = node_injection['status']
        if not status:
            status = self.defaultParticlesValues()['continuous_injection']
            self.setContinuousInjection(status)
        return status


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


    def getParticlesModel(self):
        """
        Return the current particles model.
        """
        return self.__nodeParticlesModel()['model']


    def setBreakUp(self, status):
        """
        Update the status for activation of evolution equation on the particle diameter.
        """
        self.isOnOff(status)
        node_model = self.__nodeParticlesModel()
        node_diameter = node_model.xmlInitChildNode('break_up', 'status')
        node_diameter['status'] = status


    def getBreakUp(self):
        """
        Return status for the activation of an evolution equation on the particle diameter.
        """
        node_model = self.__nodeParticlesModel()
        node_diameter = node_model.xmlInitChildNode('break_up', 'status')
        status = node_diameter['status']
        if not status:
            status = self.defaultParticlesValues()['break_up']
            self.setBreakUp(status)
        return status


    def setHeating(self, status):
        """
        Update the status for the activation of an evolution equation on the particle temperature.
        """
        self.isOnOff(status)
        node_model = self.__nodeParticlesModel()
        node_thermal = node_model.xmlInitChildNode('thermal', 'status')
        node_thermal['status'] = status
        #if status == "on":
            #node_temp = node_thermal.xmlInitChildNode('particle_temperature')
            #temp = self.getParticlesTemperatureValue()
            #node_temp.xmlSetTextNode(str(temp))
            #node_cp = node_thermal.xmlInitChildNode('particle_specific_heat')
            #cp = self.getParticlesSpecificHeatValue()
            #node_cp.xmlSetTextNode(str(cp))


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


    #def setParticlesTemperatureValue(self, value):
        #"""
        #"""
        #self.isFloat(value)
        #node_model = self.__nodeParticlesModel()
        #node_thermal = node_model.xmlInitChildNode('thermal', 'status')
        ##check if node_thermal status == "on" ?
        #node_temp = node_thermal.xmlInitChildNode('particle_temperature')
        #node_temp.xmlSetTextNode(str(value))


    #def getParticlesTemperatureValue(self):
        #"""
        #"""
        #node_model = self.__nodeParticlesModel()
        #node_thermal = node_model.xmlInitChildNode('thermal', 'status')
        #value = node_thermal.xmlGetDouble('particle_temperature')
        #if not value:
            #value = self.defaultParticlesValues()['particle_temperature']
            #self.setParticlesTemperatureValue(value)
        #return value


    #def setParticlesSpecificHeatValue(self, value):
        #"""
        #"""
        #self.isFloat(value)
        #node_model = self.__nodeParticlesModel()
        #node_thermal = node_model.xmlInitChildNode('thermal', 'status')
        ##check if node_thermal status == "on" ?
        #node_cp = node_thermal.xmlInitChildNode('particle_specific_heat')
        #node_cp.xmlSetTextNode(str(value))


    #def getParticlesSpecificHeatValue(self):
        #"""
        #"""
        #node_model = self.__nodeParticlesModel()
        #node_thermal = node_model.xmlInitChildNode('thermal', 'status')
        #value = node_thermal.xmlGetDouble('particle_specific_heat')
        #if not value:
            #value = self.defaultParticlesValues()['particle_specific_heat']
            #self.setParticlesSpecificHeatValue(value)
        #return value


    def setEvaporation(self, status):
        """
        Update the status for the activation of an evolution equation on the particle temperature.
        """
        self.isOnOff(status)
        node_model = self.__nodeParticlesModel()
        node_mass = node_model.xmlInitChildNode('evaporation', 'status')
        node_mass['status'] = status


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


    def setCoalFouling(self, status):
        """
        Update the status for coal particle fouling.
        """
        self.isOnOff(status)
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        node_coal['status'] = status
        if status == "on":
            self.coalThermoChModel = CoalThermoChemistry.CoalThermoChemistryModel("dp_FCP", self.case)
            coals = self.coalThermoChModel.getCoals()
            for icoal in range(coals.getNumber()):
                self.getThresholdTemperatureOfFouling(icoal+1)
                self.getCriticalViscosityOfFouling(icoal+1)
                self.getCoef1OfFouling(icoal+1)
                self.getCoef2OfFouling(icoal+1)
        elif status == "off":
            node_coal.xmlRemoveChild('threshold_temperature')
            node_coal.xmlRemoveChild('critical_viscosity')
            node_coal.xmlRemoveChild('fouling_coefficient_1')
            node_coal.xmlRemoveChild('fouling_coefficient_2')


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


    def setCoef1OfFouling(self, icoal, value):
        """
        Update the value for the coefficient for the coal specified.
        """
        self.isFloat(value)
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        node_visc = node_coal.xmlInitChildNode('fouling_coefficient_1', coal=icoal)
        node_visc.xmlSetTextNode(str(value))


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


    def setCoef2OfFouling(self, icoal, value):
        """
        Update the value for the coefficient for the coal specified.
        """
        self.isFloat(value)
        node_model = self.__nodeParticlesModel()
        node_coal = node_model.xmlInitChildNode('coal_fouling', 'status')
        node_visc = node_coal.xmlInitChildNode('fouling_coefficient_2', coal=icoal)
        node_visc.xmlSetTextNode(str(value))


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


    def set2WayCouplingStartIteration(self, value):
        """
        """
        self.isInt(value)
        self.isGreaterOrEqual(value, 0)
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_2way.xmlSetData('iteration_start', value)


    def get2WayCouplingStartIteration(self):
        """
        """
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        niter = node_2way.xmlGetInt('iteration_start')
        if niter == None:
            niter = self.defaultParticlesValues()['iteration_start']
            self.set2WayCouplingStartIteration(niter)
        return niter


    def set2WayCouplingDynamic(self, status):
        """
        Update the status for 2 way coupling on continuous phase dynamic markup from the XML document.
        """
        self.isOnOff(status)
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_dyn = node_2way.xmlInitChildNode('dynamic')
        node_dyn['status'] = status


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


    def set2WayCouplingMass(self, status):
        """
        Update the status markup for 2 way coupling on mass from the XML document.
        """
        self.isOnOff(status)
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_mass = node_2way.xmlInitChildNode('mass')
        node_mass['status'] = status


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


    def set2WayCouplingTemperature(self, status):
        """
        Update the status markup for 2 way coupling on temperature from the XML document.
        """
        self.isOnOff(status)
        node_2way = self.node_lagr.xmlInitChildNode('two_way_coupling')
        node_temp = node_2way.xmlInitChildNode('thermal')
        node_temp['status'] = status


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


    def setSchemeOrder(self, value):
        """
        Update value for scheme order.
        """
        self.isInt(value)
        self.isInList(value, (1,2))
        node_order = self.node_lagr.xmlInitChildNode('scheme_order', 'choice')
        node_order['choice'] = value


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


    def setTurbulentDispersion(self, status):
        """
        Update the status markup for turbulent dispersion status from the XML document.
        """
        self.isOnOff(status)
        node_turb = self.node_lagr.xmlInitChildNode('turbulent_dispersion', 'status')
        node_turb['status'] = status


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


    def setTurbulentDiffusion(self, status):
        """
        Update the status markup for turbulent diffusion status from the XML document.
        """
        self.isOnOff(status)
        node_turb = self.node_lagr.xmlInitChildNode('fluid_particles_turbulent_diffusion', 'status')
        node_turb['status'] = status


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


    def setCompleteModelStartIteration(self, iteration):
        """
        Set value for complete model start iteration.
        """
        self.isInt(iteration)
        self.isGreaterOrEqual(iteration, 0)
        self.node_lagr.xmlSetData('complete_model', iteration)


    def getCompleteModelStartIteration(self):
        """
        Return value for complete model iteration.
        """
        iteration = self.node_lagr.xmlGetInt('complete_model')
        if iteration == None:
            iteration = self.defaultParticlesValues()['complete_model_iteration']
            self.setCompleteModelStartIteration(iteration)
        return iteration


    def setCompleteModelDirection(self, value):
        """
        Set value for complete model direction.
        """
        self.isInt(value)
        self.isInList(value, (1,2,3))
        node_direction = self.node_lagr.xmlInitChildNode('complete_model_direction', 'choice')
        node_direction['choice'] = value


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


    def checkLagrangianStatus(self):
        """Check whether the Lagrangian status could be set and get."""
        mdl = LagrangianModel(self.case)
        mdl.setLagrangianStatus("on")

        assert mdl.node_lagr == self.xmlNodeFromString("""<lagrangian model="on"/>"""), \
               'Could not get lagrangian status.'

        assert mdl.getLagrangianStatus() == "on", \
               'Could not get lagrangian status.'


    def checklagrangianStatus(self):
        """Check whether the lagrangianStatus could be get."""
        from TurbulenceModel import TurbulenceModel
        mdl = LagrangianModel(self.case)
        TurbulenceModel(self.case).setTurbulenceModel('LES_Smagorinsky')

        assert mdl.lagrangianStatus() == ('off',), \
               'Could not use the lagrangianStatus method'


    def checkLagrangianModel(self):
        """Check whether the LagrangianModel could be set and get."""
        mdl = LagrangianModel(self.case)

        assert mdl.lagrangianStatus() == ('off', 'on'), \
               'Could not use the lagrangianStatus method'

        mdl.setCouplingMode("frozen")
        doc = """<lagrangian model="off">
                     <coupling_mode model="frozen"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
               'Could not get lagrangian model.'

        for mode in mdl.lagrangianCouplingMode():
            mdl.setCouplingMode(mode)
            name = mdl.getCouplingMode()
            assert mode == mdl.getCouplingMode(), \
                   'Could not use the get/setCouplingMode method for %s'% mode


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


    def checkMaxNumber(self):
        """Check whether the max number of particles method could be set and get."""
        mdl = LagrangianModel(self.case)
        n = mdl.getMaxNumber()

        assert n == mdl.defaultParticlesValues()['particles_max_number'], \
            'Could not set default values for max number of particles'

        mdl.setMaxNumber(123456)
        doc = """<lagrangian model="off">
                      <particles_max_number>123456</particles_max_number>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
             'Could not get default values for max number of particles'


    def checkContinuousInjection(self):
        """Check whether the continuous injection could be set and get."""
        mdl = LagrangianModel(self.case)
        status = mdl.getContinuousInjection()

        assert status == mdl.defaultParticlesValues()['continuous_injection'], \
            'Could not get default values for continuous injection status'

        mdl.setContinuousInjection('on')
        doc = """<lagrangian model="off">
                      <continuous_injection status="on"/>
                 </lagrangian>"""

        assert mdl.node_lagr == self.xmlNodeFromString(doc), \
            'Could not set values for continuous injection status'


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
    print os.path.basename(__file__)
    runner = unittest.TextTestRunner()
    runner.run(suite())

#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
