#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2019 EDF S.A.
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

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import sys

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# Tests suite
#-------------------------------------------------------------------------------


def starttest1():
    from Base.XMLengine import runTest
    runTest()

def starttest2():
    from Base.XMLvariables import runTest
    runTest()

def starttest3():
    from Base.XMLmodel import runTest
    runTest()

def starttest4():
    from Pages.IdentityAndPathesModel import runTest
    runTest()

def starttest5():
    from Pages.SolutionDomainModel import runTest
    runTest()

def starttest8():
    from Pages.MobileMeshModel import runTest
    runTest()

def starttest9():
    from Pages.TurbulenceModel import runTest
    runTest()

def starttest10():
    from Pages.CoalCombustionModel import runTest
    runTest()

def starttest11():
    from Pages.CurrentSpeciesModel import runTest
    runTest()

def starttest12():
    from Pages.ThermalScalarModel import runTest
    runTest()

def starttest13():
    from Pages.ThermalRadiationModel import runTest
    runTest()

def starttest14():
    from Pages.InitializationModel import runTest
    runTest()

def starttest15():
    from Pages.ReferenceValuesModel import runTest
    runTest()

def starttest16():
    from Pages.FluidCharacteristicsModel import runTest
    runTest()

def starttest17():
    from Pages.BodyForcesModel import runTest
    runTest()

def starttest18():
    from Pages.DefineUserScalarsModel import runTest
    runTest()

def starttest19():
    from Pages.Boundary import runTest, runTest2, runTest3, runTest4, runTest5
    from Pages.Boundary import runTest6, runTest7
    runTest()
    runTest2()
    runTest3()
    runTest4()
    runTest5()
    runTest6()
    runTest7()

def starttest23():
    from Pages.TimeAveragesModel import runTest
    runTest()

def starttest24():
    from Pages.SteadyManagementModel import runTest
    runTest()

def starttest25():
    from Pages.TimeStepModel import runTest
    from Pages.TimeStepModel import runTest2
    runTest()
    runTest2()

def starttest26():
    from Pages.OutputControlModel import runTest
    runTest()

def starttest27():
    from Pages.OutputVolumicVariablesModel import runTest
    runTest()

def starttest28():
    from Pages.OutputSurfacicVariablesModel import runTest
    runTest()

def starttest29():
    from Pages.ProfilesModel import runTest
    runTest()

def starttest30():
    from Pages.NumericalParamEquationModel import runTest
    runTest()

def starttest31():
    from Pages.NumericalParamGlobalModel import runTest
    runTest()

def starttest33():
    from Pages.StartRestartModel import runTest
    runTest()

def starttest34():
    from Pages.BatchRunningModel import runTest
    runTest()

def starttest35():
    from Pages.LagrangianModel import runTest
    runTest()

def starttest36():
    from Pages.ElectricalModel import runTest
    runTest()

def starttest37():
    from Pages.GasCombustionModel import runTest
    runTest()

def starttest45():
    from Pages.LocalizationModel import runTest
    runTest()

def starttest46():
    from Pages.HeadLossesModel import runTest
    runTest()

def starttest47():
    from Pages.FluidStructureInteractionModel import runTest
    runTest()

def starttest48():
    from Pages.AtmosphericFlowsModel import runTest
    runTest()

if __name__ == '__main__':

    print('STARTING GUI UNIT TESTS')
    #sys.exit(0)
    starttest1()
    starttest2()
    starttest3()
    starttest4()
    starttest5()
    starttest8()
    starttest9()
    starttest10()
    starttest11()
    starttest12()
    starttest13()
    starttest14()
    starttest15()
    starttest16()
    starttest17()
    starttest18()
    starttest19()
    starttest23()
    starttest24()
    starttest25()
    starttest26()
    starttest27()
    starttest28()
    starttest29()
    starttest30()
    starttest31()
    starttest33()
    starttest34()
    starttest35()
    #starttest36()
    #starttest37()
    #starttest38()
    #starttest39()
    #starttest40()
    #starttest41()
    #starttest42()
    #starttest43()
    #starttest44()
    starttest45()
##    starttest46()
    starttest47()
    starttest48()


#-------------------------------------------------------------------------------
# End
#-------------------------------------------------------------------------------
