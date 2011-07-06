#============================================================================
#
#     This file is part of CFDSTUDY the plug-in for Salome
#     of Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2010 EDF S.A., France
#
#     contact: saturne-support@edf.fr
#
#     CFDSTUDY is free software; you can redistribute it
#     and/or modify it under the terms of the GNU General Public License
#     as published by the Free Software Foundation; either version 2 of
#     the License, or (at your option) any later version.
#
#     CFDSTUDY is distributed in the hope that it will be
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
#============================================================================

import CFDSTUDY_ORB__POA
import SALOME_ComponentPy
import SALOME_DriverPy

class CFDSTUDY(CFDSTUDY_ORB__POA.CFDSTUDY_Gen,
            SALOME_ComponentPy.SALOME_ComponentPy_i,
            SALOME_DriverPy.SALOME_DriverPy_i):
    """
        Pour etre un composant SALOME cette classe Python
        doit avoir le nom du composant et heriter de la
        classe CFDSTUDY_Gen issue de la compilation de l'idl
        par omniidl et de la classe SALOME_ComponentPy_i
        qui porte les services generaux d'un composant SALOME
    """
    def __init__ ( self, orb, poa, contID, containerName, instanceName,
                   interfaceName ):
        #print "CFDSTUDY.__init__: ", containerName, ';', instanceName
        SALOME_ComponentPy.SALOME_ComponentPy_i.__init__(self, orb, poa,
                    contID, containerName, instanceName, interfaceName, 0)
        SALOME_DriverPy.SALOME_DriverPy_i.__init__(self, interfaceName)
        # On stocke dans l'attribut _naming_service, une reference sur
        # le Naming Service CORBA
        self._naming_service = SALOME_ComponentPy.SALOME_NamingServicePy_i( self._orb )
