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
This module defines the following classes and functions:
- displaySelectedPage
"""

#-------------------------------------------------------------------------------
# Library modules import
#-------------------------------------------------------------------------------

import os, sys, logging

#-------------------------------------------------------------------------------
# Application modules import
#-------------------------------------------------------------------------------

from code_saturne.gui.base.QtCore  import QCoreApplication
from code_saturne.model.Common import *

#-------------------------------------------------------------------------------
# displaySelectedPage direct to the good page with its name
#-------------------------------------------------------------------------------

def displaySelectedPage(page_name, root, case, stbar=None, tree=None):
    """
    This function enables to display a new page when the TreeNavigator
    send the order.
    """
    # 'win' is the frame-support of the Pages
    # 'thisPage' is the instance of classes which create thePages
    # 'page_name' is the name of the page
    #
    if page_name == tr("Paths"):
        import code_saturne.gui.studymanager_gui.PathesView as Page
        thisPage = Page.PathesView(root, case)

    elif page_name == tr("Manage cases"):
        import code_saturne.gui.studymanager_gui.ManageCasesView as Page
        thisPage = Page.ManageCasesView(root, case)

    elif page_name == tr("Define plotter"):
        import code_saturne.gui.studymanager_gui.ManagePlotterView as Page
        thisPage = Page.ManagePlotterView(root, case)

    else:
        msg = tr("Warning: the corresponding Page %s doesn't exist!") % page_name
        print(msg)
        import code_saturne.gui.case.WelcomeView as Page
        thisPage = Page.WelcomeView()

    case['current_page'] = str(page_name)

    return thisPage


def tr(text):
    """
    Translation
    """

    # Note that the matching tree entries are declared in BrowserView,
    # so the translation context requires this.
    # Merging this into BrowserView (and making the former more dynamic)
    # and using function pointers would go a long way towards
    # making the code more modular an easier to navigate thna this mess.

    return QCoreApplication.translate('BrowserView', text)


#-------------------------------------------------------------------------------
# End of Toolbox
#-------------------------------------------------------------------------------
