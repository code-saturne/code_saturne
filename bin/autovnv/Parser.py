#-------------------------------------------------------------------------------
#
#     This file is part of the Code_Saturne Scripts, element of the
#     Code_Saturne CFD tool.
#
#     Copyright (C) 1998-2011 EDF S.A., France
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

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, sys, logging
from xml.dom import minidom

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger(__file__)
log.setLevel(logging.NOTSET)
#log.setLevel(logging.DEBUG)

#-------------------------------------------------------------------------------
# Reader of the XML file
#-------------------------------------------------------------------------------

class Parser(object):
    """ Parser -- class to parse XML file."""

    def __init__ (self, XMLFileName):
        """
        Constructor of the XML reader.
        @type XMLFileName: C{String}
        @param XMLFileName: name of the xml file
        """
        self.filename = XMLFileName
        try:
            self.doc =  minidom.parse(XMLFileName)
        except:
            print "No file or syntax error"
            sys.exit(1)

        self.root = self.doc.firstChild

        if self.root.nodeName != "autoverif":
            print XMLFileName + ": Wrong XML file."
            sys.exit(1)


    def write(self):
        return self.doc.toxml()


    def __checkDirectory(self, d):
        """
        Test if the given directory exists.
        @type d: C{String}
        @param d: name of the directory
        @rtype: C{String}
        @return: the directory
        """
        log.debug("__checkDirectory(): %s" % d)
        dir = os.path.abspath(d)
        if os.path.isdir(dir):
            log.debug("__checkDirectory(): %s" % dir)
            return dir
        else:
            print "Directory: %s does not exists." % dir
            sys.exit(1)


    def getDataFromNode(self, father, childName):
        """
        @type father: C{DOM element}
        @param father: father node
        @type childName: C{String}
        @param childName: label of a markup, child of the father node
        @rtype: C{String}
        @return: data from the markup I{childName}
        """
        data = None
        l = father.getElementsByTagName(childName)

        if len(l) == 1:
            current = l.item(0)
            data = current.firstChild.data
        else:
            print "Error: in getDataFromNode several markup %s found." %  childName
            sys.exit(1)

        return data


    def addChild(self, father, childname):
        """
        Add a new node as a child of the father node.
        @type father: C{DOM element}
        @param father: father node
        @type childName: C{String}
        @param childName: label of a markup, child of the father node
        @rtype: C{DOM element}
        @return: new node I{childName}
        """
        return father.appendChild(self.doc.createElement(childname))


    def getChild(self, father, childname):
        """
        Return a node as a child of the father node.
        @type father: C{DOM element}
        @param father: father node
        @type childName: C{String}
        @param childName: label of a markup, child of the father node
        @rtype: C{DOM element}
        @return: new node I{childName}
        """
        l = father.getElementsByTagName(childname)

        if len(l) == 0:
            el = self.doc.createElement(childname)
            return father.appendChild(el)
        elif len(l) == 1:
            return l.item(0)
        else:
            print "Error: in getChild several markup %s found." %  childname
            sys.exit(1)


    def getChilds(self, father, childname):
        """
        Return a list of child of the father node.
        @type father: C{DOM element}
        @param father: father node
        @type childName: C{String}
        @param childName: label of a markup, childs of the father node
        @rtype: C{List} of C{DOM element}
        @return: list of nodes I{childName}
        """
        return father.getElementsByTagName(childname)


    def getRepository(self):
        """
        @rtype: C{String}
        @return: repository directory of all studies.
        """
        return self.__checkDirectory(self.getDataFromNode(self.root, "repository"))


    def getDestination(self):
        """
        @rtype: C{String}
        @return: destination directory of all studies.
        """
        return self.getDataFromNode(self.root, "destination")


    def getStudiesLabel(self):
        """
        Read:
            <study label='STUDY' status='on'>
                <case label='CASE' status='on' compute="on" post="on"/>
            </study>

        @rtype: C{List} of C{String}
        @return: list of the label of all studies.
        """
        labels = []

        for node in self.root.getElementsByTagName("study"):
            label  = str(node.attributes["label"].value)
            status = str(node.attributes["status"].value)
            if status == 'on':
                if label not in labels:
                    labels.append(label)
                else:
                    print "Study: %s is repeted in the xml file of paramaters" % label
                    sys.exit(1)
        return labels


    def getStudyNode(self, l):
        """
        Read:
            <study label='STUDY' status='on'>
                <case label='CASE' status='on' compute="on" post="on"/>
            </study>

        @type l: C{String}
        @param l: label of a study
        @rtype: C{DOM element}
        @return: node of the xml document for study I{l}
        """
        node = None

        for n in self.root.getElementsByTagName("study"):
            label = str(n.attributes["label"].value)
            if n.attributes["status"].value != "on":
                raise ValueError, "Error: the getStudyNode method is used with the study %s turned off " % l

            if label == l:
                node = n
                break

            if not n:
                raise ValueError, "Error: the getStudyNode does not found the node of the study %s" % l

        return node


    def getCasesLabel(self, l, attr=None):
        """
        Read:
            <study label='STUDY' status='on'>
                <case label='CASE' status='on' compute="on" post="on"/>
            </study>

        @type l: C{String}
        @param l: label of a study
        @type attr: C{String}
        @param attr: attribute I{compute} or I{post}
        @rtype: C{List} of C{String}
        @return: list of cases for study I{l}
        """
        labels = []

        for node in self.getStudyNode(l).getElementsByTagName("case"):
            status = str(node.attributes["status"].value)
            if status == 'on':
                labels.append(str(node.attributes["label"].value))

                if attr and str(node.attributes[attr].value) != 'on':
                    labels.pop()

        return labels


    def getCasesKeywords(self, l):
        """
        Read N_PROCS and USER_INPUT_FILES in:
            <study label='STUDY' status='on'>
                <case label='CASE1' status='on' compute="on" post="on">
                    <N_PROCS>2</N_PROCS>
                </case>
                <case label='CASE2' status='on' compute="on" post="on">
                    <USER_INPUT_FILES>['data1', 'data2']</USER_INPUT_FILES>
                </case>
                <case label='CASE2' status='on' compute="on" post="on">
                    <USER_INPUT_FILES>['data3', 'data4']</USER_INPUT_FILES>
                </case>
            </study>
        @type l: C{String}
        @param l: label of a study
        @rtype: C{Dictionary}
        @return: keywords and value.
        """
        data = []

        for node in self.getStudyNode(l).getElementsByTagName("case"):
            if str(node.attributes["status"].value) == 'on':
                d = {}
                d['node']  = node
                d['label'] = str(node.attributes["label"].value)
                d['compute'] = str(node.attributes["compute"].value)
                d['post']    = str(node.attributes["post"].value)
                for n in node.childNodes:
                    if n.nodeType == minidom.Node.ELEMENT_NODE and n.childNodes:
                        if n.tagName not in ("compare", "script", "data"):
                            d[n.tagName] = n.childNodes[0].data

                data.append(d)
        return data


    def setAttribute(self, node, attr, v):
        """
        Change:
            <study label='STUDY' status='on'>
                <case label='CASE1' status='on' compute="on" post="on"/>
            </study>
        To:
            <study label='STUDY' status='on'>
                <case label='CASE1' status='on' compute="off" post="on"/>
            </study>

        @type l: C{DOM Element}
        @param l: node of the I{attribute} to change
        @type attr: C{String}
        @param attr: attribute I{compute} or I{post}
        @type v: C{String}
        @param v: value of the attribute
        """
        node.attributes[attr].value = v

        f = os.path.join(self.getDestination(), self.filename)
        writer = open(f, mode="w" )
        self.doc.writexml(writer)
        writer.close()


    def getCompare(self, caseNode):
        """
        Read:
            <study label='STUDY' status='on'>
                <case label='CASE1' status='on' compute="on" post="on">
                    <compare repo="20110216-2047" dest="20110216-2147" threshold='1e-4'/>
                </case>
            </study>
        @type caseNode: C{DOM Element}
        @param caseNode: node of the current case
        @rtype: C{True} or C{False}, C{String}, C{String}, C{Float}
        @return: if the cs_io_dump/compare markup exists, and value of the threshold
        """
        nodes = caseNode.getElementsByTagName("compare")

        compare    = False
        threshold  = None
        repo       = None
        dest       = None

        if nodes:
            compare = True
            try:
                repo      = str(nodes[0].attributes["repo"].value)
                dest      = str(nodes[0].attributes["dest"].value)
                threshold = str(nodes[0].attributes["threshold"].value)
            except:
                repo      = None
                dest      = None
                threshold = None

        return compare, repo, dest, threshold


    def getScript(self, caseNode):
        """
        Read:
            <study label='STUDY' status='on'>
                <case label='CASE1' status='on' compute="on" post="on">
                    <script label="script_post.py" args="" repo="20110216-2047" dest="20110216-2147" status="on"/>
                </case>
            </study>
        @type caseNode: C{DOM Element}
        @param caseNode: node of the current case
        """
        script = False
        label, args, repo, dest = [], [], [], []

        for node in caseNode.getElementsByTagName("script"):
            if str(node.attributes["status"].value) == 'on':
                script = True
                label.append(str(node.attributes["label"].value))
                try:
                    args.append(str(node.attributes["args"].value))
                except:
                    args.append("")
                try:
                    repo.append(str(node.attributes["repo"].value))
                except:
                    repo.append("")
                try:
                    dest.append(str(node.attributes["dest"].value))
                except:
                    dest.append("")

        return script, label, args, repo, dest


    def getResult(self, node):
        """
        Read:
            <data file='profil1.dat'>
                <plot fig='1' xcol='1' ycol='2' legend='U'/>
                <plot fig='2' xcol='1' ycol='3' legend='V'/>
            </data>

        @type node: C{DOM Element}
        @param node: node of the current case
        @rtype: C{List}, C{List}
        @return: C{List} of nodes <plot>, and C{List} of file names
        """
        plots = node.getElementsByTagName("plot")
        dest  = str(node.attributes["dest"].value)
        repo  = str(node.attributes["repo"].value)
        file  = str(node.attributes["file"].value)

        return plots, file, dest, repo


    def getMeasurement(self, l):
        """
        Return the list of files of measurements.
        @type l: C{String}
        @param l: label of a study
        @rtype: C{List}
        @return: C{List} of list of nodes <plot>, and C{List} of files of measurements
        """
        nodes = []
        files = []

        for node in self.getStudyNode(l).getElementsByTagName("measurement"):
            nodes.append(node.getElementsByTagName("plot"))
            fileName = node.attributes["file"].value
            try:
                filePath = node.attributes["path"].value
                files.append(os.path.join(filePath, fileName))
            except:
                files.append(fileName)

        return nodes, files


    def getPlots(self, studyLabel):
        Idfig = []
        nodes = []
        nbfigure = 0
        for node in self.getStudyNode(studyLabel).getElementsByTagName("subplot"):
            figvalue = int(node.attributes["fig"].value)
            nbfigure += 1
            if figvalue not in Idfig:
                Idfig.append(figvalue)
                nodes.append(node)
        return nodes


    def getFigures(self, studyLabel):
        return self.getStudyNode(studyLabel).getElementsByTagName("figure")


    def getPltCommands(self, node):
        """
        """
        cmd = []
        for n in node.getElementsByTagName("plt_command"):
            if n.nodeType == minidom.Node.ELEMENT_NODE and n.childNodes:
                if n.tagName == "plt_command":
                    cmd.append(n.childNodes[0].data)
        return cmd

#-------------------------------------------------------------------------------
