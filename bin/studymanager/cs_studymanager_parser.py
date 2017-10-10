# -*- coding: utf-8 -*-

#-------------------------------------------------------------------------------

# This file is part of Code_Saturne, a general-purpose CFD tool.
#
# Copyright (C) 1998-2017 EDF S.A.
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
# Standard modules import
#-------------------------------------------------------------------------------

import os, sys, logging, re
from xml.dom import minidom
from xml.sax.handler import ContentHandler
from xml.sax import make_parser

#-------------------------------------------------------------------------------
# log config
#-------------------------------------------------------------------------------

logging.basicConfig()
log = logging.getLogger(__file__)
log.setLevel(logging.NOTSET)
#log.setLevel(logging.DEBUG)

#-------------------------------------------------------------------------------
# Checker of XML file syntax
#-------------------------------------------------------------------------------

def xmlChecker(filename):
    """Try to open the xml file, and return a message if an error occurs.

    @param filename name of the file of parameters ith its absolute path
    @return m error message
    """
    m = ""

    try:
        p = make_parser()
        p.setContentHandler(ContentHandler())
        p.parse(filename)
    except Exception as e:
        f = os.path.basename(filename)
        m = "%s file reading error. \n\n"\
            "This file is not in accordance with XML specifications.\n\n"\
            "The parsing syntax error is:\n\n%s" % (f, e)

    return m

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
        self.__repo = None
        self.__dest = None

        try:
            self.doc =  minidom.parse(XMLFileName)
        except:
            print("Error in the syntax of the xml.\n")
            msg =  xmlChecker(self.filename)
            if msg:
                print(msg)
            sys.exit(1)

        self.root = self.doc.firstChild

        if self.root.nodeName != "studymanager":
            print(XMLFileName + ": Wrong XML file. The Root markup is not studymanager.\n")
            sys.exit(1)

    #---------------------------------------------------------------------------

    def write(self):
        return self.doc.toxml()

    #---------------------------------------------------------------------------

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
            print("Error: in getDataFromNode several markup %s found." %  childName)
            sys.exit(1)

        return data

    #---------------------------------------------------------------------------

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

    #---------------------------------------------------------------------------

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
            return None
        elif len(l) == 1:
            return l.item(0)
        else:
            print("Error: in getChild several markup %s found." %  childname)
            sys.exit(1)

    #---------------------------------------------------------------------------

    def getChildren(self, father, childname):
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

    #---------------------------------------------------------------------------

    def getRepository(self):
        """
        @rtype: C{String}
        @return: repository directory of all studies.
        """
        if self.__repo == None: # Set the member
            elt = self.root.getElementsByTagName("repository");
            if len(elt) == 1:
                try:
                    self.__repo=elt.item(0).firstChild.data;
                except AttributeError:
                    msg="Parser.getRepository() >> No repository set.\n";
                    msg+="Add a path to the *.xml file or use the command "
                    msg+="line argument.\n"
                    sys.exit(msg);

        return self.__repo;

    #---------------------------------------------------------------------------

    def setRepository(self, repo_path):
        """
        @rtype: C{String}
        @param: repository directory of all studies.
        """

        self.__repo = os.path.expanduser(repo_path);
        self.__repo = os.path.abspath(self.__repo);

    #---------------------------------------------------------------------------

    def getDestination(self):
        """
        @rtype: C{String}
        @return: destination directory of all studies.
        """
        if self.__dest == None:
            self.__dest = self.getDataFromNode(self.root, "destination")

        return self.__dest

    #---------------------------------------------------------------------------

    def setDestination(self, dest_path):
        """
        @rtype: C{String}
        @param: destination directory of all studies.
        """
        self.__dest = os.path.expanduser(dest_path);
        self.__dest = os.path.abspath(self.__dest);

    #---------------------------------------------------------------------------

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
                    print("Several occurences of Study %s in xml file of parameters" % label)
                    sys.exit(1)
        return labels

    #---------------------------------------------------------------------------

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
            if label == l:
                if str(n.attributes["status"].value) != "on":
                    raise ValueError("Error: the getStudyNode method is used with the study %s turned off " % l)
                node = n
                break

            if not n:
                raise ValueError("Error: the getStudyNode does not found the node of the study %s" % l)

        return node

    #---------------------------------------------------------------------------

    def getStatusOnCasesLabels(self, l, attr=None):
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

    #---------------------------------------------------------------------------

    def getStatusOnCasesKeywords(self, l):
        """
        Read N_PROCS and USER_INPUT_FILES in:
            <study label='STUDY' status='on'>
                <case label='CASE1' run_id ="Grid 1" status='on' compute="on" post="on"/>
                <case label='CASE2' status='on' compute="on" post="on" compare="on"/>
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
                d['node']    = node
                d['label']   = str(node.attributes["label"].value)
                d['compute'] = str(node.attributes["compute"].value)
                d['post'] = str(node.attributes["post"].value)

                try:
                    d['compare'] = str(node.attributes["compare"].value)
                except:
                    d['compare'] = "on"

                try:
                    d['run_id'] = str(node.attributes["run_id"].value)
                except:
                    d['run_id'] = ""

                try:
                    d['n_procs'] = str(node.attributes["n_procs"].value)
                except:
                    d['n_procs'] = None

                try:
                    tags = str(node.attributes["tags"].value)
                    d['tags'] = [tag.strip() for tag in re.split(',', tags)]
                except:
                    d['tags'] = None

                for n in node.childNodes:
                    if n.nodeType == minidom.Node.ELEMENT_NODE and n.childNodes:
                        if n.tagName not in ("compare", "prepro", "script", "data"):
                            d[n.tagName] = n.childNodes[0].data
                data.append(d)

        return data

    #---------------------------------------------------------------------------

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

    #---------------------------------------------------------------------------

    def getCompare(self, caseNode):
        """
        Read:
            <study label='STUDY' status='on'>
                <case label='CASE1' status='on' compute="on" post="on">
                    <compare repo="" dest="" args='--section Pressure --threshold 1e-2' status="on"/>
                </case>
            </study>
        @type caseNode: C{DOM Element}
        @param caseNode: node of the current case
        @rtype: C{True} or C{False}, C{String}, C{String}, C{Float}, C{String}
        @return: if the cs_io_dump/compare markup exists, and value of the threshold
        """
        compare, nodes, threshold, args, repo, dest = [], [], [], [], [], []

        for node in caseNode.getElementsByTagName("compare"):
            try:
                if str(node.attributes["status"].value) == 'on':
                    compare.append(True)
                else:
                    compare.append(False)
            except:
                compare.append(True)
            nodes.append(node)
            repo.append(str(node.attributes["repo"].value))
            dest.append(str(node.attributes["dest"].value))
            try:
                args.append(str(node.attributes["args"].value))
            except:
                args.append(None)
            try:
                threshold.append(str(node.attributes["threshold"].value))
            except:
                threshold.append(None)

        return compare, nodes, repo, dest, threshold, args

    #---------------------------------------------------------------------------

    def getPrepro(self, caseNode):
        """
        Read:
            <study label='STUDY' status='on'>
                <case label='CASE1' status='on' compute="on" post="on">
                    <prepro label="script_pre.py" args="" status="on"/>
                </case>
            </study>
        @type caseNode: C{DOM Element}
        @param caseNode: node of the current case
        """
        prepro, label, nodes, args = [], [], [], []

        for node in caseNode.getElementsByTagName("prepro"):
            if str(node.attributes["status"].value) == 'on':
                prepro.append(True)
            else:
                prepro.append(False)

            label.append(str(node.attributes["label"].value))
            nodes.append(node)
            try:
                args.append(str(node.attributes["args"].value))
            except:
                args.append("")

        return prepro, label, nodes, args

    #---------------------------------------------------------------------------

    def getScript(self, caseNode):
        """
        Read:
            <study label='STUDY' status='on'>
                <case label='CASE1' status='on' compute="on" post="on">
                    <script label="script_post.py" args="" dest="20110216-2147" status="on"/>
                </case>
            </study>
        @type caseNode: C{DOM Element}
        @param caseNode: node of the current case
        """
        script, label, nodes, args, repo, dest = [], [], [], [], [], []

        for node in caseNode.getElementsByTagName("script"):
            if str(node.attributes["status"].value) == 'on':
                script.append(True)
            else:
                script.append(False)

            label.append(str(node.attributes["label"].value))
            nodes.append(node)
            try:
                args.append(str(node.attributes["args"].value))
            except:
                args.append("")
            try:
                repo.append(str(node.attributes["repo"].value))
            except:
                repo.append(None)
            try:
                dest.append(str(node.attributes["dest"].value))
            except:
                dest.append(None)

        return script, label, nodes, args, repo, dest

    #---------------------------------------------------------------------------

    def getResult(self, node):
        """
        @type node: C{DOM Element}
        @param node: node of the current case
        @rtype: C{List}, C{List}
        @return: C{List} of nodes <resu>, and C{List} of file names
        """
        f  = str(node.attributes["file"].value)

        if node.tagName == "data":
            plots = node.getElementsByTagName("plot")
        else:
            plots = []

        try:
            dest  = str(node.attributes["dest"].value)
        except:
            dest = None
        try:
            repo = str(node.attributes["repo"].value)
        except:
            repo = None

        return plots, f, dest, repo

    #---------------------------------------------------------------------------

    def getInput(self, node):
        """
        @type node: C{DOM Element}
        @param node: node of the current case
        @rtype: C{List}, C{List}
        @return: C{List} of nodes <input>, and C{List} of file names
        """
        f  = str(node.attributes["file"].value)

        try:
            dest  = str(node.attributes["dest"].value)
        except:
            dest = None
        try:
            repo = str(node.attributes["repo"].value)
        except:
            repo = None
        try:
            tex = str(node.attributes["tex"].value)
        except:
            tex = None

        return  f, dest, repo, tex

    #---------------------------------------------------------------------------

    def getProbes(self, node):
        """
        Read:
            <probes file='monitoring_pressure.dat' dest="" fig="3"/>

        @type node: C{DOM Element}
        @param node: node of the current case
        """
        f  = str(node.attributes["file"].value)
        try:
            dest  = str(node.attributes["dest"].value)
        except:
            dest = None
        try:
            fig  = str(node.attributes["fig"].value)
        except:
            fig = None

        return f, dest, fig

    #---------------------------------------------------------------------------

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
            filePath = node.attributes["path"].value

            if filePath == "":
              for root, dirs, fs in os.walk(os.path.join(self.getRepository(), l)):
                  if fileName in fs:
                      filePath = root
                      break
            else: # for Code_Saturne exp data are supposed to be in POST
              for root, dirs, fs in os.walk(os.path.join(self.getRepository(), l, 'POST', filePath)):
                  if fileName in fs:
                      filePath = root
                      break

            files.append(os.path.join(filePath, fileName))

        return nodes, files

    #---------------------------------------------------------------------------

    def getPostPro(self, l):
        """
        Read:
            <study label='STUDY' status='on'>
                <case label='CASE1' status='on' compute="on" post="on"/>
                <postpro label="script_post.py" args="" status="on">
                    <data file="profile.dat">
                        <plot fig="1" xcol="1" ycol="2" legend="Grid 1"/>
                    </data>
                </postpro>
            </study>

        Return the list of files of postpro.
        @type l: C{String}
        @param l: label of a study
        @rtype: C{List}
        @return: C{List} of list of nodes <postpro>
        """
        scripts, labels, nodes, args = [], [], [], []

        for node in self.getStudyNode(l).getElementsByTagName("postpro"):
            if str(node.attributes["status"].value) == 'on':
                scripts.append(True)
            else:
                scripts.append(False)

            labels.append(str(node.attributes["label"].value))
            nodes.append(node)
            try:
                args.append(str(node.attributes["args"].value))
            except:
                args.append("")

        return scripts, labels, nodes, args


    #---------------------------------------------------------------------------

    def getSubplots(self, studyLabel):
        return self.getStudyNode(studyLabel).getElementsByTagName("subplot")

    #---------------------------------------------------------------------------

    def getFigures(self, studyLabel):
        return self.getStudyNode(studyLabel).getElementsByTagName("figure")

    #---------------------------------------------------------------------------

    def getPltCommands(self, node):
        """
        """
        cmd = []
        for n in node.getElementsByTagName("plt_command"):
            if n.nodeType == minidom.Node.ELEMENT_NODE and n.childNodes:
                if n.tagName == "plt_command":
                    cmd.append(n.childNodes[0].data)
        return cmd

    #---------------------------------------------------------------------------

    def getAttributes(self, node):
        """
        Return a dictionary with attributes and value of a node.
        """
        d = {}
        for k in node.attributes.keys():
            d[k] = node.attributes[k].value
        return d

    #---------------------------------------------------------------------------

    def getAttribute(self, node, k, default = None):
        """
        Return a value of an attribute.
        """
        if k in node.attributes.keys():
            return node.attributes[k].value
        else:
            if default == None:
                raise ValueError("Error: attribute %s is mandatory!" % k)
            else:
                return default

    #---------------------------------------------------------------------------

    def getAttributeTuple(self, node, k, default = None):
        """
        Return a value of an attribute.
        """
        if k in node.attributes.keys():
            n = node.attributes[k].value
            return tuple(float(s) for s in n[1:-1].split(','))
        else:
            if default == None:
                raise ValueError("Error: attribute %s is mandatory!" % k)
            else:
                return default


#-------------------------------------------------------------------------------
