#!/bin/python3

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

import matplotlib.pyplot as plt
import os
import csv
import re

## class matching that in cs_user_profile.h but aggregating all dumped time steps

class Profile():

    def __init__(self,dirname):
        self.__dirname=dirname #directory storing the profiles results dumped by cs
        self.__name=os.path.basename(dirname) #name of the profile (last directory name)
        self.__field="" #cs field name
        self.__selection_criteria="" #cell selection criteria
        self.__dir_v="" #direction vector of the profile (char for graph legend
        self.__csv_file=dirname+os.sep+"results_profile.csv" #csv file storing the cs results
        self.__log_file=dirname+os.sep+"profile.log" #log file
        self.__n_layers=0 #number of layer of the profile
        self.__pos=[[]] #position 1D vector of the center of each layer
        self.__pos_n=[[]] #position 1D vector of the center of each layer normalized
        self.__mean_f=[[]] # list of vector of mean value per layer (1 per time step extracted)
        self.__sd_f=[[]] #list of vector of standard deviation for the field per layer (1 per time step extracted)
        self.__mean_f_n=[[]] #list of vector of mean normalized value of the field per layer (1 per time step extracted)
        self.__sd_f_n=[[]] #list of vector of standard devieation for the normalized field per layer (1 per time step extracted)
        self.__weigth=[[]] #list of vector of the weigth of each layer (1 per time step extracted)
        self.__time_step=[] #list of extracted time step
        self.__time=[] # list of extracted simulation time


    #############################################################
    #
    # Purpose of this function is to reset following structure member
    #       *  field, n_layers, dir_v, selection_criteria
    #
    #############################################

    def read_profil_setup(self):

        with open(self.__log_file,'r') as f:
            StringFile = f.read()
            find = re.findall('field : .+',StringFile)
            if find:
                find=find[0].replace("field : ",'')
            else:
                find="not found"

            self.__field = find #update field in structure

            find = re.findall('normal profile direction : .+',StringFile)
            if find:
                find=find[0].replace("normal profile direction : ",'')
            else:
                find="not found"

            self.__dir_v=find

            find = re.findall('number of layers : \d+',StringFile)
            if find:
                find=int(find[0].replace("number of layers : ",''))
            else:
                find="not found"

            self.__n_layers=find

            find = re.findall('cells selection : .+',StringFile)
            if find:
                find=find[0].replace("cells selection : ",'')
            else:
                find="not found"

            self.__selection_criteria=find

    #############################################################
    #
    # Purpose of this function is to reset following structure member
    #       *  pos, weigth,mean_f, sd_f, pos_n, mean_f_n, sd_f_n
    #
    #############################################

    def reset_results(self):

        # create a list of profile properties in the same order than
        # _dump_profile_values_csv of cs_user_profile.c
        profile_var = [self.__pos,
                       self.__weigth,
                       self.__mean_f,
                       self.__sd_f,
                       self.__pos_n,
                       self.__mean_f_n,
                       self.__sd_f_n]

        n_var=len(profile_var)

        for var_id in range(len(profile_var)): #remove last empyt list added
            profile_var[var_id]=[[]]

        self.__time_step=[]
        self.__time=[]

    #######################################################################
    #
    # Purpose of this function is to populate following parameters of the structure
    #   With the values dumped by Code Saturne in the csv file
    #       *  pos, weigth,mean_f, sd_f, pos_n, mean_f_n, sd_f_n
    #
    #####################################################################

    def read_csv_results(self):

        #first reset structure value (function should be called once per profile)
        self.reset_results()

        with open(self.__csv_file,'r',newline='') as f_in:
            sheet_csv = csv.reader(f_in)
            sheet=[]
            for row in sheet_csv:
                sheet.append(row) #each row is stored as a list

            del sheet_csv

            header=sheet[0]

            # create a list of profile properties in the same order than
            # _dump_profile_values_csv of cs_user_profile.c
            profile_var=[self.__pos,
                         self.__weigth,
                         self.__mean_f,
                         self.__sd_f,
                         self.__pos_n,
                         self.__mean_f_n,
                         self.__sd_f_n]

            n_var=len(profile_var)

            for t_id in range(1,len(sheet)):
                self.__time_step.append(int(sheet[t_id][0]))
                self.__time.append(float(sheet[t_id][1]))
                for var_id in range(len(profile_var)):
                    for layer_id in range(self.__n_layers):
                        profile_var[var_id][t_id-1].append(float(sheet[t_id][2+layer_id*n_var+var_id]))

                    profile_var[var_id].append([])

            for var_id in range(len(profile_var)): #remove last empyt list added
                del profile_var[var_id][-1]

    #######################################################################
    #
    # Purpose of this function is plot graph of profiles with uncertainties
    #
    #      Graphs generated will be dumped in figures directory
    #
    #####################################################################

    def plot_graphs(self, dpi=150):

        for time_step in range(len(self.__time_step)):
            fig=plt.figure(num=self.__name+' '+str(self.__time_step[time_step]))
            ax = fig.add_subplot(1, 1, 1)
            ax.errorbar(self.__pos[time_step],
                        self.__mean_f[time_step],
                        yerr=self.__sd_f[time_step],
                        fmt='--o',
                        capsize=5)  #errorbar layout
            ax.set_title(self.__name +'\n' + 'dir : '+self.__dir_v+'\n'+self.__selection_criteria)
            ax.set_xlabel('distance [m]')
            ax.set_ylabel(self.__field)

        path_out=self.__dirname+os.sep+'figures'

        try:
            os.makedirs(path_out) #create output directory if not exist
        except:
            pass

        #dump the graph in a dedicated directory
        figs = [plt.figure(num=name) for name in plt.get_figlabels()]
        figlabels = plt.get_figlabels()
        i = 0
        for fig in figs:
            fig.savefig(path_out+os.sep+figlabels[i]+'.png',
                        dpi=dpi,bbox_inches='tight')
            i+=1

#-------------------------------------------------------------------------------

if __name__=='__main__':

    dir_profiles = "../profiles"
    # List all directories in profiles (each is a profile)
    for profile in os.listdir(dir_profiles):
        # If the item is a directory, it is a profile
        if os.path.isdir(dir_profiles+os.sep+profile):
            Profile_test=Profile(dir_profiles+os.sep+profile) # Create profile class
            Profile_test.read_profil_setup()
            Profile_test.read_csv_results()
            Profile_test.plot_graphs()
