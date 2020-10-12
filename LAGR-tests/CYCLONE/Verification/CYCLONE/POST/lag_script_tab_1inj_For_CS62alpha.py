#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Standard modules import
#-------------------------------------------------------------------------------

import os, sys
import string
from optparse import OptionParser
from math import sqrt
from pylab import *
import glob
from numpy import *


#-------------------------------------------------------------------------------

def process_cmd_line(argv):
    """
    Processes the passed command line arguments.
    """
    parser = OptionParser(usage="usage: %prog [options]")

    parser.add_option("-r", "--repo", dest="repo", type="string",
                      help="Directory of the result in the repository")

    parser.add_option("-d", "--dest", dest="dest", type="string",
                      help="Directory of the result in the destination")

    parser.add_option("-p", "--post", dest="post", type="string",
                      help="Directory of the post in the destination")

    (options, args) = parser.parse_args(argv)

    return options


#-------------------------------------------------------------------------------
  
def main(options):

    case = options.dest.split("/")[-3]
    print(case)
    post = "../../../POST"
    pathpath = options.dest
  
    #name of file (listing or run_solver.log)
    filepath = os.path.join(options.dest, "run_solver.log" )
    #total number of particles class
    number_of_class = 7
    #type name in CS of the cyclone inlet
    inlet_BC = "BC_1"
    #type name in CS of the cyclone outlet up
    outlet_up = "BC_2"
    #type name in CS of the cyclone outlet bottom
    outlet_bot = "BC_4"
    #time to start the average to estimate the proportion of particles going through the up or the   bottom outlet
    start_average = 0.79982
    
    # Get a list of all the file paths that have the type class_*_type_BC_*.dat" from in specified directory
    fileList = glob.glob(pathpath+"/class_*_type_BC_*.dat")
    # Iterate over the list of filepaths & remove each file.
    for filePath in fileList:
        try:
            os.remove(filePath)
        except:
            print("Error while deleting file : ", filePath)
  
    FIRST_INJECTION = ones(number_of_class, bool)
    #initialization of the average flow rate for all the classes of particles
    average_inlet_flow_rate = zeros(number_of_class, float)
  
    #initialization of the average relative flow rate from start_average to the end
    average_relative_flow_rate = zeros([2,number_of_class], float)
  
    #number of lines in the listing where there is a particle of the current class going through the   outlet up
    cnt_outlet_up = zeros(number_of_class, int)
    #number of lines in the listing where there is a particle of the current class going through the outlet bot
    cnt_outlet_bot = zeros(number_of_class, int)

    time_step = start_average

    FLAG_MAIN_CALCULATION = False

    with open(filepath) as fp:  
        line = fp.readline()
        cnt = 1
        while line:
            cnt += 1
            if ("MAIN CALCULATION" in line): FLAG_MAIN_CALCULATION=True
            if ("INSTANT" in line) and FLAG_MAIN_CALCULATION:
                line_split = line.split()
                time = line_split[1]
                iteration = line_split[5]
  
            line = fp.readline()
  
            if ("Zone  Class  Mass flow rate(kg/s)      Name (type)" in line) and FLAG_MAIN_CALCULATION:
                while not("--") in line:
                    line = fp.readline()
                    line_split = line.split()
  
                    if "BC_" in line:
                        type_bnd = line_split[2]
                        for j in range(number_of_class):
                            F = open(pathpath+"/class_"+str(j+1)+"_type_"+type_bnd+".dat","a")
                            F.close()
  
                    elif not("BC_") or not("--") in line:
                        classe = line_split[0]
                        debit_classe = line_split[1]
                        debit_relative_inlet = 0.0
                        if not(FIRST_INJECTION[int(classe)-1]) : debit_relative_inlet = -float(debit_classe)*(float(time)-time_step)/average_inlet_flow_rate[int(classe)-1]
                        if float(time) > start_average :
                            if type_bnd == outlet_up :
                                cnt_outlet_up[int(classe)-1] += 1
                                average_relative_flow_rate[0,int(classe)-1] = debit_relative_inlet + average_relative_flow_rate[0,int(classe)-1]
                            if type_bnd == outlet_bot :
                                cnt_outlet_bot[int(classe)-1] += 1
                                average_relative_flow_rate[1,int(classe)-1] = debit_relative_inlet + average_relative_flow_rate[1,int(classe)-1]
                        F = open(pathpath+"/class_"+classe+"_type_"+type_bnd+".dat","a")
                        s = str("time= "+time+" iteration= "+iteration+" classe: "+ classe+" debit_classe= "+ debit_classe +" relative_debit= "+ str(debit_relative_inlet) + " type_bnd= "+type_bnd+"\n")
                        F.write(s)
                        F.close()
                        if FIRST_INJECTION[int(classe)-1] :
                            FIRST_INJECTION[int(classe)-1] = False
                            average_inlet_flow_rate[int(classe)-1] = float(debit_classe)*(float(time)-time_step)
                            print(average_inlet_flow_rate)
                time_step = float(time)

    relative_flow_rate = zeros([2,number_of_class], float)
    relative_flow_rate[:,:] = average_relative_flow_rate[:,:]*100.0
    print("total time computed for the average=")
    print(time_step-start_average)
    print("%, outlet_up")
    print(relative_flow_rate[0,:])
    print("%, outlet_bot")
    print(relative_flow_rate[1,:])
    print("%, outlet_up+outlet_bot")
    print(relative_flow_rate[0,:]+relative_flow_rate[1,:])

    Tex_output_file = os.path.join(options.dest, post, case+"_particles_through_outlets.tex" )
    Output_File = open(Tex_output_file,"w")
    Output_File.write("Particles of the class 1 to "+ str(number_of_class) +" going through the Outlet up, % \n")
    Output_File.write(" \\\\\n".join([" & ".join(map(str,relative_flow_rate[0,:]))]))
    Output_File.write("\n")
    Output_File.write("\n")
    Output_File.write("Particles of the class 1 to "+ str(number_of_class) +" going through the Outlet bottom, % \n")
    Output_File.write(" \\\\\n".join([" & ".join(map(str,relative_flow_rate[1,:]))]))
    Output_File.write("\n")
    Output_File.write("\n")
    Output_File.write("Outlet up + Outlet bottom for the particles of the class 1 to "+ str(number_of_class) +" , % \n")
    Output_File.write(" \\\\\n".join([" & ".join(map(str,relative_flow_rate[0,:]+relative_flow_rate[1,:]))]))
    Output_File.write("\n")
    Output_File.close()

    Output_Tab_File = open(pathpath+"/particles_through_bottom_outlet.dat","w")
    for j in range(number_of_class):
        Output_Tab_File.write(str(j+1) + " " + str(relative_flow_rate[1,j])+"\n")
    
    Output_Tab_File.close()
#-------------------------------------------------------------------------------

if __name__ == '__main__':
    options = process_cmd_line(sys.argv[1:])
    main(options)


