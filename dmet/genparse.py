from dataparse import DataParse
import argparse
import numpy as np
import pandas as pd
import sys
import os
import re

"""
From first filename in the list from sys create a template for parsing other filenames
--> use this template to create DataFrame for inserting information from file
--> parse file, given several possible keywords and associated functions that fill into DataFrame

"""

class GenData(object):

#read data folder, use first file to create DataFrame format
    def __init__(self,data_folder,parse_info,target_ext = None):
        self.tar_ext = target_ext
        self.data_folder = data_folder
        print self.data_folder
        self.parse_info = parse_info
        self.parse_dict = parse_info
        self.parse_keys = self.parse_dict.keys()
        self.finish = False
        self.file_list = os.listdir(data_folder)
        print "Ignoring directories:"
        for item in self.file_list:
            if os.path.isdir(self.data_folder+'/'+item):
                print item
                self.file_list.remove(item)
        if self.tar_ext != None:
            deleted_extensions = []
            print "Target extension: " + self.tar_ext
            self.file_list = [item for item in self.file_list if item.split('.')[-1] == self.tar_ext]
        model_file = self.file_list[0]
        self.dfModelGen(model_file)
        self.dirParse()

    def dfModelGen(self, model_file):
        mod_name = model_file.split('/')[-1]
        print "Model file: " + model_file
        oterms, ndata = self.nameParse(mod_name)
        if oterms == []:
            self.outfile = self.data_folder.split('/')[1]+".csv"
        else:
            self.outfile = "_".join(oterms)+".csv"
        id_terms = [tup[0] for tup in ndata]
        vnames = []
        for key in self.parse_keys:
            vnames.append(self.parse_dict[key][0])
        id_terms += [item for sublist in vnames for item in sublist] #flatten vnames, then append
        self.df = pd.DataFrame(index=range(len(self.file_list)), columns=id_terms)
#        print self.df

    @staticmethod
    def nameParse(filename):
        out_terms = []  #For non-numeric information within the filename that may become part of datafile name
        name_data = []  #A list of 2-tuples to be filled; will specify file within df
        fname = filename.split('.')[0] #remove .out
        fname_split = fname.split('_')
        for sect in fname_split: #parse remaining sections individually
            num = re.search('\d', sect)
            if num == None:
                out_terms.append(sect)
            else:
                sect_numspl = re.findall('\d+|\D+',sect) #split up strings and numbers
                if len(sect_numspl) & 1:
                    out_terms.append(sect_numspl.pop())
                while sect_numspl != []:
                    alph = sect_numspl.pop(0)
                    num = sect_numspl.pop(0)
                    try:
                        num = int(num)
                        name_data.append((alph,num))
                    except ValueError: #Some files are named 4box, ect, and this should probs be fixed
                        out_terms.append(alph+num)
        return out_terms, name_data

#Will call parse for every file in directory, and fill the df with scraped data
    def dirParse(self):
        for i,dfile in enumerate(self.file_list):
            oterms, ndata = self.nameParse(dfile)
            extracted_data = self.parse(self.data_folder+'/'+dfile)
            #fill df row with ndata, then extracted_data
            for item in ndata:
                self.df[item[0]][i] = item[1]
            for item in extracted_data:
                self.df[item[0]][i] = item[1]

#main workhorse: iterates through file and parses
    def parse(self, infile):
        extracted_data = []
        with open(infile, 'r') as f:
            for line in f:
                #if match with parse list, call function
                for item in self.parse_keys:
                    if item in line:
                        vnames = self.parse_dict[item][0]
                        func_point = self.parse_dict[item][1]
                        assert type(vnames) == list #vnames should be a list, even if it's length 1
                        extracted_data += func_point(f,line,vnames) #append new data tuples to list
        return extracted_data


#inpParam will be specialized for the children of this class
    def inpParam(self,infile):
        pass

    def isComplete(self, infile = None,line=None):
        self.finish = True

    def splitNGrab(self,line,write_var,grab_index):
        spline = line.split()
        write_var = spline[grab_index]


def grabLast(f,line,vnames):
    data = line.split()[-1]
    #Expected return is a list of 2-tuples
    return [(vnames[0],data)]

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='Folder filled with data', type=str, required=True)
    parser.add_argument('-o', help='Output filename', type=str, required=True)
    parser.add_argument('-t', help='Target extension, if specific one desired', type=str, default=None)
    args = parser.parse_args()
    data_folder = args.f
    outfile = args.o

#    parse_info = [(["E"],"INFO - Final energy:",GenData.grabLast),(["Ecen"],"INFO - Center energy:",GenData.grabLast)]
    parse_info = {"INFO - Final energy:":(["E"],grabLast), "INFO - Center energy:":(["Ecen"],grabLast)}
    gendata = GenData(data_folder,parse_info,args.t)
    gendata.df.to_csv(outfile)












