import sys


class DataParse(object):
    def __init__(self,infile):
        self.finish = False
        self.parse(infile)

#main workhorse: iterates through file and parses
    def parse(self, infile):
        parse_list = self.parse_dict.keys()
        with open(infile, 'r') as f:
            for line in f:
                #if match with parse list, call function
                for item in parse_list:
                    if item in line:
                        func_point = self.parse_dict[item]
                        func_point(f,line)

#inpParam will be specialized for the children of this class
    def inpParam(self,infile):
        pass

    def isComplete(self, infile = None,line=None):
        self.finish = True

    def splitNGrab(self,line,write_var,grab_index):
        spline = line.split()
        write_var = spline[grab_index]
