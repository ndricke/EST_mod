from dataparse import DataParse
import re
import sys


class DmetParse(DataParse):
    def __init__(self,infile):
        self.m = None
        self.l = None
        self.E = None
        self.E_cen = None

        self.parse_dict = {}
        self.parse_dict["INFO - Final energy:"] = self.getE
        self.parse_dict["INFO - Center energy:"] = self.getEc
        self.nameParse(infile)
        super(DmetParse,self).__init__(infile)

    def getE(self,f,line):
        spline = line.split()
        self.E = spline[-1]

    def getEc(self,f,line):
        spline = line.split()
        self.E_cen = spline[-1]

    def nameParseGeneric(self,infile):
        spl = infile.split('/')[-1]
        spl = spl.split('.')[0]
        data_str = spl.split('_')

        

    def nameParse(self,infile):
        spl = infile.split('/')[-1]
        spl = spl.split('.')[0]
        data_str = spl.split('_')
        self.lattype = data_str[0]
        self.l = data_str[1]
        self.p = data_str[2]
        self.m = self.grabTrailNum(data_str,'m')

#This algorithm will run into trouble if you repeat search keys or have overlap in their names
    def grabTrailNum(self,alphnum_list,search_key):
        sk_len = len(search_key)
        for item in alphnum_list:
            if search_key in item:
                return item[sk_len:]

    def grabConstr(self,data_string):
        for item in data_string:
            if 'c' in item:
                return item

    def dataWrite(self,outfile):
        dmet_data = ','.join([self.lattype,self.l,self.p,self.m,self.E,self.E_cen])
        with open(outfile,'a') as f:
            f.write(dmet_data)
            f.write('\n')


if __name__ == "__main__":
    dm_data = DmetParse('n8_m4_U4_LR2.out')
    dm_data.dataWrite('test_out.csv')
