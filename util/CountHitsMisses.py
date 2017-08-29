#!/usr/bin/env python
# encoding: utf-8
'''

@author:     cjustin
@contact:    cjustin@bcgsc.ca
'''

# from optparse import OptionParser
import sys

class HitCounter:
    
    def __init__(self, path):
        """Constructor"""
        self._paths = path
    
    def count(self):
        """Counts matches a match is correct if the filename"""
        #get IDs
        
        matchCount = 0
        total = 0
        
        for i in self._paths:
            #read files
            print i
            with open(i, 'r') as tempFH:
                for line in tempFH:
                    dat = line.split("\t")
                    id = dat[0]
                    if id == "0":
                        id = "rand"
                    readID = dat[1][0:len(id)]
                    if id == readID:
                        matchCount += 1
                    total += 1
        
        print("matches: " + str(matchCount))
        print("total: " + str(total))
        

    
if __name__ == '__main__':

#     parser = OptionParser()
#     parser.add_option("-i", "--input", dest="input", metavar="INPUT",
#                       help="input files to process")
#                 
#     (options, args) = parser.parse_args()
        
#     if options.input:
        runner = HitCounter(sys.argv[1:])
        runner.count()
#     else:
#         print ('ERROR: Missing Required Options. Use -h for help')
