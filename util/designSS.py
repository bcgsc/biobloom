#!/usr/bin/python
import sys, getopt, random
#  designSS.py
#  SSDesign
#
#  Created by Birol on 2016-03-14.
#  Copyright (c) 2016 Birol. All rights reserved.

# usage:
# designSS.py -n number_of_seeds -m allowed_number_of_misses -l seed_length
# rules
# m < n
# l > n
def usage():
    scriptName=sys.argv[0]
    print("SYNOPSIS")
    print("\t%s [-n number_of_seeds] [-m allowed_number_of_misses] [-l seed_length]" %scriptName)
    print("DESCRIPTION")
    print("\tDesigns n spaced seeds of length l. Each bit location would have n-m seeds")
    print("\trequiring a match. The set of seeds will be symetric to accommodate reverse")
    print("\tcomplements.")
    print("OPTIONS")
    print("\t-n|--seeds")
    print("\t\tnumber of spaced seeds to be designed (4)")
    print("\t-m|--misses")
    print("\t\tallowed number of missed hits (2). should be less than n")
    print("\t-l|--seedlen")
    print("\t\tseed length (100). should be more than n")
    print("\t-e|--minlen")
    print("\t\tminimum seed length (0) of at least one seed.")


def main():
    # get run time parameters
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hn:m:l:e:", ["help", "nseeds=", "misses=", "seedlen="])
    except getopt.GetoptError as err:
        # print help information and exit:
        print str("option not recognized")
        usage()
        sys.exit(2)
    nseeds = 4 # default number of seeds
    misses = 2 # default number of allowed number of misses
    seedlen = 100 # default seed length
    minlen = 0
    for opt, arg in opts:
        if opt in ("-n", "--nseeds"):
            nseeds = int(arg)
        elif opt in ("-m", "--missees"):
            misses = int(arg)
        elif opt in ("-l", "--seedlen"):
            seedlen = int(arg)
        elif opt in ("-e", "--minlen"):
            minlen = int(arg)
        elif opt in ("-h", "--help"):
            usage()
            sys.exit()
        else:
            assert False, "unhandled option"
    if (misses >= nseeds or seedlen <= nseeds):
        # print help information and exit:
        print str("### parameters violate assumptions!")
        usage()
        sys.exit(3)
    
    # initialize seeds
    seed=[[1]*((seedlen+1)//2)+[0]*(seedlen//2) for i in range(nseeds)]
    
    currMinLen = -1
    while minlen >= currMinLen:
        seed=[[1]*((seedlen+1)//2)+[0]*(seedlen//2) for i in range(nseeds)]
        seedList=range(0,nseeds)
        # scan through the seed length/2
        for i in range(0,((seedlen+1)//2)):
            # select seeds to introduce don't care positions
            dontcare=random.sample(seedList,nseeds-misses)
            # insert don't care positions
            for j in dontcare:
                seed[j][i]=0
        # make seed set symmetric
        for j in range(nseeds):
            for i in range(seedlen//2):
                seed[nseeds-j-1][seedlen-i-1]=seed[j][i]
                
        for j in range(nseeds):
            weight = 0
            for i in range(seedlen):
                if seed[j][i] != 0:
                    weight += 1
            if weight > currMinLen:
                currMinLen = weight
        
    # output seeds
    for j in range(nseeds):
        weight = 0
        for i in range(seedlen):
            sys.stdout.write(str(seed[j][i]))
            if seed[j][i] != 0:
                weight += 1
        sys.stdout.write("\t" + str(weight) + "\t" + str(seedlen-weight))
        if weight > currMinLen:
            currMinLen = weight
        print
        
    for j in range(nseeds):
        for i in range(seedlen):
            sys.stdout.write(str(seed[j][i]))
        sys.stdout.write(" ")
    print

if __name__ == "__main__":
    main()
