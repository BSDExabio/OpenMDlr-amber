#!/usr/bin/env python

import sys
from subprocess import Popen, PIPE

seqlen = [129, 76, 132, 20, 40]
names = ["1aki", "1ubq", "1ifc", "1l2y", "2erl"] #1gci
angRST = [15.0] #15.0, 30.0] #20.0, 30.0] #, 30.0, 60.0]
distRST = [1.0, 2.0, 5.0]
angforceRST = [40.0, 70.0] #0.0, 20.0, 45.0, 70.0]
distforceRST = [40.0, 70.0] #, 20.0, 45.0, 70.0]

def tup2str(angtup, disttup):
	angfst = str(angtup[0])
	angsnd = str(angtup[1])
	distfst = str(disttup[0])
	distsnd = str(disttup[1])
	strform = "f"+angfst+"_"+distfst+"_s"+angsnd+"_"+distsnd
	return(strform)


for i, n in enumerate(names):
    for a in angRST:
        for d in distRST:
            for af in angforceRST:
                for df in distforceRST:
                    run_name= n+'_'+str(a)+'_'+str(d)+'_'+str(af)+'_'+str(df)+'_tmr'
                    an_name = 'annealing_'+run_name

                    h = open(n+"/helperscript", "w+")
                    h.write('#!/bin/bash\n\nsed "s/NAME/'+n+'/g" < annealing_gen.pbs > '+ an_name+'\nsed -i "s/ANGRST/'+str(a)+'/g" '+an_name+'\nsed -i "s/DISTRST/'+str(d)+'/g" '+an_name+'\nsed -i "s/SUB/'+run_name+'/g" '+an_name+'\nsed -i "s/FINAL/'+run_name+'/g" '+ an_name+'\nsed -i "s/ANGFORCERST1/'+str(af)+'/g" '+ an_name+'\nsed -i "s/DISTFORCERST1/'+str(df)+'/g" '+an_name+'\nqsub '+ an_name)
                    h.close()


                    process = Popen('chmod u+x '+n+'/helperscript', shell=True, stdout=PIPE, stderr=PIPE)
                    stdout, stderr = process.communicate()

                    process = Popen('source '+n+'/helperscript', shell=True, stdout=PIPE, stderr=PIPE)
                    stdout, stderr = process.communicate()

