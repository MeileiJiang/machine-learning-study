import sys
from os import system, makedirs
from subprocess import call


f = open(sys.argv[1], "r")
condor = "".join(f.readlines())
f.close()


for n in range(100, 501, 50):
    open("submit_condor", "w").write("%s" % condor % {'n' : n})
    retcode = call(["condor_submit", "submit_condor"]) 
    

