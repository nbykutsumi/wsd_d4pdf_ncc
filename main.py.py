import subprocess
import os, sys

prj    = "HAPPI"
model  = "MIROC5"
expr   = "C20"
#lscen  = ["ALL","P15","P20"]
lscen  = ["P15"]

#lens   = [2]
#lens   = [3,4,5,6]
#lens   = [4,5,6]
#lens   = [7,8,9,10]
#lens   = range(12,15+1)
#lens   = range(16,20+1)
#lens   = range(22,30+1)
#lens   = range(32,40+1)
#lens   = range(1,15+1)
#lens   = range(16,40+1)
lens   = range(41,50+1)

dYM  = {
           "ALL":[[2006,1],[2015,12]]
          ,"P15":[[2106,1],[2115,12]]
          ,"P20":[[2106,1],[2115,12]]
         }

noleap       = True
tstp_runmean = "day"


lkeys = [[scen,ens] for scen in lscen for ens in lens]
for scen, ens in lkeys:
    run    = "%s-%s-%03d"%(expr,scen,ens)
    res    = "128x256"

    iYear,iMon = dYM[scen][0]
    eYear,eMon = dYM[scen][1]
    iYear_data = iYear
    eYear_data = eYear
    iMon_data  = 1
    iYearMinMax= iYear 
    eYearMinMax= eYear 

    cmd  = ["python","./main.py", prj, model, run, res, noleap, tstp_runmean
            ,iYear, iMon, eYear, eMon
            ,iYear_data, eYear_data
            ,iMon_data
            ,iYearMinMax, eYearMinMax
            ]

    cmd  = " ".join(map(str,cmd))
    p= subprocess.call(cmd, shell=True)

