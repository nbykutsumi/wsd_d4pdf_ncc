import subprocess, sys

#-- read input file name -----
fortname = sys.argv[1]
if fortname[-4:] == ".f90":
  binname  =  fortname[:-4]
else:
  print "check file name", fortname
#-- compile ------------------
cmd = "f2py --f90exec=/usr/bin/gfortran -c -m %s %s"%(binname, fortname)
#cmd = "f2py --f90exec=/opt/intel/composerxe/bin/ifort -c -m %s %s"%(binname, fortname)
subprocess.call(cmd, shell=True)





