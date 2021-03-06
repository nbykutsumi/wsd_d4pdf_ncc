from numpy import *
import os, sys

def IO_Master(prj, model, run, res, dbbaseDir):
    if prj=="JRA55":
        import IO_JRA55
        #iom  = IO_JRA55.IO_Jra55(model, run, res)
        iom  = IO_JRA55.IO_Jra55(dbbaseDir)
   
    elif prj=="HAPPI":
        import IO_HAPPI
        iom  = IO_HAPPI.IO_Happi(model, run, res)

    elif prj=="d4PDF":
        import IO_d4PDF
        iom  = IO_d4PDF.IO_d4pdf(run, dbbaseDir) 

    else:
        print("check prj, model, run, res=")
        print((prj, model, run, res))
        sys.exit()
    return iom


