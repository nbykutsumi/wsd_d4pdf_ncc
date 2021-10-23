#! /usr/bin/python
#from JRA55 import Jra55
import JRA55
from JRA55 import anl_p125, anl_surf125, fcst_phy2m125, fcst_surf125_mon, LL125 
from datetime import datetime, timedelta

#class IO_Jra55(Jra55):
class IO_Jra55(anl_p125, anl_surf125, fcst_phy2m125, fcst_surf125_mon, LL125):
    def __init__(self, dbbaseDir):
        #Jra55.__init__(self, res)

        self.dvar = {
                "ta"   :"tmp"
               ,"ua"   :"ugrd"
               ,"va"   :"vgrd"
               ,"slp"  :"Mean sea level pressure"
               ,"spfh" :"spfh"
               ,"prcp" :"Mean total precipitation"
               ,"sst"  :"Brightness temperature"
               ,"topo" :"Surface height"
               ,"land" :"Land-sea mask"
               ,"pwat" :"PWAT"
               }

        self.anl_p125         = anl_p125(dbbaseDir)
        self.anl_surf125      = anl_surf125(dbbaseDir)
        self.fcst_phy2m125    = fcst_phy2m125(dbbaseDir)
        self.fcst_surf125_mon = fcst_surf125_mon(dbbaseDir)
        self.LL125            = LL125(dbbaseDir)

        self.Lat = JRA55.Lat125(crd='sa')
        self.Lon = JRA55.Lon125(crd='sa')
        self.ny  = len(self.Lat)
        self.nx  = len(self.Lon)
        self.miss= JRA55.miss()   # -9999.


    def Load_6hrPlev(self, var, DTime, plev):
        Var  = self.dvar[var]
        return self.anl_p125.load_6hr(Var, DTime, plev, crd='sa')

    def Load_6hrSfc(self, var, DTime):
        Var  = self.dvar[var]
        return self.anl_surf125.load_6hr(Var, DTime, crd='sa')

    def Load_monSfc(self, var, Year, Mon):
        Var  = self.dvar[var]
        return self.fcst_surf125_mon.load_mon(Var, Year, Mon, crd='sa')

    def Load_monPrcp_mmd(self, Year, Mon):
        return self.load_mon_prcp_mmd(Year, Mon, crd='sa')

    #def Load_day_spfh(self, DTime, plev, verbose=False):
    #    return self.time_ave(self.dvar["spfh"],DTime, DTime+timedelta(hours=23), timedelta(hours=6), lev=plev, verbose=False)


    def Load_const(self, var):
        Var  = self.dvar[var]
        return self.LL125.load_const(Var, crd='sa')
