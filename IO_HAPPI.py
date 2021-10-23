#! /usr/bin/python
from HAPPI import Happi
from datetime import datetime, timedelta

class IO_Happi(Happi):
    def __init__(self, model, run, res):
        expr, scen, ens = run.split("-")[:3]
        ens  = int(ens)
        Happi.__init__(self)
        Happi.__call__(self, model, expr, scen, ens)

        self.dvar = {
                "ta"   :"T"
               ,"ua"   :"u"
               ,"va"   :"v"
               ,"slp"  :"slp"
               ,"spfh" :"q"
               ,"prcp" :"prcp"
               ,"sst"  :"Ts"
               ,"topo" :"topo"
               ,"land" :"lndfrc"
               }

    def Load_6hrPlev(self, var, DTime, plev):
        Var  = self.dvar[var] + "%03d"%(plev)
        return self.load_6hr(Var, DTime)

    def Load_6hrSfc(self, var, DTime):
        Var  = self.dvar[var]
        return self.load_6hr(Var, DTime)

    def Load_dayPlev(self, var, DTime, plev):
        Var  = self.dvar[var] + "%03d"%(plev)
        return self.load_day(Var, DTime)

    def Load_daySfc(self, var, DTime):
        Var  = self.dvar[var]
        return self.load_day(Var, DTime)

    def Load_monSfc(self, var, Year, Mon):
        Var  = self.dvar[var]
        return self.load_mon(Var, Year, Mon)

    def Load_monPrcp_mms(self, Year, Mon):
        return self.load_mon_prcp_mms(Year, Mon)

    def Load_monPrcp_mmh(self, Year, Mon):
        return self.load_mon_prcp_mmh(Year, Mon)

    def Load_monPrcp_mmd(self, Year, Mon):
        return self.load_mon_prcp_mmd(Year, Mon)

    def Load_day_spfh(self, DTime, plev):
        return self.Load_dayPlev("spfh", DTime, plev)

    def Load_const(self, var):
        Var  = self.dvar[var]
        return self.load_const(Var)


#    def load_monSfc(self, var, Year, Mon):
#        #return self.load_mon(var, Year, Mon) 
