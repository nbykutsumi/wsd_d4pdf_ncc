#! /usr/bin/python
#from d4PDF import H
from d4PDF import snp_6hr_2byte, avr_mon_320x640
from datetime import datetime, timedelta
import sys

class IO_d4pdf(snp_6hr_2byte):
    def __init__(self, run, dbbaseDir):
        expr, scen, ens = run.split("-")[:3]
        self.expr = expr
        self.scen = scen
        self.ens  = int(ens)
        self.dbatm   = snp_6hr_2byte(vtype='atm', dbbaseDir=dbbaseDir)
        self.dbsfc   = snp_6hr_2byte(vtype='sfc', dbbaseDir=dbbaseDir)
        self.dbsfcmon= avr_mon_320x640(vtype='sfc', dbbaseDir=dbbaseDir)

        self.Lat   = self.dbatm.Lat
        self.Lon   = self.dbatm.Lon
        self.ny    = self.dbatm.ny
        self.nx    = self.dbatm.nx
        self.miss  = self.dbatm.miss_out

    #    Happi.__init__(self)
    #    Happi.__call__(self, model, expr, scen, ens)

        self.dvar = {
                "ta"   :"T"
               ,"ua"   :"U"
               ,"va"   :"V"
               ,"slp"  :"SLP"
               ,"spfh" :""
               ,"prcp" :"PRECIPI"
               ,"sst"  :"TGEF"
               ,"topo" :"height"
               ,"land" :"ratiol"
               }

    def Load_6hrPlev(self, var, DTime, plev):
        Var  = self.dvar[var] + "%03d"%(plev)
        return self.dbatm.load_6hr(Var, self.scen, self.ens, DTime, miss_fill=self.miss)

    def Load_6hrSfc(self, var, DTime):
        Var  = self.dvar[var]
        return self.dbsfc.load_6hr(Var, self.scen, self.ens, DTime, miss_fill=self.miss)


    #def Load_dayPlev(self, var, DTime, plev):
    #    Var  = self.dvar[var] + "%03d"%(plev)
    #    return self.load_day(Var, DTime)

    #def Load_daySfc(self, var, DTime):
    #    Var  = self.dvar[var]
    #    return self.load_day(Var, DTime)

    def Load_monSfc(self, var, Year, Mon):
        Var  = self.dvar[var]
        return self.dbsfcmon.load_ave_mon(Var, self.scen, self.ens, Year, Mon, miss_fill=self.miss)

    #def Load_monPrcp_mms(self, Year, Mon):
    #    return self.load_mon_prcp_mms(Year, Mon)

    #def Load_monPrcp_mmh(self, Year, Mon):
    #    return self.load_mon_prcp_mmh(Year, Mon)

    #def Load_monPrcp_mmd(self, Year, Mon):
    #    return self.load_mon_prcp_mmd(Year, Mon)

    #def Load_day_spfh(self, DTime, plev):
    #    return self.Load_dayPlev("spfh", DTime, plev)

    def Load_const(self, var):
        if var=='topo':
            data = self.dbatm.load_topo(vname='height', miss_fill=self.miss)
        elif var=='land':
            data = self.dbatm.load_topo(vname='ratiol', miss_fill=self.miss)
        else:
            print 'check var',var
            sys.exit()
        return data

