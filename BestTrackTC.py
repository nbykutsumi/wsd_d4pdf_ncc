import sys
import IBTrACS

def BestTrack(prjName="IBTrACS"):
  if prjName == "IBTrACS":
    self = IBTrACS.IBTrACS()
  return self

def BestTrack_2D(prjName, Year, a1lon, a1lat, miss):
  if prjName == "IBTrACS":
    #self = IBTrACS.IBTrACS_2D(Year,a1lon, a1lat, miss, ver="v03r06")
    self = IBTrACS.IBTrACS_2D(Year,a1lon, a1lat, miss, ver="v03r10")
  return self


