Weather System Detector for d4PDF
====

Objective detection of weather systems.
d4PDF (Global) version.

## Requirement
- python2.7
- f2py
- gfortran
- numpy
- matplotlib (for visualization)
## Setup

### Compile fortran code
```
python f2py.make.py detect_fsub.f90
```
### Edit Configuration file (config.py)
```
cfg['databaseDir'] = '/home/utsumi/mnt/lab_work/hk02/BAMS.EEE.2020/test.HPB_NAT.m100' # Input database directory
cfg['baseDir'] = '/home/utsumi/mnt/lab_tank/utsumi/bams2020'  # output directory
```

## Usage
### Edit main.py
```
prj     = "d4PDF"
model   = "__"
run     = "XX-HPB_NAT-100"   # {expr}-{scen}-{ens}
res     = "320x640"
noleap  = False
tstp_runmean = "6hr"
logDir = "/home/utsumi/log"

iYear, iMon = [2010,1] # Start of the detection period
eYear, eMon = [2010,1] # End of the detection period

iYear_data  = 2010 # First year of the available data
eYear_data  = 2010 # Last year of the available data
iMon_data   = 1 # First month of the available data
```
### Run
```
python main.py
```

### Draw the cyclone tracks
Edit draw.tclines.py and
```
python draw.tclines.py
```

## Author
(https://github.com/nbykutsumi)
