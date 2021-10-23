import glob
import subprocess
#scmd = 'ls -l /home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM/XX-HPB-*/6hr/clist/2010/12'
##scmd = 'ls -l /home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM/XX-HPB-*/6hr/clist/1980/12'
#subprocess.call(scmd, shell=True)


#print '*******************'
#scmd = 'ls -l /home/utsumi/mnt/lab_tank/utsumi/WS/d4PDF_GCM/XX-HPB_NAT-*/6hr/clist/2010/12'
scmd = 'ls -l /media/disk2/out/WS/d4PDF_GCM/XX-HPB_NAT-*/6hr/clist/2010/12'
subprocess.call(scmd, shell=True)
