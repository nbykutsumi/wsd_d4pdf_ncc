import glob
import subprocess
#scmd = 'ls -l /home/utsumi/mnt/lab_work_hk03/d4PDF_GCM/HPB_NAT/m0*/201012/*'
#scmd = 'ls -l /work/hk03/d4PDF_GCM/HPB_NAT/m0*/201012/*'
scmd = 'ls -l /home/utsumi/mnt/lab_work_hk03/d4PDF_GCM/HPB_NAT/m0*/201012/*'
subprocess.call(scmd, shell=True)

#print '\n\n*******************'
#scmd = 'ls /home/utsumi/mnt/lab_work_hk03/d4PDF_GCM/HPB/m0*/201012/*'
#subprocess.call(scmd, shell=True)
