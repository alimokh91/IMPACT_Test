import subprocess
import re

file = './MRT_Daten_Bern/127/10005583/10005584/10009705/1000A2D1'

venc_stringXandY='sWiPMemBlock\.adFree\[8\]                   = [0-9]*\.?[0-9]*'
venc_stringZ='sWiPMemBlock\.adFree\[9\]                   = [0-9]*\.?[0-9]*'
reg_ex_number='[0-9]*\.?[0-9]*'
output = subprocess.run(['gdcmdump','--csa',file],check=True, stdout=subprocess.PIPE, universal_newlines=True)

resXandY = re.findall(venc_stringXandY, output.stdout)
resZ = re.findall(venc_stringZ, output.stdout)

vencXandY = float(re.findall(reg_ex_number,resXandY[0])[-2])
vencZ = float(re.findall(reg_ex_number,resZ[0])[-2])

print('venc along X and Y = ', vencXandY)
print('venc along Z = ', vencZ)
