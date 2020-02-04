#!/usr/local/python/bin/python

import os

pirepFile = "/home/dadriaan/projects/sae2019/data/match/hrrr/PIREPShrrrCIPPOSNEG.out"
f = open(pirepFile,"r")

data = f.readlines()

for l in data:
  fout = open("/home/dadriaan/projects/sae2019/icing-vx/tmp.in","w")
  fout.write('%s' % (l))
  fout.close()

  cmd = "python /home/dadriaan/projects/sae2019/icing-vx/extract_grid_data_CIP.py %s" % ("/home/dadriaan/projects/sae2019/icing-vx/tmp.in")
  print(cmd)
  os.system(cmd)

  if os.path.exists("/home/dadriaan/projects/sae2019/icing-vx/all.csv"):
    cmd = "tail -1 /home/dadriaan/projects/sae2019/icing-vx/tmp.csv >> /home/dadriaan/projects/sae2019/icing-vx/all.csv"
  else:
    cmd = "cat /home/dadriaan/projects/sae2019/icing-vx/tmp.csv >> /home/dadriaan/projects/sae2019/icing-vx/all.csv"
  print(cmd)
  os.system(cmd)
  
  cmd = "rm -rfv /home/dadriaan/projects/sae2019/icing-vx/tmp.csv"
  print(cmd)
  os.system(cmd)

cmd = "rm -rfv /home/dadriaan/projects/sae2019/icing-vx/tmp.in"
print(cmd)
os.system(cmd)
