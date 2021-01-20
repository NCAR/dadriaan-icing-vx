import datetime,os
st = datetime.datetime(2019,2,1,0,0,0)
et = datetime.datetime(2019,2,28,0,0,0)
incr = 86400

while st <= et:

  cmd = 'python AB_diffmap.py %s %s' % (datetime.datetime.strftime(st,'%Y%m%d'),'*.nc')
  print(cmd)
  os.system(cmd)
  
  st = st + datetime.timedelta(seconds=incr)
