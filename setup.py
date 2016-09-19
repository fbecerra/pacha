import sys, os
from subprocess import call

try:
  option = str(sys.argv[1])
except:
  print 'Usage:'
  print '1. Installation: python setup.py install'
  print '2. Clean: python setup.py clean'
  sys.exit()

if option == 'install':
  os.system('OPT="-O3 -ffast-math" python compile.py build_ext -i')
elif option == 'clean':
  os.system('rm -rf src/*.c src/*.so build/*')
else:
  print 'Invalid option!'
  sys.exit()
