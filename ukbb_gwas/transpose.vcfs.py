import sys
import re

fid = re.compile(r'^(\d*)_\1$')
a00 = re.compile('0/0')
a10 = re.compile('1/0')
a01 = re.compile('0/1')
a11 = re.compile('1/1')
an0 = re.compile('./0')
an1 = re.compile('./1')
a0n = re.compile('0/.')
a1n = re.compile('1/.')
ann = re.compile('./.')

infile = sys.argv[1]
outfile = sys.argv[2]

with open(infile, 'r') as f:
for i in range(6):
  next(f)
lis = [x.split()[9:] for x in f]

with open(outfile, 'w') as f:
for x in zip(*lis):
  str=''
  for y in x:
    y = fid.sub(r'\1\t\1', y)
    y = a00.sub('0', y)
    y = a10.sub('1', y)
    y = a01.sub('1', y)
    y = a11.sub('2', y)
    y = an0.sub('NA', y)
    y = an1.sub('NA', y)
    y = a0n.sub('NA', y)
    y = a1n.sub('NA', y)
    y = ann.sub('NA', y)
    str = str+y+'\t'
  str = str+'\n'
  f.write(str)
