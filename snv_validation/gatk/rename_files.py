import os
import re

tmp = os.listdir('.')

files = []
for f in range(len(tmp)):
    if not tmp[f].startswith('U'):
        continue
    files.append(tmp[f])

old_new = []
tmp = []
for f in range(len(files)):
    fil = files[f].split('.')
    if len(fil) == 3:
        fil[0] = fil[0] + '_000'
        tmp = fil
    elif len(fil) == 4:
        fil[0] = fil[0] + '_00' + fil[3]
        tmp = fil[:-1]
    else:
        print('Something is up')
        break

    old_new.append( [files[f], '.'.join(tmp)] )

for o in range(len(old_new)):
    os.rename(old_new[o][0], old_new[o][1])
