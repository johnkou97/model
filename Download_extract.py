import requests
import sys
pyversion = sys.version_info.major
import urllib.request
import shutil
import tarfile as tar
import os
def download(url):
    filename = url.split('/')[-1]
    print('Downloading ' + url )
    if pyversion == 2:
        r = urllib2.urlopen(url).read()
        f = open(filename, 'w')   # write it to the right filename
        f.write(r)
        f.close()
    else:
        urllib.request.urlretrieve(url, filename)
    print("File download complete")
    return filename

#download BAM data
BAM=['002','003','004','009','010','022','035','036','046','048','053','057','058','059','061','065','070','080','089','090','091','092','093','098','107','121','122','123','124','126','128']
for k in BAM:
    url='http://141.35.26.215/data/Public_tar/BAM:0'+k+':R01.tar'
    request = requests.get(url)
    if request.status_code == 200:
        download(url)

#download THC data
THC=['01','02','03','04','05','06','07','10','11','12','13','14','15','16','17','18','19','20','21','29','31','32','36']
for k in THC:
    url='http://141.35.26.215/data/Public_tar/THC:00'+k+':R01.tar' 
    request = requests.get(url)
    if request.status_code == 200:
        download(url)

os.mkdir('data')
os.mkdir('metadata')
for m in range(0,2):
    for k in range(0,10):
        for j in range(0,10):
            name = 'BAM:0%s%s%s:R01.tar' %(m,k,j)
            try:
                file=tar.open(name)
               # print('bam')
                file.extractall()
               # print('yeah')
                shutil.copy('/home/jannis/Documents/model/Public/BAM:0{0}{1}{2}/R01/metadata.txt' .format(m,k,j),'/home/jannis/Documents/model/metadata/BAM:0{0}{1}{2}.txt' .format(m,k,j))
               # print('lol')
            except OSError:
                pass

for m in range(0,2):
    for k in range(0,10):
        for j in range(0,10):
            name = 'BAM:0%s%s%s:R01.tar' %(m,k,j)
            try:

                shutil.copy('/home/jannis/Documents/model/Public/BAM:0{0}{1}{2}/R01/data.h5' .format(m,k,j),'/home/jannis/Documents/model/data/BAM:0{0}{1}{2}.h5' .format(m,k,j))
               # print('lol')
            except OSError:
                pass


for m in range(0,4):
    for k in range(0,10):
            name = 'THC:00%s%s:R01.tar' %(m,k)
            try:
                file=tar.open(name)
               # print('bam')
                file.extractall()
               # print('yeah')
                shutil.copy('/home/jannis/Documents/model/Public/THC:00{0}{1}/R01/metadata.txt' .format(m,k),'/home/jannis/Documents/model/metadata/THC:00{0}{1}.txt' .format(m,k))
               # print('lol')
            except OSError:
                pass


for m in range(0,4):
    for k in range(0,10):
            name = 'THC:00%s%s:R01.tar' %(m,k)
            try:
                shutil.copy('/home/jannis/Documents/model/Public/THC:00{0}{1}/R01/data.h5' .format(m,k),'/home/jannis/Documents/model/data/THC:00{0}{1}.h5' .format(m,k))
               # print('lol')
            except OSError:
                pass
