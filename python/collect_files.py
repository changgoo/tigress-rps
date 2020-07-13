import sys,os,glob
import shutil
from distutils.file_util import copy_file
import tigress_models

model_list=tigress_models.get_model_list(return_type='list')

base='/tigress/changgoo/'
tobase='/tigress/changgoo/public_html/TIGRESS_example_data/wind-paper/'
for pid in model_list:
    hstfname='{}{}/hst/{}.hst.p'.format(base,pid,pid)
    snfname='{}{}/hst/{}.sn.p'.format(base,pid,pid)
    hstzpfname='{}{}/hst/{}.hst_zp.p'.format(base,pid,pid)
    parfname='{}{}/{}.par'.format(base,pid,pid)

    todir='{}{}'.format(tobase,pid)
    if not os.path.isdir(todir): os.mkdir(todir)
    if not os.path.isdir(todir+'/hst'): os.mkdir(todir+'/hst')
    for f in [hstfname,snfname,hstzpfname,parfname]:
        #shutil.copy(f,f.replace(base,tobase))
        dstname,copied=copy_file(f,f.replace(base,tobase),update=1,verbose=1)
        if copied: print(dstname)

    if not os.path.isdir(todir+'/zprof_merged'): os.mkdir(todir+'/zprof_merged')
    for ph in ['phase1','phase2','phase3','phase4','phase5']:
        zpfile='{}{}/zprof_merged/{}.{}.nc'.format(base,pid,pid,ph)
        dstname,copied=copy_file(fluxfile,fluxfile.replace(base,tobase),update=1,verbose=1)
        if copied: print(dstname)
        zpfile='{}{}/zprof_merged/{}.{}-icm.nc'.format(base,pid,pid,ph)
        dstname,copied=copy_file(fluxfile,fluxfile.replace(base,tobase),update=1,verbose=1)
        if copied: print(dstname)
