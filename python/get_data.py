import sys,glob,os
import urllib.request
import tigress_models
import requests
import subprocess

def exists(path):
    r = requests.head(path)
    return r.status_code == requests.codes.ok

def web_download(file_list,dirname,overwrite=False):
    # check file
    data_url = 'http://tigress-web.princeton.edu/~changgoo/TIGRESS_example_data/rps-paper'
    for file_name in file_list:
        print("checking file: {}".format(dirname+file_name),end=' -- ')
        if os.path.isfile(dirname+file_name) and (not overwrite):
            print("found!")
        else:
            url = dirname.replace('../data',data_url)+file_name
            if exists(url):
                print("downloading",end='... ')
                urllib.request.urlretrieve(url, dirname+file_name)
                print("complete!")
            else:
                print("no such file")

def rsync_download():
    # check file
    data_url = 'changgoo@tigressdata.princeton.edu:/tigress/changgoo/public_html/TIGRESS_example_data/rps-paper/*'
    process = subprocess.run (['rsync','-Cav',data_url,'../data/'],universal_newlines=True)
    return process

# check directory
if not os.path.isdir('../data'): os.mkdir('../data')

model_list=tigress_models.get_model_list(return_type='list')
model_list=['RPS_8pc_ICM2_zmax6']

rsync=False
if len(sys.argv) > 1:
    rsync = eval(sys.argv[1])

if rsync:
    rsync_download()
else:
    kind='all'
    if len(sys.argv) > 2:
        kind = sys.argv[2]
    for pid in model_list:
        if not os.path.isdir('../data/'+pid): os.mkdir('../data/'+pid)
        if not os.path.isdir('../data/'+pid+'/hst'): os.mkdir('../data/'+pid+'/hst')
        if not os.path.isdir('../data/'+pid+'/zprof_merged'): os.mkdir('../data/'+pid+'/zprof_merged')
        if not os.path.isdir('../data/'+pid+'/vz_pdf'): os.mkdir('../data/'+pid+'/vz_pdf')
        file_list=['hst/{}.hst.p'.format(pid),
                   'hst/{}.sn.p'.format(pid),
                   'hst/{}.hst_zp.p'.format(pid),
                   '{}.par'.format(pid),]
        zp_list=['zprof_merged/{}.phase1.zprof.nc'.format(pid),
                 'zprof_merged/{}.phase2.zprof.nc'.format(pid),
                 'zprof_merged/{}.phase3.zprof.nc'.format(pid),
                 'zprof_merged/{}.phase4.zprof.nc'.format(pid),
                 'zprof_merged/{}.phase5.zprof.nc'.format(pid)]
        zpicm_list=['zprof_merged/{}.phase1-icm.zprof.nc'.format(pid),
                    'zprof_merged/{}.phase2-icm.zprof.nc'.format(pid),
                    'zprof_merged/{}.phase3-icm.zprof.nc'.format(pid),
                    'zprof_merged/{}.phase4-icm.zprof.nc'.format(pid),
                    'zprof_merged/{}.phase5-icm.zprof.nc'.format(pid)]


        if kind == 'all':
            file_list += zp_list + zpicm_list
        elif kind == 'zp':
            file_list = zp_list
        elif kind == 'zpicm':
            file_list = zpicm_list
        web_download(file_list,'../data/{}/'.format(pid),overwrite=True)

