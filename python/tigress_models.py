from data_container import data_container
import numpy as np
import sys,glob,os
import xarray as xr

import pyathena.preprocessing as pre

def get_model_list(return_type='dict'):

    ism_models=dict()
    ism_models['R8']=dict(pid='R8_4pc_newacc')
    ism_models['R8-8pc']=dict(pid='R8_8pc_newacc')
    ism_models['ISM']=dict(pid='R8_8pc_metal')
    ism_models['ISM-hydro']=dict(pid='R8_8pc_metal_hydro')

    icm_models=dict()
    icm_models['ICM-P1']=dict(pid='RPS_8pc_ICM0_newacc')
    icm_models['ICM-P3']=dict(pid='RPS_8pc_ICM1_newacc')
    icm_models['ICM-P3h']=dict(pid='RPS_4pc_ICM1_newacc')
    icm_models['ICM-P7']=dict(pid='RPS_8pc_ICM2_newacc')
    icm_models['ICM-P7h']=dict(pid='RPS_4pc_ICM2_newacc')
    icm_models['ICM-P14']=dict(pid='RPS_8pc_ICM3_newacc')

    all_models={**ism_models,**icm_models}

    Models = all_models

    model_list=[]
    for key in Models.keys():
        model_list.append(Models[key]['pid'])
    if return_type == 'list':
        return model_list
    else:
        return Models

def recal(base='/tigress/changgoo/',flux_recal=True):
    Models=init(read=False)
    for key in Models.keys():
        pid=Models[key]['pid']
        dc=data_container(pid,create=flux_recal)

if __name__ == '__main__':
    flux=False
    if narg>1: flux=eval(sys.argv[1])
    recal(flux_recal=flux)
