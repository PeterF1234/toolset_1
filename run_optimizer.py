#!/usr/bin/env python3

import os
import bench_opt as BO
import add_labels_optimized_defonly as labels
# import matplotlib.pyplot as plt
import numpy as np
# import sys
import concurrent.futures
import time

def opt_spect(error_function):
    lambdarange = np.arange(200,700)
    SIGMA = 0.2
    lambda_shift = 1
    ref_pad = 50
    files = os.listdir()
    ref_csv = "PTZ_1.csv"
    header = "reference,molecule,functional,basis,error_function,error_value,bandwidth,shift_factor,Wavelength [nm],Intensity [a.u.]\n"

    
    errfn_name = error_function.__name__
    os.mkdir(errfn_name)
    dumpfile_name = 'optimized_params_'+errfn_name+'.csv'
    # optimize the spectra
    print(f'process for {errfn_name} has started')
    for logfile in files:
        if logfile[-4:] == '.log':
            print(f'Working on file: {logfile}')
            BO.print_result(ref_csv,logfile,lambdarange,SIGMA,lambda_shift,ref_pad,dumpfile=dumpfile_name,writeopt=True,errorfunc=error_function,workdir=errfn_name,logscale_hmap=True)
            # sys.exit(1)
    # build the database
    opt_csv_files = []
    outfile = ref_csv[:-4] + "_database_" + errfn_name + ".csv"
    for i in os.listdir(errfn_name):
        errorfile = 'optimized_params_' + errfn_name + '.csv'
        if i[-4:] == ".csv" and i != errorfile:
            opt_csv_files.append(i)
        else:
            print(f'Skipping file {i} ...')
    labels.build_filedata(opt_csv_files,outfile,dumpfile_name,header,workdir=errfn_name)
    # databases.append(outfile)
    return f'Done with {errfn_name}!'

if __name__ ==  '__main__':
    start = time.perf_counter()

    # error_function = BO.r_square # r_square is a function defined in bench_opt.py
    # error_functions = [BO.rmsle, BO.mae, BO.mse, BO.r_square, BO.lg_mse, BO.rmsd]
    error_functions = [BO.rmsle, BO.mae, BO.mse, BO.r_square]
    # error_functions = [BO.mse, BO.r_square]
    ref_csv = "PTZ_1.csv"
    
    BO.save_expt(ref_csv,ref_pad=50) # reformat the reference spectrum for the final database
    #sys.exit(1)
    
    databases = []
    
    with concurrent.futures.ProcessPoolExecutor() as executor:
        e = executor.map(opt_spect, error_functions)
        for i in e:
            print(i)

    for i in error_functions:
        errfn_name = i.__name__
        dbname = ref_csv[:-4] + "_database_" + errfn_name + ".csv"
        databases.append(dbname)
    # merge the databases and add the reference spectrum
    finaldb_name = "database_merged.csv"
    labels.merge_db(databases,merged_name=finaldb_name)
    ref_filename = ref_csv[:-4] + "_integerx.csv" # BO.save_expt creates a file with this name
    labels.add_reference(ref_filename, finaldb_name)
    
    finish = time.perf_counter()
    print(f'Finished in {round(finish-start, 2)} second(s)')

