#!/usr/bin/env python3

#import numpy as np
import pandas
import matplotlib.pyplot as plt
import seaborn as sns
import time


def do_wide_errors(errors):
    # compute the occcurrences of a each functional having an optimized bandwidth that is considered wide
    wide_errors = errors[errors['bandwidth'] > 0.35]
    wide_errors_functionals = wide_errors.groupby("functional", as_index=False).size()
    wide_errors_functionals = wide_errors_functionals.rename(columns={"size": "counts"})
    # wide_errors_functionals.plot(kind='bar', ylabel='Counts')
    sns.barplot(x="counts", y="functional", data=wide_errors_functionals, color="orangered").set(ylabel=None)
    sns.despine()
    plt.tight_layout()
    plt.savefig('wide_functionals.png', dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.savefig('wide_functionals.svg', bbox_inches='tight', pad_inches=0)
    plt.close()

def do_wide_errors_molecule(errors):
    # compute the occcurrences of a each molecule having an optimized bandwidth that is considered wide
    wide_errors = errors[errors['bandwidth'] > 0.35]
    wide_errors_molecules = wide_errors.groupby("molecule", as_index=False).size()
    wide_errors_molecules = wide_errors_molecules.rename(columns={"size": "counts"})
    wide_errors_molecules_withzero = pandas.merge(errors,wide_errors_molecules,how='left',on='molecule').fillna(0)
    wide_errors_molecules_withzero = wide_errors_molecules_withzero.astype({'counts': 'int64'})
    wide_errors_to_plot = wide_errors_molecules_withzero.groupby("molecule", as_index=False).median()
    wide_errors_to_plot["maxcounts"]=40 # 40 = len(errors.functional.unique())*len(errors.error_function.unique())
    # plot
    e_min = 5 # 0-4 = within 10% of the total 40
    colors = ['royalblue' if (x < e_min) else "orangered" for x in wide_errors_to_plot["counts"]]
    # e_min = wide_errors_molecules["counts"].min()
    # colors = ['royalblue' if (x < e_min*1.2) else "orangered" for x in wide_errors_molecules["counts"]]
    sns.set_context("paper", font_scale=1.0, rc={"lines.linewidth": 1.75,"axes.linewidth": 1.75,'xtick.major.width': 1.6,'ytick.major.width': 1.6})
    sns.barplot(x="maxcounts", y="molecule", data=wide_errors_to_plot, color='lightgrey', zorder=0).set(ylabel=None)
    sns.barplot(x="counts", y="molecule", data=wide_errors_to_plot, palette=colors, zorder=1).set(ylabel=None)
    sns.despine()
    plt.tight_layout()
    plt.savefig('wide_molecules.png', dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.savefig('wide_functionals.svg', bbox_inches='tight', pad_inches=0)
    plt.close()

def do_optimal_fitparams(errors):
    # compute the occcurrences of a each functional having an optimized bandwidth and shift that are considered optimal
    optimals = errors[errors['bandwidth'] < 0.30]
    optimals = optimals[optimals['shift_factor'] < 1.2]
    optimals = optimals[optimals['shift_factor'] > 0.8]
    optimal_functionals = optimals.groupby("functional", as_index=False).size()
    optimal_functionals = optimal_functionals.rename(columns={"size": "counts"})
    # plot
    sns.barplot(x="counts", y="functional", data=optimal_functionals, color="seagreen").set(ylabel=None)
    sns.despine()
    plt.tight_layout()
    plt.savefig('optimal_functionals.png', dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.savefig('optimal_functionals.svg', bbox_inches='tight', pad_inches=0)
    plt.close()

def do_optimal_fitparams_molecules(errors):
    # compute the occcurrences of a each molecule having an optimized bandwidth and shift that are considered optimal
    optimals = errors[errors['bandwidth'] < 0.30]
    optimals = optimals[optimals['shift_factor'] < 1.2]
    optimals = optimals[optimals['shift_factor'] > 0.8]
    optimal_molecules = optimals.groupby("molecule", as_index=False).size()
    optimal_molecules = optimal_molecules.rename(columns={"size": "counts"})
    optimal_molecules_withzero = pandas.merge(errors,optimal_molecules,how='left',on='molecule').fillna(0)
    optimal_molecules_withzero = optimal_molecules_withzero.astype({'counts': 'int64'})
    optimals_to_plot = optimal_molecules_withzero.groupby("molecule", as_index=False).median()
    optimals_to_plot["maxcounts"]=40 # 40 = len(errors.functional.unique())*len(errors.error_function.unique())
    # plot
    e_max = optimals_to_plot["counts"].max()
    colors = ['royalblue' if (x > e_max*0.95) else "green" for x in optimals_to_plot["counts"]]
    sns.set_context("paper", font_scale=1.0, rc={"lines.linewidth": 1.75,"axes.linewidth": 1.75,'xtick.major.width': 1.6,'ytick.major.width': 1.6})
    # sns.barplot(x="counts", y="molecule", data=all_molecules, color='grey').set(ylabel=None)
    sns.barplot(x="maxcounts", y="molecule", data=optimals_to_plot, color='lightgrey', zorder=0).set(ylabel=None)
    sns.barplot(x="counts", y="molecule", data=optimals_to_plot, palette=colors, zorder=1).set(ylabel=None)
    sns.despine()
    plt.tight_layout()
    plt.savefig('optimal_molecules.png', dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.savefig('optimal_functionals.svg', bbox_inches='tight', pad_inches=0)
    plt.close()

def within_10_percent(errors):
    optimals = errors[errors['shift_factor'] < 1.1]
    optimals = optimals[optimals['shift_factor'] > 0.9]
    optimal_functionals = optimals.groupby("functional", as_index=False).size()
    optimal_functionals = optimal_functionals.rename(columns={"size": "counts"})
    # plot
    sns.barplot(x="counts", y="functional", data=optimal_functionals, color="royalblue").set(ylabel=None)
    sns.despine()
    plt.tight_layout()
    plt.savefig('optimal_functionals_within_10_percent.png', dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.savefig('optimal_functionals_within_10_percent.svg', bbox_inches='tight', pad_inches=0)
    plt.close()
    
def do_spectrum_error_function(expt, pred, molecule, functional):
    pred['error_function'].replace({'r_square': '$r^2$'},inplace=True)
    molecule_expt = expt[expt['molecule'] == molecule]
    molecule_pred = pred[pred['molecule'] == molecule]
    molecule_pred_functional = molecule_pred[molecule_pred['functional'] == functional]
    sns.lineplot(x="Wavelength [nm]", y="Intensity [a.u.]", data=molecule_expt, color="black", linestyle='dashed')
    sns.lineplot(x="Wavelength [nm]", y="Intensity [a.u.]", data=molecule_pred_functional, hue="error_function").legend(frameon=False,fontsize='12',ncol=2)
    sns.despine()
    plt.tight_layout()
    plt.xlim(300,600)
    # plt.ylim(0,0.1)
    if functional == '$\\omega$-B97XD':
        plt.savefig(f'{molecule}_error_functions_fit_wB97XD.png', dpi=300, bbox_inches='tight', pad_inches=0)
        print(f'{molecule}_error_functions_fit_wB97XD.png' + ' has been created!')
    else:
        plt.savefig(f'{molecule}_error_functions_fit_{functional}.png', dpi=300, bbox_inches='tight', pad_inches=0)
        print(f'{molecule}_error_functions_fit_{functional}.png' + ' has been created!')
    # plt.savefig(f'{molecule}_error_functions_fit_{molecule}.svg', bbox_inches='tight', pad_inches=0)
    plt.close()

def do_all_spectrum_error_function(expt, pred):
    molecules = list(expt.molecule.unique())
    functionals = list(pred.functional.unique())
    for i in molecules:
        for j in functionals:
            do_spectrum_error_function(expt, pred, i, j)
    
def do_avg_errors(errors, error_function):
    sub_errors = errors[errors['error_function'] == error_function]
    error_means = sub_errors.groupby("functional", as_index=False).mean()
    error_means = error_means.rename(columns={"error_value": f"mean {error_function} error"})
    # plot
    sns.barplot(x=f"mean {error_function} error", y="functional", data=error_means, color="royalblue").set(ylabel=None)
    sns.despine()
    plt.tight_layout()
    plt.xlim(error_means[f"mean {error_function} error"].min()*0.8,error_means[f"mean {error_function} error"].max()*1.05)
    plt.savefig(f'error_means_{error_function}.png', dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.savefig(f'error_means_{error_function}.svg', bbox_inches='tight', pad_inches=0)
    plt.close()

def do_avg_errors_diffcolor(errors, error_function):
    sub_errors = errors[errors['error_function'] == error_function]
    error_means = sub_errors.groupby("functional", as_index=False).mean()
    error_means = error_means.rename(columns={"error_value": f"mean {error_function} error"})
    # plot
    e_min = error_means[f"mean {error_function} error"].min()
    e_max = error_means[f"mean {error_function} error"].max()
    if error_function == "r_square" or error_function == '$r^2$':
        colors = ['green' if (x > e_max*0.95) else "royalblue" for x in error_means[f"mean {error_function} error"]]
    else:
        colors = ['green' if (x < e_min*1.15) else "royalblue" for x in error_means[f"mean {error_function} error"]]
    sns.barplot(x=f"mean {error_function} error", y="functional", data=error_means, palette=colors).set(ylabel=None)
    sns.despine()
    plt.tight_layout()
    plt.xlim(e_min*0.8,e_max*1.05)
    plt.ticklabel_format(style='scientific',scilimits=(-2,3),axis='x')
    plt.savefig(f'error_means_{error_function}.png', dpi=300, bbox_inches='tight', pad_inches=0)
    # plt.savefig(f'error_means_{error_function}.svg', bbox_inches='tight', pad_inches=0)
    plt.close()

def do_avg_errors_all(errors):
    errors['error_function'].replace({'r_square': '$r^2$'},inplace=True)
    error_functions = list(errors.error_function.unique())
    for i in error_functions:
        do_avg_errors(errors, i)

def do_avg_errors_all_diffcolor(errors):
    errors['error_function'].replace({'r_square': '$r^2$'},inplace=True)
    error_functions = list(errors.error_function.unique())
    for i in error_functions:
        do_avg_errors_diffcolor(errors, i)
    
if __name__ ==  '__main__':
    start = time.perf_counter()
    
    plt.style.use("seaborn-paper")
    sns.set_context("paper", font_scale=1.75, rc={"lines.linewidth": 1.75,"axes.linewidth": 1.75,'xtick.major.width': 1.6,'ytick.major.width': 1.6})
    
    fulldb = pandas.read_csv("OPC_database.csv", header=0)
    molecule_list = list(fulldb.molecule.unique())
    # default header looks like this: "reference,molecule,functional,basis,error_function,error_value,bandwidth,shift_factor,Wavelength [nm],Intensity [a.u.]\n"
    #fulldb=fulldb.astype({'Wavelength [nm]': 'float64'})
    #fulldb=fulldb.astype({'Intensity [a.u.]': 'float64'})
    #to debug: fulldb.dtypes
    expt = fulldb[fulldb['reference'] == 'yes']
    expt = expt.drop(columns=["reference","functional","basis","error_function","error_value","bandwidth","shift_factor"])
    
    pred = fulldb[fulldb['reference'] == 'no']
    pred = pred[pred["error_function"].isin(["mae","mse","r_square","rmsle"])]
    errors_list = list(pred.error_function.unique())
    errors = pred[pred['Wavelength [nm]'] == 333] # this is the database without the spectral data
    errors = errors.drop(columns=["reference","Wavelength [nm]","Intensity [a.u.]"])
    errors = errors.astype({'error_value': 'float64','bandwidth': 'float64','shift_factor': 'float64'})
    
    # do_optimal_fitparams(errors)
    # do_wide_errors(errors)
    # within_10_percent(errors)
    # do_avg_errors(errors, "mse")
    # do_avg_errors_diffcolor(errors, "mse")
    # do_avg_errors_all_diffcolor(errors)
    # do_avg_errors_all(errors)
    # do_spectrum_error_function(expt, pred, "Me2-Acr_1", "M06")
    # do_all_spectrum_error_function(expt, pred)
    # do_wide_errors_molecule(errors)
    do_optimal_fitparams_molecules(errors)
    # do_wide_errors_molecule2(errors)

    # to show a subset of these databases: errors[['functional','molecule']]
    # errors_r2 = errors[errors['error_function'] == 'r_square']
    # errors_r2_M06 = errors_r2[errors_r2['functional'] == 'M06']
    # example: errors_r2_M06[['functional','molecule','error_function','error_value']]
    
    # sns.scatterplot(x="error_value", y="molecule", data=errors_r2, hue="functional").legend(frameon=False,fontsize='12',ncol=1,bbox_to_anchor=(1.0,1.0))
    # sns.scatterplot(x="error_value", y="functional", data=errors_r2, hue="molecule",legend=False)
    
    
    
    # sns.displot(data=mse_errors, x="shift_factor", y="bandwidth", kind="kde")
    # plt.xlim(errors.shift_factor.min(),errors.shift_factor.max()) 
    # plt.xlim(0.575,1.475)
    # plt.ylim(0,0.55)
    
    # plt.setp(s.get_legend().get_texts(), fontsize='22')
    finish = time.perf_counter()
    print(f'Finished in {round(finish-start, 2)} second(s)')
