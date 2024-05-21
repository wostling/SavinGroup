#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  mass_scan.py
#  
#  Copyright 2020 Daniel Schury <daniel.schury@columbia.edu>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  
import sys                                                                                # obligatory sys unit
import argparse                                                                            # commandline parameters
import matplotlib.pyplot as plt                                                        # separate plots
from cycler import cycler                                                                # own color cycler
from lmfit.models import LorentzianModel, GaussianModel, LinearModel, ConstantModel    # fitting
import numpy as np                                                                        # math
from decimal import Decimal                                                            # get exponent base 10
import PySimpleGUI as sg                                                                # simple gui
import threading                                                                        # threading

own_color_cycler = cycler('color', ['#e66101', '#fdb863', '#b2abd2', '#5e3c99', '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c'])
params = {
    
    'axes.prop_cycle': own_color_cycler,
}
plt.rcParams.update(params)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

#search for highest point where to put label text
def text_index(x, y, index):
    index = x.index(index)
    return max(y[index],y[index+1],y[index-1])

def prepare_gauss_model(idx, x, y, fid, fit_model):
    Prefix = "g"+str(idx)+"_"
    if fit_model == "Gauss":
        gauss  = GaussianModel(prefix=Prefix)
    else:
        gauss  = LorentzianModel(prefix=Prefix)
    if idx == 1:
        center = -3
    elif idx == 2:
        center = 3
    pars = gauss.make_params()
    pars[Prefix+"center"].set(center)
    if fid:
        pars[Prefix+"sigma"].set(0.1)
        pars[Prefix+"amplitude"].set(0.2)
    else:
        pars[Prefix+"sigma"].set(0.1)
        pars[Prefix+"amplitude"].set(0.5)
    return gauss, pars

def DoFitting(x,y, fid, fit_model):
    #initialize parameter
    gmod = ConstantModel()
    pars = gmod.make_params()
    pars["c"].set(0.01)
    for i in [1,2]:
        gmod_dummy, pars_dummy = prepare_gauss_model(i, x, y, fid, fit_model)
        gmod += gmod_dummy
        pars += pars_dummy
    #do the fitting
    fgauss = gmod.fit(y, pars, x=x)
    #prepare the fitted gauss
    #plt.plot(x, fgauss.best_fit, label="Fiducial Fit")
    if fit_model == "Lorentz":
        label = "Lorentzian Fit"
    else:
        label = "Gaussian Fit"
    if not fid:
        plt.plot(x, fgauss.best_fit, label=label)
    return fgauss.best_values

def fexp(number):
    (sign, digits, exponent) = Decimal(number).as_tuple()
    return len(digits) + exponent - 1

def fman(number):
    return Decimal(number).scaleb(-fexp(number)).normalize()

def plot_spectrum(data, filename, second, fit_model, baseline):
    path = filename
    # initiliase the figure
    fig, ax = plt.subplots(1, 1)
    # plot the measured values
    if filename.find("leg2") >= 0:
        leg = 2
    else:
        leg = 1
    fit_data = ["Leg {}".format(leg)]
    ax.plot(data[0], data[1], label="BPM {}".format(leg))
    if second == "fiducial":
        ax.plot(data[2], data[3], label="Fid {}".format(leg))
    elif second == "both":
        ax.plot(data[2], data[3], label="BPM 2")
    # fit fiducials
    if second != "both":
        params = DoFitting(data[2], data[3], True, fit_model)
        center1 = params["g1_center"]
        center2 = params["g2_center"] 
        # fit peaks
        params = DoFitting(data[0], data[1], False, fit_model)
        if baseline:
            y = [params["c"], params["c"]]
            x = [data[0][0], data[0][-1]]
            ax.plot(x, y, label="Baseline")
            fit_data.append("Baseline: {:.2e} A".format(params["c"]))
            print(fit_data[-1])
        if fit_model == "Gauss":
            sigma1 = params["g1_sigma"]*10
            sigma2 = params["g2_sigma"]*10
            fit_data.append("Sigma Horizontal: {:.3f} mm\nSigma Vertical: {:.3f} mm".format(sigma1, sigma2))
            print(fit_data[-1])
        if fit_model == "Lorentz":
            fwhm1 = 2*params["g1_sigma"]*10
            fwhm2 = 2*params["g2_sigma"]*10
        else:
            fwhm1 = 2*np.sqrt(2*np.log(2))*params["g1_sigma"]*10
            fwhm2 = 2*np.sqrt(2*np.log(2))*params["g2_sigma"]*10
        fit_data.append("FWHM Horizontal: {:.3f} mm\nFWHM Vertical: {:.3f} mm".format(fwhm1, fwhm2))
        print(fit_data[-1])
        offset1 = center1 - params["g1_center"]
        offset1 *= 10
        offset2 = center2 - params["g2_center"]
        offset2 *= 10
        fit_data.append("Offset Horizontal: {:.3f} mm\nOffset Vertical: {:.3f} mm".format(offset1, offset2))
        title = "Horizontal: FWHM = {:.3f} mm    Offset = {:.3f} mm".format(fwhm1, offset1)
        title += "\nVertical: FWHM = {:.3f} mm    Offset = {:.3f} mm".format(fwhm2, offset2)
        ax.set_title(title)
        print(fit_data[-1])
        if fit_model == "Lorentz":
            path = path.replace(".pdf", "_lorentz.pdf")
        else:
            path = path.replace(".pdf", "_gauss.pdf")
        # get indices for offset and fwhm
    # configure plot
    plt.legend()
    ax.grid(linestyle=':', linewidth='0.5', color='black')
    ax.set_ylabel("Intensity (a.u.)")
    ax.set_xlabel("Position (cm)")
    ax.set_xlim(left=-6, right=6)
    plt.xticks(ticks=[-6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6])
    plt.tight_layout()
    plt.savefig(path)
    plt.close()
    with open(str.replace(path, ".pdf", ".dat"), "w") as filehandle:
        filehandle.writelines("%s\n" % place for place in fit_data)
    #return [offset1, offset2, fwhm1, fwhm2]

def analyze_spectrum_thread(window, datafile, fiducial, fit_model, baseline, commandline=False):
    print("Processing file {}".format(datafile))
    with open(datafile, "r") as file:
        fileinput = file.read().splitlines()
    # Convert list to matrix
    data = list()
    for i in fileinput:
        data.append(i.split())
    # Divide in x and y values
    bpm1_x = []
    bpm1_y = []
    fid1_x = []
    fid1_y = []
    bpm2_x = []
    bpm2_y = []
    fid2_x = []
    fid2_y = []
    for line in data:
        bpm1_x.append(float(line[0]))
        bpm1_y.append(float(line[1]))
        fid1_x.append(float(line[2]))
        fid1_y.append(float(line[3]))
        bpm2_x.append(float(line[4]))
        bpm2_y.append(float(line[5]))
        fid2_x.append(float(line[6]))
        fid2_y.append(float(line[7]))
        if bpm1_y[-1] == 0:
            bpm1_x.pop(-1)
            bpm1_y.pop(-1)
            fid1_x.pop(-1)
            fid1_y.pop(-1)
        if bpm2_y[-1] == 0:
            bpm2_x.pop(-1)
            bpm2_y.pop(-1)
            fid2_x.pop(-1)
            fid2_y.pop(-1)
    print("Start plotting...")
    # leg 1
    if len(bpm1_x) != 0:
        data = [bpm1_x, bpm1_y, fid1_x, fid1_y]
        # Get base path
        path = datafile.replace(".dat", "_leg1.pdf")
        print("Leg 1:")
        if fiducial:
            second = "fiducial"
        else:
            second = "none"
        plot_spectrum(data, path, second, fit_model, baseline)
    # leg 2
    if len(bpm2_x) != 0:
        data = [bpm2_x, bpm2_y, fid2_x, fid2_y]
        # Get base path
        path = datafile.replace(".dat", "_leg2.pdf")
        print("Leg 2:")
        plot_spectrum(data, path, second, fit_model, baseline)
    # both legs
    if (len(bpm1_x) != 0) and (len(bpm2_x) != 0):
        data = [bpm1_x, bpm1_y, bpm2_x, bpm2_y]
        path = datafile.replace(".dat", "_both.pdf")
        plot_spectrum(data, path, "both", fit_model, baseline)
    if not commandline:
        window.write_event_value('-THREAD DONE-', '')

def analyze_spectrum(window, filename, fiducial, fitmodel, baseline, commandline=False):
    threading.Thread(target=analyze_spectrum_thread, args=(window, filename, fiducial, fitmodel, baseline), daemon=True).start()

def main(args):
    if not args.commandline:
        sg.theme("SystemDefaultForReal")
        layout = [
            [sg.Text('Select BPM profile data file to analyze', size=(30, 1), font=("Sans", 25))],
            [sg.Text('Profile Data'),
             sg.InputText(args.scan),
             sg.FilesBrowse(key="-FILENAME-", file_types=(("BPM Profiles", "*.dat"),))
            ],
            [sg.Frame(layout=[
                [sg.Checkbox('Fiducial', default=args.fiducial, key="-FIDUCIAL-"),  sg.Checkbox('Baseline', default=args.baseline, key="-BASELINE-")],
                [sg.Text('Choose Fit Model')],
                [sg.InputOptionMenu(('Gauss', 'Lorentz'), key="-FITMODEL-")],
                ], title='Options', relief=sg.RELIEF_SUNKEN, tooltip='Use these to set flags'),
                sg.Ok(tooltip='Click to analyze the data', key="-BOK-"), sg.Text("", size=(30,1), key="-STATUS-"),
            ],
            [sg.Text('Leave GUI: ', size=(45,1)), sg.Button('LEAVE', size=(15,1)), sg.Text('\n')],
            [sg.Text('Stop programme: ', size=(45,1)), sg.Button('CANCEL', size=(15,1)), sg.Text('\n')]
        ]
        window = sg.Window('BPM Profiles', layout, default_element_size=(40, 1), grab_anywhere=False).Finalize()
        if args.lorentz:
            window["-FITMODEL-"].update(values=("Lorentz", "Gauss"))
        while True:
            event, values = window.read()
            
            files = str(values["-FILENAME-"])
            file_list= files.split(";")
            num = len(file_list)
        
            if event in ('LEAVE'):
                break
            if event in (sg.WIN_CLOSED, 'CANCEL'):
                exit()
            if event == "-BOK-":
                window["-STATUS-"].update("Analyzing…")
                window["-BOK-"].update(disabled=True)
                window.refresh()
                for i in range(num):
                    filepath = file_list[i]
                    analyze_spectrum(window, filepath, values["-FIDUCIAL-"], values["-FITMODEL-"], values["-BASELINE-"])
            elif event == '-THREAD DONE-':
                window["-STATUS-"].update("Done…")
                window["-BOK-"].update(disabled=False)
                window.refresh()
        window.close()
    else:
        if args.lorentz:
            fitmodel = "Lorentz"
        else:
            fitmodel = "Gauss"
        analyze_spectrum_thread(None,  args.scan,  args.fiducial, fitmodel,  args.baseline, args.commandline)
    return 0

if __name__ == '__main__':
    #set commandline parameters
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--scan", type=str, default="", help="BPM scan file")
    parser.add_argument("-f", "--fiducial", action="store_true", help="plot fiducials")
    parser.add_argument("-b", "--baseline", action="store_true", help="plot baseline")
    parser.add_argument("-l", "--lorentz", action="store_true", help="fit Lorentz curve")
    parser.add_argument("-c", "--commandline", action="store_true", help="no gui commandline only")
    args = parser.parse_args()
    sys.exit(main(args))
