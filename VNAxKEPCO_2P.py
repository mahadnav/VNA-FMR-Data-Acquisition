import os
import time
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import data_container

from RS_VNA_ZVA40 import RS_VNA_Z, zva40_VNA
from bop50_8d import KEPCO_BOP


def Plot_dPdH(Data):
    f = plt.figure('VNA-FMR 2P dP/dH', (6, 4))
    
    minfield = Data['h'].min()
    maxfield = Data['h'].max()
    freq = Data['f']

    if not(f.axes):
        plt.subplot()
    ax = f.axes[0]

    if not(ax.lines):
        ax.plot([],[],'go-', markersize=5)
        ax.set_xlim([minfield, maxfield])
        ax.set_ylim([-1E-10, 1E-10])
        
    line = ax.lines[-1]
    line.set_data(Data['h'], Data['dP/dH']*1000)
    ax.set_xlabel('Field (Oe)')
    ax.set_ylabel('dP/dH')
    ax.grid(True)

    #check Y scale
    ymax = np.nan_to_num(Data['dP/dH']).max()*1000
    ymin = np.nan_to_num(Data['dP/dH']).min()*1000
    dy  = ymax - ymin
    yc = (ymax + ymin)/2
    ymin, ymax = ax.get_ylim()
    ymax = np.max([yc + dy*1.1/2, ymax])
    ymin = np.min([yc - dy*1.1/2, ymin])
    ax.set_ylim([ymin, ymax])

    f.tight_layout()
    f.canvas.draw()
    
    plotname = f'dPdH Plot {round(freq/1e9)} GHz @ {minfield} - {maxfield} Oe'
    f.savefig(rf'C:\Users\plyslab\Desktop\VNA Automation\Figs\{plotname}.png')

def Plot_ColorMap(Data):
        f = plt.figure('VNA-FMR-2P', (5, 4))

        minfield = Data['h'].min()
        maxfield = Data['h'].max()
        minfreq = Data['f'].min()/1E9
        maxfreq = Data['f'].max()/1E9

        extent = np.array([minfield,
                            maxfield,
                            minfreq,
                            maxfreq])

        if not(f.axes):
            plt.subplot()
        ax = f.axes[0]
        ax.clear()
        ax.imshow(Data['ColorMap'].T,
                aspect='auto',
                origin='lower',
                cmap='coolwarm',
                extent=extent)
        
        ax.set_title(f'{minfreq}-{maxfreq}GHz @ {minfield}-{maxfield} Oe')
        ax.set_xlabel('Field (Oe)')
        ax.set_ylabel('Freq (GHz)')
        plt.grid(False)
    
        f.tight_layout()
        f.canvas.draw()
        
        plotname = f'Frequency Sweep {minfreq} - {maxfreq} GHz @ {minfield} - {maxfield} Oe'
        f.savefig(rf'C:\Users\plyslab\Desktop\VNA Automation\Figs\{plotname}.png')

class Experiment2P():
    def __init__(self, logFilePath=None):
        if logFilePath is None:
            if not os.path.isdir(os.path.abspath('./Experiment_Logs')):
                os.mkdir(os.path.abspath('./Experiment_Logs'))
            logFilePath = './Experiment_Logs/FMR_log_{}.log'.format(self._get_timestring())
        with open(logFilePath, 'w') as log:
            log.write('SpinLab Instruments LogFile @ {}'.format(datetime.utcnow()) + '\n')
        self._logFile = os.path.abspath(logFilePath)
        self._logWrite('OPEN_')

        self.vna = RS_VNA_Z(logFile=self._logFile)
        self.zva = zva40_VNA(self.vna)
        self.PS = KEPCO_BOP(logFile=self._logFile)
        
        # Experimental data
        self.Data = data_container.DataContainer()
        self.Data.file_id = '.VNA_2P_Raw' #S11 vs hs vs fs

        self.ColorMapData = data_container.DataContainer()
        self.ColorMapData.file_id = '.VNA_ColorMap' #PAbs vs hs vs fs
        
        self.Data_Osc = data_container.DataContainer()
        self.Data_Osc.file_id = '.VNA_Osc_2P_Raw'  # S11 vs h  vs hosc
        self.Data_Osc_Info = 'VNA Osc Data'

        self.Data_dPdH = data_container.DataContainer()
        self.Data_dPdH.file_id = '.VNA_2P_dPdH'  # dP/dH vs h (fixed freq)
        self.dPdH_Info = 'VNA 2P dP/dH Data'
        
        self.PS.CurrentMode()
        self.PS.voltage = 40
        self.PS.current = 0
        self.HperOut = 663
        
        self.channel1 = self.zva.VNA.Ch1

        self._welcome()
        
    def __del__(self):
        self._logWrite('CLOSE')
        del self.PS
        del self.vna
        del self.zva

    def __str__(self):
        return 'VNA 2-Port FMR Experiment @ ' + datetime.utcnow()
    
    def _logWrite(self, action, value=''):
            if self._logFile is not None:
                with open(self._logFile, 'a') as log:
                    timestamp = datetime.utcnow()
                    log.write('%s %s : %s \n' % (timestamp, action, repr(value)))
    _log = _logWrite
    
    def _welcome(self):
        print("Welcome to the VNA 2-Port FMR Experiment!")
        print("Here are some default experiment parameters.\n")
        self._print_parameters()
        
    def _print_parameters(self):
        parameters = {
            'PS Output Current (A)': self.PS.MeasuredCurrent,
            'PS Output Voltage (V)': self.PS.MeasuredVoltage,
            'PS Output Mode (Current/Voltage)': self.PS.OperationMode,
            'Log File': self._logFile}
        for key, val in parameters.items():
            print(key, ':\t', val)
        
    def _get_timestring(self):
        now = datetime.now()
        return '{}-{}-{}_{}-{}-{}'.format(now.year, now.month, now.day, now.hour, now.minute, now.second)
        
    def current2field(self, current):
            return current * self.HperOut
        
    def field2current(self, field):
            return field / self.HperOut
        
    def PlotColorMap(self, i=None):
        Pabs_ref = 1 \
                   - np.abs(self.Data['S11_Ref'])**2 \
                   - np.abs(self.Data['S21_Ref'])**2

        if i is not None:
            # Update only i column
            Pabs = 1 \
                   - np.abs(self.Data['S11'][i])**2 \
                   - np.abs(self.Data['S21'][i])**2
            if self.Data['h'][0] > self.Data['h'][-1]:
                i = -1 - i
            self.ColorMapData['ColorMap'][i] = Pabs - Pabs_ref
        else:
            Pabs = 1 \
                   - np.abs(self.Data['S11'])**2 \
                   - np.abs(self.Data['S21'])**2
            self.ColorMapData['ColorMap'] = Pabs - Pabs_ref[None,:]
            if self.Data['h'][0] > self.Data['h'][-1]:
                self.ColorMapData['ColorMap'] = Pabs[::-1] - Pabs_ref[None,:]
        Plot_ColorMap(self.ColorMapData)
        
    def MeasureRef(self, ref_field, start_freq, stop_freq, num_points, file_name):
        current = self.field2current(ref_field)
        self.PS.set_current(current)
        time.sleep(0.5)
        
        channel1_instance = self.zva.VNA.Ch1
        channel1_instance.SetSweep(start_freq, stop_freq, num_points)
        self.zva.set_traces_SParameters_2P()
        
        self.Data['S11_Ref'] = self.zva.getSData(0, True)
        self.Data['S21_Ref'] = self.zva.getSData(1, False)
        self.Data['S22_Ref'] = self.zva.getSData(2, False)
        self.Data['S12_Ref'] = self.zva.getSData(3, False)
        time.sleep(5)
        self.PS.current = 0
        time.sleep(0)
        
        
        if file_name is not None:
            self.Data.save(file_name)

    def Measure(self, start_freq, stop_freq, num_points, fields, file_name):
        
        self.channel1.SetSweep(start_freq, stop_freq, num_points)
        self.zva.set_traces_SParameters_2P()
        
        frequencies = self.zva.frequencies
        
        self.Data['h'] = fields
        self.Data['f'] = frequencies
        data_shape = (len(fields), len(frequencies))
        
        self.Data['S11'] = np.zeros(data_shape, dtype=complex)
        self.Data['S21'] = np.zeros(data_shape, dtype=complex)
        self.Data['S22'] = np.zeros(data_shape, dtype=complex)
        self.Data['S12'] = np.zeros(data_shape, dtype=complex)
        
        self.ColorMapData['h'] = self.Data['h']
        self.ColorMapData['f'] = self.Data['f']
        self.ColorMapData['ColorMap'] = np.zeros(data_shape, dtype=float)
        self.ColorMapData['ColorMap'] += np.nan
    
        
        # Loop for each field
        start = time.time()
        for i, h in enumerate(fields):
            current = self.field2current(h)
            self.PS.set_current(current)
            time.sleep(0.02)
            self.Data['S11'][i] = self.zva.getSData(0)
            self.Data['S21'][i] = self.zva.getSData(1)
            self.Data['S22'][i] = self.zva.getSData(2)
            self.Data['S12'][i] = self.zva.getSData(3)
            self.PlotColorMap(i)
            
            if i == 0:
                end = time.time()
                total_time = end - start
                est_time = len(fields) * total_time
                print(f'Estimated total time for data acquisition: {round(est_time/60, 1)} min')

        if file_name is not None:
            self.Data.save(file_name)
        self.PS.current = 0
        self.PS.BEEP()
        
        
    def PlotdPdH(self, i=None):
        import pandas as pd
        from scipy.signal import savgol_filter
        delta_H = self.Data_Osc['delta_H']
        Pabs_upper = 1 - np.abs(self.Data_Osc['S11_upper'])**2 - np.abs(self.Data_Osc['S21_upper'])**2
        Pabs_lower = 1 - np.abs(self.Data_Osc['S11_lower'])**2 - np.abs(self.Data_Osc['S21_lower'])**2
        A_Pabs = (Pabs_upper - Pabs_lower) / 2*delta_H
        
        A_Pabs_flat = A_Pabs.flatten()
        window_size = 10
        A_Pabs_smoothed = pd.Series(A_Pabs_flat).rolling(window=window_size, center=True).mean().values
        
        if i is not None:
            self.Data_dPdH['dP/dH'][i] = A_Pabs_smoothed[i]
        else:
            self.Data_dPdH['dP/dH'] = A_Pabs_smoothed
        Plot_dPdH(self.Data_dPdH)
        
    
    def Measure_dPdH(self, freq, fields, file_name, delta_H=1, read_reps=3, returnData=None):

        self.zva.set_traces_SParameters_2P()
        self.channel1.SetSweep(freq, freq, num_points=1)
        
        self.Data_Osc['h'] = fields
        self.Data_Osc['f'] = freq
        self.Data_Osc['delta_H'] = delta_H
        
        data_shape = (len(fields), 1)
        self.Data_Osc['S11_upper'] = np.zeros(data_shape, dtype=complex)
        self.Data_Osc['S11_lower'] = np.zeros(data_shape, dtype=complex)
        self.Data_Osc['S21_upper'] = np.zeros(data_shape, dtype=complex)
        self.Data_Osc['S21_lower'] = np.zeros(data_shape, dtype=complex)
        self.Data_Osc.info = self.Data_Osc_Info

        self.Data_dPdH['h'] = fields
        self.Data_dPdH['f'] = freq
        self.Data_dPdH['dP/dH'] = np.zeros_like(fields) + np.nan
        extra_info = ['',
                      'Frequency : %(f)0.6f GHz' % {'f':freq/1E9},
                      'delta H :%(delta_H)d' % {'delta_H':delta_H},
                      '']
        self.Data_dPdH.info = self.dPdH_Info + '\n'.join(extra_info)
        
        start = time.time()
        for i, h in enumerate(fields):
            S11repArray_lower = []
            S21repArray_lower = []
            S11repArray_upper = []
            S21repArray_upper = []
            
            
            time.sleep(0.1) 
            for reps in range(read_reps):
                current_upper = self.field2current(h) + self.field2current(delta_H)
                current_lower = self.field2current(h) - self.field2current(delta_H)
                
                # H0 - delta H
                self.PS.set_current(current_lower)
                S11_lower = self.zva.getSData(0)
                S21_lower = self.zva.getSData(1)
                S11repArray_lower.append(S11_lower)
                S21repArray_lower.append(S21_lower)
                
                # H0 + delta H
                self.PS.set_current(current_upper)
                S11_upper = self.zva.getSData(0)
                S21_upper = self.zva.getSData(1)
                S11repArray_upper.append(S11_upper)
                S21repArray_upper.append(S21_upper)
                
            self.Data_Osc['S11_lower'][i] = np.mean(S11repArray_lower)
            self.Data_Osc['S11_upper'][i] = np.mean(S11repArray_upper)
            self.Data_Osc['S21_lower'][i] = np.mean(S21repArray_lower)
            self.Data_Osc['S21_upper'][i] = np.mean(S21repArray_upper)
                
            self.PlotdPdH(i)
            
            if i == 0:
                end = time.time()
                total_time = end - start
                est_time = len(fields) * total_time
                print(f'Estimated total time for data acquisition: {round(est_time/60, 1)} min')
                
        plt.show()

        if file_name is not None:
            self.Data_Osc.save(file_name)
            self.Data_dPdH.savetxt(file_name + '.txt', keys=['h', 'dP/dH'])

        self.PS.current = 0
        self.PS.BEEP()
        
        if returnData is not None:
            return self.Data_dPdH['dP/dH']