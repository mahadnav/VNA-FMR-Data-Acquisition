import time as _time
import numpy as _np
import pandas as pd
from instrument_base import InstrumentBase as _InstrumentBase
from instrument_base import InstrumentChild as _InstrumentChild


class RS_VNA_Z(_InstrumentBase):
    def __init__(self, GPIB_Address=20, GPIB_Device=0, ResourceName=None, logFile=None):
        if ResourceName is None:
            ResourceName = 'GPIB%d::%d::INSTR' % (GPIB_Device, GPIB_Address)
        super().__init__(ResourceName, logFile)
        self._IDN = 'R&S VNA'
        self.write('*CLS')
        self.write('SYST:COMM:GPIB:RTER EOI')
        self.VI.write_termination = None
        self.VI.read_termination = None
        self.write('FORM:DATA REAL, 64')
        self.values_format.is_binary = True
        self.values_format.datatype = 'd'  # float 64 bits
        self.values_format.is_big_endian = False
        self.Ch1 = Channel(self, 1)
        #  self.Ch2 = Channel(self, 2)

    def query(self, command):
        # Remove \n terminations by software
        returnQ = super().query(command)
        return returnQ.strip('\n')



class Channel(_InstrumentChild):
    def __init__(self, parent, ChanNum=1, ID='Auto'):
        super().__init__(parent)
        self._number = ChanNum
        if ID == 'Auto':
            self.ID = 'Ch%d' % self.number
        else:
            self.ID = ID
        self.write('CONF:CHAN%d:STAT ON' % self.number)

    @property
    def number(self):
        return self._number

    @property
    def traces(self):
        TrcNames = self.query('CONF:TRAC:CAT?').strip('\'').split(',')[1::2]
        vTraces = []
        for trN in TrcNames:
            tr = Trace(self, trN)
            if tr.channel_number == self.number:
                vTraces.append(tr)
        return vTraces

    @property
    def bandwidth(self):
        return self.query_float('SENS%(Ch)d:BWID?' % {'Ch': self.Number})

    @bandwidth.setter
    def bandwidth(self, newBW):
        self.write('SENS%(Ch)d:BWID %(BW)0.9E' %
                   {'Ch':self.Number, 'BW':newBW})

    def SetSweep(self, start, stop, num_points, na=None):
        self.write('SENS%(Ch)d:SWE:TYPE LIN' %
                   {'Ch':self.number})

        self.write('SENS%(Ch)d:SWE:POIN %(n)d' %
                   {'Ch':self.number, 'n':num_points})
        self.write('SENS%(Ch)d:FREQ:STAR %(f)0.9E' %
                   {'Ch':self.number, 'f':start})
        self.write('SENS%(Ch)d:FREQ:STOP %(f)0.9E' %
                   {'Ch':self.number, 'f':stop})

        if na is not None:
            self.write('SENS%(Ch)d:AVER:COUN %(n)d' %
               {'Ch':self.number, 'n':na})
            if na > 1:
                self.write('SENS%(Ch)d:AVER:STAT ON' % {'Ch':self.number})
            else:
                self.write('SENS%(Ch)d:AVER:STAT OFF' % {'Ch':self.number})

    def SetFrequencies(self, fs):
        self.write('SENS%(Ch)d:SEGM1:CLE' % {'Ch': self.number})

        fs = _np.atleast_1d(fs).copy()
        fs.sort()
        if len(fs)%2 != 0:
            fs = _np.r_[fs, [fs[-1]]]
        for sgn in range(len(fs)//2):
            self.write('SENS%(Ch)d:SEGM%(sg)d:ADD' %
                       {'Ch': self.number, 'sg': sgn+1})
            self.write('SENS%(Ch)d:SEGM%(sg)d:FREQ:START %(f)0.9E' %
                       {'Ch': self.number, 'sg': sgn+1, 'f':fs[2*sgn]})
            self.write('SENS%(Ch)d:SEGM%(sg)d:FREQ:STOP %(f)0.9E' %
                       {'Ch': self.number, 'sg': sgn+1, 'f':fs[2*sgn+1]})
            npts = 2
            if fs[2*sgn+1] == fs[2*sgn]:
                npts = 1
            self.write('SENS%(Ch)d:SEGM%(sg)d:SWE:POIN %(npts)d' %
                       {'Ch':self.number, 'sg':sgn+1, 'npts':npts})
        self.write('SENS%(Ch)d:SWE:TYPE SEGM' %
                   {'Ch':self.number})


    def getSTIM(self):
        return self.query_values('CALC%(Ch)d:DATA:STIM?' % {'Ch': self.number})

    def INIT(self):
        while True:
            self.write('INIT%(Ch)d:IMM' % {'Ch': self.number})
            if self.query('SYST:ERR:ALL?') == '0,"No error"':
                break


class Trace(_InstrumentChild):
    def __init__(self, parent, Name='Auto'):
        super().__init__(parent)
        self.name = Name
        self._ChNumber = self.query_int('CONF:TRAC:CHAN:NAME:ID? \'%s\''
                                        % Name)

    @property
    def channel_number(self):
        '''Channel Number'''
        return self._ChNumber

    def getNewData(self):
        ChN = self.channel_number
        self.write('INIT:CONT OFF')
        averStat = self.query('SENS%(Ch)d:AVER:STAT?' % {'Ch': ChN})
        if averStat == '1':
            averCount = self.query_int('SENS%(Ch)d:AVER:COUN?' % {'Ch': ChN})
            self.write('SENS%(Ch)d:AVER:CLE' % {'Ch': ChN})
        else:
            averCount = 1 
        for i in range(averCount):
            self.parent.INIT()
            while self.query('CALC%(Ch)d:DATA:NSW:COUN?' % {'Ch': ChN}) == '0':
                _time.sleep(0.01)

    def getFDAT(self, New=False):
        '''
        Return formatted trace data,
        accordingly to the selected trace format
        '''
        ChN = self.channel_number
        self.write('CALC%(Ch)d:PAR:SEL \'%(N)s\'' %
                   {'Ch': ChN, 'N': self.name})
        if New:
            self.getNewData()
        return self.query_values('CALC%(Ch)d:DATA? FDAT' % {'Ch': ChN})

    def getSDAT(self, New=False):
        '''
        Returns unformatted trace data :
        Real and imaginary part of each measurement point
        For wave quantities the unit is Volts
        '''
        ChN = self.channel_number
        self.write('CALC%(Ch)d:PAR:SEL \'%(N)s\'' %
                   {'Ch': ChN, 'N': self.name})
        if New:
            self.getNewData()
        SDAT = self.query_values('CALC%(Ch)d:DATA? SDAT' % {'Ch': ChN})
        return SDAT[::2] + 1.0j*SDAT[1::2]

    def SaveSData(self, fileName, New=False):
        f = self.parent.getSTIM()
        S = self.getSDAT(New)
        df = pd.DataFrame(data={'Frequency': f/1E9, 'S_Real': S.real, 'S_img': S.imag})
        df.to_csv(fileName, index=False, header=True)
        
    def SaveFData(self, fileName, New=False):
        f = self.parent.getSTIM()
        F = self.getFDAT(New)
        formatted = _np.array(F).flatten()
        df = pd.DataFrame(data={'Frequency': f/1E9, 'Formatted Data': formatted})
        df.to_csv(fileName, index=False, header=True)
        
        
        
        
class zva40_VNA(object):

    def __init__(self, VNA_instrument):
        self.VNA = VNA_instrument
    
    @property
    def frequencies(self):
        return self.VNA.Ch1.getSTIM()

    def cleanDisplayArea(self):
        self.VNA.write('SYST:DISP:UPDATE ON')
        areas = self.VNA.query('DISP:CAT?').strip('\'').split(',')[::2]
        if not '1' in areas:
            self.VNA.write('DISP:WIND1:STAT ON')
        for area in areas:
            if area != '1':
                self.VNA.write('DISP:WIND%s:STAT OFF' % area)
        self.VNA.write('CALC1:PAR:DEL:ALL')

    def set_traces_SParameters_1P(self):
        self.cleanDisplayArea()

        self.VNA.write('CALC1:PAR:DEF \'Trc1\', S11')
        self.VNA.write('CALC1:PAR:SEL \'Trc1\'')
        self.VNA.write('CALC1:FORM MLIN')
        self.VNA.write('DISP:WIND1:TRAC1:FEED \'Trc1\'')
        self.VNA.write('DISP:WIND1:TRAC1:Y:PDIV 0.11')
        self.VNA.write('DISP:WIND1:TRAC1:Y:RPOS 50')
        self.VNA.write('DISP:WIND1:TRAC1:Y:RLEV 0.5')

    def set_traces_SParameters_2P(self):
        self.cleanDisplayArea()

        self.VNA.write('CALC1:PAR:DEF \'Trc1\', S11')
        self.VNA.write('CALC1:PAR:SEL \'Trc1\'')
        self.VNA.write('CALC1:FORM MLIN')
        self.VNA.write('CALC1:PAR:DEF \'Trc2\', S21')
        self.VNA.write('CALC1:PAR:SEL \'Trc2\'')
        self.VNA.write('CALC1:FORM MLIN')
        self.VNA.write('CALC1:PAR:DEF \'Trc3\', S22')
        self.VNA.write('CALC1:PAR:SEL \'Trc3\'')
        self.VNA.write('CALC1:FORM MLIN')
        self.VNA.write('CALC1:PAR:DEF \'Trc4\', S12')
        self.VNA.write('CALC1:PAR:SEL \'Trc4\'')
        self.VNA.write('CALC1:FORM MLIN')

        self.VNA.write('DISP:WIND1:TRAC1:FEED \'Trc1\'')
        self.VNA.write('DISP:WIND1:TRAC2:FEED \'Trc2\'')
        self.VNA.write('DISP:WIND1:TRAC3:FEED \'Trc3\'')
        self.VNA.write('DISP:WIND1:TRAC4:FEED \'Trc4\'')

        self.VNA.write('DISP:WIND1:TRAC1:Y:PDIV 0.11')
        self.VNA.write('DISP:WIND1:TRAC2:Y:PDIV 0.11')
        self.VNA.write('DISP:WIND1:TRAC3:Y:PDIV 0.11')
        self.VNA.write('DISP:WIND1:TRAC4:Y:PDIV 0.11')

        self.VNA.write('DISP:WIND1:TRAC1:Y:RPOS 50')
        self.VNA.write('DISP:WIND1:TRAC2:Y:RPOS 50')
        self.VNA.write('DISP:WIND1:TRAC3:Y:RPOS 50')
        self.VNA.write('DISP:WIND1:TRAC4:Y:RPOS 50')
        
        if self.VNA.query_float('DISP:WIND1:TRAC1:Y:RPOS?') != 50:
            self.VNA.write('DISP:WIND1:TRAC1:Y:RPOS 5')
            self.VNA.write('DISP:WIND1:TRAC2:Y:RPOS 5')
            self.VNA.write('DISP:WIND1:TRAC3:Y:RPOS 5')
            self.VNA.write('DISP:WIND1:TRAC4:Y:RPOS 5')

        self.VNA.write('DISP:WIND1:TRAC1:Y:RLEV 0.5')
        self.VNA.write('DISP:WIND1:TRAC2:Y:RLEV 0.5')
        self.VNA.write('DISP:WIND1:TRAC3:Y:RLEV 0.5')
        self.VNA.write('DISP:WIND1:TRAC4:Y:RLEV 0.5')

    def set_traces_WaveQuantities(self):
        self.cleanDisplayArea()

        self.VNA.write('CALC1:PAR:DEF \'Trc5\', R1')
        self.VNA.write('CALC1:FORM MLIN')
        self.VNA.write('CALC1:PAR:DEF \'Trc6\', R2')
        self.VNA.write('CALC1:FORM MLIN')
        self.VNA.write('CALC1:PAR:DEF \'Trc7\', A')
        self.VNA.write('CALC1:FORM MLIN')
        self.VNA.write('CALC1:PAR:DEF \'Trc8\', B')
        self.VNA.write('CALC1:FORM MLIN')

        self.VNA.write('DISP:WIND1:TRAC5:FEED \'Trc5\'')
        self.VNA.write('DISP:WIND1:TRAC6:FEED \'Trc6\'')
        self.VNA.write('DISP:WIND1:TRAC7:FEED \'Trc7\'')
        self.VNA.write('DISP:WIND1:TRAC9:FEED \'Trc8\'')

    def backup_sweep(self):
        freqStart = self.VNA.query('SENS1:FREQ:STAR?')
        freqStop = self.VNA.query('SENS1:FREQ:STOP?')
        sweepPoints = self.VNA.query('SENS1:SWE:POIN?')
        self._sweep_data = [freqStart, freqStop, sweepPoints]

    def restore_sweep(self):
        freqStart, freqStop, sweepPoints = self._sweep_data
        self.VNA.write('SENS1:SWE:POIN %s' %sweepPoints)
        self.VNA.write('SENS1:FREQ:STAR %s' %freqStart)
        self.VNA.write('SENS1:FREQ:STOP %s' %freqStop)
        
    def getSData(self, TrN, new=True):
        '''
        Get unformatted trace data:
        Real and imaginary part of each measurement point.
        2 values per trace point irrespective of the selected trace format
        '''
        X = self.VNA.Ch1.traces[TrN].getSDAT(new)
        return X.copy()

    def getFData(self, TrN, new=True):
        '''
        Get formatted trace data:
        1 values per trace point as the selected trace format
        '''
        X = self.VNA.Ch1.traces[TrN].getFDAT(new)
        return X.copy()
        
    def _MeasureSpectra(self):
        Ss = []
        S = self.VNA.Ch1.traces[0].getFDAT(True)
        Ss.append(S)
        return Ss