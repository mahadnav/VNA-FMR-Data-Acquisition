import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
import os
from scipy.ndimage import gaussian_filter, convolve
from scipy.signal import spline_filter, medfilt2d
from scipy.optimize import curve_fit
plt.style.use('seaborn-v0_8')

class filters:
    @staticmethod
    def Filter_Spline(data, a=0.1):
        """
        Apply a spline filter to the data.

        Parameters:
        - data: Input data array
        - a: Spline filter parameter (default is 0.1)

        Returns:
        Spline-filtered data.
        """
        import scipy.signal
        data = scipy.signal.spline_filter(data, a)
        return data

    @staticmethod
    def FilterBG_MedianH(data):
        """
        Subtract median background along the horizontal axis.

        Parameters:
        - data: Input data array

        Returns:
        Data with horizontal median background subtracted.
        """
        m = np.median(data, axis=0)
        bg = np.zeros_like(data)
        for i, x in enumerate(m):
            bg[:, i] = np.mean(data[:, i][data[:, i] < x])
        data = data - bg
        return data

    @staticmethod
    def FilterBG_Median(data):
        """
        Subtract median background along the vertical axis.

        Parameters:
        - data: Input data array

        Returns:
        Data with vertical median background subtracted.
        """
        m = np.median(data, axis=1)
        bg = np.zeros_like(data)
        for i, x in enumerate(m):
            bg[i] = np.mean(data[i][data[i] < x])
        data = data - bg
        return data

    @staticmethod
    def filter(fil, data):
        """
        Apply various filters to the data.

        Parameters:
        - fil: Filter type ('gaussian', 'median', 'spline', 'fft', 'None')
        - data: Input data array

        Returns:
        Filtered data.
        """
        if fil == 'gaussian':
            data = gaussian_filter(data, sigma=0.5)
        elif fil == 'median':
            data = medfilt2d(data, kernel_size=[1, 3])
        elif fil == 'spline':
            data = spline_filter(data, 10)
        elif fil == 'fft':
            fft_result = np.fft.fft2(data)
            fft_result[np.abs(fft_result) < np.percentile(np.abs(fft_result), 96)] = 0
            data = np.fft.ifft2(fft_result).real
        elif fil == 'None':
            data = data
        return data

    @staticmethod
    def laplacian_reg(data, reg_strength):
        """
        Apply Laplacian regularization to the data.

        Parameters:
        - data: Input data array
        - reg_strength: Regularization strength

        Returns:
        Regularized data.
        """
        laplacian = np.array([[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]])
        regularized_data = data + reg_strength * convolve(data, laplacian, mode='nearest')
        return regularized_data
    
class Crop:
    @staticmethod
    def cut_freq(freq, lower_freq, upper_freq):
        """
        Crop frequency array within specified range.

        Parameters:
        - freq: Original frequency array
        - lower_freq: Lower bound of the cropped range
        - upper_freq: Upper bound of the cropped range

        Returns:
        Cropped frequency array, lower index, and upper index.
        """
        lower_index = max(0, math.floor(len(freq) * (lower_freq - min(freq)) / (max(freq) - min(freq))))
        upper_index = min(math.ceil(len(freq) * (upper_freq - min(freq)) / (max(freq) - min(freq))), len(freq))
        
        if lower_freq == upper_freq:
            lower_index = np.clip(math.ceil(len(freq) * (lower_freq - min(freq)) / (max(freq) - min(freq))), 0, len(freq)-1)
            upper_index = lower_index

        return freq[lower_index:upper_index], lower_index, upper_index

    @staticmethod
    def cut_field(field, lower_field, upper_field):
        """
        Crop field array within specified range.

        Parameters:
        - field: Original field array
        - lower_field: Lower bound of the cropped range
        - upper_field: Upper bound of the cropped range

        Returns:
        Cropped field array, lower index, and upper index.
        """
        lower_index = math.floor(len(field) / (max(field) - min(field)) * (lower_field - min(field)))
        if lower_field != upper_field:
            upper_index = math.ceil(len(field) / (max(field) - min(field)) * (upper_field - min(field)))
        else:
            upper_index = lower_index + 1  # Set upper_index to the next index to avoid an empty array
        
        return field[lower_index:upper_index], lower_index, upper_index

class Plot:
    @staticmethod
    def find(condition):
        """
        Find indices where the condition is True.

        Parameters:
        - condition: Boolean array

        Returns:
        Array of indices where the condition is True.
        """
        res, = np.nonzero(np.ravel(condition))
        return res
    
    @staticmethod
    def plotH(data, H, f, field):
        """
        Plot a horizontal cut of data at a specific magnetic field.

        Parameters:
        - data: Data array
        - H: Magnetic field value for the cut
        - f: Frequency array
        - field: Magnetic field array

        Returns:
        Horizontal cut of data at the specified magnetic field.
        """
        hi = Plot.find(field >= H)[0]
        plt.plot(f, data[:, hi].real, linewidth=2, label=f'{H} Oe')
        plt.xlabel('Frequency (GHz)')

    @staticmethod
    def plotF(data, F, h, freq):
        """
        Plot a vertical cut of data at a specific frequency.

        Parameters:
        - data: Data array
        - F: Frequency value for the cut
        - h: Magnetic field array
        - freq: Frequency array

        Returns:
        Vertical cut of data at the specified frequency.
        """
        fi = Plot.find(freq >= F)[0]
        plt.plot(h, data[fi].real, linewidth=2, label=f'{F} GHz')
        plt.xlabel('Field (Oe)')
    
    @staticmethod
    def Hdata(data, H, field):
        """
        Extract the data of a horizontal cut at a specific magnetic field.

        Parameters:
        - data: Data array
        - H: Magnetic field value for the cut
        - f: Frequency array
        - field: Magnetic field array

        Returns:
        Horizontal cut of data at the specified magnetic field.
        """
        hi = Plot.find(field >= H)[0]
        return data[:, hi].real

    @staticmethod
    def Fdata(data, F, freq):
        """
        Extract the data of a vertical cut at a specific frequency.

        Parameters:
        - data: Data array
        - F: Frequency value for the cut
        - h: Magnetic field array
        - freq: Frequency array

        Returns:
        Vertical cut of data at the specified frequency.
        """
        fi = Plot.find(freq >= F)[0]
        return data[fi].real
    
    def shift_dPdH(directory_path):
        """
        Plot dP/dH data for multiple frequencies from CSV files in a directory.

        Parameters:
        - directory_path: Path to the directory containing CSV files.

        Returns:
        None.
        """
        csv_files = [file for file in os.listdir(directory_path) if file.endswith('.csv')]
        for csv_file in csv_files:
            csv_file_path = os.path.join(directory_path, csv_file)
            df = pd.read_csv(csv_file_path)[5:-4]
            frequency = csv_file.split('_')[5]
            plt.plot(df.iloc[:, 0], df.iloc[:, 1], 'o-', linewidth=3, label=fr'{frequency} GHz')
        plt.xlabel('Field (Oe)')
        plt.ylabel('dP/dH (a.u.)')
        plt.legend()
        
    def S_parameter(freq, field, specified_freq, S11, S21=False):
        ind = Crop.cut_freq(freq, specified_freq, specified_freq)[1]

        plt.figure(figsize=(8,3))
        plt.subplot(1,2,1)
        plt.plot(field, S11.imag[:, ind],'ro-', markersize=3, linewidth=1)
        plt.title('Imaginary')
        plt.xlabel('Field (Oe)')
        plt.ylabel('Im S11 (au)')
        plt.grid(False)

        plt.subplot(1,2,2)
        plt.plot(field, S11.real[:, ind], 'bo-', markersize=3, linewidth=0.5)
        plt.title('Real')
        plt.xlabel('Field (Oe)')
        plt.ylabel('Real S11 (au)')
        plt.grid(False)

        plt.subplots_adjust(left=0.2, right=0.9, bottom=1, top=1.5, wspace=0.6, hspace=1)
        
        if S21 is not False:
            plt.figure(figsize=(8,3))
            plt.subplot(1,2,1)
            plt.plot(field, S21.imag[:, ind],'ro-', markersize=3, linewidth=1)
            plt.title('Imaginary')
            plt.xlabel('Field (Oe)')
            plt.ylabel('Im S21 (a.u.)')
            plt.grid(False)

            plt.subplot(1,2,2)
            plt.plot(field, S21.real[:, ind], 'bo-', markersize=3, linewidth=0.5)
            plt.title('Real')
            plt.xlabel('Field (Oe)')
            plt.ylabel('Real S21 (a.u.)')
            plt.grid(False)

            plt.subplots_adjust(left=0.2, right=0.9, bottom=1, top=1.5, wspace=0.6, hspace=1)
        
class FWHM:
    @staticmethod
    def lin_interp(x, y, i, half):
        """
        Linear interpolation for finding half-max points.

        Parameters:
        - x, y: Input arrays
        - i: Index
        - half: Half-max value

        Returns:
        Interpolated x-coordinate for the half-max value.
        """
        return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

    @staticmethod
    def half_max_x(x, y):
        """
        Find x-coordinates of half-max points.

        Parameters:
        - x, y: Input arrays

        Returns:
        List of x-coordinates of half-max points and indices where crossing occurs.
        """
        arr = []
        half = max(y)/2.0
        signs = np.sign(np.add(y, -half))
        zero_crossings = (signs[0:-2] != signs[1:-1])
        zero_crossings_i = np.where(zero_crossings)[0]
        
        for i in range(len(zero_crossings_i)):
            arr.append(FWHM.lin_interp(x, y, zero_crossings_i[i], half))
        
        return arr, zero_crossings_i

    @staticmethod
    def fwhm(data, field, frequency, freq_used):
        """
        Calculate and plot the Full Width at Half Maximum (FWHM) for a given frequency.

        Parameters:
        - data: Data array
        - field: Magnetic field array
        - frequency: Frequency value
        - freq_used: Frequencies used in the data

        Returns:
        None.
        """
        x = field
        y = Plot.Fdata(data, F=frequency, freq=freq_used)

        # find the two crossing points
        hmx = FWHM.half_max_x(x,y)[0]
        index = FWHM.half_max_x(x,y)[1]
        plt.figure(figsize=(6,4))
        plt.plot(x, y, 'k', linewidth=3)
        plt.title(f"Full Width Half Maximum at {frequency} GHz")
        for i in range(0, len(hmx), 2):
            fwhm = abs(hmx[i] - hmx[i + 1])
            print(f"FWHM: {fwhm}")
            plt.plot([hmx[i], hmx[i+1]], [(max(y))/2, (max(y))/2], 'ro-', markersize=4)
        plt.ylabel('$P_{abs}$ (a.u.)')
        plt.xlabel('Field (Oe)')
        plt.show()
        
class Calc:
    @staticmethod
    def AbsPower1P(S11, S11_ref, bg='Ref'):
        """
        Calculate absorbed power for 1-port measurements.

        Parameters:
        - S11: S11 data array
        - S11_ref: Reference S11 data array
        - bg: Background subtraction method ('Ref' or 'None', default is 'Ref')

        Returns:
        Array of absorbed power values.
        """
        AbsPower = 1 - np.abs(S11)**2
        m = np.median(AbsPower, axis=1)
        if bg == 'Ref':
            bg = 1 - np.abs(S11_ref)**2
            outArray = AbsPower - bg[:, None].T
        else:
            bg = np.zeros_like(m)
            for i, x in enumerate(m):
                bg[i] = np.mean(AbsPower[i][AbsPower[i] < x])
                outArray = AbsPower - bg[:, np.newaxis]    
        return outArray.T
    
    @staticmethod
    def AbsPower2P(S11, S21, S11_Ref, S21_Ref, bg='Ref'):
        """
        Calculate absorbed power for 2-port measurements.

        Parameters:
        - S11: S11 data array
        - S21: S21 data array
        - S11_Ref: Reference S11 data array
        - S21_Ref: Reference S21 data array
        - bg: Background subtraction method ('Ref' or 'None', default is 'Ref')

        Returns:
        Array of absorbed power values.
        """
        AbsPower = 1 - np.abs(S11)**2 - np.abs(S21)**2
        m = np.median(AbsPower, axis=1)
        if bg == 'Ref':
            bg = 1 - np.abs(S11_Ref)**2 - np.abs(S21_Ref)**2
            outArray = AbsPower - bg[:, None].T
        else:
            bg = np.zeros_like(m)
            for i, x in enumerate(m):
                bg[i] = np.mean(AbsPower[i][AbsPower[i] < x])
                outArray = AbsPower - bg[:, np.newaxis]    
        return outArray.T
    
    def concentrated_dpdh(data, frequency, field_used, freq_used, field_below, field_above):
        """
        Extract a concentrated range of magnetic fields for a given frequency.

        Parameters:
        - data: Data array
        - frequency: Frequency value
        - field_used: Magnetic field array
        - freq_used: Frequencies used in the data
        - field_below: Range below the maximum field
        - field_above: Range above the maximum field

        Returns:
        Lower and upper bounds for the concentrated range of magnetic fields.
        """
        d = Plot.Fdata(data, F=frequency, freq=freq_used)
        max_field = field_used[np.argmax(d)]
        lower = max_field - field_below
        if lower < 0:
            lower = 0
        upper = max_field + field_above
        return lower, upper

class Load:
    @staticmethod
    def NPZ_1P(fileName):
        """
        Load data from a 1-port NPZ file.

        Parameters:
        - fileName: Name of the NPZ file without extension

        Returns:
        Tuple containing information string, frequency array, magnetic field array,
        reference S11 data array, and S11 data array.
        """
        npz_file = np.load(fileName + '.VNA_1P_Raw.npz')
        Info = str(npz_file['Info'])
        f = npz_file['f']
        h = npz_file['h']
        S11_Ref = npz_file['S11_Ref']
        S11 = npz_file['S11']
        return Info, f/1e9, h, S11_Ref, S11
    
    @staticmethod
    def NPZ_2P(fileName):
        """
        Load data from a 2-port NPZ file.

        Parameters:
        - fileName: Name of the NPZ file without extension

        Returns:
        Tuple containing information string, frequency array, magnetic field array,
        reference S11, S21, S22, and S12 data arrays, and measured S11, S21, S22, and S12 data arrays.
        """
        npz_file = np.load(fileName + '.VNA_2P_Raw.npz')
        Info = str(npz_file['Info'])
        f = npz_file['f']
        h = npz_file['h']
        S11_Ref = npz_file['S11_Ref']
        S21_Ref = npz_file['S21_Ref']
        S22_Ref = npz_file['S22_Ref']
        S12_Ref = npz_file['S12_Ref']
        S11 = npz_file['S11']
        S12 = npz_file['S12']
        S21 = npz_file['S21']
        S22 = npz_file['S22']
        return Info, f/1e9, h, S11_Ref, S21_Ref, S22_Ref, S12_Ref, S11, S12, S21, S22
    
    def NPZ_Fit(fileName):
        """
        Load fitted data from an NPZ file.

        Parameters:
        - fileName: Name of the NPZ file without extension

        Returns:
        Tuple containing frequency, field array, raw data, fitted data.
        """
        npz_file = np.load(fileName)
        f = npz_file['f']
        h = npz_file['field']
        raw = npz_file['Raw']
        fitted = npz_file['Fitted']
        return f, h, raw, fitted

class Fitting:
    @staticmethod
    def X(A, G, x0, x):
        """
        Lorentzian component function.

        Parameters:
        - A: Amplitude
        - G: DeltaH (Full Width at Half Maximum)
        - x0: Center position
        - x: Input array

        Returns:
        Calculated Lorentzian component values.
        """
        return A * (G**2 / (4*(x - x0)**2 + G**2)) + 1E-35
    
    @staticmethod
    def lorentzian_function(x, A1, x01, G1, A2, x02, G2, offset):
        """
        Combined Lorentzian function with two components.

        Parameters:
        - x: Input array
        - A1, x01, G1: Parameters for the first Lorentzian component
        - A2, x02, G2: Parameters for the second Lorentzian component
        - offset: Offset parameter

        Returns:
        Calculated combined Lorentzian values.
        """
        return Fitting.X(A1, G1, x01, x) + Fitting.X(A2, G2, x02, x) + offset
    
    @staticmethod
    def fit_lorentzian(x, y, initial_guess):
        """
        Fit the Lorentzian function to data.

        Parameters:
        - x, y: Input data
        - initial_guess: Initial guess for the fitting parameters

        Returns:
        Fitted parameters if successful, None otherwise.
        """
        try:
            params, _ = curve_fit(Fitting.lorentzian_function, x, y, p0=initial_guess)
            return params
        except RuntimeError as e:
            print(f"Error fitting Lorentzian: {e}")
            return None
        
    def lorentzian(freq_used, field_used, cutData, filepath, plotting=True):
        """
        Perform Lorentzian fitting on the given data.

        Parameters:
        - freq_used: List of frequencies
        - field_used: Magnetic field array
        - cutData: Data to be fitted
        - filepath: Filepath to save the results
        - plotting: Flag to enable/disable plotting

        Returns:
        None
        """
        rows = []
        results_df = pd.DataFrame(columns=['Frequency (GHz)', 'Amplitude 1', 'H_fmr 1', 'DeltaH 1',
                                        'Amplitude 2', 'H_fmr 2', 'DeltaH 2'])
        for fr in freq_used:
            x = field_used
            y = Plot.Fdata(cutData, F=fr, freq=freq_used)
            x01 = x[np.argsort(y)[-1]]
            x02 = -x01 + 20  # + 20 is added as a small correction for the location of the second peak
            
            # find the two crossing points
            hmx = FWHM.half_max_x(x, y)[0]
            index = FWHM.half_max_x(x, y)[1]
            
            # Handle the IndexError and break the loop
            try:
                fwhm = [abs(hmx[i] - hmx[i + 1]) for i in range(0, len(hmx), 2)]
            except IndexError:
                print("Error: list index out of range. Skipping this iteration.")
                break
            
            peak1 = max(y[:round(len(y)/2)])
            peak2 = max(y[round(len(y)/2):])
            initial_guess = [peak1, x01, fwhm[0], peak2, x02, fwhm[-1], min(y)]
            fwhm.clear()
            
            params = Fitting.fit_lorentzian(x, y, initial_guess=initial_guess)

            # Extract the fitted parameters
            A1_fit, x01_fit, G1_fit, A2_fit, x02_fit, G2_fit, offset_fit = params

            if plotting:
                plt.figure(figsize=(6, 4))
                plt.plot(x, y, 'ko', markersize=5, label='Experimental')
                plt.plot(x, Fitting.lorentzian_function(x, A1_fit, x01_fit, G1_fit, A2_fit, x02_fit, G2_fit, offset_fit), 'r-', label='Fitted Curve')
                plt.xlabel('Field (Oe)')
                plt.title(f'{fr} GHz')
                plt.legend()
                plt.show()

            if x01_fit > 0:
                temp = x01_fit
                x01_fit = x02_fit
                x02_fit = temp
                
            # Append the fitted parameters to the DataFrame
            current_row = {'Frequency (GHz)': round(fr, 2),
                        'Amplitude 1': A1_fit,
                        'H_fmr 1': x01_fit,
                        'DeltaH 1': G1_fit,
                        'Amplitude 2': A2_fit,
                        'H_fmr 2': x02_fit,
                        'DeltaH 2': G2_fit}
            rows.append(current_row)
        results_df = pd.DataFrame(rows)
        results_df.to_csv(filepath, index=False, header=True)
        
    def fit_plot(freq, frequency, field, raw_data, fitted_data):
        """
        Plot the raw and fitted data for a specific frequency.

        Parameters:
        - freq: Frequency for which data is to be plotted
        - frequency: Array of frequencies
        - field: Magnetic field array
        - raw_data: Raw data array
        - fitted_data: Fitted data array

        Returns:
        None
        """
        index = np.where(frequency == freq)[0]
        if index.size > 0:
            index = index[0]
            field_row = field[index, :]
            raw_data_row = raw_data[index, :]
            fitted_data_row = fitted_data[index, :]

            # Plotting raw and fitted data
            plt.plot(field_row, raw_data_row, 'ko')
            plt.plot(field_row, fitted_data_row, linewidth=2, label=f'{freq} GHz')
            plt.xlabel('Field (Oe)')
            plt.ylabel('$dP/dH (a.u.)$')
            plt.legend()