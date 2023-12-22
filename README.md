# VNA Data Acquisition Notebook

This Jupyter Notebook is designed for the VNA FMR experiment, providing tools for efficient data acquisition and analysis. The notebook is tailored to work with the Rohde & Schwarz VNA ZVA 40 and KEPCO Bipolar Operational Power Supply (BOP) in the Spin Physics Lab @ LUMS, a project of PhysLab.

## Prerequisites

Before running the notebook, ensure that you have the required Python libraries installed. You can install them using the following command:

```bash
pip install numpy matplotlib pyvisa lmfit scipy
```

## Usage

1. **Initialize VNA and KEPCO**: Configure the VNA and KEPCO instruments by setting up the experiment according to your requirements. Ensure that the VISA resource manager is initialized.

```python
import pyvisa
rm = pyvisa.ResourceManager()
# Initialize other instruments as needed
```

2. **1 Port Experiment**:
   - Customize the experiment parameters such as frequency range, fields, and sample name.
   - Run the code to measure and visualize S-parameters, power absorbed spectra, filtered data, and dP/dH.

```python
from VNAxKEPCO_1P import Experiment1P
E1 = Experiment1P()
# Customize parameters and run measurements
```

3. **2 Port Experiment**:
   - Similar to the 1 Port experiment, customize parameters for the 2 Port experiment.
   - Run the code to measure and visualize S-parameters, power absorbed spectra, filtered data, and dP/dH.

```python
from VNAxKEPCO_2P import Experiment2P
E2 = Experiment2P()
# Customize parameters and run measurements
```

4. **Data Analysis**:
   - Explore the provided functions for filtering, plotting, and analyzing the acquired data.
   - Customize the code to suit your specific analysis needs.


## Outputs

The notebook generates various plots, including raw and filtered S-parameter plots, power absorbed spectra, and dP/dH plots, depending on the specific analysis performed. Data is saved in designated directories.

## Credits

This notebook was developed for PhysLab @ LUMS.

- **Developer:** Mahad Naveed (BSc Physics @ LUMS 2023)
- **Supervisor:** Dr. Sabieh Anwar
- **Mentor:** Dr. Adnan Raza

Explore PhysLab: www.physlab.org

Feel free to adapt and customize the notebook for your experiment. For inquiries, contact Mahad Naveed at mahadn20@gmail.com.

Happy experimenting!
