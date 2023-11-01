# Get Phase Tools
code to extract the continuous instantaneous phase of experimental signals according to [1]

## content:

**get_phase** folder and its content: Matlab scripts implementing the methods proposed in [1].  
**gen_utils** folder: Matlab scripts required to reproduce the analyses reported in [1].  
**data** folder: files analysed in [1] and used to exemplify the method proposed there.  
**demo.mlx** file: demo file illustrating the functioning of the matlab functions.  
**paper_analyses.m** file: Matlab script to launch the analyses reported in [1].  
**getPhaseToolsM.py** file: Python functions implementing the methods proposed in [1].  
**requirements.txt** file: requirement file to install Python dependencies via PIP.  
**requirements.yml** file: requirement file to install Python dependencies via conda.  


# Installation Instructions:

## Matlab:

In order to use the proposed methods in you own workflow you need to add the folder **get_phase** and its subfolders to the Matlab path. This can be done by going to the folder containing the folder **get_phase** and run:

addpath(genpath('./get_phase'))

This is not required when launching the demo.mlx file from the folder containing **get_phase** and **data** or the script **paper_analyses.m** from the folder containing **get_phase**, **gen_utils** and **data**.
 

## Python:  
To install the modules required for the functioning og the methods proposed in [1]
go in the directory containitng the files **requirements.txt** and **requirements.yml**

type 

pip install -r requirements.txt

or

conda env update --name [envname] --file requirements.yml

where **[envname]** should be substituted by the name of the environment where you want to use Get Phase Tools.

# Usage: 
The scripts are accessible via two main functions having the same names and arguments in both Matlab and Python implementations: **getPhaseMask** and **mEMDdenoise**.

## getPhaseMask
This function computes the instantaneous phase of an input signal given its sampling frequency.
Its first argument is the input signal and its second argument the sampling frequency. The following optional arguments are:  
**m** (positive integer; optional, default: 16): number of filtered points for Savitzky-Golay differentiator.  
**n** (positive integer; optional, default: 5): polinomial order of the differentiator.  
**nMasks** (positive integer; optional, default 22): number of masks used to extract the independent mode function via masked EMD.  
**ampCoeff** (positive real; optional, default: 2): coefficient determining the amplitude of the masks as a proportion of 4*std(signal) (as this is a rough estimate of the signal's range from the std if the signal's values are normally distributed).  
**quadMethod** (a string or a cell of two strings in Matlab and a string or a list of two strings in Python, default: {'h','h'} in Matlab and ['h','h'] in Python): Method to be use in the computation of the quadrature signal 'h' stands for Hilbert and 'q'. If two strings are provided a different method will be adopted in in the first or the second part of the algorithm.  
**threshs** (scalar or vector of two positive real values close to zero, default: [1E-10,1E-10]  ): threshold for refined amplitude normalization. If two values, different thresholds will be used in the two parts of the algorithm.  

The function's outputs are:  
**PHI**: sthe signal's phase.  
**IMF1**: output of masked sifting.  
**PHI**: initial phase estimte.  
**centredSig**: centred signal.  
**mask**: masking signal with inital phase set to 0


## mEMDdeNoise
The function removes high frequency noise from an input signal given its sampling frequency.
Its first argument is the input signal and its second argument the sampling frequency. The other (optional) arguments are:  

**nMasks** (positive integer; optional, default: 22): number of mask signal for masked EMD application .  
**ampCoeff**(positive real; optional, default: 2*estimated range of first IMF): amplitude of the mask signals .  
**nReps**(positive integer; optional, default: 100): number of repetitions of simulated random processes .  
**m** (positive integer; optional, default: 16): length of the Savitzky-Golay differentiator.  
**n** (positive integer; optional, default: 5): order of the Savitzky-Golay differentiator.

the function's outputs are:  

**filtered**: filtered signal
**imf**: all signal's IMFs.  
**imfF**: frequencies of extracted IMFs.  
**filteredidx**: indexes of IMFs composing the filtered signal.  
**noiseStd**: standard deviation of the estimated noise (sum of the random components).
