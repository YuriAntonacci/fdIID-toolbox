# fdIID-toolbox
This Matlab toolbox provides a framework for frequency-specic quantification of the information shared between a target and two source time series (even multivariate). 
The interaction information decomposition is provided in the frequency domain and related to the time domain through the block-coherence
estimator. The decomposition can be computed directly from the estimation
of Power Spectral Density matrix through parametric and non-parametric approaches.
The framework is illustrated in simulations of linearly interacting stochastic
processes, showing how it allows to retrieve amounts of information shared by
the processes within specifc frequency bands which are otherwise not detectable
by time-domain information measures, as well as coupling features which are
not detectable by spectral measures [1]. The performance of different parametric and non-parametric 
estimators are compared in simulative setting with a
four-variate VAR process of interacting processes with a pre defined coupling
structure [2]. Its application is also illustrated on three different datasets:

1)  Electroencephalographic signals recorded from 109 healthy subjects performing a motor execution task [3]

2) Time series recorded from a unidirectionally-coupled ring of 32 electronic non-linear chaotic oscillator [4]

3) Three time series representative of the dynamic activity of three environmental variables: global land temperature, global temperature of the
ocean and carbon dioxide volume [5]

[1]-Faes, Luca, et al. "Information Decomposition in the Frequency Domain:
a New Framework to Study Cardiovascular and Cardiorespiratory Oscillations."
bioRxiv (2020).

[2]-Antonacci, Yuri, et al. "Measuring High-Order Interactions in Rhythmic
Processes: A Framework for the Spectral Information Decomposition of Multi-
variate Time Series.", IEEE Access (2021, sub)

[3]- https://physionet.org/content/eegmmidb/1.0.0/

[4]- Minati, Ludovico, et al. "Apparent remote synchronization of amplitudes:
A demodulation and interference effect." Chaos: An Interdisciplinary Journal
of Nonlinear Science 28.6 (2018): 063124.
[5]- E. system research laboratory, \Earth system research laboratory," https:
//www.esrl.noaa.gov/. NASA, \Nasa - goddard institute for space studies,"
https://data.giss.nasa.gov/

The code is provided free of charge. It is neither exhaustively tested nor
particularly well documented. The authors accept no liability for its use. Use,
modification and redistribution of the code is allowed in any way users see fit.
Authors ask only that authorship is acknowledged and ref. [1]-[2] is cited upon
utilization of the code in integral or partial form.To get started, we recommend
that you run and work through the demonstration scripts (README.pdf for more information)
