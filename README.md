# Digital-filter-in-python-and-C
Program to design digital filters and anyze their magnitude characteristics after performing filtering on signal. Because of
pandemic I simulate work of DSP with program in C (processing of signal) and osciloscope, generator etc with program in Python.

Filter Design folder:
FilterDesign.py - its purpose is compute coefficients of desired filter it can also display magnitude, phase characteristics
and zeros and poles of filter

filter_bench.py - its program which merges FilterDesign and filtering program in C. In filter_bench sinusoidal signal is generated
in range of frequencies of user choice and it displays magnitude characteristic, filtering can be performed for various types
of arithmetic (floating point, fixed point arithmetic in Qn.m code) and user can analyse characteristics diffrences.

levinson is part of spectrum library for some reason when i was writing that program it must be in folder where was FilterDesign.py

Filtering folder:
filtr.c - there are implemented structures of filters
for FIR - direct in floating point and fixed point
for IIR direct, biquad and lattice in floating and fied point. Some of functions are exported to .dll file and used in python part
of project.

Tests folder contains images of few filter designs and diferences in characteristics when various structures are used.



