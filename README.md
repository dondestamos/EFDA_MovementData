# EFDA_MovementData
A lightweight mex-precompiled Matlab libraries demonstrating and extending two "fdasrvf" (https://github.com/jdtuck/fdasrvf_MATLAB) functions for elastic alignment of biological movement signals. Ready to run on most Windows 10/11 machines.
Submitted with manuscript Krotov, S. Razavian, Sadeghi, and Sternad
"Time-Warping to Extract Salient Features and Decouple Spatial and Temporal Variability in Biological Signals – A Tutorial Review with Application to Human Movement Data"

![Figure 1 - Problem, TW, Variabilities](https://github.com/user-attachments/assets/a8cfe9a6-15b8-482e-99b9-cc123c5e8b25)

![Figure 3 - results synthetic estimates demo](https://github.com/user-attachments/assets/98b31501-1f7c-47af-8a6b-ff7f9d08710d)

![Figure 4 - Synthetic Variabilities Horiz](https://github.com/user-attachments/assets/d3773c36-f3a4-4c68-b45d-612e8213ce74)

![Figure 5 - whip task methods](https://github.com/user-attachments/assets/64e03dd5-5c30-4e67-b998-bdafa5abce52)


Features:
- demonstrations of several "elementary" warping functions and their effects on an exemplary signal (see Warp shape examples)
- demonstration of EFDA alignment three-way using either of two signals as a template or to their estimated Kärcher mean (EFDA_conceptDemos)
- generating synthetic 2-gaussian-peak signals with random noise injected in various (configurable) ways; their following alignment and averaging aiming to estimate signal parameters from that mean. Via Time-padding, time-normalization, and EFDA time-warping; extraction of variabilities.
- alignment via the same 3 options (plus one, time-padding with rightward alignment) of the experimental data of hand speed during whip manipulation task (see Krotov & Russo 2022, Royal Society Open Science). Extraction of variabilities.
- EFDA_alignmentMain.m as a wrapper of fdawarp.time_warping with additional features:
   - quick 1-second alignment on downsampled data using a subset of original signals
   - using derivative or integral for alignment (with following integration or derivation) to keep them invariant - if needed for the analysis.
   - estimating execution time - becomes relevant for large or long datasets
   - optional filtering
   - elastic distances and a few variability estimates in the result object
   - result object fields named in a user-friendlier way

See the functions from folders (1.) and (2.) for more description.

Dependencies are included here from the full library "fdasrvf" (https://github.com/jdtuck/fdasrvf_MATLAB), with edits only in the fdawarp.time_warping class/method adding history of iterations and additional variability definitions to the output structure. 

Mex-files were precompiled on 64-bit Windows 10 and are expected to run straight ahead on most of the modern Windows platforms. For other platforms, use the original library and implement changes from here manually, if needed.


