DIVAIntAnalysis
===============

Scripts for analyzing the results from the DIVA Intensity Shift Experiment (CadLab)

Main program: intShift_analysis_1.m 

You will need to modify lines 15 and 24 for 
  1) where the two data folders (DN and UP) are and 
  2) where the individual-subject plots will be saved

Usage: intShift_analysis_1(bCD, nStresses, reverseOpt, otherOpts)
Inputs: 
  bCD - a 0/1 variable that indicates whether contrast distance (1) or asbolute values (0) are shown. 
  nStresse - vector for specifying the stress positions to include in the analysis.
            e.g., [1, 2] - both sentences with stress on word 1 and sentences with stress on word 2
                  [1] - only sentences with stress on word 1
                  
  reverseOpt - option to reverse stressed and unstressed, for looking at the absolute values from the unstressed words. 
          NOTE: this option only works under bCD = 0. Do not use this option with bCD = 1.
          For example, if you want to look at the absolute values from the unstressed words (words 1 or 2) of all sentences 
          Do: 
            intShift_analysis_1(0, [1, 2], 'reverse')
          Or, if you want to look at the absolute values from the unstressed words on position 1, do:
            intShift_analysis_1(0, [2], 'reverse')
            Note that we put [2], instead of [1] here, because when word 2 is stressed, word 1 is unstressed. 
            
  otherOpts - additional options, such as "showByEpoch" and "showIndS" (see examples below).
            
Usage exapmles:
  1. Contrast distance, from stressed words as both position 1 and 2: intShift_analysis_1(1, [1, 2])
  2. Same as above, but show epoch-by-epoch data: intShift_analysis_1(1, [1, 2], 'showByEpoch')
  3. Same as above, but show data from individual subjects: intShift_analysis_1(1, [1, 2], 'showIndS')
  4. Absolute values, from stressed words at positions 1 and 2: intShift_analysis_1(0, [1, 2])
  5. Absolute values, from only stressed words at position 2: intShift_analysis_1(0, [2])
  6. Absolute values, from unstressed words at positions 1 and 2: intShift_analysis_1(0, [1, 2], 'reverse')
  7. Show composite prosody adaptation (CPA) scores, by phase: intShift_analysis_1(1, [1, 2], 'cpa')
  8. Show CPA scores, by phase and by epoch: intShift_analysis_1(1, [1, 2], 'cpa', 'showByEpoch')

  9. Show CPA scores and perfor permutation test (10000 times): intShift_analysis_1(1, [1, 2], 'cpa', 'permute', 10000)
  
  10. Analyze the correlation between the adaptation measures and the pertContr:
	intShift_analysis_1(1, [1, 2], 'cpa', 'corrPert')
	
	Note: you may need to revise the "L3_DATA_DIR" variable to set the correct path to the level-3 (in Kevin's parlance) csv files.
	
  11. Use the rank-sum and signed rank tests, instead of the default t-tests; perform random permutation for 10,000 iterations, and do not show the uncorrected p-values in the figures:
         intShift_analysis_1(1, [1, 2], 'cpa', 'loadCache', 'rs', 'permute', 10000, 'noUnc')
