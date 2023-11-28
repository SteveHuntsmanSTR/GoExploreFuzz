# GoExploreFuzz
Copyright 2023 STR, except for get_Ahp.m, which is adapted from https://github.com/zboyd2/hitting_probabilities_metric and available under the MIT license at https://github.com/zboyd2/hitting_probabilities_metric/blob/master/LICENSE, herein included as hpLicense.txt).

The creatively named goExploreFuzzScript.m MATLAB script takes weeks to run as-is to reproduce results from the paper, but is also useful as an example. It calls various MATLAB functions that we list in the same order as the TeX comments or appendices (either one or the other is present) in whichever preprint version is at hand:

*  goExploreFuzzVersion6.m is the main fuzzer function
*  dynamicExecution.m dynamically executes a program specified by an underlying control flow graph produced by fuzzableCFG.m
*  fuzzableCFG.m produces an annotated control flow graph and associated arrays suitable for dynamic execution by dynamicExecution.m
*  functionalEditDistance.m implements a weighted edit distance that accepts function handles for costs
*  generateInitialFuzz.m generates initial inputs with distinct corresponding traces/paths (differs from generateLandmarkFuzz.m in that only distinctness is sought, not maximal diversity)
*  generateLandmarksFuzz.m generates diverse landmark inputs as gauged by traces/paths (differs from generateInitialFuzz.m in that diversity is sought, not merely distinctness)
*  generateLandmarksFuzzVerbose.m is a variant of the above used for benchmarking
*  get_Ahp.m is a variant of the version at https://github.com/zboyd2/hitting_probabilities_metric; see above
*  partialCouponCollection.m calculates the expected time for the general coupon collection problem
*  stateCell.m maps states to cells
*  strongCutoff.m produces the "strong cutoff" of a symmetric dissimilarity space
*  temp20230508_update is the example update function used in the script
*  temp20230511 is the example local generator used in the script
