% This script creates a DoseTable for the AH Loewe workflow and saves it to
% DoseTable.mat using the default baseline dose and overdose multiplier.

clear

formulation = 'R';
dose_original = 120;
overdose_multiplier = 1;

DoseTable = BuildDoseTable_AH(formulation, dose_original, overdose_multiplier);
save('DoseTable.mat', 'DoseTable');
