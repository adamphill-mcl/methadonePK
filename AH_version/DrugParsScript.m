% This script creates DrugPars for the AH Loewe workflow and saves it to
% DrugPars.mat. Set hill_mode to 'ideal' or 'real' before running.

clear

hill_mode = 'ideal';

DrugPars = BuildDrugPars_AH(hill_mode);
save('DrugPars.mat', 'DrugPars');
