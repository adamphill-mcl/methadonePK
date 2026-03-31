% This script creates a DoseTable which contains the times and doses 
% of the dosing reginme to be appliedd. The results are saved in DoseTable.mat.

clear

% Define the regular dose
dose_original = 180;                         % Avg dose of 0.087 mg/kg multiplid by average weight of 53.6kg = 4.6632mg
F1_R_bioavailability = 0.718;                 % Percentage of dose to pass FPM, 0.718 for R-methadone, 0.606 for S-methadone
F1_S_bioavailability = 0.606;                 % Percentage of dose to pass FPM, 0.718 for R-methadone, 0.606 for S-methadone
R_dose = (dose_original / 2)*F1_R_bioavailability;    % Concentration of the drug to be added
S_dose = (dose_original / 2)*F1_S_bioavailability;    % Concentration of the drug to be added
FD = 4;                                          %Final dose, select either 1x for normal dose, or 2 or 4 for OD/overdose

% Choose a Methadone flag (comment one out)
%flag = 'R';     % for R-Methadone
%flag = 'S';     % for S-Methadone

% Define the dosing regime here (edit as required)
DoseArray1 = [
    struct( 'Conc',0,            'dur',  24,  'MethadoneFlag','R')  % nil dose, followed by 24 hour interval.
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')  % 1st dose, followed by 24 hour interval.
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')  % 2nd dose, followed by 24 hour interval.
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose,       'dur',  24,  'MethadoneFlag','R')
    struct( 'Conc',R_dose*FD,       'dur',  24,  'MethadoneFlag','R')  % Final dose, followed by 24 hour interval.
    ];

% Convert the array of structs into a table
DoseTable_R = struct2table(DoseArray1);
save('DoseTable_R.mat','DoseTable_R');

DoseArray2 = [
    struct( 'Conc',0,        'dur',  24,  'MethadoneFlag','S')  % nil dose, followed by 24 hour interval.
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')  % 1st dose, followed by 24 hour interval.
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')  % 2nd dose, followed by 24 hour interval.
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose,   'dur',  24,  'MethadoneFlag','S')
    struct( 'Conc',S_dose*FD,   'dur',  24,  'MethadoneFlag','S')  % Final dose, followed by 24 hour interval.
    ];

% Convert the array of structs into a table
DoseTable_S = struct2table(DoseArray2);
save('DoseTable_S.mat','DoseTable_S');
