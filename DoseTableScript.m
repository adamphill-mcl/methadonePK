% This script creates a DoseTable which contains the times and doses 
% of the dosing reginme to be appliedd. The results are saved in DoseTable.mat.

clear

% Define the regular dose
dose_original = 4.6632;                     % Avg dose of 0.087 mg/kg multiplid by average weight of 53.6kg
F1_bioavailability = 0.718;                 % Percentage of dose to pass FPM
dose = dose_original*F1_bioavailability;    % Concentration of the drug to be added

% Choose a Methadone flag (comment one out)
%flag = 'R';     % for R-Methadone
flag = 'S';     % for S-Methadone

% Define the dosing regime here (edit as required)
DoseArray = [
    struct( 'Conc',0,     'dur',  24,  'MethadoneFlag',flag)  % nil dose, followed by 24 hour interval.
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)  % 1st dose, followed by 24 hour interval.
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)  % 2nd dose, followed by 24 hour interval.
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur',  24,  'MethadoneFlag',flag)
    struct( 'Conc',dose,  'dur', 24,  'MethadoneFlag',flag)  % Final dose, followed by 24 hour interval.
    ];

% Convert the array of structs into a table
DoseTable = struct2table(DoseArray);
save('DoseTable.mat','DoseTable');
