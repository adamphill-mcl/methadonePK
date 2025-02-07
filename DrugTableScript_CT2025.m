% This script reconstructs the drug-response curve of Methadone
% using the IC50 from Supplementary Table S2 of Llopis-Lorente, J.,
% Gomis-Tena, J., Cano, J., Romero, L., Saiz, J., & Trenor, B. (2020)
% using the method described in Heitmann, Vandenbberg, Hill (2023).
% The Hill coefficient are assumed to be h=1 in all cases.

clear

% EFTPC, IC50 and Hill coefficient for Methadone from Llopis-Lorente.
EFTPC    = 507;        % Effective Free Therapeutic Dose
IC50IKr  = 7000;       % IC50 for IKr
hIKr     = 0.82;       % Hill coefficent for IKr (not used)
IC50INaL  = 31800;      % IC50 for INa (not used) %AH edited
hINaL     = 1.37;       % Hill coefficient for INa (not used) %AH edited
% IC50INa  = 31800;      % IC50 for INa (not used) 
% hINa     = 1.37;       % Hill coefficient for INa (not used)
IC50ICaL = 37400;      % IC50 for ICaL
hICaL    = 1.67;       % Hill coefficient for ICaL (not use)

% construct the output table
DrugTable = table();

% for each Cmax ....
for Cmax = 0:0.1:100
    % Concentration of the drug at the current Cmax 
    Conc = Cmax * EFTPC;

    % Calculate the scaling factor for conductance, assuming a Hill
    % coefficient of h=1. See Eqn (7) in Heitmann et al (2023).
    GKrScale  = 1./(1 + (Conc./IC50IKr)); 
    GCaLScale = 1./(1 + (Conc./IC50ICaL)); 
    GNaLScale = 1./(1 + (Conc./IC50INaL));
   
    % collate the output into a single data structure
    DS = [];
        DS.Compound  = "Methadone";         % Name of Compound
        DS.Class     = 1;                   % Methadone has Class 1 clinical risk label (Credible Meds)
        DS.EFTPC     = EFTPC;               % Effective Free Therapeutic Dose
        DS.Cmax      = Cmax;                % Dose relative to EFTPC
        DS.Conc      = Conc;                % Concentration of drug in nM (based on)
        DS.GKrScale  = GKrScale;            % Scaling factor for IKr
        DS.GCaLScale = GCaLScale;           % Scaling factor for GCaL
        DS.GNaLScale = GNaLScale;           % Scaling factor for GNaL
        DS.LogGKrScale  = log(GKrScale);    % Natural logarithm
        DS.LogGCaLScale = log(GCaLScale);   % Natural logarithm
        DS.LogGNaLScale = log(GNaLScale);   % Natural logarithm
        
        % Append the data to the output table as a new row
        DrugTable = [DrugTable ; struct2table(DS)];  
    end    

% Compute the drug risk scores for Methadone using the two-channel
% axis-of-arrhythmia by Heitmann et al (unpublished).
BCaL =  0.81316;    % Coefficient of the axis of arrhythmia in LogGCaLScale
BKr  = -0.58204;    % Coefficient of the axis of arrhythmia in LogGKrScale
BNaL = 0.377;       % Coefficient of the axis of arrhythmia in LogGNaLScale
DrugTable.RiskScore = BCaL * DrugTable.LogGCaLScale + BKr * DrugTable.LogGKrScale + BNaL * DrugTable.LogGNaLScale;

save('DrugTable.mat','DrugTable');
