
% Ensure the Brain Dynamics Toolbox is in the matlab PATH.
addpath ../../../bdtoolbox-2023a/

% Load the Methadone dose-risk response curve
load('DrugTable.mat','DrugTable');

% Load the relevant dose schedule (DoseTable) from file
load('DoseTable','DoseTable'); %change this every time

% Run the Aruldhas model using the given dose schedule
flag = 'S';     % Aruldhas flag must be 'R' or 'S' methadone
dt = 0.1;       % Simulation time step (hours)
RunTable = RunModel(DoseTable,flag,dt);

% Time-dependent concentration of the drug in the central compartment
tsim = RunTable.t;       % time points of the simulation
%Conc = RunTable.A3;      % central compartment. WARNING: this probably needs rescaling for different units
Conc = RunTable.A3/1000/309.445*1000000000; %convert to nm. mw for metadone %AH edit

% Interpolate the RiskScore for Methadone from the dose-risk curve in DrugTable
RunTable.RiskScore = interp1(DrugTable.Conc, DrugTable.RiskScore, Conc);

% Plot the simulation results
stackedplot(RunTable);
plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;

save('RunModel_S-meth_test','RunTable','DoseTable') %change this every time