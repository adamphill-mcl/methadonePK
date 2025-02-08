function [outputArg1,outputArg2] = ModelMain(formulation,CypScore)
    %% Combined script for racaemic methadone. We can run R and S individually, 
    % but need to combine before calculating the combined risk score
    
    % Ensure the Brain Dynamics Toolbox is in the matlab PATH.
    addpath ../../../bdtoolbox-2023a/
    
    % Load the Methadone dose-risk response curve
    load('DrugTable.mat','DrugTable');
    
    % Load the relevant dose schedule (DoseTable) from file
    load('DoseTable','DoseTable'); %change this every time

    load('DrugPars.mat');
    
    %IC50s - update from Cliffs data
    IC50IKr_meth  = 7e-6;       % IC50 for IKr
    hIKr_meth     = 0.82;       % Hill coefficent for IKr
    IC50IKr_metab  = 7e-6;       % IC50 for IKr
    hIKr_metab     = 0.82;       % Hill coefficent for IKr
    IC50INaL_meth  = 31.8e-6;      % IC50 for INa 
    hINaL_meth     = 1.37;       % Hill coefficient for INa
    IC50INaL_metab  = 31.8e-6;      % IC50 for INa 
    hINaL_metab     = 1.37;       % Hill coefficient for INa
    IC50ICaL_meth = 37.4e-6;      % IC50 for ICaL
    hICaL_meth    = 1.67;       % Hill coefficient for ICaL
    IC50ICaL_metab = 37.4e-6;      % IC50 for ICaL
    hICaL_metab    = 1.67;       % Hill coefficient for ICaL
    
    % Run the Aruldhas model using the given dose schedule

    switch formulation
        case 'R'
            mkdir(strcat(formulation, '_', num2str(CypScore)))
            flag = 'R';     % Aruldhas flag must be 'R' or 'S' methadone
            dt = 0.1;       % Simulation time step (hours)
            RunTable = RunModel(DoseTable,flag,dt);
    
            % Time-dependent concentration of the drug in the central compartment
            tsim = RunTable.t;       % time points of the simulation
            %Conc = RunTable.A3;      % central compartment. WARNING: this probably needs rescaling for different units
            RunTable.Conc_meth = log10(RunTable.A2/1000/309.445); %convert to molar. mw for metadone. log units
            RunTable.Conc_metab = log10(RunTable.A3/1000/309.445); %convert to molar. confirm mw for metabolite. log units
        
            %calculate the combined block for the parent and metbolite for
            %each channel
            
            RunTable.IKrBlock = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_metab')) - RunTable.Conc_metab) .* DrugPars.h('IKr_metab'))));
            RunTable.ICaLBlock = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('ICaL_methR'))));
            RunTable.INaLBlock = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('INaL_methR'))));

            %Calculate the risk score for parent and metabolite combined using the two-channel
% axis-of-arrhythmia by Heitmann et al (unpublished).
           
            BCaL =  0.81316;    % Coefficient of the axis of arrhythmia in LogGCaLScale
            BKr  = -0.58204;    % Coefficient of the axis of arrhythmia in LogGKrScale
            BNaL = 0.377;       % Coefficient of the axis of arrhythmia in LogGNaLScale
            RunTable.RiskScore = BCaL .* log(1-RunTable.ICaLBlock) + BKr .* log(1-RunTable.IKrBlock) + BNaL .* log(1-RunTable.INaLBlock);
            
%             % Interpolate the RiskScore for Methadone from the dose-risk curve in DrugTable
%             RunTable.RiskScore = interp1(DrugTable.Conc, DrugTable.RiskScore, Conc);
            
            % Plot the simulation results
            figure(1)
            stackedplot(RunTable, ["A1'", "A2", "A3", "A4", "RiskScore"]);
            figure(2)
            plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
            
%             save(strcat('./', name, '/', sheets(s), '.mat'),'mdl','tbl', 'violin_stats_summary', 'em_means')
            save(strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore)),'RunTable','DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'overlaid','.fig'))
        
        case 'S'
            mkdir(strcat(formulation, '_', num2str(CypScore)))
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
            figure(1)
            stackedplot(RunTable);
            figure(2)
            plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
            
%             save(strcat('./', name, '/', sheets(s), '.mat'),'mdl','tbl', 'violin_stats_summary', 'em_means')
            save(strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore)),'RunTable','DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'overlaid','.fig'))
        
        case 'racaemic'
            mkdir(strcat(formulation, '_', num2str(CypScore)))
            flag = 'R';     % Aruldhas flag must be 'R' or 'S' methadone
            dt = 0.1;       % Simulation time step (hours)
            RunTable = RunModel(DoseTable,flag,dt);
    
            % Time-dependent concentration of the drug in the central compartment
            tsim = RunTable.t;       % time points of the simulation
            %Conc = RunTable.A3;      % central compartment. WARNING: this probably needs rescaling for different units
            Conc = RunTable.A3/1000/309.445*1000000000; %convert to nm. mw for metadone %AH edit
            
            % Interpolate the RiskScore for Methadone from the dose-risk curve in DrugTable
            RunTable.RiskScore = interp1(DrugTable.Conc, DrugTable.RiskScore, Conc);
            
            % Plot the simulation results
            figure(1)
            stackedplot(RunTable);
            figure(2)
            plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
            
%             save(strcat('./', name, '/', sheets(s), '.mat'),'mdl','tbl', 'violin_stats_summary', 'em_means')
            save(strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore)),'RunTable','DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'overlaid','.fig'))
        









        otherwise
            error('Input formulation must be ''R'', ''S'', or "racaemic" only.');
    end
    
%     flag = 'S';     % Aruldhas flag must be 'R' or 'S' methadone
%     dt = 0.1;       % Simulation time step (hours)
%     RunTable = RunModel(DoseTable,flag,dt);
%     
%     % Time-dependent concentration of the drug in the central compartment
%     tsim = RunTable.t;       % time points of the simulation
%     %Conc = RunTable.A3;      % central compartment. WARNING: this probably needs rescaling for different units
%     Conc = RunTable.A3/1000/309.445*1000000000; %convert to nm. mw for metadone %AH edit
%     
%     % Interpolate the RiskScore for Methadone from the dose-risk curve in DrugTable
%     RunTable.RiskScore = interp1(DrugTable.Conc, DrugTable.RiskScore, Conc);
%     
%     % Plot the simulation results
%     stackedplot(RunTable);
%     plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
%     
%     save('RunModel_S-meth_test','RunTable','DoseTable') %change this every time
%     outputArg1 = inputArg1;
%     outputArg2 = inputArg2;
end