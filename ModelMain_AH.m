function [outputArg1,outputArg2] = ModelMain_AH(formulation,CypScore)
    %% Combined script for racaemic methadone. We can run R and S individually, 
    % but need to combine before calculating the combined risk score
    
    % Ensure the Brain Dynamics Toolbox is in the matlab PATH.
    addpath ../../../bdtoolbox-2023a/
    
    % Load the Methadone dose-risk response curve
    load('DrugTable.mat','DrugTable');
    
    % Load the relevant dose schedule (DoseTable) from file
    load('DoseTable','DoseTable'); %change this every time

     % Load the dose response parameters from Cliffs data from file
     load('DrugPars.mat');
    
    V_R = 176; %apparent volumes of central and peripheral from model
    V_S = 98.3;
    % Run the Aruldhas model using the given dose schedule

    switch formulation
        case 'R'
            mkdir(strcat(formulation, '_', num2str(CypScore)))
            flag = 'R';     % Aruldhas flag must be 'R' or 'S' methadone
            dt = 0.1;       % Simulation time step (hours)
            RunTable = RunModel_AH(DoseTable,flag,dt,CypScore);
    
            % Time-dependent concentration of the drug in the central compartment
            tsim = RunTable.t;       % time points of the simulation
%             RunTable.Conc_meth = log10(RunTable.A2/1000/309.445); %convert to molar. mw for metadone. log units
%             RunTable.Conc_metab = log10(RunTable.A3/1000/309.445); %convert to molar. confirm mw for metabolite. log units
            RunTable.Conc_meth = log10(RunTable.A2/1000/V_R/309.445); %above doesnt take into account the apparent volume of compartment
            RunTable.Conc_metab = log10(RunTable.A3/1000/V_R/309.445); %convert to molar. confirm mw for metabolite. log units
        
            %calculate the combined block for the parent and metbolite for
            %each channel
            RunTable.IKrBlock = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_metab')) - RunTable.Conc_metab) .* DrugPars.h('IKr_metab'))));
            RunTable.ICaLBlock = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('ICaL_metab'))));
            RunTable.INaLBlock = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('INaL_metab'))));

            %Calculate the risk score for parent and metabolite combined using the two-channel
            % axis-of-arrhythmia by Heitmann et al (unpublished).
           
            BCaL =  0.81316;    % Coefficient of the axis of arrhythmia in LogGCaLScale
            BKr  = -0.58204;    % Coefficient of the axis of arrhythmia in LogGKrScale
            BNaL = 0.377;       % Coefficient of the axis of arrhythmia in LogGNaLScale
            RunTable.RiskScore = BCaL .* log(1-RunTable.ICaLBlock) + BKr .* log(1-RunTable.IKrBlock) + BNaL .* log(1-RunTable.INaLBlock);
            
            
            % Plot the simulation results
            figure(1)
%             stackedplot(RunTable, ["A1'", "A2", "A3", "A4", "RiskScore"]);
            stackedplot(RunTable);
            figure(2)
            plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
            
%             save the figures and data
            save(strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore)),'RunTable','DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'overlaid','.fig'))
        
        case 'S'
            mkdir(strcat(formulation, '_', num2str(CypScore)))
            flag = 'S';     % Aruldhas flag must be 'R' or 'S' methadone
            dt = 0.1;       % Simulation time step (hours)
            RunTable = RunModel_AH(DoseTable,flag,dt);
    
            % Time-dependent concentration of the drug in the central compartment
            tsim = RunTable.t;       % time points of the simulation
            RunTable.Conc_meth = log10(RunTable.A2/1000/V_S/309.445); %convert to molar. mw for metadone. log units
            RunTable.Conc_metab = log10(RunTable.A3/1000/V_S/309.445); %convert to molar. confirm mw for metabolite. log units
        
            %calculate the combined block for the parent and metbolite for
            %each channel
            RunTable.IKrBlock = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methS')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methS')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methS')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methS'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_metab')) - RunTable.Conc_metab) .* DrugPars.h('IKr_metab'))));
            RunTable.ICaLBlock = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methS')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methS')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methS')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methS'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('ICaL_metab'))));
            RunTable.INaLBlock = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methS')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methS')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methS')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methS'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('INaL_metab'))));

            %Calculate the risk score for parent and metabolite combined using the two-channel
            % axis-of-arrhythmia by Heitmann et al (unpublished).
           
            BCaL =  0.81316;    % Coefficient of the axis of arrhythmia in LogGCaLScale
            BKr  = -0.58204;    % Coefficient of the axis of arrhythmia in LogGKrScale
            BNaL = 0.377;       % Coefficient of the axis of arrhythmia in LogGNaLScale
            RunTable.RiskScore = BCaL .* log(1-RunTable.ICaLBlock) + BKr .* log(1-RunTable.IKrBlock) + BNaL .* log(1-RunTable.INaLBlock);
            
            
            % Plot the simulation results
            figure(1)
%             stackedplot(RunTable, ["A1'", "A2", "A3", "A4", "RiskScore"]);
            stackedplot(RunTable);
            figure(2)
            plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
            
%             save the figures and data
            save(strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore)),'RunTable','DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'overlaid','.fig'))
        
        case 'racaemic'
            mkdir(strcat(formulation, '_', num2str(CypScore)))
            DoseTable.Conc = DoseTable.Conc./2; %we are going to run this twice, with half conc each of R and S enantiomers
            for ii = 1:2
                if ii == 1
                    flag = 'R';     % Aruldhas flag must be 'R' or 'S' methadone
                else
                    flag = 'S';

                    dt = 0.1;       % Simulation time step (hours)
                    RunTable = RunModel(DoseTable,flag,dt);
            
                    % Time-dependent concentration of the drug in the central compartment
                    
                    
                     %calculate the combined block for the parent and metbolite for
                    %each channel
                    if ii == 1 %not going o twork as need to acocunt for A cmopartments too - move runmodel ino tloop and write two tabel whihc are then compbined
                        RunTable.Conc_meth_R = log10(RunTable.A2/1000/309.445); %convert to molar. mw for metadone. log units
                        RunTable.Conc_metab_R = log10(RunTable.A3/1000/309.445); %convert to molar. confirm mw for metabolite. log units
                        RunTable.IKrBlock_R = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_metab')) - RunTable.Conc_metab) .* DrugPars.h('IKr_metab'))));
                        RunTable.ICaLBlock_R = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('ICaL_metab'))));
                        RunTable.INaLBlock_R = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methR')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('INaL_metab'))));
                    else
                        RunTable.Conc_meth_S = log10(RunTable.A2/1000/309.445); %convert to molar. mw for metadone. log units
                        RunTable.Conc_metab_S = log10(RunTable.A3/1000/309.445); %convert to molar. confirm mw for metabolite. log units
                        RunTable.IKrBlock_S = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methS')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methS')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methS')) - RunTable.Conc_meth) .* DrugPars.h('IKr_methS'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_metab')) - RunTable.Conc_metab) .* DrugPars.h('IKr_metab'))));
                        RunTable.ICaLBlock_S = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methS')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methS')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methS')) - RunTable.Conc_meth) .* DrugPars.h('ICaL_methS'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('ICaL_metab'))));
                        RunTable.INaLBlock_S = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methS')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methS')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methS')) - RunTable.Conc_meth) .* DrugPars.h('INaL_methS'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_metab')) - RunTable.Conc_metab) .* DrugPars.h('INaL_metab'))));
                    end
                        %Calculate the risk score for parent and metabolite combined using the two-channel
                        % axis-of-arrhythmia by Heitmann et al (unpublished).
                   
                    BCaL =  0.81316;    % Coefficient of the axis of arrhythmia in LogGCaLScale
                    BKr  = -0.58204;    % Coefficient of the axis of arrhythmia in LogGKrScale
                    BNaL = 0.377;       % Coefficient of the axis of arrhythmia in LogGNaLScale
                    RunTable.RiskScore = BCaL .* log(1-RunTable.ICaLBlock) + BKr .* log(1-RunTable.IKrBlock) + BNaL .* log(1-RunTable.INaLBlock);
                    
                    
                    % Plot the simulation results
                    figure(1)
        %             stackedplot(RunTable, ["A1'", "A2", "A3", "A4", "RiskScore"]);
                    stackedplot(RunTable);
                    figure(2)
                    plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
                    
        %             save the figures and data
                    save(strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore)),'RunTable','DoseTable') %change this every time
                    savefig(figure(1),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'stacked','.fig'))
                    savefig(figure(2),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'overlaid','.fig'))
                end
            end

    end

end