function [outputArg1,outputArg2] = ModelMain_Loewe(formulation,CypScore)
    %% Combined script for racaemic methadone. We can run R and S individually, 
    % but need to combine before calculating the combined risk score
    
    % Ensure the Brain Dynamics Toolbox is in the matlab PATH.
%     addpath ../bdtoolbox-2023a/
    
    % Load the Methadone dose-risk response curve
%     load('DrugTable.mat','DrugTable');
    
    % Load the relevant dose schedule (DoseTable) from file
    load('DoseTable','DoseTable'); %change this every time

     % Load the dose response parameters from Cliffs data from file
     load('DrugPars2.mat');
    
    V_R = 176; %apparent volumes of central and peripheral from model
    V_S = 98.3;
    % Run the Aruldhas model using the given dose schedule

    switch formulation
        case 'R'
            mkdir(strcat(formulation, '_', strrep(num2str(CypScore),'.','_')))
            flag = 'R';     % Aruldhas flag must be 'R' or 'S' methadone
            dt = 0.1;       % Simulation time step (hours)
            RunTable = RunModel_AH(DoseTable,flag,dt,CypScore);
    
            % Time-dependent concentration of the drug in the central compartment
            tsim = RunTable.t;       % time points of the simulation
            RunTable.Conc_meth = RunTable.A2/1000/V_R/309.445; %above doesnt take into account the apparent volume of compartment
            RunTable.Conc_metab = RunTable.A3/1000/V_R/309.445; %convert to molar. confirm mw for metabolite. log units
        
            %calculate the combined block for the parent and metbolite for
            %each channel
            Loewe_methR_IKr = (RunTable.Conc_meth/DrugPars.IC50s('IKr_methR')) .^ DrugPars.h('IKr_methR')
            Loewe_metab_IKr = (RunTable.Conc_meth/DrugPars.IC50s('IKr_metab')) .^ DrugPars.h('IKr_metab')
            RunTable.IKrBlock = (Loewe_methR_IKr + Loewe_metab_IKr) ./ (Loewe_methR_IKr + Loewe_metab_IKr +1);
            
            Loewe_methR_ICaL = (RunTable.Conc_meth/DrugPars.IC50s('ICaL_methR')) .^ DrugPars.h('ICaL_methR')
            Loewe_metab_ICaL = (RunTable.Conc_meth/DrugPars.IC50s('ICaL_metab')) .^ DrugPars.h('ICaL_metab')
            RunTable.ICaLBlock = (Loewe_methR_ICaL + Loewe_metab_ICaL) ./ (Loewe_methR_ICaL + Loewe_metab_ICaL +1);

            Loewe_methR_INaL = (RunTable.Conc_meth/DrugPars.IC50s('INaL_methR')) .^ DrugPars.h('INaL_methR')
            Loewe_metab_INaL = (RunTable.Conc_meth/DrugPars.IC50s('INaL_metab')) .^ DrugPars.h('INaL_metab')
            RunTable.INaLBlock = (Loewe_methR_INaL + Loewe_metab_INaL) ./ (Loewe_methR_INaL + Loewe_metab_INaL +1);

            %Calculate the risk score for parent and metabolite combined using the two-channel
            % axis-of-arrhythmia by Heitmann et al (unpublished).
           
            BCaL =  0.81316;    % Coefficient of the axis of arrhythmia in LogGCaLScale
            BKr  = -0.58204;    % Coefficient of the axis of arrhythmia in LogGKrScale
            BNaL = 0.377;       % Coefficient of the axis of arrhythmia in LogGNaLScale
            RunTable.RiskScore = BCaL .* log(1-RunTable.ICaLBlock) + BKr .* log(1-RunTable.IKrBlock) + BNaL .* log(1-RunTable.INaLBlock);
            
            
            % Plot the simulation results
            figure(1)
            %stackedplot(RunTable, ["A1'", "A2", "A3", "A4", "RiskScore"]);
            stackedplot(RunTable);
            figure(2)
            plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
            
%             save the figures and data
            save(strcat('./',formulation, '_', strrep(num2str(CypScore),'.','_'),'/',formulation, '_', strrep(num2str(CypScore),'.','_')),'RunTable','DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', strrep(num2str(CypScore),'.','_'),'/',formulation, '_', strrep(num2str(CypScore),'.','_'),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', strrep(num2str(CypScore),'.','_'),'/',formulation, '_', strrep(num2str(CypScore),'.','_'),'_', 'overlaid','.fig'))
        
        case 'S'
            mkdir(strcat(formulation, '_', strrep(num2str(CypScore),'.','_')))
            flag = 'S';     % Aruldhas flag must be 'R' or 'S' methadone
            dt = 0.1;       % Simulation time step (hours)
            RunTable = RunModel_AH(DoseTable,flag,dt,CypScore);
    
            % Time-dependent concentration of the drug in the central compartment
            tsim = RunTable.t;       % time points of the simulation
            RunTable.Conc_meth = RunTable.A2/1000/V_S/309.445; %convert to molar. mw for metadone. log units
            RunTable.Conc_metab = RunTable.A3/1000/V_S/309.445; %convert to molar. confirm mw for metabolite. log units
            
             %calculate the combined block for the parent and metbolite for
            %each channel
            Loewe_methS_IKr = (RunTable.Conc_meth/DrugPars.IC50s('IKr_methS')) .^ DrugPars.h('IKr_methS')
            Loewe_metab_IKr = (RunTable.Conc_meth/DrugPars.IC50s('IKr_metab')) .^ DrugPars.h('IKr_metab')
            RunTable.IKrBlock = (Loewe_methS_IKr + Loewe_metab_IKr) ./ (Loewe_methS_IKr + Loewe_metab_IKr +1);
            
            Loewe_methS_ICaL = (RunTable.Conc_meth/DrugPars.IC50s('ICaL_methS')) .^ DrugPars.h('ICaL_methS')
            Loewe_metab_ICaL = (RunTable.Conc_meth/DrugPars.IC50s('ICaL_metab')) .^ DrugPars.h('ICaL_metab')
            RunTable.ICaLBlock = (Loewe_methS_ICaL + Loewe_metab_ICaL) ./ (Loewe_methS_ICaL + Loewe_metab_ICaL +1);

            Loewe_methS_INaL = (RunTable.Conc_meth/DrugPars.IC50s('INaL_methS')) .^ DrugPars.h('INaL_methS')
            Loewe_metab_INaL = (RunTable.Conc_meth/DrugPars.IC50s('INaL_metab')) .^ DrugPars.h('INaL_metab')
            RunTable.INaLBlock = (Loewe_methS_INaL + Loewe_metab_INaL) ./ (Loewe_methS_INaL + Loewe_metab_INaL +1);

            %Calculate the risk score for parent and metabolite combined using the two-channel
            % axis-of-arrhythmia by Heitmann et al (unpublished).
           
            BCaL =  0.81316;    % Coefficient of the axis of arrhythmia in LogGCaLScale
            BKr  = -0.58204;    % Coefficient of the axis of arrhythmia in LogGKrScale
            BNaL = 0.377;       % Coefficient of the axis of arrhythmia in LogGNaLScale
            RunTable.RiskScore = BCaL .* log(1-RunTable.ICaLBlock) + BKr .* log(1-RunTable.IKrBlock) + BNaL .* log(1-RunTable.INaLBlock);
            
            
            % Plot the simulation results
            figure(1)
            %stackedplot(RunTable, ["A1'", "A2", "A3", "A4", "RiskScore"]);
            stackedplot(RunTable);
            figure(2)
            plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
            
%             save the figures and data
            save(strcat('./',formulation, '_', strrep(num2str(CypScore),'.','_'),'/',formulation, '_', strrep(num2str(CypScore),'.','_')),'RunTable','DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', strrep(num2str(CypScore),'.','_'),'/',formulation, '_', strrep(num2str(CypScore),'.','_'),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', strrep(num2str(CypScore),'.','_'),'/',formulation, '_', strrep(num2str(CypScore),'.','_'),'_', 'overlaid','.fig'))
        
        case 'racemic'
            mkdir(strcat(formulation, '_', strrep(num2str(CypScore),'.','_')))

            %we are going to run this twice, with half conc each of R and S enantiomers
            DoseTableRAC = DoseTable;
            DoseTableRAC.Conc = DoseTableRAC.Conc./2;
                                 
            % R-Methadone -----------
            flag = 'R';     % Aruldhas flag must be 'R' or 'S' methadone
            dt = 0.1;       % Simulation time step (hours)
            RunTableR = RunModel_AH(DoseTableRAC,flag,dt,CypScore);
    
            % Time-dependent concentration of the drug in the central compartment
            %tsim = RunTableR.t;       % time points of the simulation
            RunTableR.Conc_meth = RunTableR.A2/1000/V_R/309.445; %above doesnt take into account the apparent volume of compartment
            RunTableR.Conc_metab = RunTableR.A3/1000/V_R/309.445; %convert to molar. confirm mw for metabolite. log units
            
           %Store overall results in new table
            MyTable = table();
            MyTable.t = RunTableR.t;
            MyTable.conc_methR = RunTableR.Conc_meth;
            MyTable.conc_metabR = RunTableR.Conc_metab;
            
            
            % S-Methadone ------------
            flag = 'S';     % Aruldhas flag must be 'R' or 'S' methadone
            %dt = 0.1;       % Simulation time step (hours)
            RunTableS = RunModel_AH(DoseTableRAC,flag,dt,CypScore);
    
            % Time-dependent concentration of the drug in the central compartment
            %tsim = RunTableS.t;       % time points of the simulation
            RunTableS.Conc_meth = RunTableS.A2/1000/V_S/309.445; %convert to molar. mw for metadone. log units
            RunTableS.Conc_metab = RunTableS.A3/1000/V_S/309.445; %convert to molar. confirm mw for metabolite. log units
            
            MyTable.conc_methS = RunTableS.Conc_meth;
            MyTable.conc_metabS = RunTableS.Conc_metab;
            
            %calculate the combined block for the parent and metbolite for each channel
            Loewe_methR_IKr = (RunTableR.Conc_meth/DrugPars.IC50s('IKr_methR')) .^ DrugPars.h('IKr_methR')
            Loewe_metabR_IKr = (RunTableR.Conc_meth/DrugPars.IC50s('IKr_metab')) .^ DrugPars.h('IKr_metab')
            Loewe_methS_IKr = (RunTableS.Conc_meth/DrugPars.IC50s('IKr_methS')) .^ DrugPars.h('IKr_methS')
            Loewe_metabS_IKr = (RunTableS.Conc_meth/DrugPars.IC50s('IKr_metab')) .^ DrugPars.h('IKr_metab')

            MyTable.IKrBlock = (Loewe_methR_IKr + Loewe_metabR_IKr + Loewe_methS_IKr + Loewe_metabS_IKr) ./ (Loewe_methR_IKr + Loewe_metabR_IKr + Loewe_methS_IKr + Loewe_metabS_IKr +1);

            Loewe_methR_ICaL = (RunTableR.Conc_meth/DrugPars.IC50s('ICaL_methR')) .^ DrugPars.h('ICaL_methR')
            Loewe_metabR_ICaL = (RunTableR.Conc_meth/DrugPars.IC50s('ICaL_metab')) .^ DrugPars.h('ICaL_metab')
            Loewe_methS_ICaL = (RunTableS.Conc_meth/DrugPars.IC50s('ICaL_methS')) .^ DrugPars.h('ICaL_methS')
            Loewe_metabS_ICaL = (RunTableS.Conc_meth/DrugPars.IC50s('ICaL_metab')) .^ DrugPars.h('ICaL_metab')

            MyTable.ICaLBlock = (Loewe_methR_ICaL + Loewe_metabR_ICaL + Loewe_methS_ICaL + Loewe_metabS_ICaL) ./ (Loewe_methR_ICaL + Loewe_metabR_ICaL + Loewe_methS_ICaL + Loewe_metabS_ICaL +1);

            Loewe_methR_INaL = (RunTableR.Conc_meth/DrugPars.IC50s('INaL_methR')) .^ DrugPars.h('INaL_methR')
            Loewe_metabR_INaL = (RunTableR.Conc_meth/DrugPars.IC50s('INaL_metab')) .^ DrugPars.h('INaL_metab')
            Loewe_methS_INaL = (RunTableS.Conc_meth/DrugPars.IC50s('INaL_methS')) .^ DrugPars.h('INaL_methS')
            Loewe_metabS_INaL = (RunTableS.Conc_meth/DrugPars.IC50s('INaL_metab')) .^ DrugPars.h('INaL_metab')

            MyTable.INaLBlock = (Loewe_methR_INaL + Loewe_metabR_INaL + Loewe_methS_INaL + Loewe_metabS_INaL) ./ (Loewe_methR_INaL + Loewe_metabR_INaL + Loewe_methS_INaL + Loewe_metabS_INaL +1);
            
            %Calculate the risk score for parent and metabolite combined using the two-channel
            % axis-of-arrhythmia by Heitmann et al (unpublished).
            BCaL =  0.81316;    % Coefficient of the axis of arrhythmia in LogGCaLScale
            BKr  = -0.58204;    % Coefficient of the axis of arrhythmia in LogGKrScale
            BNaL = 0.377;       % Coefficient of the axis of arrhythmia in LogGNaLScale
            
            MyTable.RiskScore = BCaL .* log(1-MyTable.ICaLBlock) + BKr .* log(1-MyTable.IKrBlock) + BNaL .* log(1-MyTable.INaLBlock);
                    
                    
            % Plot the simulation results
            figure(1)
            % stackedplot(RunTable, ["A1'", "A2", "A3", "A4", "RiskScore"]);
            %stackedplot(RunTableR);
            
            figure(2)
            %stackedplot(RunTableS);
            
            figure(3)
            MyTable = table2timetable(MyTable);
            
            stackedplot(MyTable);
            
            %plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
                    
            % save the figures and data
            save(strcat('./',formulation, '_', strrep(num2str(CypScore),'.','_'),'/',formulation, '_', strrep(num2str(CypScore),'.','_')),'RunTableR', 'RunTableS', 'MyTable', 'DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', strrep(num2str(CypScore),'.','_'),'/',formulation, '_', strrep(num2str(CypScore),'.','_'),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', strrep(num2str(CypScore),'.','_'),'/',formulation, '_', strrep(num2str(CypScore),'.','_'),'_', 'overlaid','.fig'))
    end  % switch

end % function