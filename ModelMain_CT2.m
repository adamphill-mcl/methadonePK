function [outputArg1,outputArg2] = ModelMain_AH(formulation,CypScore)
    %% Combined script for racaemic methadone. We can run R and S individually, 
    % but need to combine before calculating the combined risk score
    
    % Ensure the Brain Dynamics Toolbox is in the matlab PATH.
%     addpath ../bdtoolbox-2023a/
    
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
            %stackedplot(RunTable);
            figure(2)
            %plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
            
%             save the figures and data
            save(strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore)),'RunTable','DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'overlaid','.fig'))
        
        case 'S'
            mkdir(strcat(formulation, '_', num2str(CypScore)))
            flag = 'S';     % Aruldhas flag must be 'R' or 'S' methadone
            dt = 0.1;       % Simulation time step (hours)
            RunTable = RunModel_AH(DoseTable,flag,dt,CypScore);
    
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
            %stackedplot(RunTable);
            figure(2)
            %plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
            
%             save the figures and data
            save(strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore)),'RunTable','DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'overlaid','.fig'))
        
        case 'racemic'
            mkdir(strcat(formulation, '_', num2str(CypScore)))

            %we are going to run this twice, with half conc each of R and S enantiomers
            DoseTableR = DoseTable;
            DoseTableR.Conc = DoseTable.Conc./2;
            DoseTableS = DoseTable;
            DoseTableS.Conc = DoseTableS.Conc./2;
                                 
            % R-Methadone -----------
            flag = 'R';     % Aruldhas flag must be 'R' or 'S' methadone
            dt = 0.1;       % Simulation time step (hours)
            RunTableR = RunModel_AH(DoseTableR,flag,dt,CypScore);
    
            % Time-dependent concentration of the drug in the central compartment
            %tsim = RunTableR.t;       % time points of the simulation
            RunTableR.Conc_meth = log10(RunTableR.A2/1000/V_R/309.445); %above doesnt take into account the apparent volume of compartment
            RunTableR.Conc_metab = log10(RunTableR.A3/1000/V_R/309.445); %convert to molar. confirm mw for metabolite. log units
                       
            % S-Methadone ------------
            flag = 'R';     % Aruldhas flag must be 'R' or 'S' methadone
            %dt = 0.1;       % Simulation time step (hours)
            RunTableS = RunModel_AH(DoseTableS,flag,dt,CypScore);
    
            % Time-dependent concentration of the drug in the central compartment
            %tsim = RunTableS.t;       % time points of the simulation
%             RunTableS.Conc_meth = log10(RunTableS.A2/1000/V_S/309.445); %convert to molar. mw for metadone. log units
%             RunTableS.Conc_metab = log10(RunTableS.A3/1000/V_S/309.445); %convert to molar. confirm mw for metabolite. log units
            RunTableS.Conc_meth = log10(RunTableR.A2/1000/V_R/309.445); %above doesnt take into account the apparent volume of compartment
            RunTableS.Conc_metab = log10(RunTableR.A3/1000/V_R/309.445); %convert to molar. confirm mw for metabolite. log units                        
            
            %calculate the combined block for the parent and metbolite for each channel
            IKrBlock_partEnanR = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTableR.Conc_meth) .* DrugPars.h('IKr_methR'))));
%             IKrBlock_partEnanS = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methS')) - RunTableS.Conc_meth) .* DrugPars.h('IKr_methS'))));
            IKrBlock_partEnanS = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTableS.Conc_meth) .* DrugPars.h('IKr_methR'))));
            IKrBlock_partMetabR =  (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_metab')) - RunTableR.Conc_metab) .* DrugPars.h('IKr_metab'))));
            IKrBlock_partMetabS =  (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_metab')) - RunTableS.Conc_metab) .* DrugPars.h('IKr_metab'))));
            IKrBlock_R = IKrBlock_partEnanR + ((1-IKrBlock_partEnanR) .* IKrBlock_partMetabR);
            IKrBlock_S = IKrBlock_partEnanS + ((1-IKrBlock_partEnanS) .* IKrBlock_partMetabS);
            IKrBlock_Total = IKrBlock_R + ((1-IKrBlock_R) .* IKrBlock_S);              
            RunTableR.IKrBlock = IKrBlock_R;
            RunTableS.IKrBlock = IKrBlock_S;


            ICaLBlock_partEnanR = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methR')) - RunTableR.Conc_meth) .* DrugPars.h('ICaL_methR'))));
            ICaLBlock_partEnanS = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methS')) - RunTableS.Conc_meth) .* DrugPars.h('ICaL_methS'))));
            ICaLBlock_partMetabR =  (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_metab')) - RunTableR.Conc_metab) .* DrugPars.h('ICaL_metab'))));
            ICaLBlock_partMetabS =  (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_metab')) - RunTableS.Conc_metab) .* DrugPars.h('ICaL_metab'))));
            ICaLBlock_R = ICaLBlock_partEnanR + ((1-ICaLBlock_partEnanR) .* ICaLBlock_partMetabR);
            ICaLBlock_S = ICaLBlock_partEnanS + ((1-ICaLBlock_partEnanS) .* ICaLBlock_partMetabS);
            ICaLBlock_Total = ICaLBlock_R + ((1-ICaLBlock_R) .* ICaLBlock_S);
            RunTableR.ICaLBlock = ICaLBlock_R;
            RunTableS.ICaLBlock = ICaLBlock_S;
            
            INaLBlock_partEnanR = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methR')) - RunTableR.Conc_meth) .* DrugPars.h('INaL_methR'))));
            INaLBlock_partEnanS = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methS')) - RunTableS.Conc_meth) .* DrugPars.h('INaL_methS'))));
            INaLBlock_partMetabR =  (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_metab')) - RunTableR.Conc_metab) .* DrugPars.h('INaL_metab'))));
            INaLBlock_partMetabS =  (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_metab')) - RunTableS.Conc_metab) .* DrugPars.h('INaL_metab'))));
            INaLBlock_R = INaLBlock_partEnanR + ((1-INaLBlock_partEnanR) .* INaLBlock_partMetabR);
            INaLBlock_S = INaLBlock_partEnanS + ((1-INaLBlock_partEnanS) .* INaLBlock_partMetabS);
            INaLBlock_Total = INaLBlock_R + ((1-INaLBlock_R) .* INaLBlock_S);
            RunTableR.INaLBlock = INaLBlock_R;
            RunTableS.INaLBlock = INaLBlock_S;
                      
            %Old AH R equations
            %RunTableR.IKrBlock_R = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTableR.Conc_meth) .* DrugPars.h('IKr_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methR')) - RunTableR.Conc_meth) .* DrugPars.h('IKr_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_metab')) - RunTableR.Conc_metab) .* DrugPars.h('IKr_metab'))));
            %RunTableR.ICaLBlock_R = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methR')) - RunTableR.Conc_meth) .* DrugPars.h('ICaL_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methR')) - RunTableR.Conc_meth) .* DrugPars.h('ICaL_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_metab')) - RunTableR.Conc_metab) .* DrugPars.h('ICaL_metab'))));
            %RunTableR.INaLBlock_R = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methR')) - RunTableR.Conc_meth) .* DrugPars.h('INaL_methR')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methR')) - RunTableR.Conc_meth) .* DrugPars.h('INaL_methR'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_metab')) - RunTableR.Conc_metab) .* DrugPars.h('INaL_metab'))));
            
            %Old AH S equations           
            %RunTableS.IKrBlock_S = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methS')) - RunTableS.Conc_meth) .* DrugPars.h('IKr_methS')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_methS')) - RunTableS.Conc_meth) .* DrugPars.h('IKr_methS'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('IKr_metab')) - RunTableS.Conc_metab) .* DrugPars.h('IKr_metab'))));
            %RunTableS.ICaLBlock_S = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methS')) - RunTableS.Conc_meth) .* DrugPars.h('ICaL_methS')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_methS')) - RunTableS.Conc_meth) .* DrugPars.h('ICaL_methS'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('ICaL_metab')) - RunTableS.Conc_metab) .* DrugPars.h('ICaL_metab'))));
            %RunTableS.INaLBlock_S = (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methS')) - RunTableS.Conc_meth) .* DrugPars.h('INaL_methS')))) + (1 - (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_methS')) - RunTableS.Conc_meth) .* DrugPars.h('INaL_methS'))))) .* (1 ./ (1 + 10.^((log10(DrugPars.IC50s('INaL_metab')) - RunTableS.Conc_metab) .* DrugPars.h('INaL_metab'))));

            
            %Calculate the risk score for parent and metabolite combined using the two-channel
            % axis-of-arrhythmia by Heitmann et al (unpublished).
            BCaL =  0.81316;    % Coefficient of the axis of arrhythmia in LogGCaLScale
            BKr  = -0.58204;    % Coefficient of the axis of arrhythmia in LogGKrScale
            BNaL = 0.377;       % Coefficient of the axis of arrhythmia in LogGNaLScale
            
            RiskScore = BCaL .* log(1-ICaLBlock_Total) + BKr .* log(1-IKrBlock_Total) + BNaL .* log(1-INaLBlock_Total);
                    
                    
            % Plot the simulation results
            figure(1)
            % stackedplot(RunTable, ["A1'", "A2", "A3", "A4", "RiskScore"]);
            %stackedplot(RunTableR);
            
            figure(2)
            %stackedplot(RunTableS);
            
            figure(3)
            MyTable = table();
            MyTable.t = RunTableR.t;
            MyTable.concR_enan = RunTableR.Conc_meth;
            MyTable.concR_metab = RunTableR.Conc_metab;
            MyTable.concS_enan = RunTableS.Conc_meth;
            MyTable.concS_metab = RunTableS.Conc_metab;
            MyTable.IKr_Block = IKrBlock_Total;
            MyTable.ICaL_Block = ICaLBlock_Total;
            MyTable.INaL_Block = INaLBlock_Total;
            %MyTable.s = RunTableS.t; %these are probably identical, but Stewart suggested there's a chance that time series can be calculated varyingly between runs
            %MyTable.RA1 = RunTableR.A1;
            %MyTable.RA2 = RunTableR.A2;
            %MyTable.RA3 = RunTableR.A3;
            %MyTable.RA4 = RunTableR.A4;
            %MyTable.SA1 = RunTableS.A1;
            %MyTable.SA2 = RunTableS.A2;
            %MyTable.SA3 = RunTableS.A3;
            %MyTable.SA4 = RunTableS.A4;
            MyTable.RiskScore =RiskScore;
            
            MyTable = table2timetable(MyTable);
            
            stackedplot(MyTable);
            
            %plot(RunTable.t,RunTable.A1,'DisplayName','RunTable.A1');hold on;plot(RunTable.t,RunTable.A2,'DisplayName','RunTable.A2');plot(RunTable.t,RunTable.A3,'DisplayName','RunTable.A3');plot(RunTable.t,RunTable.A4,'DisplayName','RunTable.A4');plot(RunTable.t,RunTable.RiskScore,'DisplayName','RunTable.RiskScore');hold off;
                    
            % save the figures and data
            save(strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore)),'RunTableR', 'RunTableS', 'MyTable', 'DoseTable') %change this every time
            savefig(figure(1),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'stacked','.fig'))
            savefig(figure(2),strcat('./',formulation, '_', num2str(CypScore),'/',formulation, '_', num2str(CypScore),'_', 'overlaid','.fig'))
    end  % switch

end % function