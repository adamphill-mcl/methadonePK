function sys = Aruldhas2021_CT_Adult(flag,CypScore,RF,BW)
    % Aruldhas2021_CT_Adult  Pharmacokinetics of R- and S-Methadone in adults
    % Four-compartment pharmacokinetic model of Methadone and its
    % metabolites. Implemented for the Brain Dynamics Toolbox.
    %
    % SYNTAX
    %   sys = Aruldhas('R');        % initialised for R-Methadone
    %   sys = Aruldhas('S');        % initialised for S-Methadone
    %
    % Example:
    %   flag = 'R';                 % flag must be 'R' or 'S'
    %   sys = Aruldhas2021_CT_Adult(flag,CypScore,RF,BW);   % construct the system struct
    %   gui = bdGUI(sys);           % open the Brain Dynamics GUI
    %
    % Authors
    %   Stewart Heitmann (2024)

    % Handle to our ODE function
    sys.odefun = @odefun;
    if nargin < 3 || isempty(RF)
        RF = 1;
    end
    if nargin < 4 || isempty(BW)
        BW = 70;
    end
        
    switch flag
        case 'R'
            % Population means for R-Methadone (Table 2)
            Ka  = 0.55;                                         % tablet (1/h), per Foster 2003
            CL  = 9.66;                                         % per Foster 2000
            %CYP_score = 2;                                     % Set 2.0 for ultra rapid, 1.5 for rapid, 1 for normal, 0.5 for intermediate and 0 for poor metaboliser
            CLF = (0.217*(1+0.745*(CypScore-1))*(BW/70)^0.75);  % (L/h), note: this equation assumes heterozygous genotype of rs2246709 (not indicated) and a body weight of 70kg (not indicated)
            CL1 = (7.5638) * RF;                                % (L/h) now set these as constant, based on CLF fraction of CL for normal metaboloser
            %CL1 = (CL-CL2)*RF;
            CL3 = 19.99 * RF;                                   % (L/h), per Shi, 2022
            Q4  = 79.8;                                         % (L/h), per Csajka, 2019
            V2  = 145;                                          % (L), per Foster, 2003
            V4  = 444;                                          % (L), per Foster, 2003
            F1  = 0.5;                                          % arbitrarily set, as per Csajka, 2016
        case 'S'
            % Population means for S-Methadone (Table 2)
            Ka  = 0.62;                                         % tablet (1/h), as per Foster, 2003
            CL  = 9.54;                                         % per Foster, 2000
            %CYP_score = 1;                                     % Set 2.0 for ultra rapid, 1.5 for rapid, 1 for normal, 0.5 for intermediate and 0 for poor metaboliser
            CLF = (0.135*(1+0.636*(CypScore-1))*(BW/70)^0.75);  % (L/h), note: this equation assumes a dominant genotype for rs11882424 (not indicated) and a body weight of 70kg (not indicated)
            CL1 = (8.2521) * RF;                                 %  (L/h)
            %CL1 = (CL-CL2)*RF;
            CL3 = 19.99 * RF;                                   % (L/h), per Shi, 2022
            Q4  = 69.1;                                         % (L/h), per Csajka, 2016
            V2  = 94;                                           % (L), per Foster, 2003
            V4  = 246;                                          % (L) Per Foster, 2003
            F1  = 0.5;                                          % arbitrarily set, as per Csajka, 2016
        otherwise
            error('Input flag must be ''R'' or ''S'' only. See help text for Aruldhas2021.');
    end

    % ODE parameter definitions
    sys.pardef = [
        struct('name','Ka',  'value',Ka,  'lim',[0 30])
        struct('name','CL',  'value',CL,  'lim',[0 20])
        struct('name','CL1', 'value',CL1, 'lim',[0 20])
        struct('name','CLF', 'value',CLF, 'lim',[0 1])
        struct('name','CL3', 'value',CL3, 'lim',[0 30])
        struct('name','Q4',  'value',Q4,  'lim',[0 100])
        struct('name','V2',  'value',V2,  'lim',[1 200])
        struct('name','VF',  'value',1,   'lim',[1e-2 2])
        struct('name','V4',  'value',V4,  'lim',[1 500])
        ];
    
    % ODE variable definitions
    sys.vardef = [
        struct('name','A1', 'value',F1, 'lim',[0 200])
        struct('name','A2', 'value',0,  'lim',[0 200])
        struct('name','A3', 'value',0,  'lim',[0 200])
        struct('name','A4', 'value',0,  'lim',[0 200])
        ];

    % Latex (Equations) panel
    sys.panels.bdLatexPanel.title = 'Equations'; 
    sys.panels.bdLatexPanel.latex = { 
        '$\textbf{Aruldhas2021}$'
        ''
        'A four-compartment pharmacokinetic model,'
        '\quad $\dot A_1(t) = -A_1(t) \, K_{12}$'
        '\quad $\dot A_2(t) =  A_1(t) \, K_{12} + A_4(t) \, K_{42} - A_2(t) \, (K_{23} + K_{24} + K_{20})$'
        '\quad $\dot A_3(t) =  A_2(t) \, K_{23} - A_3(t) \, K_{30}$'
        '\quad $\dot A_4(t) =  A_2(t) \, K_{24} - A_4(t) \, K_{42}$'
        'where' 
        '\quad $\dot A_1(t)$ is the concentration of the drug in the gastro-intestinal compartment,'
        '\quad $\dot A_2(t)$ is the concentration of the drug in the central compartment,'
        '\quad $\dot A_3(t)$ is the concentration of the metabolite in the central compartment,'
        '\quad $\dot A_4(t)$ is the concentration of the drug in the peripheral compartment,'
        'and'
        '\quad $K_{12} = K_a$ is the absorbtion rate of the drug into the central compartment from the gut,'
        '\quad $K_{20} = CL1/V2$ is the clearance rate of the drug from the central compartment by other routes,'
        '\quad $K_{23} = CL2/V2$ is the clearance rate of the drug from the central compartment by metabolism,'
        '\quad $K_{30} = CL3/V3$ is the clearance rate of the metabolite from the central compartment by renal elimination,'
        '\quad $K_{24} = Q4/V2$ is the rate of drug transfer from the central compartment into the peripheral compartment,'
        '\quad $K_{42} = Q4/V4$ is the rate of drug transfer from the peripheral compartment into the central compartment,'
        '\quad $CL2 = CL * CLF$,'
        '\quad $CL1$ is an independently specified non-metabolic clearance term,'
        '\quad $V2$ is the volume of the central drug compartment,'
        '\quad $V3 = V2 * VF$ is the volume of the central metabolite compartment,'
        '\quad $VF = 1$ is a scaling factor for defining V3 in relation to V2,'
        '\quad $V4$ is the volume of the peripheral drug compartment.'
        ''
        '$\textbf{Reference}$'
        'Aruldhas BW, Quinney SK, Overholser BR, et al. (2021) Pharmacokinetic modeling of R and'
        '\quad S-Methadone and their metabolites to study the effects of various covariates in post-operative children.'
        '\quad CPT Pharmacometrics Syst Pharmacol (10) 1183-1194.'
        };

    % Time Portrait panel 
    sys.panels.bdTimePortrait(1).selector1={1 1 1};
    sys.panels.bdTimePortrait(1).selector2={2 1 1};
    sys.panels.bdTimePortrait(2).selector1={3 1 1};
    sys.panels.bdTimePortrait(2).selector2={4 1 1};

    % Solver panel
    sys.panels.bdSolverPanel = [];
    
    % Default time span (optional)
    sys.tspan = [0 100];
    
    % ODE solver options (optional)
    sys.odeoption.RelTol = 1e-6;        % Relative Tolerance
    sys.odeoption.InitialStep = 1;      % Required by odeEul solver
end

% The ODE function.
% The variables A and dAdt are both (4x1) vectors.
% All other parameters are scalars.
function dAdt = odefun(t,A,Ka,CL,CL1,CLF, CL3,Q4,V2,VF,V4) 
    A1 = A(1);      % conc of drug in gastro-intestinal compartment 
    A2 = A(2);      % conc of drug in central compartment
    A3 = A(3);      % conc of metabolte in central compartment
    A4 = A(4);      % conc of drug in peripheral compartment

    % Constants
    CL2 = CL*CLF;
    %CL1 = (CL-CL2)*RF;
    V3  = V2*VF;
    K12 = Ka;
    K20 = CL1/V2;
    K23 = CL2/V2;
    K30 = CL3/V3;
    K24 = Q4/V2;
    K42 = Q4/V4;

    % Differential equations
    dA1dt = -A1*K12;
    dA2dt =  A1*K12 + A4*K42 - A2*(K23 + K24 + K20);
    dA3dt =  A2*K23 - A3*K30;
    dA4dt =  A2*K24 - A4*K42;

    dAdt = [dA1dt; dA2dt; dA3dt; dA4dt];
end

