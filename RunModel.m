% RunModel runs the Aruldhas pharmacokinetic model using the dosing
% schedule in DoseTable (see DoseTableScript). The flag ('R' or 'S')
% denotes the sub-type of methadone simulated by the Aruldhas model.
% The time step of the simulation output is defined by dt.
% The results are returned in table format.
function RunTable = RunModel(DoseTable,flag,dt)
    % Construct an instance of the model
    sys = Aruldhas2021_CT2025(flag);       % construct the system
    sys = bdSetVar(sys,'A1',0);     % set initial conditions to zero
    sys = bdSetVar(sys,'A2',0);
    sys = bdSetVar(sys,'A3',0);
    sys = bdSetVar(sys,'A4',0);
    sys.tspan = [-1 0];             % set time span to end at zero
    %bdGUI(sys);

    % init result 
    RunTable = table();

    % for each row in DoseTable ...
    ndoses = size(DoseTable,1);
    for row=1:ndoses
        % Get the current dose and duration
        dose = DoseTable.Conc(row);
        dur  = DoseTable.dur(row);
    
        % Add the dose to the A1 compartment 
        A1 = bdGetVar(sys,'A1');        % Get the final value of A1 from the previous run
        A1 = A1 + dose;                 % Add the drug
        sys = bdSetVar(sys,'A1',A1);    % Apply the new value of A1
    
        % Extend the model for the specified duration
        sys.tspan = sys.tspan(2) + [0 dur];  % Advance the time span from the previous run
        [sys,sol] = bdEvolve(sys,1);                 % Run and evolve.
    
        % Extract the results of the simulation
        tdomain = sys.tspan(1) : dt : sys.tspan(2)-dt;
        A1 = bdEval(sol,tdomain,sys,'A1');            % time signal of A1 from sol
        A2 = bdEval(sol,tdomain,sys,'A2');            % time signal of A2 from sol
        A3 = bdEval(sol,tdomain,sys,'A3');            % time signal of A3 from sol
        A4 = bdEval(sol,tdomain,sys,'A4');            % time signal of A4 from sol
    
        % Reformat time
        tpoints = duration(tdomain',0,0);

        % Append the results
        RunTable = [RunTable ; struct2table(struct('t',tpoints, 'A1',A1', 'A2',A2','A3',A3', 'A4',A4'))];
    end

    % Convert the results to a timetable
    RunTable = table2timetable(RunTable);
end