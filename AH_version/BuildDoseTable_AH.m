function DoseTable = BuildDoseTable_AH(formulation, doseOriginal, overdoseMultiplier)
    % Build the fixed daily dosing schedule used by the Loewe workflow.

    if nargin < 2 || isempty(doseOriginal)
        doseOriginal = 120;
    end
    if nargin < 3 || isempty(overdoseMultiplier)
        overdoseMultiplier = 1;
    end

    validateattributes(doseOriginal, {'numeric'}, {'scalar', 'real', 'nonnegative'}, mfilename, 'doseOriginal');
    validateattributes(overdoseMultiplier, {'numeric'}, {'scalar', 'real', 'nonnegative'}, mfilename, 'overdoseMultiplier');

    F1_bioavailability = 0.718;
    dose = doseOriginal * F1_bioavailability;

    ndailyDoses = 26;
    doseSequence = [0; repmat(dose, ndailyDoses, 1)];
    doseSequence(end) = doseSequence(end) * overdoseMultiplier;
    durationHours = repmat(24, numel(doseSequence), 1);
    formulationColumn = repmat(string(formulation), numel(doseSequence), 1);

    DoseTable = table(doseSequence, durationHours, formulationColumn, ...
        'VariableNames', {'Conc', 'dur', 'MethadoneFlag'});
end
