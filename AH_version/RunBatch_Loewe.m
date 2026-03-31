function BatchSummary = RunBatch_Loewe(batchInput)
    % Run a batch of Loewe-model scenarios from a table or table-containing file.

    if nargin < 1 || isempty(batchInput)
        error('Provide a batch table or a file path to one.');
    end

    batchTable = loadBatchTable(batchInput);
    requiredVars = ["formulation", "baseline_dose", "overdose", "RF", "BW", "CYPScore"];
    missingVars = requiredVars(~ismember(requiredVars, string(batchTable.Properties.VariableNames)));
    if ~isempty(missingVars)
        error('Batch table is missing required variables: %s', strjoin(missingVars, ', '));
    end

    nRuns = height(batchTable);
    OutputDir = strings(nRuns,1);
    Status = strings(nRuns,1);
    Message = strings(nRuns,1);
    ResultHeight = nan(nRuns,1);
    PeakRiskScore = nan(nRuns,1);
    FinalRiskScore = nan(nRuns,1);

    for row = 1:nRuns
        try
            [result, details] = ModelMain_Loewe( ...
                char(string(batchTable.formulation(row))), ...
                batchTable.CYPScore(row), ...
                batchTable.BW(row), ...
                batchTable.RF(row), ...
                batchTable.baseline_dose(row), ...
                batchTable.overdose(row));

            OutputDir(row) = string(details.OutputDir);
            Status(row) = "ok";
            Message(row) = "";
            ResultHeight(row) = height(result);
            PeakRiskScore(row) = max(result.RiskScore);
            FinalRiskScore(row) = result.RiskScore(end);
        catch ME
            Status(row) = "error";
            Message(row) = string(ME.message);
        end
    end

    BatchSummary = batchTable;
    BatchSummary.OutputDir = OutputDir;
    BatchSummary.Status = Status;
    BatchSummary.Message = Message;
    BatchSummary.ResultHeight = ResultHeight;
    BatchSummary.PeakRiskScore = PeakRiskScore;
    BatchSummary.FinalRiskScore = FinalRiskScore;

    save('BatchSummary.mat', 'BatchSummary');
end

function batchTable = loadBatchTable(batchInput)
    if istable(batchInput)
        batchTable = batchInput;
        return
    end

    if ~(ischar(batchInput) || (isstring(batchInput) && isscalar(batchInput)))
        error('batchInput must be a table or a file path.');
    end

    filePath = char(batchInput);
    [~, ~, ext] = fileparts(filePath);

    switch lower(ext)
        case {'.csv', '.txt'}
            batchTable = readtable(filePath);
        case {'.xlsx', '.xls'}
            batchTable = readtable(filePath);
        case '.mat'
            loaded = load(filePath);
            varNames = fieldnames(loaded);
            tableMask = structfun(@istable, loaded);
            if ~any(tableMask)
                error('No MATLAB table found in %s.', filePath);
            end
            batchTable = loaded.(varNames{find(tableMask, 1, 'first')});
        otherwise
            error('Unsupported batch file type: %s', ext);
    end
end
