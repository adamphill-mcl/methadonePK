function [result, details] = ModelMain_Hand(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, hillMode)
    % Combined driver for R-, S-, and racemic methadone simulations using the Hand model.

    % Apply defaults so callers can omit trailing arguments during ad hoc runs.
    if nargin < 2 || isempty(CypScore)
        CypScore = 1;
    end
    if nargin < 3 || isempty(BW)
        BW = 70;
    end
    if nargin < 4 || isempty(RF)
        RF = 1;
    end
    if nargin < 5 || isempty(doseOriginal)
        doseOriginal = 120;
    end
    if nargin < 6 || isempty(overdoseMultiplier)
        overdoseMultiplier = 1;
    end
    if nargin < 7 || isempty(hillMode)
        hillMode = 'ideal';
    end

    % Load the pharmacodynamic parameters (IC50 and Hill slope) for the
    % requested Hill-mode assumption, then define model-wide constants.
    DrugPars = BuildDrugPars_AH(hillMode);
    volumes.R = 176;
    volumes.S = 98.3;
    dt = 0.1;
    mwMethadone = 309.445;
    riskCoeff = struct('BCaL', 0.691, 'BKr', -0.617, 'BNaL', 0.377);

    % Build names used both for output folders and figure titles so each run
    % is self-describing and does not overwrite previous outputs.
    runId = composeRunId(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, hillMode);
    runLabel = composeRunLabel(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, hillMode);
    outputDir = fullfile(pwd, runId);

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    switch formulation
        case {'R', 'S'}
            % For single-enantiomer runs, build one dose schedule, run the PK
            % model, convert amounts to molar concentrations, then combine the
            % parent and metabolite effects with the Hand model on each channel.
            DoseTable = BuildDoseTable_AH(formulation, doseOriginal, overdoseMultiplier);
            runTable = RunModel_AH(DoseTable, formulation, dt, CypScore, RF, BW);
            runTable = addConcentrations(runTable, volumes.(formulation), mwMethadone);
            runTable = addBlocksAndRiskHand(runTable, DrugPars, formulation, riskCoeff);

            % Produce the same two plot styles used by the earlier workflow so
            % Hand outputs can be inspected in the same way as Loewe outputs.
            figStacked = figure();
            plotStackedRun(runTable, runLabel);

            figOverlay = figure();
            plotSingleRun(runTable, runLabel);

            saveRunOutputs(outputDir, runId, runTable, DoseTable, figStacked, figOverlay);

            result = runTable;
            details = buildDetailsStruct(outputDir, formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, hillMode, DoseTable, DrugPars);

        case 'racemic'
            % For racemic dosing, split the original oral dose evenly before
            % applying the enantiomer-specific bioavailability factors. This
            % yields distinct R and S dose tables that are then simulated
            % independently through the PK model.
            DoseTableR = BuildDoseTable_AH('R', doseOriginal ./ 2, overdoseMultiplier);
            DoseTableS = BuildDoseTable_AH('S', doseOriginal ./ 2, overdoseMultiplier);

            runTableR = RunModel_AH(DoseTableR, 'R', dt, CypScore, RF, BW);
            runTableR = addConcentrations(runTableR, volumes.R, mwMethadone);

            runTableS = RunModel_AH(DoseTableS, 'S', dt, CypScore, RF, BW);
            runTableS = addConcentrations(runTableS, volumes.S, mwMethadone);

            % Combine the four racemic contributors (R parent, R metabolite,
            % S parent, S metabolite) with the Hand model channel by channel.
            result = buildRacemicTableHand(runTableR, runTableS, DrugPars, riskCoeff);

            figStacked = figure();
            plotStackedRun(result, runLabel);

            figOverlay = figure();
            plotRacemicRun(result, runLabel);

            saveRacemicOutputs(outputDir, runId, runTableR, runTableS, result, DoseTableR, DoseTableS, figStacked, figOverlay);

            details = buildDetailsStruct(outputDir, formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, hillMode, [], DrugPars);
            details.DoseTableR = DoseTableR;
            details.DoseTableS = DoseTableS;
            details.RunTableR = runTableR;
            details.RunTableS = runTableS;

        otherwise
            error('Unsupported formulation "%s". Use ''R'', ''S'', or ''racemic''.', formulation);
    end
end

function details = buildDetailsStruct(outputDir, formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, hillMode, DoseTable, DrugPars)
    % Package the simulation settings and key inputs alongside the outputs so
    % downstream analysis can recover exactly how a result was generated.
    details = struct( ...
        'OutputDir', outputDir, ...
        'Formulation', formulation, ...
        'CypScore', CypScore, ...
        'BW', BW, ...
        'RF', RF, ...
        'DoseOriginal', doseOriginal, ...
        'OverdoseMultiplier', overdoseMultiplier, ...
        'HillMode', char(string(hillMode)), ...
        'DrugPars', DrugPars, ...
        'DoseTable', DoseTable);
end

function runId = composeRunId(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, hillMode)
    % Folder-safe identifier used for output files.
    runId = sprintf('%s_hand_hill_%s_cyp_%s_bw_%s_rf_%s_dose_%s_od_%s', ...
        formulation, ...
        char(string(hillMode)), ...
        numToTag(CypScore), ...
        numToTag(BW), ...
        numToTag(RF), ...
        numToTag(doseOriginal), ...
        numToTag(overdoseMultiplier));
end

function tag = numToTag(value)
    tag = strrep(num2str(value, '%.6g'), '.', '_');
end

function label = composeRunLabel(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, hillMode)
    % Human-readable figure title summarising the run conditions.
    label = sprintf('%s | Hand | Hill %s | CYP %s | BW %s kg | RF %s | Dose %s | OD x%s', ...
        formulation, ...
        char(string(hillMode)), ...
        num2str(CypScore, '%.6g'), ...
        num2str(BW, '%.6g'), ...
        num2str(RF, '%.6g'), ...
        num2str(doseOriginal, '%.6g'), ...
        num2str(overdoseMultiplier, '%.6g'));
end

function runTable = addConcentrations(runTable, apparentVolume, molecularWeight)
    % Convert central-compartment amounts to molar concentrations for the
    % parent drug (A2) and metabolite (A3). These concentrations drive the
    % downstream channel block calculations.
    runTable.Conc_meth = runTable.A2 ./ 1000 ./ apparentVolume ./ molecularWeight;
    runTable.Conc_metab = runTable.A3 ./ 1000 ./ apparentVolume ./ molecularWeight;
end

function runTable = addBlocksAndRiskHand(runTable, drugPars, enantiomer, riskCoeff)
    % Combine parent and metabolite channel block using the Hand model. Each
    % channel is treated independently, then the three channel blocks are
    % collapsed into the scalar risk score.
    runTable.IKrBlock = handBlockFromContributors([runTable.Conc_meth, runTable.Conc_metab], ...
        {['IKr_meth' enantiomer], 'IKr_metab'}, drugPars);
    runTable.ICaLBlock = handBlockFromContributors([runTable.Conc_meth, runTable.Conc_metab], ...
        {['ICaL_meth' enantiomer], 'ICaL_metab'}, drugPars);
    runTable.INaLBlock = handBlockFromContributors([runTable.Conc_meth, runTable.Conc_metab], ...
        {['INaL_meth' enantiomer], 'INaL_metab'}, drugPars);
    runTable.RiskScore = riskScoreFromBlocks(runTable.ICaLBlock, runTable.IKrBlock, runTable.INaLBlock, riskCoeff);
end

function myTable = buildRacemicTableHand(runTableR, runTableS, drugPars, riskCoeff)
    % Assemble a racemic output table so the final object mirrors the older
    % racemic workflow, but with Hand-model channel blocks.
    myTable = table();
    myTable.t = runTableR.t;
    myTable.conc_methR = runTableR.Conc_meth;
    myTable.conc_metabR = runTableR.Conc_metab;
    myTable.conc_methS = runTableS.Conc_meth;
    myTable.conc_metabS = runTableS.Conc_metab;

    % Compute each channel block from the full set of racemic contributors.
    myTable.IKrBlock = handBlockFromContributors( ...
        [runTableR.Conc_meth, runTableR.Conc_metab, runTableS.Conc_meth, runTableS.Conc_metab], ...
        {'IKr_methR', 'IKr_metab', 'IKr_methS', 'IKr_metab'}, drugPars);

    myTable.ICaLBlock = handBlockFromContributors( ...
        [runTableR.Conc_meth, runTableR.Conc_metab, runTableS.Conc_meth, runTableS.Conc_metab], ...
        {'ICaL_methR', 'ICaL_metab', 'ICaL_methS', 'ICaL_metab'}, drugPars);

    myTable.INaLBlock = handBlockFromContributors( ...
        [runTableR.Conc_meth, runTableR.Conc_metab, runTableS.Conc_meth, runTableS.Conc_metab], ...
        {'INaL_methR', 'INaL_metab', 'INaL_methS', 'INaL_metab'}, drugPars);

    myTable.RiskScore = riskScoreFromBlocks(myTable.ICaLBlock, myTable.IKrBlock, myTable.INaLBlock, riskCoeff);
    myTable = table2timetable(myTable);
end

function block = handBlockFromContributors(concentrationMatrix, keys, drugPars)
    % Evaluate Hand-model block for one channel.
    %
    % Each row of concentrationMatrix contains the contributor
    % concentrations present at a particular simulation timepoint. The keys
    % array maps those concentrations to the corresponding IC50 and Hill
    % slope entries in DrugPars.
    nRows = size(concentrationMatrix, 1);
    block = zeros(nRows, 1);

    % Precompute the effect grid used to numerically reconstruct the Hand
    % mixture curve. The grid is dense near 0 and 1 where Hill curves are
    % most nonlinear.
    xGrid = handEffectGrid();

    for row = 1:nRows
        cRow = concentrationMatrix(row, :);
        keep = cRow > 0;

        % If all contributors are absent, the channel block is zero.
        if ~any(keep)
            continue
        end

        cRow = cRow(keep);
        keyRow = keys(keep);

        % If only one contributor is present, the Hand model reduces to the
        % original mono-therapy Hill curve for that contributor.
        if numel(cRow) == 1
            block(row) = monoHillBlock(cRow, drugPars.IC50s(keyRow{1}), drugPars.h(keyRow{1}));
            continue
        end

        % Hand combines drugs by mixture ratio, so use the relative
        % contribution of each concentration to the total concentration at
        % this timepoint.
        totalDose = sum(cRow);
        weights = cRow ./ totalDose;
        sMix = zeros(size(xGrid));

        for idx = 1:numel(cRow)
            ic50 = drugPars.IC50s(keyRow{idx});
            hill = drugPars.h(keyRow{idx});

            % The Hand effect-sensitivity curve is the derivative of effect
            % with respect to dose, expressed as a function of effect rather
            % than concentration. A weighted sum of these sensitivities
            % defines the additive mixture.
            sMix = sMix + weights(idx) .* hillEffectSensitivity(xGrid, ic50, hill);
        end

        % Recover the mixture dose-effect curve by integrating dc/dx = 1/s(x)
        % over the effect grid. Then interpolate the effect reached at the
        % total contributor concentration observed at this timepoint.
        doseCurve = [0; cumtrapz(xGrid, 1 ./ sMix)];
        effectCurve = [0; xGrid];
        [doseCurve, uniqueIdx] = unique(doseCurve, 'stable');
        effectCurve = effectCurve(uniqueIdx);

        block(row) = interp1(doseCurve, effectCurve, totalDose, 'linear', 1 - 1e-10);
        block(row) = min(max(block(row), 0), 1 - 1e-10);
    end
end

function xGrid = handEffectGrid()
    % Cached nonuniform effect grid used for Hand-model numerical
    % integration. The dense sampling near the edges improves stability when
    % Hill slopes create very sharp curvature close to 0 or 1 block.
    persistent cachedGrid
    if isempty(cachedGrid)
        low = logspace(-10, -2, 300);
        mid = linspace(0.01, 0.99, 1200);
        high = 1 - fliplr(logspace(-10, -2, 300));
        cachedGrid = unique([low, mid, high]');
    end
    xGrid = cachedGrid;
end

function sensitivity = hillEffectSensitivity(effect, ic50, hill)
    % Closed-form Hand effect-sensitivity for a Hill curve with Emax = 1.
    % This is s(x) = f'(f^{-1}(x)), which is the central quantity the Hand
    % model averages across mixture components.
    effect = min(max(effect, 1e-12), 1 - 1e-12);
    sensitivity = (hill ./ ic50) .* effect .^ ((hill - 1) ./ hill) .* (1 - effect) .^ ((hill + 1) ./ hill);
    sensitivity = max(sensitivity, realmin);
end

function block = monoHillBlock(concentration, ic50, hill)
    % Standard single-agent Hill block curve used when only one contributor
    % is present at a timepoint.
    ratio = (concentration ./ ic50) .^ hill;
    block = ratio ./ (1 + ratio);
end

function score = riskScoreFromBlocks(ICaLBlock, IKrBlock, INaLBlock, coeff)
    % Collapse channel-specific block into the scalar risk score used
    % elsewhere in the workflow.
    score = coeff.BCaL .* log(1 - ICaLBlock) + coeff.BKr .* log(1 - IKrBlock) + coeff.BNaL .* log(1 - INaLBlock);
end

function plotStackedRun(runTable, runLabel)
    % Stacked overview of all tracked outputs for quick qualitative review.
    stackedplot(runTable);
    sgtitle(runLabel, 'Interpreter', 'none');
end

function plotSingleRun(runTable, runLabel)
    % Overlay the main PK states and the derived risk score in one panel.
    plot(runTable.t, runTable.A1, 'DisplayName', 'A1');
    hold on;
    plot(runTable.t, runTable.A2, 'DisplayName', 'A2');
    plot(runTable.t, runTable.A3, 'DisplayName', 'A3');
    plot(runTable.t, runTable.A4, 'DisplayName', 'A4');
    plot(runTable.t, runTable.RiskScore, 'DisplayName', 'RiskScore');
    hold off;
    title(runLabel, 'Interpreter', 'none');
    legend('Location', 'best');
end

function plotRacemicRun(myTable, runLabel)
    % Overlay racemic component concentrations and the resulting risk score.
    plot(myTable.t, myTable.conc_methR, 'DisplayName', 'conc_methR');
    hold on;
    plot(myTable.t, myTable.conc_metabR, 'DisplayName', 'conc_metabR');
    plot(myTable.t, myTable.conc_methS, 'DisplayName', 'conc_methS');
    plot(myTable.t, myTable.conc_metabS, 'DisplayName', 'conc_metabS');
    plot(myTable.t, myTable.RiskScore, 'DisplayName', 'RiskScore');
    hold off;
    title(runLabel, 'Interpreter', 'none');
    legend('Location', 'best');
end

function saveRunOutputs(outputDir, runId, runTable, doseTable, figStacked, figOverlay)
    % Persist single-enantiomer outputs using the same file conventions as
    % the earlier model for easier side-by-side comparison.
    baseName = fullfile(outputDir, runId);
    RunTable = runTable;
    DoseTable = doseTable;
    save(baseName, 'RunTable', 'DoseTable');
    savefig(figStacked, [baseName '_stacked.fig']);
    savefig(figOverlay, [baseName '_overlaid.fig']);
end

function saveRacemicOutputs(outputDir, runId, runTableR, runTableS, myTable, doseTableR, doseTableS, figStacked, figOverlay)
    % Persist racemic outputs, including the separate R and S dose tables
    % used to generate the run.
    baseName = fullfile(outputDir, runId);
    RunTableR = runTableR;
    RunTableS = runTableS;
    MyTable = myTable;
    DoseTableR = doseTableR;
    DoseTableS = doseTableS;
    save(baseName, 'RunTableR', 'RunTableS', 'MyTable', 'DoseTableR', 'DoseTableS');
    savefig(figStacked, [baseName '_stacked.fig']);
    savefig(figOverlay, [baseName '_overlaid.fig']);
end
