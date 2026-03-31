function [result, details] = ModelMain_Loewe(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier)
    % Combined driver for R-, S-, and racemic methadone simulations.

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

    pars = load('DrugPars.mat', 'DrugPars');
    DrugPars = pars.DrugPars;
    DoseTable = BuildDoseTable_AH(formulation, doseOriginal, overdoseMultiplier);

    volumes.R = 176;
    volumes.S = 98.3;
    dt = 0.1;
    mwMethadone = 309.445;
    riskCoeff = struct('BCaL', 0.691, 'BKr', -0.617, 'BNaL', 0.377);
    runId = composeRunId(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier);
    runLabel = composeRunLabel(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier);
    outputDir = fullfile(pwd, runId);

    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end

    switch formulation
        case {'R', 'S'}
            runTable = RunModel_AH(DoseTable, formulation, dt, CypScore, RF, BW);
            runTable = addConcentrations(runTable, volumes.(formulation), mwMethadone);
            runTable = addBlocksAndRisk(runTable, DrugPars, formulation, riskCoeff);

            figStacked = figure();
            plotStackedRun(runTable, runLabel);

            figOverlay = figure();
            plotSingleRun(runTable, runLabel);

            saveRunOutputs(outputDir, runId, runTable, DoseTable, figStacked, figOverlay);

            result = runTable;
            details = buildDetailsStruct(outputDir, formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, DoseTable);

        case 'racemic'
            DoseTableRac = DoseTable;
            DoseTableRac.Conc = DoseTableRac.Conc ./ 2;

            runTableR = RunModel_AH(DoseTableRac, 'R', dt, CypScore, RF, BW);
            runTableR = addConcentrations(runTableR, volumes.R, mwMethadone);

            runTableS = RunModel_AH(DoseTableRac, 'S', dt, CypScore, RF, BW);
            runTableS = addConcentrations(runTableS, volumes.S, mwMethadone);

            result = buildRacemicTable(runTableR, runTableS, DrugPars, riskCoeff);

            figStacked = figure();
            plotStackedRun(result, runLabel);

            figOverlay = figure();
            plotRacemicRun(result, runLabel);

            saveRacemicOutputs(outputDir, runId, runTableR, runTableS, result, DoseTable, figStacked, figOverlay);

            details = buildDetailsStruct(outputDir, formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, DoseTable);
            details.DoseTableRacemic = DoseTableRac;
            details.RunTableR = runTableR;
            details.RunTableS = runTableS;

        otherwise
            error('Unsupported formulation "%s". Use ''R'', ''S'', or ''racemic''.', formulation);
    end
end

function details = buildDetailsStruct(outputDir, formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier, DoseTable)
    details = struct( ...
        'OutputDir', outputDir, ...
        'Formulation', formulation, ...
        'CypScore', CypScore, ...
        'BW', BW, ...
        'RF', RF, ...
        'DoseOriginal', doseOriginal, ...
        'OverdoseMultiplier', overdoseMultiplier, ...
        'DoseTable', DoseTable);
end

function runId = composeRunId(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier)
    runId = sprintf('%s_cyp_%s_bw_%s_rf_%s_dose_%s_od_%s', ...
        formulation, ...
        numToTag(CypScore), ...
        numToTag(BW), ...
        numToTag(RF), ...
        numToTag(doseOriginal), ...
        numToTag(overdoseMultiplier));
end

function tag = numToTag(value)
    tag = strrep(num2str(value, '%.6g'), '.', '_');
end

function label = composeRunLabel(formulation, CypScore, BW, RF, doseOriginal, overdoseMultiplier)
    label = sprintf('%s | CYP %s | BW %s kg | RF %s | Dose %s | OD x%s', ...
        formulation, ...
        num2str(CypScore, '%.6g'), ...
        num2str(BW, '%.6g'), ...
        num2str(RF, '%.6g'), ...
        num2str(doseOriginal, '%.6g'), ...
        num2str(overdoseMultiplier, '%.6g'));
end

function runTable = addConcentrations(runTable, apparentVolume, molecularWeight)
    runTable.Conc_meth = runTable.A2 ./ 1000 ./ apparentVolume ./ molecularWeight;
    runTable.Conc_metab = runTable.A3 ./ 1000 ./ apparentVolume ./ molecularWeight;
end

function runTable = addBlocksAndRisk(runTable, drugPars, enantiomer, riskCoeff)
    runTable.IKrBlock = combinedLoeweBlock(runTable.Conc_meth, runTable.Conc_metab, drugPars, ['IKr_meth' enantiomer], 'IKr_metab');
    runTable.ICaLBlock = combinedLoeweBlock(runTable.Conc_meth, runTable.Conc_metab, drugPars, ['ICaL_meth' enantiomer], 'ICaL_metab');
    runTable.INaLBlock = combinedLoeweBlock(runTable.Conc_meth, runTable.Conc_metab, drugPars, ['INaL_meth' enantiomer], 'INaL_metab');
    runTable.RiskScore = riskScoreFromBlocks(runTable.ICaLBlock, runTable.IKrBlock, runTable.INaLBlock, riskCoeff);
end

function myTable = buildRacemicTable(runTableR, runTableS, drugPars, riskCoeff)
    myTable = table();
    myTable.t = runTableR.t;
    myTable.conc_methR = runTableR.Conc_meth;
    myTable.conc_metabR = runTableR.Conc_metab;
    myTable.conc_methS = runTableS.Conc_meth;
    myTable.conc_metabS = runTableS.Conc_metab;

    totalIKr = loeweTerm(runTableR.Conc_meth, drugPars, 'IKr_methR') + ...
        loeweTerm(runTableR.Conc_metab, drugPars, 'IKr_metab') + ...
        loeweTerm(runTableS.Conc_meth, drugPars, 'IKr_methS') + ...
        loeweTerm(runTableS.Conc_metab, drugPars, 'IKr_metab');
    myTable.IKrBlock = totalIKr ./ (totalIKr + 1);

    totalICaL = loeweTerm(runTableR.Conc_meth, drugPars, 'ICaL_methR') + ...
        loeweTerm(runTableR.Conc_metab, drugPars, 'ICaL_metab') + ...
        loeweTerm(runTableS.Conc_meth, drugPars, 'ICaL_methS') + ...
        loeweTerm(runTableS.Conc_metab, drugPars, 'ICaL_metab');
    myTable.ICaLBlock = totalICaL ./ (totalICaL + 1);

    totalINaL = loeweTerm(runTableR.Conc_meth, drugPars, 'INaL_methR') + ...
        loeweTerm(runTableR.Conc_metab, drugPars, 'INaL_metab') + ...
        loeweTerm(runTableS.Conc_meth, drugPars, 'INaL_methS') + ...
        loeweTerm(runTableS.Conc_metab, drugPars, 'INaL_metab');
    myTable.INaLBlock = totalINaL ./ (totalINaL + 1);

    myTable.RiskScore = riskScoreFromBlocks(myTable.ICaLBlock, myTable.IKrBlock, myTable.INaLBlock, riskCoeff);
    myTable = table2timetable(myTable);
end

function block = combinedLoeweBlock(parentConc, metabConc, drugPars, parentKey, metabKey)
    total = loeweTerm(parentConc, drugPars, parentKey) + loeweTerm(metabConc, drugPars, metabKey);
    block = total ./ (total + 1);
end

function term = loeweTerm(concentration, drugPars, key)
    term = (concentration ./ drugPars.IC50s(key)) .^ drugPars.h(key);
end

function score = riskScoreFromBlocks(ICaLBlock, IKrBlock, INaLBlock, coeff)
    score = coeff.BCaL .* log(1 - ICaLBlock) + coeff.BKr .* log(1 - IKrBlock) + coeff.BNaL .* log(1 - INaLBlock);
end

function plotStackedRun(runTable, runLabel)
    stackedplot(runTable);
    sgtitle(runLabel, 'Interpreter', 'none');
end

function plotSingleRun(runTable, runLabel)
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
    plot(myTable.t, myTable.conc_methR, 'DisplayName', 'conc\_methR');
    hold on;
    plot(myTable.t, myTable.conc_metabR, 'DisplayName', 'conc\_metabR');
    plot(myTable.t, myTable.conc_methS, 'DisplayName', 'conc\_methS');
    plot(myTable.t, myTable.conc_metabS, 'DisplayName', 'conc\_metabS');
    plot(myTable.t, myTable.RiskScore, 'DisplayName', 'RiskScore');
    hold off;
    title(runLabel, 'Interpreter', 'none');
    legend('Location', 'best');
end

function saveRunOutputs(outputDir, runId, runTable, doseTable, figStacked, figOverlay)
    baseName = fullfile(outputDir, runId);
    RunTable = runTable;
    DoseTable = doseTable;
    save(baseName, 'RunTable', 'DoseTable');
    savefig(figStacked, [baseName '_stacked.fig']);
    savefig(figOverlay, [baseName '_overlaid.fig']);
end

function saveRacemicOutputs(outputDir, runId, runTableR, runTableS, myTable, doseTable, figStacked, figOverlay)
    baseName = fullfile(outputDir, runId);
    RunTableR = runTableR;
    RunTableS = runTableS;
    MyTable = myTable;
    DoseTable = doseTable;
    save(baseName, 'RunTableR', 'RunTableS', 'MyTable', 'DoseTable');
    savefig(figStacked, [baseName '_stacked.fig']);
    savefig(figOverlay, [baseName '_overlaid.fig']);
end
