function DrugPars = BuildDrugPars_AH(hillMode)
    % Build the Loewe dose-response parameter table using ideal or real Hill coefficients.

    if nargin < 1 || isempty(hillMode)
        hillMode = 'ideal';
    end

    parameters = ["IKr_methR";
        "IKr_methS";
        "IKr_metab";
        "INaL_methR";
        "INaL_methS";
        "INaL_metab";
        "ICaL_methR";
        "ICaL_methS";
        "ICaL_metab"];

    IC50s = [1.89e-6;
        2.09e-6;
        3.29e-5;
        4.5e-6;
        9.34e-6;
        4.13e-5;
        3.08e-5;
        2.07e-5;
        1.19e-4];

    switch lower(string(hillMode))
        case "ideal"
            h = ones(size(IC50s));
        case "real"
            h = [0.7631;
                0.8279;
                0.7941;
                0.9880;
                1.4110;
                1.5620;
                1.492;
                0.9227;
                0.9776];
        otherwise
            error('Unsupported hillMode "%s". Use ''ideal'' or ''real''.', hillMode);
    end

    DrugPars = table(IC50s, h, RowNames=parameters);
end
