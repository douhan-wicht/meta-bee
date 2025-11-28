function update_exchange_file(exchangesMat, savePath, varargin)
    % update_exchange_file(exchangesMat, savePath, mediaFiles...)
    %
    % exchangesMat  - path to original Exchanges.mat
    % savePath      - path to updated Exchanges.mat
    % varargin      - ANY number of XLSX files containing a compound ID column

    %% Load the existing Exchanges.mat (preserve ALL fields)
    data = load(exchangesMat);

    % Check for required fields
    if isfield(data, 'EXC') && isfield(data, 'EXCFormula')
        EXC = data.EXC;
        EXCFormula = data.EXCFormula;
    else
        error('EXC or EXCFormula not found in %s', exchangesMat);
    end

    %% Collect compound IDs from all provided Excel files
    all_compounds = {};

    for i = 1:length(varargin)
        mediaFile = varargin{i};
        fprintf("Loading media file: %s\n", mediaFile);

        % Load using detectImportOptions to mitigate messy headers
        opts = detectImportOptions(mediaFile, 'VariableNamingRule', 'preserve');
        tbl = readtable(mediaFile, opts);

        % Debug output
        disp("Detected columns:");
        disp(tbl.Properties.VariableNames);

        % Find candidate ID columns (flexible matching)
        vars = string(tbl.Properties.VariableNames);

        match = vars(contains(lower(vars), "compound") & contains(lower(vars), "id"));

        % Error if no matching column found
        if isempty(match)
            error('File %s must contain a column with both "compound" and "id" in its header.', mediaFile);
        end

        % Use the first detected good column
        colName = char(match(1));
        fprintf("Using detected CompoundID column: %s\n", colName);

        compounds = tbl.(colName);
        all_compounds = [all_compounds; compounds];
    end

    %% Deduplicate compound list
    compounds = unique(all_compounds);
    fprintf("Total unique compounds loaded: %d\n", numel(compounds));

    %% Generate exchange reactions
    [exc_rxns, exc_names] = create_exchange_reactions(compounds);

    %% Generate transporter reactions (_e to _c)
    [tr_e_c, tr_e_c_names] = create_transporter_reactions(compounds, 3);

    %% Append new reactions to existing EXC and EXCFormula
    EXC = [
        EXC;
        exc_names(:);   % EXC_BOTH_%s_e
    ];

    EXCFormula = [
        EXCFormula;
        exc_rxns(:);    % %s_e ->
    ];

    ExtraT = [
        ExtraT;
        tr_e_c_names(:) % T_e_to_c_%s_e
    ];

    ExtraTFormulas = [
        ExtraTFormulas;
        tr_e_c(:)       % %s_e <=> %s_c
    ];

    %% Update the data structure while keeping ALL original fields
    data.EXC = EXC;
    data.EXCFormula = EXCFormula;

    %% Save updated structure back to .mat, preserving all fields
    fprintf("Saving updated Exchanges file to %s\n", savePath);
    save(savePath, '-struct', 'data');
end
