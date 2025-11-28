function update_exchange_file(exchangesMat, savePath, varargin)
    % update_exchange_file(exchangesMat, savePath, mediaFiles...)
    %
    % exchangesMat  - path to original GenericExchanges.mat
    % savePath      - path to updated Exchanges.mat
    % varargin      - ANY number of XLSX files containing a compound ID column

    fprintf("\n=== Updating Exchanges File ===\n");

    %% --------------------------------------------------------
    % Load the existing GenericExchanges.mat
    % ---------------------------------------------------------
    data = load(exchangesMat);

    requiredFields = ["EXC","EXCFormula","ExtraT","ExtraTFormulas"];
    for f = requiredFields
        if ~isfield(data, f)
            error("Field '%s' missing from %s", f, exchangesMat);
        end
    end

    EXC            = data.EXC;
    EXCFormula     = data.EXCFormula;
    ExtraT         = data.ExtraT;
    ExtraTFormulas = data.ExtraTFormulas;

    fprintf("Loaded input exchange dataset: %s\n\n", exchangesMat);

    %% --------------------------------------------------------
    % Collect compound IDs from all media Excel files
    % ---------------------------------------------------------
    all_compounds = {};

    for i = 1:length(varargin)
        mediaFile = varargin{i};
        fprintf("Loading media file: %s\n", mediaFile);

        opts = detectImportOptions(mediaFile, 'VariableNamingRule', 'preserve');
        tbl = readtable(mediaFile, opts);

        vars = string(tbl.Properties.VariableNames);

        % Flexible detection: look for column containing both "compound" AND "id"
        match = vars(contains(lower(vars), "compound") & contains(lower(vars), "id"));

        if isempty(match)
            error('File %s must contain a compound ID column (e.g., "CompoundID").', mediaFile);
        end

        colName = char(match(1));
        fprintf("  → Using column: %s\n", colName);

        compounds = tbl.(colName);
        all_compounds = [all_compounds; compounds];
    end

    compounds = unique(all_compounds);
    fprintf("\nTotal unique compounds discovered: %d\n\n", numel(compounds));

    %% --------------------------------------------------------
    % Generate NEW exchange reactions
    % ---------------------------------------------------------
    [exc_rxns, exc_names] = create_exchange_reactions(compounds);

    % Remove duplicates: only keep reactions NOT already present
    is_new_exc = ~ismember(exc_names, EXC);
    exc_names_new  = exc_names(is_new_exc);
    exc_rxns_new   = exc_rxns(is_new_exc);

    fprintf("Exchange reactions to add: %d (skipped %d duplicates)\n", ...
        numel(exc_names_new), numel(exc_names) - numel(exc_names_new));

    %% --------------------------------------------------------
    % Generate NEW transporter reactions (e → c)
    % ---------------------------------------------------------
    [tr_e_c, tr_e_c_names] = create_transporter_reactions(compounds, 3);

    % Remove duplicates
    is_new_tr = ~ismember(tr_e_c_names, ExtraT);
    tr_e_c_names_new = tr_e_c_names(is_new_tr);
    tr_e_c_new       = tr_e_c(is_new_tr);

    fprintf("Transporter reactions to add: %d (skipped %d duplicates)\n", ...
        numel(tr_e_c_names_new), numel(tr_e_c_names) - numel(tr_e_c_names_new));

    %% --------------------------------------------------------
    % Append ONLY the new reactions
    % ---------------------------------------------------------
    EXC            = [EXC;            exc_names_new(:)];
    EXCFormula     = [EXCFormula;     exc_rxns_new(:)];

    ExtraT         = [ExtraT;         tr_e_c_names_new(:)];
    ExtraTFormulas = [ExtraTFormulas; tr_e_c_new(:)];

    %% --------------------------------------------------------
    % Save updated data back to a MAT file
    % ---------------------------------------------------------
    data.EXC            = EXC;
    data.EXCFormula     = EXCFormula;
    data.ExtraT         = ExtraT;
    data.ExtraTFormulas = ExtraTFormulas;

    fprintf("\nSaving updated Exchanges file to: %s\n", savePath);
    save(savePath, '-struct', 'data');
    fprintf("✔ Saved successfully.\n");

end
