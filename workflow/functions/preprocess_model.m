function preprocess_model = preprocess_model(cobraPath, modelPath, importsFolder, savePath)
% preprocessModel  Preprocess a metabolic model using multiple imported files.
%
%   preprocessedModel = preprocessModel(cobraPath, modelPath, importsFolder, savePath)
%
%   Inputs:
%       cobraPath      – path to COBRA toolbox
%       modelPath      – path to the .mat file containing lgu50
%       importsFolder  – path to folder containing all auxiliary .mat files
%       savePath       – folder where preprocessed model is saved
%
%   Output:
%       preprocessedModel – the processed COBRA model structure

% =========================================================
% Add COBRA Toolbox
% =========================================================
addpath(genpath(cobraPath));

% =========================================================
% Load Base Model (name agnostic)
% =========================================================
vars = load(modelPath);
varNames = fieldnames(vars);

modelVar = '';
for k = 1:numel(varNames)
    candidate = vars.(varNames{k});
    if isstruct(candidate) && isfield(candidate,'rxns') && isfield(candidate,'mets')
        modelVar = varNames{k};
        break;
    end
end

if isempty(modelVar)
    error('No COBRA model found in: %s', modelPath);
end

model = vars.(modelVar);
len = length(model.rxns);

% =========================================================
% Load Compartment Data
% =========================================================
load(fullfile(importsFolder,"CompartmentData.mat"), "CompartmentData");

% =========================================================
% Add Extracellular metabolites
% =========================================================
f = 1:length(model.mets);
model.mets(f) = strcat(model.mets(f), '_c');

% =========================================================
% Adding Exchanges & Transporters
% =========================================================
load(fullfile(importsFolder,"MyExchanges.mat"));
stoich = extract_stoichiometries(EXCFormula);

for i = 1:length(EXC)
    model = addReaction(model, EXC{i}, 'reactionFormula', EXCFormula{i}, ...
        'reversible', 1, 'lowerBound', 0, 'upperBound', 50, ...
        'subSystem', 'Exchanges', 'stoichCoeffList', stoich{i});
end

% --- Transport from extracellular → cytoplasm
stoich = extract_stoichiometries(ExtraTFormulas);
for i = 1:length(ExtraT)
    model = addReaction(model, ExtraT{i}, 'reactionFormula', ExtraTFormulas{i}, ...
        'reversible', 1, 'lowerBound', -50, 'upperBound', 50, ...
        'subSystem', 'ExtracellularTransportToCytoplasm', 'stoichCoeffList', stoich{i});
end

% =========================================================
% Add ATP Maintenance Reaction
% =========================================================
load(fullfile(importsFolder,"ATPM.mat"));
stoich = extract_stoichiometries(atpmformula);

model = addReaction(model, char(atpmID), 'reactionFormula', char(atpmformula), ...
    'reversible', 1, 'lowerBound', 8.39, 'upperBound', 1000, ...
    'subSystem', 'non-growth Associated ATP');

% =========================================================
% Modify reaction directionality (Oxygen, Carbon, ATP, blocked)
% =========================================================

% Oxygen
load(fullfile(importsFolder,"OxyCulprit.mat"));
idx = ismember(model.rxns, Oxylb);   model.lb(idx) = 0;
idx = ismember(model.rxns, OxyUb);   model.ub(idx) = 0;

% Carbon
load(fullfile(importsFolder,"CarbonCulprit.mat"));
idx = ismember(model.rxns, Carbonlb); model.lb(idx) = 0;
idx = ismember(model.rxns, Carbonub); model.ub(idx) = 0;

% ATP
load(fullfile(importsFolder,"ATPCulprit.mat"));
idx = ismember(model.rxns, ATPlb);    model.lb(idx) = 0;
idx = ismember(model.rxns, ATPub);    model.ub(idx) = 0;

% Block unbalanced reactions
load(fullfile(importsFolder,"UnbalancedBlocked.mat"));
idx = ismember(model.rxns, UnbalancedBlocked);
model.lb(idx) = 0;
model.ub(idx) = 0;

% =========================================================
% Remove generic reactions using regex patterns only
% =========================================================
AllFormulas = printRxnFormula(model, model.rxns, 1, 1, 1);

% Regex patterns for generic / placeholder reaction terms
regexPatterns = {
    'acceptor'                 % acceptor
    'donor'                    % donor
    'alcohol'                  % alcohol
    'substrate'                % substrate
    'product'                  % product
    'protein'                  % protein
    'peptide'                  % peptide / polypeptide
    'polypeptide'
    'enzyme'
    '\bR\b'                    % isolated "R"
    'R[-–]'                    % R- group or R– group
    'R[-–]OH'                  % R–OH, R-OH
    '\bX\b'                    % X placeholder
    '\bY\b'                    % Y placeholder
    'R\d+'                     % R1, R2, R3...
};

% Initialize mask
Generic = false(size(AllFormulas));

% Evaluate all regex patterns
for k = 1:numel(regexPatterns)
    pat = regexPatterns{k};
    Generic = Generic | ~cellfun('isempty', regexp(AllFormulas, pat, 'ignorecase'));
end

% Remove reactions flagged as generic
model = removeRxns(model, model.rxns(Generic));

% =========================================================
% Add biomass reaction
% =========================================================

BOFfile = fullfile(importsFolder, "BOF.mat");
tmp = load(BOFfile);

% --- Check required variables ---
if ~isfield(tmp, "biomassEQN")
    error("BOF.mat does not contain variable 'biomassEQN'. Cannot add biomass reaction.");
end

% If stoichiometry is missing, auto-generate it
if isfield(tmp, "stoich")
    stoich = tmp.stoich;
else
    fprintf("No 'stoich' field found in BOF.mat — generating stoichiometries automatically...\n");
    
    % Using the biomass equation (string) as input
    stoich = extract_stoichiometries({tmp.biomassEQN});
end

% --- Convert stoichiometry to numeric vector ---
try
    stoichVec = cell2mat(stoich);
catch
    error("Stoichiometries could not be converted into numeric vector. Check BOF.mat formatting.");
end

% --- Add biomass reaction ---
model = addReaction(model, 'Biomass_rxn', ...
    'reactionFormula', char(tmp.biomassEQN), ...
    'reversible', 0, ...
    'lowerBound', 0, ...
    'upperBound', 50, ...
    'subSystem', 'Biomass', ...
    'stoichCoeffList', stoichVec);

% --- Set objective ---
model = changeObjective(model, {'Biomass_rxn'});

% =========================================================
% Add manually curated quinone reactions
% =========================================================
% load(fullfile(importsFolder,'QinoneRxns.mat'));
% for i = 1:length(ids)
%     model = addReaction(model, ids{i}, 'reactionFormula', Formulas{i}, ...
%         'reversible',1,'lowerBound',-100,'upperBound',100, ...
%         'subSystem','GapFilledManually','stoichCoeffList', stoich{i});
% end

% =========================================================
% Assign compartments
% =========================================================
model.comps = {'c'; 'e'};
model.compNames = {'Cytoplasm'; 'Extracellular'};

f = contains(model.mets,'_c'); model.metComps(f) = 1; model.metCompSymbol(f) = {'c'};
f = contains(model.mets,'_e'); model.metComps(f) = 2; model.metCompSymbol(f) = {'e'};

% =========================================================
% Cleanup fields & bounds
% =========================================================
model.rxnECNumbers = model.eccodes;
model = rmfield(model,"eccodes");
model.rxnECNumbers(len:length(model.rxns)) = {'NA'};

model.CompartmentData = CompartmentData;

model.metSEEDID = erase(model.mets, {'_c','_e'});

model.lb(model.lb==-1000) = -100;
model.ub(model.ub==1000) = 100;

% =========================================================
% FINAL MODEL OUTPUT
% =========================================================
% Ensure the directory exists
saveDir = fileparts(savePath);
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Ensure the file ends with .mat
[~,~,ext] = fileparts(savePath);
if isempty(ext)
    savePath = strcat(savePath, '.mat');
elseif ~strcmp(ext, '.mat')
    error('The output file must have a .mat extension.');
end

% Save final model
preprocessedModel = model;
save(savePath, 'preprocessedModel');

end

