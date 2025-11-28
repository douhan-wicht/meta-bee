function [exchange_reactions, reaction_names] = create_exchange_reactions(compoundIDs)
    % Generates exchange reactions and their names for a list of compound IDs:
    % compoundIDs - Cell array containing compound IDs as strings
    % exchange_reactions - Cell array containing formatted exchange reactions
    % reaction_names - Cell array containing names for each exchange reaction

    % Initialize the output cell arrays
    exchange_reactions = cell(size(compoundIDs));
    reaction_names = cell(size(compoundIDs));

    % Loop through each compound ID in the input cell array
    for i = 1:length(compoundIDs)
        % Get the current compound ID
        compoundID = compoundIDs{i};

        % Create the exchange reaction string
        reaction = sprintf('%s_e -> ', compoundID);
        
        % Create the reaction name string
        reaction_name = sprintf('EXC_BOTH_%s_e', compoundID);

        % Store the reaction and its name in the output cell arrays
        exchange_reactions{i} = reaction;
        reaction_names{i} = reaction_name;
    end
end
