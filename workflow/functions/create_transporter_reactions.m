function [transporter_reactions, reaction_names] = create_transporter_reactions(compoundIDs, reactionType)
    % Generates transporter reactions along with their names based on the specified type:
    % reactionType = 1: From extracellular (_e) to periplasmic (_p)
    % reactionType = 2: From periplasmic (_p) to cytoplasmic (_c)
    % reactionType = 3: From extracellular (_e) to cytoplasmic (_c)
    % compoundIDs - Cell array containing compound IDs as strings
    % transporter_reactions - Cell array containing formatted transporter reactions
    % reaction_names - Cell array containing names for each transporter reaction

    % Initialize the output cell arrays
    transporter_reactions = cell(size(compoundIDs));
    reaction_names = cell(size(compoundIDs));

    % Loop through each compound ID in the input cell array
    for i = 1:length(compoundIDs)
        % Get the current compound ID
        compoundID = compoundIDs{i};

        % Determine the reaction based on the type and generate corresponding names
        if reactionType == 1
            % Create the reaction string (extracellular to periplasmic)
            reaction = sprintf('%s_e <=> %s_p', compoundID, compoundID);
            reaction_name = sprintf('T_e_to_p_%s_e', compoundID);
        elseif reactionType == 2
            % Create the reaction string (periplasmic to cytoplasmic)
            reaction = sprintf('%s_p <=> %s_c', compoundID, compoundID);
            reaction_name = sprintf('T_p_to_c_%s_p', compoundID);
        elseif reactionType == 3
            % Create the reaction string (extracellular to cytoplasmic)
            reaction = sprintf('%s_e <=> %s_c', compoundID, compoundID);
            reaction_name = sprintf('T_e_to_c_%s_e', compoundID);
        else
            error('Invalid reaction type. Use 1 for _e to _p and 2 for _p to _c.');
        end

        % Store the reaction and its name in the output cell arrays
        transporter_reactions{i} = reaction;
        reaction_names{i} = reaction_name;
    end
end
