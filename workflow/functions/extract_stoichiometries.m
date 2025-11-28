function stoichiometriesList = extract_stoichiometries(reactionsList)
    stoichiometriesList = cell(1, numel(reactionsList));

    for r = 1:numel(reactionsList)
        reaction = reactionsList{r};
        
        % Split reaction into substrates and products
        reaction_parts = split(reaction, ' <=> ');
        if numel(reaction_parts) < 2
            reaction_parts = split(reaction, ' -> ');
        end

        substrates = split(reaction_parts{1}, ' + ');
        products = split(reaction_parts{2}, ' + ');

        stoichiometries = [];
        
        % Parse substrates
        for i = 1:numel(substrates)
            substrate_parts = split(substrates{i}, ' ');
            if numel(substrate_parts) > 1 && ~isnan(str2double(substrate_parts{1}))
                stoichiometry = -str2double(substrate_parts{1});
            else
                stoichiometry = -1;
            end
            stoichiometries = [stoichiometries, stoichiometry];
        end
        
        % Parse products
        for i = 1:numel(products)
            product_parts = split(products{i}, ' ');
            if numel(product_parts) > 1 && ~isnan(str2double(product_parts{1}))
                stoichiometry = str2double(product_parts{1});
            else
                stoichiometry = 1;
            end
            stoichiometries = [stoichiometries, stoichiometry];
        end

        stoichiometriesList{r} = stoichiometries;

    end
    stoichiometriesList=transpose(stoichiometriesList);
end
