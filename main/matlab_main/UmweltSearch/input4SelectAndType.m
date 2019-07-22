function inputs = input4SelectAndType(selection, t_back, inputDataFull, type)

	if islogical(selection)
		selection = find(selection);
	end
	inputs = inputDataFull((end-t_back+1):end, selection, 1, :);
	if strcmp(type, 'image')
		% do nothing... job already done above
    elseif strcmp(type, 'vector')          
		noInps = length(selection);
		inputs = reshape(inputs, t_back*noInps, size(inputDataFull, 4));
	else 
		error('Undefined type.')
	end



end