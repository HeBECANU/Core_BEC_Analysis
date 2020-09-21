function state_out = term2vec(state)
    % accepts n^SL_J_mJ string and returns [n l s j mj] vector - note the
    % ordering
    % Also turns a vector into term strings if needed?
    % outputs in form 
    term_symbols = 'SPDFGHIJKL';
    if isstring(state) || ischar(state)
        state_out = [str2num(state(1)),strfind(term_symbols,state(4))-1,0.5*(str2num(state(3))-1),str2num(state(6)),str2num(state(8))];
    else
        % it's a vector, so make it into a string
        state_out = '';

        if len(state) == 3
            state_out = sprintf('%u_%u%s',state(1),2*state(2)+1,state(3));
        elseif len(state) == 4    
            state_out = sprintf('%u_%u%s_%u',state(1),2*state(2)+1,state(3),state(4));
        elseif len(state) == 5
            state_out = sprintf('%u_%u%s_%u_%u',state(1),2*state(2)+1,state(3),state(4),state(5));
        end

    end
end