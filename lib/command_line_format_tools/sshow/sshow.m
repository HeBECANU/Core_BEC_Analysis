function sshow(s)

% Displays the tree structure of a struct variable and the types of its
% terminal nodes (leaves). 
% Input:  s      - struct
% Ourput: nothing - for now, will eventually return a tree
% Improvements: Show variable dimensions as in MATLAB's default display
        % Return a representation of the branch structure
        % 
sshow_core(s,0);

end

function sshow_core(s,depth)

if depth==0
    sname = inputname(1);
    fprintf('Structure %s:\n',sname)
end
    if isstruct(s)
        fnames = fields(s);
        nfields = numel(fnames);
        buffer = '-';
        for k=1:depth
            buffer = [buffer,'--'];
        end
        for fieldnum = 1:nfields
           fprintf([buffer,' ',fnames{fieldnum},'\n']) 
           sshow_core(s.(fnames{fieldnum}),depth+1)
        end
    else
        fprintf('\b - %s\n',class(s))
    end
end