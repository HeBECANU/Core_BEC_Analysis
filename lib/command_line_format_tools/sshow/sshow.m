function sshow(s)

% Displays the tree structure of a struct variable and the types of its
% terminal nodes (leaves). 
% Input:  s      - struct
% Ourput: nothing - for now, will eventually return a tree
% Improvements: Show variable dimensions as in MATLAB's default display
        % Return a representation of the branch structure
        % 
sname = inputname(1);
fprintf('Structure %s:\n',sname)
sshow_core(s,0);

end

function sshow_core(s,depth) 
    if isstruct(s)
        fnames = fields(s);
        nfields = numel(fnames);
        if depth >0
            buffer = ' |- ';
        else
            buffer = '--';
        end
        for k=1:depth
            buffer = ['  ',buffer];
        end
        for fieldnum = 1:nfields
           this_name = fnames{fieldnum};
           txt = '                ';
           txt(1:length(this_name)) = this_name;
           
           datlabel = [' - ',class(s.(fnames{fieldnum}))];
           if isstruct(s.(fnames{fieldnum}))
               datlabel = '';
           end
           
%            fprintf([buffer,' ',datlabel,' - ',fnames{fieldnum},'\n']) 
           fprintf([buffer,' ',txt,datlabel,'\n']) 
           sshow_core(s.(fnames{fieldnum}),depth+1)
        end
    else
%         fprintf('\b - %s\n',class(s))
    end
end