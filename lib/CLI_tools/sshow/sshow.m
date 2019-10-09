function sshow(varargin)

% Displays the tree structure of a struct variable and the types of its
% terminal nodes (leaves). 
% Input:  s      - struct
% Ourput: nothing - for now, may eventually return a tree
% Improvements: 
        % Show variable dimensions as in MATLAB's default display
        % Return a representation of the branch structure
        % Align the name/type columns more cleanly
        % Output field contents (if singleton) or first few (if
            % array/cell?)
        % Display warning if non-struct passed as input
        
        

if nargin == 1
   depth =0;
   sname = inputname(1);
   fprintf('Contents of structure %s:\n',sname)
   sshow(varargin{1},0);
else
    depth = varargin{2};
    s = varargin{1};
    if isstruct(s)
        fnames = fields(s);
        nfields = numel(fnames);
        buffer = ' >';
        for k=1:depth
            buffer = ['   ',buffer];
        end
        for fieldnum = 1:nfields
           fprintf([buffer,' ',fnames{fieldnum},'\n']) 
           sshow(s.(fnames{fieldnum}),depth+1)
        end
    else
        fprintf('\b - %s\n',class(s))
    end
end


end