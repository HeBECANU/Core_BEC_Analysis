function msg_out = cli_header(varargin)
    % A function that makes nice sub headers in the console
    % Works OK with numeric substitution, but needs a length checking
    % function.
    
    % Examples:
%         header('Header examples')
%         header(2,'Some options:')
%         header(3,'%u %u,...',[1,2])
    
    % Improvements; 
%             switch to varargin
%             Tests
    args = varargin;
    if iscell(varargin{1})
        fprintf("cli_header accepts varargin, cell args will be removed next version.\n")
        args = varargin{1};
    end
    v = 1;
    if numel(args) == 1 %text only
        msg = args{1};
        lvl = 0;%If no header level is passed, assumed to be top-level
        
    else
        if isfloat(args{1}) %Passing header level
            lvl = args{1};
            msg = args{2};
            if numel(args)>2
                v = 3;
            end
        else % No header level passed, but other values present
            lvl = 0;
            msg = args{1};
            v = 2;
        end
    end
    msg_out = sprintf(msg,args{v:end});
%     blank = '------------------------------------------------------------';
    blank = '                                                            ';
    marker = blank(1:2*lvl+1);
    fprintf([marker,' - ',msg_out,'\n'])

end