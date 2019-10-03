function msg_out = cli_header(varargin)
    % A function that makes nice sub headers in the console
    % Works OK with numeric substitution, but needs a length checking
    % function.
    % syntax
    % arg1= header level
    
    
    % Examples:
%         header('Header examples')
%         header(2,'Some options:')
%         header(3,'%u %u,...',[1,2])
    
    % Improvements; 
%             switch to varargin
%             Tests
%             standard documentation
%             explain syntax
%             commenting

% Done
% removed single letter var names
% simplified

%recomendations
% remove multi arguments in favor of cli_header(string,level, opp args (inc type))

    args = varargin;
    if iscell(varargin{1})
        fprintf("cli_header accepts varargin, cell args will be removed next version.\n")
        args = varargin{1};
    end
    
    if numel(args) == 1 %text only
        msg = args{1};
        msg_level = 0;%If no header level is passed, assumed to be top-level
        argument_index=2;
    else
        if isnumeric(args{1}) %Passing header level
            msg_level = args{1};
            msg = args{2};
            argument_index=3;
        else % No header level passed, but other values present
            msg_level = 0;
            msg = args{1};
            argument_index = 2;
        end
    end
    
    msg_out = sprintf(msg,args{argument_index:end});
%     blank = '------------------------------------------------------------';
    blank = '                                                            ';
    marker = blank(1:2*msg_level+1);
    fprintf([marker,' - ',msg_out,'\n'])

end