function cli_header(varargin)
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
    argument_index=1;
    msg_level = 0;%If no header level is passed, assumed to be top-level
    if iscell(varargin{1})
        fprintf("cli_header accepts varargin, cell args will be removed next version.\n")
    end
    if nargin == 1 %text only
        msg = varargin{1};
        msg_level = 0;%If no header level is passed, assumed to be top-level
        argument_index=2;
    else %passed level, or args, or both
        if isnumeric(varargin{1}) % Header level passed as argument
            msg_level = varargin{1};
            msg = varargin{2};
            argument_index = 3;
        else % No header level passed, but other values present
            msg_level = 0;
            msg = varargin{1};
            argument_index = 2;
        end
    end
    
    msg_out = sprintf(msg,varargin{argument_index:end});
%     blank = '------------------------------------------------------------';
    blank = '                                                            ';
    marker = blank(1:2*msg_level+1);
    fprintf([marker,' - ',msg_out,'\n'])

end