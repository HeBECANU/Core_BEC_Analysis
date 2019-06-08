function cli_header(varargin)
    % A function that makes nice sub headers in the console
    % Works OK with numeric substitution, but needs a length checking
    % function.
    
    % Examples:
%         header('Header examples')
%         header({2,'Some options:'})
%         header({3,'%u %u,...',[1,2]})
    
    % Improvements; 
%             switch to varargin
%             Tests
    vals = [];
    if nargin == 1 %text only
        msg = varargin{1};
        lvl = 0;%If no header level is passed, assumed to be top-level
    else
        if isfloat(varargin{1}) %Passing header level
            lvl = varargin{1};
            msg = varargin{2};
            vals = varargin{3:end};
        else % No header level passed, but other values present
            lvl = 0;
            msg = varargin{1};
            vals = varargin{2:end};
        end
    end
    msg_out = sprintf(msg,vals);
    blank = '------------------------------------------------------------';
    marker = blank(1:10-3*lvl);
    fprintf([marker,'   ',msg_out,'\n'])

end