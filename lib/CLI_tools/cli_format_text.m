function cli_format_text(in_string,type,level,varargin)
% format some text either as centered text or as an element of a
% hierarchical outline
% Example
% cli_format_text('456','h',1)
% cli_format_text('456','h',2)
% cli_format_text('456','h',3)
% returns
% ----  456
% --------  456
% ------------  456

% cli_format_text('4567','c',3)
% cli_format_text('','c',3)
% ====================================  4567  ====================================
% ================================================================================


% TODO
% standard function format
% usage examples
% syntax and goals


cli_width=80;
space_padding=2;
sep_types={'-','+','=','_','*'};
hierarchical_spacing=4;
hierarchical_sep='-';

if isstring(in_string)
    in_string=char(in_string);
end

if ~ischar(in_string)
    error('first agument must be a string or char vector')
end


format_types={'c','center','centre','h','hierarchical'};


if sum(strcmp(format_types,type))==0
    error('format must be supported type %s',sprintf('%s ,',format_types{:}))
end


if nargin==2 || isempty(level)
    level=1;
end

if mod(level,1)~=0
    error('level must be integer')
end


if sum(strcmp(type,{'c','center','centre'}))
   
    if level>numel(sep_types)
        error('level not valid')
    end
    sep_choice=sep_types{level};
    
    %do not pad the string with spaces if it is empty, for making will
    %widht headers
    if ~isempty(in_string) || ~strcmp(in_string,'')
        padded_str=cat(2,repmat(' ',1,space_padding),in_string,repmat(' ',1,space_padding));
    else
        padded_str=in_string;
    end
    
    remaining_width=cli_width-numel(padded_str);
    if remaining_width < 2
        error('cant make it this wide')
    end
    
    if mod(remaining_width,2)==0
        add_one_extra_to_left=0;
    else
        add_one_extra_to_left=1;
    end
    each_side_padd=floor(remaining_width/2);
    cent_str=cat(2,repmat(sep_choice,1,each_side_padd+add_one_extra_to_left),...
        padded_str,...
        repmat(sep_choice,1,each_side_padd));
    out_str=cent_str;
    
elseif sum(strcmp(type,{'h','hierarchical'}))
    
	padded_str=cat(2,repmat(' ',1,space_padding),in_string);
     
    h_str=cat(2,repmat(hierarchical_sep,1,level*hierarchical_spacing),...
        padded_str);
    out_str=h_str;
    
end
  
fprintf(cat(2,out_str,'\n'))



end