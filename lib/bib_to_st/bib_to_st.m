function cell_of_st=bib_to_st(fname_or_dir,varargin)
% if pass dir then will import all the .bib files
% can cat them together
% otherwise will just return the single file bib entries

p = inputParser;
is_c_logical=@(in) isequal(in,true) || isequal(in,false); %can x be cast as a logical
addOptional(p,'auth_st',true,is_c_logical);
parse(p,varargin{:});
parsed_input= p.Results;



if isfolder(fname_or_dir)
    files=dir(fullfile(fname_or_dir,'*.bib'));
    fnames={files.name};
    fnames(strcmp(fnames,'.'))=[];
    fnames(strcmp(fnames,'..'))=[];
    cell_of_st=cell(1,numel(fnames));
    for ii=1:numel(fnames)
        file_path=fullfile(fname_or_dir,fnames{ii});
        tmp=single_file_to_cell_st(file_path);
        cell_of_st{ii}=tmp;
    end
    cell_of_st=cat(2,cell_of_st{:});
else
    cell_of_st=single_file_to_cell_st(fname_or_dir);
end

% change the author format to a structure
if parsed_input.auth_st
    for ii=1:numel(cell_of_st)
        this_entry=cell_of_st{ii};
        if isfield(this_entry,'author')
            this_entry.author=parse_auth_list_to_st(this_entry.author);
        end
        cell_of_st{ii}=this_entry;
    end
end

% change the year to an integer
for ii=1:numel(cell_of_st)
    this_entry=cell_of_st{ii};
    if isfield(this_entry,'year')
        this_entry.year=str2double(this_entry.year);
    end
    cell_of_st{ii}=this_entry;
end


end


function st_out=parse_auth_list_to_st(auth_str)
    auths=strsplit(auth_str,'and');
    st_out=cell(1,numel(auths));
    for ii=1:numel(auths)
        this_auth=auths{ii};
        this_auth=strip(this_auth,'both',' ');
        last_first=strsplit(this_auth,',');
        if numel(last_first)>2
            error('expected only one comma per and in auth list')
        elseif numel(last_first)==2
            this_st=[];
            this_st.first=strip(last_first{2},'both',' ');
            this_st.last=strip(last_first{1},'both',' ');
        else
            % if we didnt have a comma then we just use the last space
            first_last=strsplit(this_auth,' ');
            this_st.first=strip(strjoin(first_last(1:end-1)),'both',' ');
            this_st.last=strip(first_last{end},'both',' ');
        end
        st_out{ii}=this_st;
    end

end

function cell_of_bib_st=single_file_to_cell_st(file_path)
    if isfile(file_path)
        rawstr = fileread(file_path);

        % remove comments
        % from https://tex.stackexchange.com/questions/451836/deleting-comments-from-tex-file
        entries_proc=regexprep(rawstr,'(?<!\\)%.*',' ');
        % remove newlines
        entries_proc=regexprep(entries_proc,'[\n\r]+',' ');
        entries_proc=strsplit(entries_proc,'@');
        % entries_proc has each bibtex entry as a element
        if numel(entries_proc)<2
            warning('no bibtex entries')
        else
            entries_proc=entries_proc(2:end); % remove the header (or empty)
            cell_of_bib_st=[];
            for jj=1:numel(entries_proc)
                entry_st=[];
                entry_str=entries_proc{jj};
                [entry_type,entry_str]=strtok(entry_str,'{');
                % scan from right to left untill the first } there is no direction argument so we flip
                % this is to remove random text between the bibtex entries (which there shouldnt be)
                [~,entry_str]=strtok(fliplr(entry_str),'}');  
                entry_str=fliplr(entry_str);  
                if numel(strfind(entry_str,'{'))~=numel(strfind(entry_str,'}'))
                    error('mismatching squiggly brackets')
                end
                entry_str=strip(entry_str,'both',' ');
                entry_str=strip(entry_str,'left','{');
                entry_str=strip(entry_str,'right','}');
                entry_str=strsplit(entry_str,',');
                entry_st.bib_key=entry_str{1};
                entry_st.type=entry_type;
                entry_str(1)=[];
                % merge the entries back together unitll num('{') = num('}') in entry_str
                kk=1;
                while kk<numel(entry_str)
                    if numel(strfind(entry_str{kk},'{'))~=numel(strfind(entry_str{kk},'}'))
                        entry_str{kk}=strjoin(entry_str(kk:kk+1),',');
                        entry_str(kk+1)=[];
                        kk=kk-1;
                    end
                    kk=kk+1;
                end
                % alt algo (this works fine but just a slightly less brute force approach)
                % find the indicies of comma_idx=',' and eq_idx='='
                % for each eq_idx find the nearest lower comma_idx, the chars between is the field name
                % the field value is given by going to the next equals going back to the prev comma
                % which is at the field we are initerested in
                
                % now we have cells with ' year = {} ' 
                for kk=1:numel(entry_str)
                    entry_this=strsplit(entry_str{kk},'=');
                    % there could be an equals in the text?
                    if numel(entry_this)>2
                        entry_this=strjoin(entry_this(2:end),'=');
                    end
                    entry_key=strip(entry_this{1});
                    entry_val=strip(entry_this{2});
                    entry_val=strip(entry_val,'left','{');
                    entry_val=strip(entry_val,'right','}');
                    entry_st.(entry_key)=entry_val;
                end
                cell_of_bib_st{jj}=entry_st;
            end
        end
    else
        warning('file %s \n does not exist',file_path)
    end


end