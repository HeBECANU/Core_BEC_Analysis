function [data]=import_data(import_opts)
%legacy code for people that dont have the phased array toolbox

import=1;
if exist('importsave.mat','file')==2 && ~import_opts.force_forc && ~import_opts.force_reimport
    load('importsave.mat','data','import_opts_old')
    if isequal(import_opts_old,import_opts)
        import=0; 
        fprintf('loaded old data\n')
    else
        clear('data')
    end
end

%import_opts.dld_xy_rot=0.61;

if import
    fprintf('%04i\n',0)
    data.txy={};
    for ii=1:size(import_opts.shot_num,2)
        if ~(exist([import_opts.dir,import_opts.file_name,num2str(import_opts.shot_num(ii)),'.txt'],'file')==2)
                    data.txy{ii}=[];
                    fprintf('\n no_file %04i \n %04i\n',ii,ii)
        else
            %if the txy_forc does not exist or the use_txy flag is low (re) make it
            if ~(exist([import_opts.dir,import_opts.file_name,'_txy_forc',num2str(import_opts.shot_num(ii)),'.txt'],'file')==2) || import_opts.force_forc
                dld_raw_to_txy([import_opts.dir,import_opts.file_name],import_opts.shot_num(ii),import_opts.shot_num(ii));
            end
            txydata=txy_importer([import_opts.dir,import_opts.file_name],num2str(import_opts.shot_num(ii)));
            txydata=masktxy(txydata,import_opts.txylim); %mask for counts in the window txylim    
            
            txy_rot=zeros(size(txydata));
            sin_theta = sin(import_opts.dld_xy_rot);
            cos_theta = cos(import_opts.dld_xy_rot);
            txy_rot(:,1) = txydata(:,1);
            txy_rot(:,2) = txydata(:,2)*cos_theta...
                - txydata(:,3)*sin_theta;
            txy_rot(:,3) = txydata(:,2)*sin_theta...
                + txydata(:,3)*cos_theta;
            data.txy{ii}=txy_rot;
            data.total_num{ii}=size(txy_rot,1);
        end %file exists condition
    fprintf('\b\b\b\b%04i',ii)
    end
    import_opts_old=import_opts;
    save('importsave.mat','data','import_opts_old')
    fprintf('\nsaved\n')
end
end





