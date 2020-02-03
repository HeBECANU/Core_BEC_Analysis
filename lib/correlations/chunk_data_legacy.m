function counts_chunked=chunk_data(counts,norm_samp_factor,sort_dir,num_chunks)
    %data chunking function for normalization
    %see calc_any_g2_type for usage
    
    %input
    %norm_samp_factor-should be about the correlation amp for equal noise contibution from in-shot & between shot
    %data.txy
    %data.data.num_counts
    %to do
    %   -documentation
    
    num_counts=cellfun(@(x)size(x,1),counts);
    
    
    
    total_counts=sum(num_counts);
    
    if nargin<4
        norm_chunk_size=round(mean(num_counts)*sqrt(norm_samp_factor));
        num_chunks=ceil(total_counts/norm_chunk_size);
    else
        norm_chunk_size=floor(total_counts/num_chunks);
    end
    
    all_counts=vertcat(counts{:});
    all_counts=all_counts(randperm(total_counts),:);
    iimax=num_chunks;
    counts_chunked=cell(1,iimax);
    for ii=1:iimax
%         if nargin>3
%             norm_chunk_size=round(randn(1)+mean(num_counts));
%         end
        min_idx=(ii-1)*norm_chunk_size+1;
        if ii==iimax
            max_idx=total_counts;
        else
            max_idx=ii*norm_chunk_size;
        end
        tmp_data=all_counts(min_idx:max_idx,:);
        if ~isnan(sort_dir)
            [~,order]=sort(tmp_data(:,sort_dir));
            tmp_data=tmp_data(order,:);
        end
        counts_chunked{ii}=tmp_data;
    end
    if size(vertcat(counts_chunked{:}),1)~=total_counts
        warning('lost counts')
    end
end

