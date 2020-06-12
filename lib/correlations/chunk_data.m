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
    
    num_shots = size(counts,2);
    
    total_counts=sum(num_counts);
    
    if nargin<4
        norm_chunk_size=round(mean(num_counts)*sqrt(norm_samp_factor));
        num_chunks=ceil(total_counts/norm_chunk_size);
    else
        norm_chunk_size=floor(total_counts/num_chunks);
    end
    
    if num_chunks<max(num_counts)
        num_chunks = max(num_counts);
    end
    
    counts_chunked=cell(1,num_chunks);
    for ii=1:num_shots
        this_counts = counts{ii};
        chunk_dist = randperm(num_chunks);
        for jj = 1:num_counts(ii)
            box = chunk_dist(jj);
            counts_chunked{box}=[counts_chunked{box};this_counts(jj,:)];
        end
    end
    if ~isnan(sort_dir)
        for ii = 1:num_chunks
            tmp_data = counts_chunked{ii};
            [~,order]=sort(tmp_data(:,sort_dir));
            tmp_data=tmp_data(order,:);
            counts_chunked{ii} = tmp_data;
        end
    end
    if size(vertcat(counts_chunked{:}),1)~=total_counts
        warning('lost counts')
    end
end

