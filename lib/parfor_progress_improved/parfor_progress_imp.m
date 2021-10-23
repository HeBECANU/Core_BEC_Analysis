function percent = parfor_progress_imp(N)
%PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
%   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
%   your working directory, and then keeping track of the parfor loop's
%   progress within that file. This workaround is necessary because parfor
%   workers cannot communicate with one another so there is no simple way
%   to know which iterations have finished and which haven't.
%
%   PARFOR_PROGRESS(N) initializes the progress monitor for a set of N
%   upcoming calculations.
%
%   PARFOR_PROGRESS updates the progress inside your parfor loop and
%   displays an updated progress bar.
%
%   PARFOR_PROGRESS(0) deletes parfor_progress.txt and finalizes progress
%   bar.
%
%   To suppress output from any of these functions, just ask for a return
%   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
%   returns the percentage of completion.
%
%   Example:
%
%      N = 100;
%      parfor_progress(N);
%      parfor i=1:N
%         pause(rand); % Replace with real code
%         parfor_progress;
%      end
%      parfor_progress(0);
%
%   See also PARFOR.

% By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/

% modified to use a binary file using the comment by philip
% https://au.mathworks.com/matlabcentral/fileexchange/32101-progress-monitor-progress-bar-that-works-with-parfor
%"Nice contribution. However, as Timo already found out, while fitting a huge amount of data 
%I noticed that the progress monitor gets extremely slow for very long loops. In the beginning
%it is not so bad but with every loop iteration the call to 'fscanf' gets slower and slower.
%I fixed this by reading & writing to a simple binary file (i am writing the expected loop length 
%to the first 4 bytes and the current progress to the next 4 bytes - if you think that uint32 is 
%not enough, feel free to use uint64). In addition, the filename is now stored in the folder 
%determined by matlab's 'tempname'-function (less filesystem clutter if the delete fails & e.g. tempdir 
%on an ssd is much faster than working dir somewhere on the network)."
% this makes this improved version about 30 times faster than stock

% bugs and impovements
% at the moment all calls to parfor_progress_imp go to the same file
% this could be improved by passing an argument which is some random number
% that could identify a particular process  

% could also be greatly improved by recording the time of the last N
% updates and extrapolating to the dinish time  

narginchk(0,1)

if nargin < 1
    N = -1;
end

percent = 0;
w = 50; % Width of progress bar


if N > 0
    f = fopen(fullfile(tempdir,'parfor_progress.bin'), 'w');
    if f<0
        error('Do you have write permissions for %s?', pwd);
    end
    fwrite(f,tic,'uint64'); %current system time in 1e-7 seconds
    fwrite(f,N,'uint64');
    fwrite(f,0,'uint64');
    fclose(f);

    if nargout == 0
        print_cli_data(0,0,0,w);
    end
elseif N == 0
    fname = fullfile(tempdir,'parfor_progress.bin');
    if (exist(fname, 'file')~=2)
        warning('parfor_progress.bin not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.bin.');
    else
        f = fopen(fname, 'r');
        A = fread(f,3,'uint64');
        time_start_pointer=uint64(A(1));
        fclose(f);
        delete(fname);
        if nargout == 0
            run_time_s=toc(time_start_pointer);
            print_cli_data(1,run_time_s,0,w);
            fprintf('done.\n')
        end
    end
else
    fname = fullfile(tempdir,'parfor_progress.bin');
    if (exist(fname, 'file')~=2)
        error('parfor_progress.bin not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.bin.');
    end
    
    f = fopen(fname, 'r+');
    A = fread(f,3,'uint64');
    time_start_pointer=uint64(A(1));
    todo = A(2);
    progress = A(3) + 1;
    % move to the progress position in the file 8*2=15 bytes
    fseek(f, 16, 'bof');
    fwrite(f,progress,'uint32');
    fclose(f);
    progress_frac=progress/todo;
    
    if nargout == 0
        %delete_chars=repmat('\b', 1, (w+10)); %char(8)
        fprintf(repmat('\b',[1,w+12+36])); %9 wide for the percent 36 for the progress
        run_time_s=toc(time_start_pointer);
        %est_time_full=runtime/progress_frac
        %est_time_from_now=est_time_full(1-progress_frac)
        %est_time_from_now=(runtime/progress_frac)(1-progress_frac)
        eta_s=(run_time_s/progress_frac)*(1-progress_frac);
        print_cli_data(progress_frac,run_time_s,eta_s,w)
        
    end
end

end


function print_cli_data(progress_frac,run_time_s,eta_s,w)
n_bars=round(progress_frac*w);
n_bars=min(n_bars,w);
progress_bar=[repmat('=', 1, n_bars) , '>', repmat(' ', 1, w - n_bars)];
fprintf('%5.1f%%[%s]\n',progress_frac*100,progress_bar);
run_time_hms = fix(mod(run_time_s, [0, 3600, 60]) ./ [3600, 60, 1]);
eta_hms = fix(mod(eta_s, [0, 3600, 60]) ./ [3600, 60, 1]);
fprintf('       runtime %02u:%02u:%02u eta %02u:%02u:%02u.\n',run_time_hms,eta_hms)
end


% % function percent = parfor_progress(N)
% % %PARFOR_PROGRESS Progress monitor (progress bar) that works with parfor.
% % %   PARFOR_PROGRESS works by creating a file called parfor_progress.txt in
% % %   your working directory, and then keeping track of the parfor loop's
% % %   progress within that file. This workaround is necessary because parfor
% % %   workers cannot communicate with one another so there is no simple way
% % %   to know which iterations have finished and which haven't.
% % %
% % %   PARFOR_PROGRESS(N) initializes the progress monitor for a set of N
% % %   upcoming calculations.
% % %
% % %   PARFOR_PROGRESS updates the progress inside your parfor loop and
% % %   displays an updated progress bar.
% % %
% % %   PARFOR_PROGRESS(0) deletes parfor_progress.txt and finalizes progress
% % %   bar.
% % %
% % %   To suppress output from any of these functions, just ask for a return
% % %   variable from the function calls, like PERCENT = PARFOR_PROGRESS which
% % %   returns the percentage of completion.
% % %
% % %   Example:
% % %
% % %      N = 100;
% % %      parfor_progress(N);
% % %      parfor i=1:N
% % %         pause(rand); % Replace with real code
% % %         parfor_progress;
% % %      end
% % %      parfor_progress(0);
% % %
% % %   See also PARFOR.
% % 
% % % By Jeremy Scheff - jdscheff@gmail.com - http://www.jeremyscheff.com/
% % 
% % error(nargchk(0, 1, nargin, 'struct'));
% % 
% % if nargin < 1 
% % N = -1; 
% % end
% % percent = 0; 
% % w = 50; % Width of progress bar
% % if N > 0 
% % f = fopen(fullfile(tempdir,'parfor_progress.bin'), 'w'); 
% % if f<0 
% % error('Do you have write permissions for %s?', pwd); 
% % end 
% % fwrite(f,N,'uint32'); 
% % fwrite(f,0,'uint32'); 
% % fclose(f); 
% % 
% % if nargout == 0 
% % disp([' 0%[>', repmat(' ', 1, w), ']']); 
% % end 
% % elseif N == 0 
% % delete(fullfile(tempdir,'parfor_progress.bin')); 
% % percent = 100; 
% % 
% % if nargout == 0 
% % disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']); 
% % end 
% % else 
% % fname = fullfile(tempdir,'parfor_progress.bin'); 
% % if ~exist(fname, 'file') 
% % error('parfor_progress.bin not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.bin.'); 
% % end
% % f = fopen(fname, 'r+'); 
% % A = fread(f,2,'uint32'); 
% % todo = A(1); 
% % progress = A(2) + 1; 
% % fseek(f, 4, 'bof'); 
% % fwrite(f,progress,'uint32'); 
% % fclose(f); 
% % 
% % percent = progress/todo * 100; 
% % 
% % if nargout == 0 
% % perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage 
% % disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']); 
% % end 
% % end
