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
    fwrite(f,N,'uint32');
    fwrite(f,0,'uint32');
    fclose(f);

    if nargout == 0
        disp([' 0%[>', repmat(' ', 1, w), ']']);
    end
elseif N == 0
    delete(fullfile(tempdir,'parfor_progress.bin'));
    percent = 100;

    if nargout == 0
        disp([repmat(char(8), 1, (w+9)), char(10), '100%[', repmat('=', 1, w+1), ']']);
    end
else
    fname = fullfile(tempdir,'parfor_progress.bin');
    if ~exist(fname, 'file')
        error('parfor_progress.bin not found. Run PARFOR_PROGRESS(N) before PARFOR_PROGRESS to initialize parfor_progress.bin.');
    end

    f = fopen(fname, 'r+');
    A = fread(f,2,'uint32');
    todo = A(1);
    progress = A(2) + 1;
    fseek(f, 4, 'bof');
    fwrite(f,progress,'uint32');
    fclose(f);

    percent = progress/todo * 100;

    if nargout == 0
        perc = sprintf('%3.0f%%', percent); % 4 characters wide, percentage
        disp([repmat(char(8), 1, (w+9)), char(10), perc, '[', repmat('=', 1, round(percent*w/100)), '>', repmat(' ', 1, w - round(percent*w/100)), ']']);
    end
end

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
