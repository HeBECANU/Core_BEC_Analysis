% this is a very clean way of showing progress in a loop
% it works by printing delete charaters to remove what has been written previously
% display is easiest to read if what you print is zero padded to give fixed width
% eg %04u, 01.1e
% bryce henson 2020-02-13

fprintf('\nrunning simulation\n')
fprintf('   ')
output_chars=fprintf('%04u tol %01.1e',0,0);



for ii=1:20
    out_tol=pi*ii;
    fprintf(repmat('\b',[1,output_chars]))
    output_chars=fprintf('%04u tol %01.1e',ii,out_tol);
    pause(0.1)
end
fprintf('\n')


%% changing line length

fprintf('\nrunning simulation\n')
fprintf('   ')
output_chars=fprintf('%4u tol %1.1e',0,0);



for ii=1:20
    out_tol=pi*ii*100;
    fprintf(repmat('\b',[1,output_chars]))
    output_chars=fprintf('abc %4u tol %.2f',ii,out_tol);
    pause(0.1)
end
fprintf('\n')

%% timing
fprintf('\nrunning simulation\n')
fprintf('   ')
output_chars=fprintf('%4u tol %1.1e',0,0);

tic;

for ii=1:20
    out_tol=pi*ii*100;
    fprintf(repmat('\b',[1,output_chars]))
    output_chars=fprintf('abc %4u tol %.2f',ii,out_tol);
end
fprintf('\n')


time_to_write=toc;
fprintf('time per update %.2f ms \n',1e3*time_to_write/ii)
