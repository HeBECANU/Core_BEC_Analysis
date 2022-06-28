clear all
filedir = 'E:\test';
files = dir(fullfile(filedir, '*.csv'));
he3_peak_difftime = {};
he4_peak_difftime = {};
he3_loc_difftime = {};
he4_loc_difftime = {};
stfig('SFP');
clf
for i=1:length(files)
    filename = fullfile(files(i).folder, files(i).name);
    M = readmatrix(filename);
    PD_voltage=M(:,1);
    ramp_voltage=M(:,2);
    stfig('SFP');
    plot(PD_voltage);
    hold on
    plot(ramp_voltage)
    xlabel('sample')
    ylabel('voltage')
    
    %find peak locations
    
    %find peak heights
    
    [pks,locs,w,p] = findpeaks(PD_voltage,'MinPeakHeight',5);
    he3_peak_difftime{i} = pks(1:2:end-1);
    he4_peak_difftime{i} = pks(2:2:end);
    he3_loc_difftime{i} = locs(1:2:end-1);
    he4_loc_difftime{i} = locs(2:2:end);

    he3_peak_avg(i) = mean(pks(1:2:end-1));
    he4_peak_avg(i) = mean(pks(2:2:end));
    %pks(1:2:end-1)./pks(2:2:end)
end
%%
stfig('average peak height')
plot(he3_peak_avg.')
hold on
plot(he4_peak_avg.')
xlabel('file')
ylabel('voltage')