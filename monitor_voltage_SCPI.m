% user params
samp_period = 10;%sample period in seconds
max_time  = 60*15;%max run time in seconeds


%find scope
t=tcpclient("150.203.181.15",5025);

%write out its name
writeread(t,"*IDN?")

%http://150.203.181.15/

%set it up so it monitors the max and min + rms of chanel 1
%:MEASure:VRMS? CYCL,DC,CHAN1
%:MEASure:VMIN? CHAN1
%:MEASure:VMAX? CHAN1
%:MEASure:VRMS? DISP,DC,CHAN1

% meas_list = {':MEASure:VRMS? CYCL,DC,CHAN1',':MEASure:VRMS? DISP,DC,CHAN1',':MEASure:VMIN? CHAN1',':MEASure:VMAX? CHAN1'};

meas_list = {':MEASure:VRMS? DISP,DC,CHAN1'};

data = cell(size(meas_list));
strt_time = posixtime(datetime);
run_time = posixtime(datetime)-strt_time;
%every given interval measure the voltage
while run_time<max_time
    
    if run_time>samp_period
        pause(samp_period)
    end
    
    for ii = 1:length(meas_list)
       meas_c =  str2double(writeread(t,meas_list{ii}));
       meas_time = posixtime(datetime);
       
       data{ii} = [data{ii};[meas_time,meas_c]];
    end
    run_time = posixtime(datetime)-strt_time;
    
end