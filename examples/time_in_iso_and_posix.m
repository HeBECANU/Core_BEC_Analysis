% get time in iso and posix format
nowdt=datetime('now','TimeZone','local','Format', 'yyyy-MM-dd HH:mm:ss.SSSxxxxx');
obj.log_time_iso=strrep(char(nowdt),' ','T');
obj.log_time_posix=posixtime(nowdt);
            
     