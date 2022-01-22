function time_posix=get_time_now_posix
    nowdt=datetime('now','TimeZone','local','Format', 'yyyy-MM-dd HH:mm:ss.SSSxxxxx');
    time_posix=posixtime(nowdt);
end