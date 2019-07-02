function sshow_test()
s_test = [];
s_test.time = 1200;
s_test.name = 'Anastasia';
s_test.data = [1 2 3;1 0 1];
s_test.meta.check = 0;
s_test.meta.records = {'A','B'};
s_test.more = 'cookie';
s_test.recur = s_test;
disp('=====')
sshow(s_test)


end