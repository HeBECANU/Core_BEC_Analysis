%analyze the oscilations during the ramp

%import_opts.dir='X:\ML_Optimal_Transport\ml_shunt_continued\insitu_14_150_best_20170926\';
import_opts.dir='\\amplpc29\Users\TDC_user\ProgramFiles\my_read_tdc_gui_v1.0.1\dld_output\'
import_opts.file_name='d';
import_opts.force_reimport=0;
import_opts.force_forc=1;

import_opts.dld_xy_rot=0.61;
import_opts.txylim=[[0,2];[-30e-3, 30e-3];[-30e-3, 30e-3]];     %tight XY lims to eliminate hot spot from destroying pulse widths


%%%%------------------------done user var-------------------------------%

addpath('Colormaps') 
import_opts.shot_num=find_data_files(import_opts);
data=import_data_legacy(import_opts);



