function movie = tbx_cfg_acid_dwi_series_movie

% reference image
movie_ref         = cfg_files;
movie_ref.tag     = 'movie_ref';
movie_ref.name    = 'Reference image';
movie_ref.help    = {'Select reference image.'};
movie_ref.filter  = 'image';
movie_ref.ufilter = '.*';
movie_ref.num     = [0 1];

% 1. dataset
movie_set1         = cfg_files;
movie_set1.tag     = 'movie_set1';
movie_set1.name    = '1. dMRI dataset';
movie_set1.help    = {'Select the first dataset.'};
movie_set1.filter  = 'image';
movie_set1.ufilter = '.*';
movie_set1.num     = [1 Inf];

% 2. dataset
movie_set2         = cfg_files;
movie_set2.tag     = 'movie_set2';
movie_set2.name    = '2. dMRI dataset (optional)';
movie_set2.help    = {'Select the second dataset (or click Done for none).'};
movie_set2.filter  = 'image';
movie_set2.ufilter = '.*';
movie_set2.num     = [0 Inf];
movie_set2.val     = {{''}};

% 3. dataset
movie_set3         = cfg_files;
movie_set3.tag     = 'movie_set3';
movie_set3.name    = '3. dMRI dataset (optional)';
movie_set3.help    = {'Select the third dataset (or click Done for none)'};
movie_set3.filter  = 'image';
movie_set3.ufilter = '.*';
movie_set3.num     = [0 Inf];
movie_set3.val     = {{''}};

% dummy whether all slices
movie_dummy_allslices        = cfg_menu;
movie_dummy_allslices.tag    = 'movie_dummy_allslices';
movie_dummy_allslices.name   = 'Display all slices?';
movie_dummy_allslices.help   = {''
                    'Here you can specify whether you want to display all slices or only a single slice.'
                                ''};
movie_dummy_allslices.labels = {'All slices','Single slice'};
movie_dummy_allslices.values = {1 0};
movie_dummy_allslices.val    = {0};

% dummy which order to display
movie_dummy_order      = cfg_menu;
movie_dummy_order.tag  = 'movie_dummy_order';
movie_dummy_order.name = 'If all slices, in which order?';
movie_dummy_order.help = {''
'When displaying all slices, in which order?.'
'0: inner loop across volumes'
'1: inner loop across slices'
''};
movie_dummy_order.labels = {'First volumes','First slices'};
movie_dummy_order.values = {0 1};
movie_dummy_order.val    = {0};

% slice position
movie_slice         = cfg_entry;
movie_slice.tag     = 'movie_slice';
movie_slice.name    = 'If single slice, which slice?';
movie_slice.help    = {'Provide an integer number that defines the slice at which the movie is shown. Note that it must be within the dimensions of the image. If no slice number is specified, the middle slice is used.'};
movie_slice.strtype = 'e';
movie_slice.num     = [0 Inf];
movie_slice.val     = {[]};

% dummy for saving the movie
movie_dummy_movie      = cfg_menu;
movie_dummy_movie.tag  = 'movie_dummy_movie';
movie_dummy_movie.name = 'Save Movie?';
movie_dummy_movie.help = {''
'Here you can choose to save the movie that is generated.'
''};
movie_dummy_movie.labels = {'Yes','No'};
movie_dummy_movie.values = {1 0};
movie_dummy_movie.val    = {1};

% dummy for saving the movie
movie_dummy_contour      = cfg_menu;
movie_dummy_contour.tag  = 'movie_dummy_contour';
movie_dummy_contour.name = 'Plot contour lines?';
movie_dummy_contour.help = {''
'Here you can choose if contour lines of a selected dataset is shown on all images.'
''};
movie_dummy_contour.labels = {'No','Of reference image','Of Dataset 1','Of Dataset 2','Of Dataset 3'};
movie_dummy_contour.values = {4 0 1 2 3};
movie_dummy_contour.val    = {4};

% pause
movie_pause         = cfg_entry;
movie_pause.tag     = 'movie_pause';
movie_pause.name    = 'Time interval';
movie_pause.help    = {'Video is paused for this period after each image'};
movie_pause.strtype = 'e';
movie_pause.num     = [1 1];
movie_pause.val     = {0.01};

% % Interval reference image
intv_ref         = cfg_entry;
intv_ref.tag     = 'intv_ref';
intv_ref.name    = 'Intensity range for reference image';
intv_ref.help    = {'Provide the intensity limits for the reference image (default is for b=0 image: [0 1000]).'};
intv_ref.strtype = 'e';
intv_ref.num     = [1 2];
intv_ref.val     = {[1 1000]};
% 
% % Interval source images
intv_src         = cfg_entry;
intv_src.tag     = 'intv_src';
intv_src.name    = 'Intensity range for source images';
intv_src.help    = {'Provide the intensity limits for the source images (default is for b=1000 image: [0 300]).'};
intv_src.strtype = 'e';
intv_src.num     = [1 2];
intv_src.val     = {[1 300]};

movie      = cfg_exbranch;
movie.tag  = 'movie';
movie.name = 'DWI series movie';
movie.val  = {movie_ref, movie_set1, movie_set2, movie_set3, intv_ref, intv_src, movie_dummy_contour, movie_dummy_allslices, movie_dummy_order, movie_slice, ...
    movie_dummy_movie, movie_pause};
movie.help = {
                    ''
};
movie.prog = @local_movie;
movie.vout = @vout_movie;

end

function out = local_movie(job)
    acid_dwi_series_movie(char(job.movie_ref), char(job.movie_set1), char(job.movie_set2), char(job.movie_set3),...
        job.movie_dummy_allslices, job.movie_dummy_order, job.movie_slice, job.intv_ref, job.intv_src, job.movie_dummy_movie, job.movie_pause, job.movie_dummy_contour);
    
    out.invol1 = job.movie_set1(:);
    out.invol2 = job.movie_set2(:);
    out.invol3 = job.movie_set3(:);
end

function dep = vout_movie(~)
    dep(1)            = cfg_dep;
    dep(1).sname      = '1st input dataset';
    dep(1).src_output = substruct('.','invol1'); 
    dep(1).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(2)            = cfg_dep;
    dep(2).sname      = '2nd input dataset';
    dep(2).src_output = substruct('.','invol2'); 
    dep(2).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
    dep(3)            = cfg_dep;
    dep(3).sname      = '3rd input dataset';
    dep(3).src_output = substruct('.','invol3'); 
    dep(3).tgt_spec   = cfg_findspec({{'filter','image','strtype','e'}});
end