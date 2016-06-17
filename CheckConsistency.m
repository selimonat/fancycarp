%% this script is used to test the mrt/main branch.
runs = 1:2;
for sub = 1:2
    s  = Subject(sub);
    s.SanityCheck(0:s.total_run,'amount','mrt');%check the number of files in the project folder.
    s.plot_log(1);%test whether the log file from the PTB sessions is compatible.
    %% test the fmri preprocessing pipeline.
    s.preprocess_pipeline(runs);%start the preprocessing pipeline.
    figure;s.plot_motionparams(runs(1));%plot the motion parameters
    s.CreateModels(runs(1:2));%create onsets.
    s.FitHRF(1:2);%make a first-level model
    s.VolumeGroupAverage(0,'mrt/w_ss_data.nii');    
end

