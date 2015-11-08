classdef Pupil < Project
    
    
    properties
        STOP       = 2000;
        START      = -1000;
    end
    
    properties
        time =[];
        y    =[];
        pupil  =[];
        subject = [];
        deltacsp= [];
        trialid = [];
        phase      = [];
        oddball    = [];
        ucs        = [];
        file       = [];
    end
    
    
    methods
        %%
        function [obj]=Pupil(subjects,runs)
            
            for run = runs(:)'
                for subject = subjects(:)'
                    %
                    dummy   = load(obj.path2data(subject,run,'eye'));
                    trials  = dummy.trials;
                    ttrials = length(trials);
                    %% detect maximum number of samples
                    sizet = zeros(1,ttrials);
                    for nt = 1:ttrials;
                        sizet(nt) = length(trials(nt).right.samples.time);
                    end
                    sizet = max(sizet);
                    %% decode the PTB message for trial finger print
                    for nt = 1:length(trials)
                        %extract relevant info first from the text message
                        T = trials(nt).TRIALID.msg(1,:);
                        %add them as msg so that .getfixmat can create the fields
                        pairs = regexp(T,'[\w-]*','match');
                        for np = 1:2:length(pairs)-1
                            trials(nt).(lower(pairs{np})).msg = pairs{np+1};
                        end
                    end
                    %% HERE REMOVE THE TRIALS WHICH ARE NOT COMING FROM THE REQUIRED PHASE i.e. ratings
                    %% read the data and corresponding time values
                    time   = nan(sizet,ttrials);
                    data   = time;
                    for nt = 1:ttrials
                        current                              = double(trials(nt).right.samples.time);
                        time(sizet-length(current)+1:end,nt) = current;
                        data(sizet-length(current)+1:end,nt) = double(trials(nt).right.samples.pupil(:));
                    end
                    %%
                    %                     data = conv2( data, ones(smooth_ks,1)./smooth_ks, 'same' );
                    toi                     = (time >= obj.START) & (time <= obj.STOP);
                    time(~toi)              = NaN;
                    data(~toi)              = NaN;
                    %
                    for nt = 1:ttrials
                        start_point(nt) = find(toi(:,nt),1);
                    end
                    start_point_smallest = min(start_point);
                    %% Shift the data points so that each entry on the matrix corresponds across trials
                    for nt = 1:ttrials
                        time(:,nt) = circshift( time(:,nt),start_point_smallest-start_point(nt));
                        data(:,nt) = circshift( data(:,nt),start_point_smallest-start_point(nt));
                    end
                    %%
                    i    = [find(~isnan(sum(time,2)),1):max(find(~isnan(sum(time,2))))];
                    time = time(i,:);
                    data = data(i,:);
                    %%
                    for nt = 1:ttrials
                        tBlink = length(trials(nt).right.blink.start);
                        for nb = 1:tBlink
                            %I need to exclude those entries where start and stop are both
                            %equal to zero.
                            fprintf('Trial: %03d, %03d blink detected\n',nt,tBlink);
                            if ~(trials(nt).right.blink.start(nb) == 0 & trials(nt).right.blink.end(nb) == 0)
                                t = time(:,nt);
                                data( (t > (trials(nt).right.blink.start(nb) - 100)) & (t < (trials(nt).right.blink.end(nb)) + 100),nt ) = NaN;
                            end
                        end
                    end
                    %% baseline correction with mean shifting (only with pupil)
                    for nt = 1:ttrials
                        data(:,nt) = data(:,nt) - nanmean(data(time(:,nt) < 0,nt));
                    end
                    
                    %%
                    data_conv     = abs(diff(conv2(data, [1 1 1 1 1]'./5 ,'valid' )));
                    threshold     = nanmean(range(data_conv)) + nanstd(range(data_conv));
                    outlier_ok    = sum( data_conv >= repmat(threshold,size(data_conv))) == 0;
                    %%
                    for nt = 1:ttrials
                        obj.subject(nt)    = subject;
                        obj.deltacsp(nt)   = str2num(trials(nt).deltacsp.msg);
                        obj.trialid(nt)    = str2num(trials(nt).trialid.msg);
                        obj.phase(nt)      = str2num(trials(nt).phase.msg);
                        obj.oddball(nt)    = str2num(trials(nt).oddball.msg);
                        obj.ucs(nt)        = str2num(trials(nt).ucs.msg);
                        obj.file(nt)       = str2num(trials(nt).file.msg);
                    end
                    obj.time = time;
                    obj.pupil = data;
                    %%
                    %                     obj.cond =
                    
                    
                end
            end
        end
    end
    
    
end
