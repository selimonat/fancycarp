classdef Pupil < Project
    
    
    properties
        STOP       =  1464;
        START      = -250;
    end
    
    properties
        time =[];
        y    =[];
        pupil  =[];
        pupils = [];
        subject = [];
        deltacsp= [];
        trialid = [];
        phase      = [];
        oddball    = [];
        ucs        = [];
        file       = [];
        blink      = [];
    end
    
    
    methods
        %%
        function [obj]=Pupil(subjects,runs)
            
            
            for run = runs(:)'
                subjectss    = [];
                deltacsp   = [];
                trialid    = [];
                phase      = [];
                oddball    = [];
                ucs        = [];
                file       = [];                 
                for subject = subjects(:)'
                   
                   
                    %%                                        
                    path_eye  = sprintf('%seye/data.mat',obj.pathfinder(subject,run));                    
                    dummy     = load(path_eye);
                    trials    = dummy.trials;
                    ttrials   = length(trials);
                    %% decode the PTB message for trial finger print
                    for nt = 1:length(trials)
                        %extract relevant info first from the text message
                        T     = trials(nt).TRIALID.msg(1,:);
                        %add them as msg so that .getfixmat can create the fields
                        pairs = regexp(T,'[\w-]*','match');
                        for np = 1:2:length(pairs)-1
                            trials(nt).(lower(pairs{np})) = str2num(pairs{np+1});
                        end
                    end
                    %% REMOVE THE TRIALS WHICH ARE NOT COMING FROM THE REQUIRED PHASE i.e. ratings
                    valid     = [trials(:).phase];
                    valid     = valid < 2;
                    trials    = trials(valid);%unfortunately the first trials of phase one is never recorded, bug reported (not yet solved) in exp_FearAmy.m
                    ttrials   = sum(valid);
                    %%
                    eye = {'right'};
                    if isfield(trials(1),'left') 
                        if isstruct(trials(1).left)
                            eye = {'left'};
                        end                    
                    end                    
                    %% detect maximum number of samples
                    sizet = zeros(1,ttrials);
                    for nt = 1:ttrials;
                        sizet(nt) = length(trials(nt).(eye{1}).samples.time);
                    end
                    sizet = max(sizet);                                        
                    %% read the data and corresponding time values
                    time   = nan(sizet,ttrials);
                    data   = time;
                    for nt = 1:ttrials
                        current                              = double(trials(nt).(eye{1}).samples.time);                        
                        time(sizet-length(current)+1:end,nt) = current;
                        data(sizet-length(current)+1:end,nt) = double(trials(nt).(eye{1}).samples.pupil(:));
                    end
                    %% find zero points and shift all
                    for nt = 1:ttrials
                        zero_point(nt) = find(time(:,nt)==0);
                    end
                    target_position = min(zero_point);
                    % Shift the data points so that each entry on the matrix corresponds across trials
                    for nt = 1:ttrials
                        time(:,nt) = circshift( time(:,nt),target_position-zero_point(nt));
                        data(:,nt) = circshift( data(:,nt),target_position-zero_point(nt));
                    end
                    %%
                    %                     data = conv2( data, ones(smooth_ks,1)./smooth_ks, 'same' );
                    toi                  = (time >= obj.START) & (time <= obj.STOP);                    
                    time(~toi)           = NaN;
                    data(~toi)           = NaN;                                        
                    %%
                    i        = [find(~isnan(sum(time,2)),1):max(find(~isnan(sum(time,2))))];
                    time     = time(i,:);                    
                    data     = data(i,:);
                    data_ori = data;                                            
                    
                    %%
                    for nt = 1:ttrials
                        start  = trials(nt).(eye{1}).blink.start;
                        stop   = trials(nt).(eye{1}).blink.end;                        
                        %don't care about blinks if they are out of TOI                       
                        i      = (start >= obj.START) & (stop <= obj.STOP);                        
                        start  = start(i);
                        stop   = stop(i);
                        if start == 0 & stop == 0
                            tBlink(nt) = 0;                            
                        else                            
                            tBlink(nt) = length(start);
%                             fprintf('Trial: %03d, %03d blink detected\n',nt,tBlink);
                        end
                        %% run through blinks and nanize the shit
                        for nb = 1:tBlink(nt)
                            %I need to exclude those entries where start and stop are both
                            %equal to zero.                                                        
                            if ~(trials(nt).(eye{1}).blink.start(nb) == 0 & trials(nt).(eye{1}).blink.end(nb) == 0)
                                t                                                                                                              = time(:,nt);                                                                
                                data( (t > (trials(nt).(eye{1}).blink.start(nb) - 100)) & (t < (trials(nt).(eye{1}).blink.end(nb)) + 100),nt ) = NaN;                                                                                                
                            end
                        end                        
                        dummy                         = data(~isnan(data(:,nt)),nt);
%                         try
%                             dummy                         = sgolayfilt(dummy,5,751,[],1);                                
%                             data_s(~isnan(data(:,nt)),nt) = dummy;
%                             data_s(isnan(data(:,nt)),nt)  = NaN;
%                         catch
%                             fprintf('smoothing didn''t work.\n');
%                             data_s(:,nt) = NaN;
%                         end
                        
%                         
%                         if tBlink > 0
%                             plot(time(:,nt),data_ori(:,nt));
%                             hold on;
%                             plot(time(:,nt),data(:,nt),'r.');                                                        
%                             plot(time(:,nt),data_s(:,nt),'m','linewidth',5);                            
%                             hold off;
%                             drawnow;
%                             pause;
%                         end
                    end
                    %%                                        
                    %% baseline correction with mean shifting (only with pupil)
                    for nt = 1:ttrials
                        data(:,nt) = data(:,nt) - nanmean(data(time(:,nt) < 0,nt));
                    end                    
%                     %%
%                     data_conv     = abs(diff(conv2(data, [1 1 1 1 1]'./5 ,'valid' )));
%                     threshold     = nanmean(range(data_conv)) + nanstd(range(data_conv));
%                     outlier_ok    = sum( data_conv >= repmat(threshold,size(data_conv))) == 0;
                    %%
                    
                    for nt = 1:ttrials
                        subjectss(nt)  = subject;
                        deltacsp(nt)   = trials(nt).deltacsp;
                        trialid(nt)    = trials(nt).trialid;
                        phase(nt)      = trials(nt).phase;
                        oddball(nt)    = trials(nt).oddball;
                        ucs(nt)        = trials(nt).ucs;
                        file(nt)       = trials(nt).file;
                    end
                    %                     obj.pupils = data_s;
                                       
                    obj.time           = [obj.time     time];
                    obj.pupil          = [obj.pupil    data];
                    obj.subject        = [obj.subject  subjectss];
                    obj.deltacsp       = [obj.deltacsp deltacsp];
                    obj.trialid        = [obj.trialid  trialid];
                    obj.phase          = [obj.phase    phase];
                    obj.oddball        = [obj.oddball  oddball];
                    obj.ucs            = [obj.ucs      ucs];
                    obj.file           = [obj.file     file];
                    obj.blink          = [obj.blink    tBlink];
                end
            end
        end
    end
    
    
end
