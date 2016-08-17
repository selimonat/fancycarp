classdef Pupil < Project
    
    
    properties
        STOP       =  1464;
        START      = -250;
        blinking_subjects = [5 15 25 28 29];
        clean      = 1;
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
        mbi        = [];
        tnan    =[];
    end
    
    
    methods
        %%
        function [obj]=Pupil(subjects,runs)                        
            for run = runs(:)'
                subjectss    = [];
                deltacsp     = [];
                trialid      = [];
                phase        = [];
                oddball      = [];
                ucs          = [];
                file         = [];
                tBlink       = [];
                for subject = subjects(:)'
                    fprintf('Processing Subject %i...\n',subject)
                    %% get the data first                                   
                    path_eye  = sprintf('%seye/data.mat',obj.pathfinder(subject,run));                    
                    dummy     = load(path_eye);
                    trials    = dummy.trials;
                    ttrials   = length(trials);%total number of trials
                    %% decode the PTB message for trial finger print
                    for nt = 1:ttrials
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
                    ttrials   = sum(valid);%update total trials
                    %% detect which eye is recorded.
                    eye = {'right'};
                    if isfield(trials(1),'left') 
                        if isstruct(trials(1).left)
                            eye = {'left'};
                        end                    
                    end                    
                    %% detect maximum number of samples recorded in a trial
                    sizet = zeros(1,ttrials);
                    for nt = 1:ttrials;
                        sizet(nt) = length(trials(nt).(eye{1}).samples.time);
                    end
                    sizet = max(sizet);                                        
                    %% read the data and corresponding time values, place to container
                    time   = nan(sizet,ttrials);
                    data   = time;%initialize variables.
                    for nt = 1:ttrials%run through valid trials and collect
                        current = double(trials(nt).(eye{1}).samples.time);                        
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
                    %% cut out only the wanted region
                    %                     data = conv2( data, ones(smooth_ks,1)./smooth_ks, 'same' );
                    toi                  = (time >= obj.START) & (time <= obj.STOP);                    
                    time(~toi)           = NaN;
                    data(~toi)           = NaN;
                    %% cut the stuff to same size
                    i        = [find(~isnan(sum(time,2)),1):max(find(~isnan(sum(time,2))))];
                    time     = time(i,:);                    
                    data     = data(i,:);
                    data_ori = data;                                            
                    
                    %% take care of the blinks in every single trial
                    for nt = 1:ttrials                        
                        %start and stop times of blinks in this trial
                        start      = trials(nt).(eye{1}).blink.start;
                        stop       = trials(nt).(eye{1}).blink.end;
                        %clean the start and stop from 0 entries
                        i          = find(sum([start;stop] == 0) == 2);
                        start(i)   = [];stop(i) = [];
                        %discard blinks if they are out of TOI.
                        i          = (start >= obj.START) | (stop <= obj.STOP);
                        start      = start(i);
                        stop       = stop(i);                        
                        tBlink(nt) = length(start);                        
                        tnan(nt)= 0;
                        fprintf('Trial: %03d, %03d blink detected\n',nt,tBlink(nt));
                        %% run through blinks and nanize the blink periods
                        for nb = 1:tBlink(nt)
                            %I need to exclude those entries where START and STOP are both
                            %equal to zero.
                            if ~(trials(nt).(eye{1}).blink.start(nb) == 0 & trials(nt).(eye{1}).blink.end(nb) == 0)
                                t  = time(:,nt);                                                                
                                data( (t > (trials(nt).(eye{1}).blink.start(nb) - 100)) & (t < (trials(nt).(eye{1}).blink.end(nb)) + 100),nt ) = NaN;                                                                                                
                                tnan(nt)= sum(isnan(data(:,nt)));
                            end
                        end
                        dummy = data(~isnan(data(:,nt)),nt);
%                         if tBlink(nt) == 1
%                             keyboard
%                         end
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
                    %% baseline correction with mean shifting (only with pupil)
                    for nt = 1:ttrials
                        data(:,nt) = (data(:,nt) - nanmean(data(time(:,nt) < 0,nt)));
                    end                                        
                    %% interpolate NaNs but no extrapolate
                    dataq = [];                    
                    for nt = 1:ttrials
                        Y = data(:,nt);
                        if sum(isnan(Y)) < length(Y)*.5;
                            i = ~isnan(Y);
                            X = time(:,nt);
                            dataq(:,nt)   = interp1(X(i),Y(i),X);
                        else
                            dataq(:,nt)   = NaN;
                        end
                    end
                    data = dataq;
                    %% smotth
                    data              = sgolayfilt(data,5,301,[],1);
                    %%                    
                    for nt = 1:ttrials
                        subjectss(nt)  = subject;
                        deltacsp(nt)   = trials(nt).deltacsp;
                        trialid(nt)    = trials(nt).trialid;
                        phase(nt)      = trials(nt).phase;
                        oddball(nt)    = trials(nt).oddball;
                        ucs(nt)        = trials(nt).ucs;
                        file(nt)       = trials(nt).file+1;
                        mbi(nt)        = trials(nt).mblock;
                    end                   
% % %                     %sort by conditions
% % %                     [~,i]    = sort(deltacsp);                    
% % %                     deltacsp = deltacsp(i);
% % %                     subjectss=subjectss(i);
% % %                     trialid  = trialid(i);
% % %                     phase    = phase(i);
% % %                     oddball  = oddball(i);
% % %                     ucs      = ucs(i);
% % %                     file     = file(i);
% % %                     tBlink   = tBlink(i);
% % %                     mbi      = mbi(i);
                    %
                    deltacsp(deltacsp == 3000) = 180+45;
                    deltacsp(deltacsp == 500)  = 180+90;
                    deltacsp(deltacsp == 1000) = 180+135;
                    %
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
                    obj.tnan           = [obj.tnan     tnan];
                    obj.mbi            = [obj.mbi      mbi];                  
                end
            end
            if obj.clean
                invalid = obj.tnan > 0;
                obj.time(:,invalid)      = [];
                obj.pupil(:,invalid)     = [];
                obj.subject(invalid)   = [];
                obj.deltacsp(invalid)  = [];
                obj.trialid(invalid)   = [];
                obj.phase(invalid)     = [];
                obj.oddball(invalid)   = [];
                obj.ucs(invalid)       = [];
                obj.file(invalid)      = [];
                obj.blink(invalid)     = [];
                obj.mbi(invalid)       = [];
                obj.tnan(invalid)       = [];
            end
        end
        
        function [mat]=plot_subject_matrix(obj,subject)
            ffigure(subject);
            i     = obj.subject == subject;
            mat   = obj.pupil(:,i);
            mat   = mat';
            [d u] = GetColorMapLimits(mat(:),2);
            subplot(1,3,1)
            imagesc(obj.time(:,1),1:size(mat,1),mat,[d u]);            
            colorbar
            subplot(1,3,2)
            imagesc(isnan(mat));
            fprintf('%i trials with detected blinks...\n',sum(sum(isnan(mat),2) > 0));
            subplot(1,3,3)
            plot(obj.tnan(:,i),1:sum(i),'r');axis ij;            
        end        
%         function [mat] = clean_outlier(obj,subject)
%             %
%             outlier = 1;
%             blink   = 1;
%             i       = obj.subject == subject;
%             mat     = obj.pupil(:,i);
%             %remove outliers
%             if outlier
%                 ave     = nanmean(mat);
%                 m       = nanmean(ave);
%                 s       = nanstd(ave);
%                 invalid = (ave > (m+s))|(ave < (m-s));
%                 mat(:,invalid) = [];            
%             end
%             if blink
%                 i = sum(isnan(mat)) == 0;
%                 mat(:,~i) = [];
%             end     
%             mat = detrend(mat);
%         end
        function plot_subject_timecourse(obj,subject)
                 
            D = [];
            for ncond = unique(obj.deltacsp)
                i  = obj.subject == subject&obj.deltacsp==ncond;
                D  = [D nanmean(obj.pupil(:,i),2)];
            end
            set(gca,'colororder',GetFearGenColors,'NextPlot', 'replacechildren');
            plot(D(:,1:8));
            hold off;            
        end
        
        function D = plot_subject_tuning(obj,subject)
                 
            D = [];
            for ncond = unique(obj.deltacsp)
                i  = obj.subject == subject&obj.deltacsp==ncond;
                D  = [D nanmean(obj.pupil(:,i),2)];
            end
            set(gca,'colororder',GetFearGenColors,'NextPlot', 'replacechildren');
            D = D(1500,1:8);
            bar(D);                      
        end
        
        function H = histogram(obj)
            subjects = unique(obj.subject);
            for nsubject = subjects;
                i = obj.subject == nsubject;
                for ncond = unique(obj.deltacsp(i));
                    index                                         = logical(i.*(obj.deltacsp==ncond));
                    H.Count(nsubject-min(subjects)+1,ncond./45+4) = sum(double(index));%number of trials
                    H.Blink(nsubject-min(subjects)+1,ncond./45+4) = sum(obj.blink(index));
                end                
            end
        end
    end
end
