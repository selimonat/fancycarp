classdef Pupil < Project
    
    
    properties
        STOP       =  1464;
        START      = -250;
        color_order =  [0.0784    0.3284    1.0000;...
                        0.5784    0.0784    1.0000;...
                        1.0000    0.0784    0.8284;...
                        1.0000    0.0784    0.0784;...
                        1.0000    0.8284    0.0784;...
                        0.5784    1.0000    0.0784;...
                        0.0784    1.0000    0.3284;...
                        0.0784    1.0000    1.0000;...
                        0.5784    0.5784    0.5784;...%null
                        0.0784    0.0784    0.0784]
        clean      = 1;
        
    end
    
    properties
        time =[];        
        pupil  =[];
        x      = [];
        y      =[]        
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
        discarded_trials=NaN;
    end
    
    
    methods        
        function [obj]=Pupil(subjects,runs)
            for run = runs(:)'                                
                for subject = subjects(:)'                    
                    fprintf('Processing Subject %i...\n',subject)
                    % get the data first
                    path_eye  = sprintf('%seye/data.mat',obj.pathfinder(subject,run));
                    dummy     = load(path_eye);
                    trials    = dummy.trials;
                    ttrials   = length(trials);%total number of trials
                    % decode the PTB message for trial finger print
                    for nt = 1:ttrials
                        %extract relevant info first from the text message
                        T     = trials(nt).TRIALID.msg(1,:);
                        %add them as msg so that .getfixmat can create the fields
                        pairs = regexp(T,'[\w-]*','match');
                        for np = 1:2:length(pairs)-1
                            trials(nt).(lower(pairs{np})) = str2num(pairs{np+1});
                        end
                    end
                    % REMOVE THE TRIALS WHICH ARE NOT COMING FROM THE REQUIRED PHASE i.e. ratings
                    valid     = [trials(:).phase];
                    valid     = valid < 2;
                    trials    = trials(valid);%unfortunately the first trials of phase one is never recorded, bug reported (not yet solved) in exp_FearAmy.m
                    ttrials   = sum(valid);%update total trials
                    % detect which eye is recorded.
                    eye = {'right'};
                    if isfield(trials(1),'left')
                        if isstruct(trials(1).left)
                            eye = {'left'};
                        end
                    end
                    % detect maximum number of samples recorded in a trial
                    sizet = zeros(1,ttrials);
                    for nt = 1:ttrials;
                        sizet(nt) = length(trials(nt).(eye{1}).samples.time);
                    end
                    sizet = max(sizet);
                    % read the data and corresponding time values, place to container
                    time   = nan(sizet,ttrials);
                    data   = time;%initialize variables.
                    datax  = time;%initialize variables.
                    datay  = time;%initialize variables.
                    tBlink = [];
                    tnan   = [];
                    for nt = 1:ttrials%run through valid trials and collect
                        current                               = double(trials(nt).(eye{1}).samples.time);
                        time(sizet-length(current)+1:end,nt)  = current;
                        data(sizet-length(current)+1:end,nt)  = double(trials(nt).(eye{1}).samples.pupil(:));
                        datax(sizet-length(current)+1:end,nt) = double(trials(nt).(eye{1}).samples.x(:));
                        datay(sizet-length(current)+1:end,nt) = double(trials(nt).(eye{1}).samples.y(:));
                    end
                    % find zero points and shift all
                    for nt = 1:ttrials
                        zero_point(nt) = find(time(:,nt)==0);
                    end
                    target_position = min(zero_point);
                    % Shift the data points so that each entry on the matrix corresponds across trials
                    for nt = 1:ttrials
                        time(:,nt)  = circshift( time(:,nt),target_position-zero_point(nt));
                        data(:,nt)  = circshift( data(:,nt),target_position-zero_point(nt));
                        datax(:,nt) = circshift( datax(:,nt),target_position-zero_point(nt));
                        datay(:,nt) = circshift( datay(:,nt),target_position-zero_point(nt));
                    end
                    % cut out only the wanted region                    
                    toi                  = (time >= obj.START) & (time <= obj.STOP);
                    time(~toi)           = NaN;
                    data(~toi)           = NaN;
                    datay(~toi)          = NaN;
                    datax(~toi)          = NaN;
                    % cut the stuff to same size
                    i                = [find(~isnan(sum(time,2)),1):max(find(~isnan(sum(time,2))))];
                    time             = time(i,:);
                    data             = data(i,:);
                    datay            = datay(i,:);
                    datax            = datax(i,:);                   
                    % remove trials where there is no data
                    invalid          = sum(datay == 100000000) == size(datay,1);
                    datay(:,invalid) = [];
                    datax(:,invalid) = [];
                    data(:,invalid)  = [];                    
                    time(:,invalid)  = [];
                    trials(invalid)  = [];
                    ttrials          = size(data,2);                    
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
                        tnan(nt)   = 0;                        
                        %% run through blinks and nanize the blink periods
                        for nb = 1:tBlink(nt)
                            %I need to exclude those entries where START and STOP are both
                            %equal to zero.
                            if ~(trials(nt).(eye{1}).blink.start(nb) == 0 & trials(nt).(eye{1}).blink.end(nb) == 0)
                                t  = time(:,nt);
                                data( (t > (trials(nt).(eye{1}).blink.start(nb) - 100)) & (t < (trials(nt).(eye{1}).blink.end(nb)) + 100),nt ) = NaN;
                                datax( (t > (trials(nt).(eye{1}).blink.start(nb) - 100)) & (t < (trials(nt).(eye{1}).blink.end(nb)) + 100),nt ) = NaN;
                                datay( (t > (trials(nt).(eye{1}).blink.start(nb) - 100)) & (t < (trials(nt).(eye{1}).blink.end(nb)) + 100),nt ) = NaN;
                                tnan(nt)= sum(isnan(data(:,nt)));
                            end
                        end
                        dummy = data(~isnan(data(:,nt)),nt);                       
                    end                                    
                    cprintf([1 1 0],'%03d blinks detected\n',sum(tBlink));
                    cprintf([0 1 0],'%i  trials no data...\n',sum(invalid));
                    cprintf([0 1 1],'%3.5g percent of NaNs samples (blink correction)...\n',sum(isnan(data(:)))./length(data(:))*100);
                    % baseline correction with mean shifting (only with pupil)
                    for nt = 1:ttrials
                        data(:,nt) = (data(:,nt) - nanmean(data(time(:,nt) < 0,nt)));
                    end
                    %% interpolate NaNs but no extrapolate, interpolate only if more than 50% of data
                    dataq     = [];
                    discarded = 0;
                    for nt = 1:ttrials
                        Y = data(:,nt);
                        if sum(isnan(Y)) < length(Y)*.5;
                            i = ~isnan(Y);
                            X = time(:,nt);
                            dataq(:,nt)   = interp1(X(i),Y(i),X);
                        else
                            discarded = discarded + 1;
                            dataq(:,nt)   = NaN;
                        end
                    end                    
                    data = dataq;
                    cprintf([1 0 1],'%i trials less than 50 percent valid samples...\n',discarded);                    
                    %% smooth
                    data              = sgolayfilt(data,5,301,[],1);                                        
                    %%                    
                    obj.time           = [obj.time     time];
                    obj.pupil          = [obj.pupil    data];
                    obj.x              = [obj.x        datax];
                    obj.y              = [obj.y        datay];
                    obj.subject        = [obj.subject  repmat(subject,1,ttrials)];
                    obj.deltacsp       = [obj.deltacsp [trials(:).deltacsp]];
                    obj.trialid        = [obj.trialid  [trials(:).trialid]];
                    obj.phase          = [obj.phase    [trials(:).phase]];
                    obj.oddball        = [obj.oddball  [trials(:).oddball]];
                    obj.ucs            = [obj.ucs      [trials(:).ucs]];
                    obj.file           = [obj.file     [trials(:).file]+1];
                    obj.blink          = [obj.blink    tBlink];
                    obj.tnan           = [obj.tnan     tnan];
                    obj.mbi            = [obj.mbi      [trials(:).mblock]];                    
                    %sanity check
                    if length(obj.tnan) ~= length(obj.ucs)
                        keyboard
                    end
                end
            end
            %at this point data can still contain NaN, as we did not have
            %extrapolation.
            if obj.clean
                figure(1);imagesc(data);drawnow;hold on;[y x] = find(isnan(data) == 1);plot(x,y,'ko');hold off;
                invalid = obj.tnan > 0;
                obj.time(:,invalid)      = [];
                obj.pupil(:,invalid)     = [];
                obj.x(:,invalid)         = [];
                obj.y(:,invalid)         = [];
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
                cprintf([1 1 .5],'Cleaned %i nan containined trials...\n',sum(invalid));                
            end
            %
            obj.deltacsp(obj.deltacsp == 3000) = 180+45;%null trials
            obj.deltacsp(obj.deltacsp == 500)  = 0;%180+90;%UCS
            obj.deltacsp(obj.deltacsp == 1000) = 0;%180+90;%oddball
            %
            obj.discarded_trials = NaN(1,max(obj.subject));
            for ns = unique(obj.subject)
                obj.discarded_trials(ns) = 585-sum(obj.subject == ns);%discarded trials
            end
        end
        
        function D = get_evoked(obj,subjects)
            %returns data of SUBJECT in a struct D. Different fields of D
            %contains data as [time,condition].                        
            
            for subject = subjects(:)'                
                for ncond = unique(obj.deltacsp)
                    col                    = ncond./45+4;
                    i                      = obj.subject == subject&obj.deltacsp==ncond;
                    for fields = {'pupil' 'y' 'x'};
                        D(subject).(fields{1}).mean(:,col)    = nanmean(obj.(fields{1})(:,i),2);
                        D(subject).(fields{1}).median(:,col)  = nanmedian(obj.(fields{1})(:,i),2);
                        D(subject).(fields{1}).std(:,col)     = nanstd(obj.(fields{1})(:,i),1,2);
                    end
                end
            end
        end
        
        function D = get_singletrials(obj,subjects)
            %returns trial amplitudes together with trial identifiers in a
            %structure. MEAN is taken across time points.            
            D = [];
            for subject = subjects(:)'              
                D(subject).mean   = nan(65,9);
                D(subject).median = nan(65,9);
                D(subject).count  = zeros(65,9);                
                for i = find(obj.subject == subject)                    
                    mbi                  = obj.mbi(i);
                    cond                 = obj.deltacsp(i)./45+4;
                    if cond < 10
                        D(subject).mean(mbi,cond)     = nanmean(nanmean(obj.pupil(:,i),2));
                        D(subject).median(mbi,cond)   = nanmedian(nanmedian(obj.pupil(:,i),2));
                        D(subject).count(mbi,cond)    = 1;
                    end
                end
            end
        end
            
        function plot_evoked(obj,subject)
            %will plot pupil responses overlaid on the average y
            %position of the eyes for subject SUBJECT.
            %for ns = 5:44;figure(10);subplot(7,6,ns-4);p.plot_timecourse_subject_pupil(ns);end
            %will do it for all subjects.
                        
            t = obj.time(:,1);
            D = obj.get_evoked(subject);
            P = [D(:).pupil];%get the pupil out from valid subjects
            P = cat(3,P(:).median);%store the evoked across 3 dimension
            P = nanmean(P,3);%average                        
            Y = [D(:).y];%get the pupil out from valid subjects
            Y = cat(3,Y(:).median);%store the evoked across 3 dimension
            Y = demean(nanmean(nanmean(Y,3),2));%average                        
            set(gca,'colororder',obj.color_order,'NextPlot', 'replacechildren','fontsize',12);
            plot(t,P(:,1:10));            
            hold on;
            plot(t,Y,'k--');
            xlim([min(t) max(t)]);                        
            hold off
            box off                                
            xlim([-250 1500]);
            set(gca,'xtick',[0],'xgrid','on','fontsize',12,'xticklabel',{'ON'})                 
            title(sprintf('id:%i',subject));
            drawnow
        end        
        function plot_group_average(obj)
            %will plot the group average
            %%
            t = obj.time(:,1);
            %get the data
            D = obj.get_evoked(unique(obj.subject));
            P = [D(:).pupil];%get the pupil out from valid subjects
            P = cat(3,P(:).mean);%store the evoked across 3 dimension
            P = nanmean(P,3);%average                        
            Y = [D(:).y];%get the pupil out from valid subjects
            Y = cat(3,Y(:).median);%store the evoked across 3 dimension
            Y = demean(nanmean(nanmean(Y,3),2));%average                        
            set(gca,'colororder',obj.color_order,'NextPlot', 'replacechildren','fontsize',12);
            plot(t,P(:,1:10),'linewidth',3);            
            hold on;
            plot(t,Y,'k--','linewidth',3);            
            axis tight;
            hold off
            box off 
            xlabel('time (ms)','fontsize',20);
            ylabel('pupil size','fontsize',20);
            ylim([-120 120])
            SetTickNumber(gca,5,'y');
            xlim([-250 1500]);
            set(gca,'xtick',[0 750 1400 1500],'xgrid','on','fontsize',20,'xticklabel',{'ON' 'Switch' 'UCS' '       OFF'})
            legend({'-135' '-90' '-45' 'CS+' '+45' '+90' '+135' '+180' 'Null' 'UCS' 'y pos'},'location','southwest','fontsize',12)
            legend boxoff          
            axis square
            title('Group Average')            
        end        
        function [out]=timecourse_difference(obj,subjects)            
            % plots the difference between relevant and control conditions.
            out1 = [];out2 = [];
            for subject = subjects
                D     =  obj.get_evoked(subject);
                out1 = [out1 D(subject).pupil.mean(:,4)-D(subject).pupil.mean(:,8)];
                out2 = [out2 D(subject).pupil.mean(:,6)-D(subject).pupil.mean(:,2)];%orthogonal conditions                
            end
            % plot
            M = mean(out1,2);
            S = std(out1,1,2);
            plot(obj.time(:,1),M,'r','linewidth',3);
            hold on
            plot(obj.time(:,1),M+S/2,'r-.','linewidth',1);
            plot(obj.time(:,1),M-S/2,'r-.','linewidth',1);
            M = mean(out2,2);
            S = std(out2,1,2);
            plot(obj.time(:,1),M,'c','linewidth',3);            
            plot(obj.time(:,1),M+S/2,'c-.','linewidth',1);
            plot(obj.time(:,1),M-S/2,'c-.','linewidth',1);
            axis tight;
            hold off;
            box off;
            xlabel('time (ms)','fontsize',20);
            ylabel('\Delta pupil','fontsize',20);                        
            xlim([-250 1500]);
            set(gca,'xtick',[0 750 1400 1500],'xgrid','on','fontsize',20,'xticklabel',{'ON' 'Switch' 'UCS' '       OFF'});
            legend({'-135' '-90' '-45' 'CS+' '+45' '+90' '+135' '+180' 'Null' 'UCS' 'y pos'},'location','southwest','fontsize',12);
            legend boxoff;
            axis square;
            legend({'0 vs. 180' '-90 vs. +90'},'location','northwest','fontsize',12);
            legend boxoff;
            title('Group Average');
            hold off;            
        end
        
        function D = plot_tuning_subject(obj,subject)
                  
            keyboard
            D = obj.get_evoked(subject);
            M = mean(D(subject).pupil.mean(:,1:8));
            Project.plot_bar(M');
            title(sprintf('id:%02d',subject),'fontsize',12);%subject and face id
            set(gca,'fontsize',12);
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
      
        
        
        
              
        function D = plot_group_tuning(obj,time)
                 
            D = [];
            for ncond = unique(obj.deltacsp)
                i  = obj.deltacsp==ncond;
                D  = [D nanmean(obj.pupil(:,i),2)];
            end
            set(gca,'colororder',GetFearGenColors,'NextPlot', 'replacechildren');
            [~,i]=min(abs(obj.time(:,1)-time));
            o = D(i,1:9);
            bar(o);                      
        end
        
        
        
        
        function H = histogram(obj)
            %computes for each subject and condition the number of recorded
            %trials in .Count. In .BLINK the number of detected blinks.
            %
            subjects = unique(obj.subject);
            for nsubject = subjects;
                i = obj.subject == nsubject;
                for ncond = unique(obj.deltacsp(i));
                    index                                         = logical(i.*(obj.deltacsp==ncond));
                    H.Count(nsubject-min(subjects)+1,ncond./45+4) = sum(double(index));%number of trials
                    H.Blink(nsubject-min(subjects)+1,ncond./45+4) = sum(obj.blink(index));
                end                
            end                        
            %% plot the count as a bar plot
            figure(1);clf
            subplot(2,1,1)
            set(gcf,'position',[275         276        1406         679]);
            bar(5:44,mean(H.Count(:,1:8),2),.9,'k');
            set(gca,'ytick',linspace(0,64,5),'xtick',subjects,'xticklabel',subjects);
            box off
            set(gca,'ytick',linspace(0,64,5),'xtick',subjects,'xticklabel',subjects);
            set(gca,'ytick',linspace(0,64,5),'xtick',subjects,'xticklabel',subjects,'ygrid','on');
            xlabel('Conditions')
            xlabel('Subjects')
            ylabel('Trials')                        
            subplot(2,1,2)
            bar(1:44,obj.discarded_trials,.9,'k');
            set(gca,'ytick',linspace(0,585,4),'xtick',subjects,'xticklabel',subjects);
            box off
            set(gca,'ytick',linspace(0,585,4),'xtick',subjects,'xticklabel',subjects);
            set(gca,'ytick',linspace(0,585,4),'xtick',subjects,'xticklabel',subjects,'ygrid','on');
            xlabel('Conditions')
            xlabel('Subjects')
            ylabel('# Discarded Trials')                        
        end
    end
end
