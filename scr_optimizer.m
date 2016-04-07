%% get the subject object
viz = 1;
addpath('/Users/onat/Documents/Code/Matlab/fancycarp/');
% options  = optimset('algorithm',{'levenberg-marquardt',.01},'display','iter','MaxFunEvals',25000,'maxiter',25000,'tolx',10^-12,'tolfun',10^-12);
    options  = optimset('algorithm',{'levenberg-marquardt',.01},'display','iter','MaxFunEvals',500000,'maxiter',500000,'tolx',10^-12,'tolfun',10^-12,'OutputFcn',@scr_optimizer_plot);
global time
global onsets
global original
filter_string = 'calibration';
for subject = 38%[Project.subjects_1500 Project.subjects_600];
    tic
%     try
        fprintf('Working on subject %02d, started at %s\n',subject,datestr(now,'hh:mm:ss'))
        s  = Subject(subject);
        %detect the calibration phases
        phase_id = cellfun(@(x) ~isempty(regexp(x,filter_string)), s.scr.BlockNames );
        s.scr.cut(phase_id);%cut the required segment in
        s.scr.smooth('sgolay');%smooth the data to get pixel noise out
        s.scr.tonic_tfals;%detect the tonic component
        s.scr.downsample(50);%makes everthing faster
        if viz
            s.scr.plot;
        end
        %% collect the events, data and time variable to build the generative model...
        event_id                           = cellfun(@(x) ~isempty(regexp(x,filter_string)), s.scr.event_name );%detect the relevant event column
        [onset_sample x]                   = find(s.scr.event(:,event_id));%in samples
        % exclude samples too far from first and last events.
        %
        onsets                             = s.scr.time(onset_sample)./1000;%onset times in s
        data                               = s.scr.phasic;%data we want to model
        original                           = data;
        time                               = s.scr.time./1000;
        data(time < (min(onsets)-10))      = 0;
        data(time > (max(onsets)+10))      = 0;
        fun                                = @(params) sum(abs((data - scr_model( time , onsets(:), params ))));
        funlsq                             = @(params) (abs(data - scr_model( time , onsets(:), params ))).^2;
        %%
        R    = range(data);
        t    = 100000;        
        err  = [];
        tau1 = zeros(1,t);tau2=tau1;latency=tau1;amp1=tau1;
        n    = 0;
        while n < t
            n                = n +1;
            tau1(n)          = randsample(linspace(1,7,100),1);
            tau2(n)          = tau1(n) - randsample(linspace(0,tau1(n),100),1);
            amp1(n)          = randsample(linspace(R*.1,R,1000),1);
            latency(n)       = randsample(linspace(1,3,10),1);
            params           = [ones(length(onsets),1)*amp1(n) ;latency(n);tau1(n); tau2(n)];
            err(n)           = fun(params);
            if 1
                plot(data);
                hold on;
                plot(scr_model( s.scr.time./1000 , onsets(:), params ),'r');
                hold off
                ylim([0 R].*1.4)
                title(mat2str(err(n)));
                drawnow;
            end
        end
        %plot the best random one and replot the data
        [~,n]   = min(err);
        params  = [ones(length(onsets),1)*amp1(n) ;latency(n);tau1(n); tau2(n)];
        %%
        LB      = [zeros(1,length(onsets))   0 0 0];
        UB      = [ones(1,length(onsets))*4  5 10 10];
        params0 = params;
        if viz
            clf;
            plot(s.scr.time./1000,data,'k');
            hold on;
            plot(s.scr.time./1000,scr_model(s.scr.time/1000,onsets,params0),'r');
            title(mat2str(fun(params0)));
            drawnow;
        end
        %% fit everything                
        % x  = fmincon(fun,params0,[],[],[],[],LB,UB,[],options)
        %     x = fminsearch(fun,params0,options);
        x = lsqnonlin(funlsq,params0,LB,[],options);
        clf;
        set(gcf,'position',[156        1334        1677         449]);
        plot(time,data,'k');
        hold on;
        plot(time,scr_model(time,onsets,x),'r');
%         SaveFigure(sprintf('/Users/onat/Pictures/SCR/scr_model/subject_%02d_calibration.png',subject),'-r250')
        %% store
        param_fit{subject} = x;
%         save('/Users/onat/Pictures/SCR/scr_model/data.mat','param_fit');
%     catch
%         fprintf('Didn''t work for subject %02d', subject);
%     end
    fprintf('Finished subject %02d in %g minutes...\n',subject,toc);
end
