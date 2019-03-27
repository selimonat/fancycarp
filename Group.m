classdef Group < Project
    properties (Hidden,Constant)
        mean_correction = 0;%decides if mean correction should be applied
        align_tunings   = 1;%should ratings be aligned to CS+ face
    end
    properties
        subject
        ids
        pmf
        tunings
        SI
        sigma_cond
        sigma_test
    end
    
    methods
        %%
        function group = Group(subjects)
            c = 0;
            for s = subjects(:)'
                fprintf('subject: %03d\n',s)
                c                = c+1;
                dummy            = Subject(s);
                group.subject{c} = dummy;
                group.ids        = subjects;
                try
                    group.getPMF;
                end
            end
        end
        function csps = getcsp(self)
            csps = [];
            for s = 1:length(self.subject)
                csps = [csps self.subject{s}.csp];
            end
        end
        %
        function T = GroupTable(self)
            
            %%
            T = table();
            sub = [];
            for ns = 1:length(self.subject)
                sub = [sub ns];
                try
                    T = [T;[self.subject{ns}.feargen_scr(:).param_table]];
                catch
                    dummy = array2table(nan(1,size(T,2)),'variablenames',T.Properties.VariableNames);
                    T = [T;dummy];
                end
            end
            T.subject = sub(:);
            %%
            %             T2 = table();sub = [];
            %             for ns = 1:length(self.subject)
            %                 sub = [sub ns]
            %                 try
            %                     T2 = [T2;[self.subject{ns}.feargen_scr(:).param_table]];
            %                 catch
            %                     dummy = array2table(nan(1,size(T2,2)),'variablenames',T2.Properties.VariableNames);
            %                     T2 = [T2;dummy];
            %                 end
            %             end
            %             %%
            %             T3 = [T T2];
        end
        function FitStan(self)
            %will fit data a VonMises function using the MLE estimates as
            %initial values;
            %% raw data
            data.x=[];data.y=[];data.m=[];c = 0;
            for phase  = 2:4
                c  = c+1;
                for ns = 1:length(self.subject)
                    data.y = [data.y self.subject{ns}.get_rating(phase).y(:)];
                    data.x = [data.x self.subject{ns}.get_rating(phase).x(:)];
                    data.m = [data.m c];
                end
            end
            data.X = size(data.x,1);
            data.T = size(data.y,2);
            data.M = length(unique(data.m));
            data.x = data.x(:,1);
            %% initial estimates
            T                = self.GroupTable;
            init             = [];
            init.amp         = [T.rating_amp_02'    T.rating_amp_03'      T.rating_amp_04'];
            init.kappa       = max(0.01,([T.rating_kappa_02'  T.rating_kappa_03'    T.rating_kappa_04']));
            init.offset      = [T.rating_offset_02' T.rating_offset_03'   T.rating_offset_04'];
            init.sigma_y     = [T.rating_sigma_y_02' T.rating_sigma_y_03' T.rating_sigma_y_04'];
            %
            for m = unique(data.m)
                init.sigma_y(m)     = median(init.sigma_y(data.m == m));
                init.mu_amp(m)      = median(init.amp(data.m == m));
                init.mu_offset(m)   = median(init.offset(data.m == m));
                init.mu_kappa(m)    = max(0.01,(median(init.kappa(data.m == m))));
            end
            %
            init.sigma_amp   = repmat(0.5,1,length(unique(data.m)));
            init.sigma_offset= repmat(0.5,1,length(unique(data.m)));
            init.sigma_kappa = repmat(0.5,1,length(unique(data.m)));
            init.sigma_y     = repmat(0.5,1,length(unique(data.m)));
            
            %%
            addpath('/home/onat/Documents/Code/C++/cmdstan-2.12.0/');
            addpath('/home/onat/Documents/Code/Matlab/MatlabProcessManager/');
            addpath('/home/onat/Documents/Code/Matlab/MatlabStan/')
            cd ~/Desktop/FitVonMises/;
            !rm FitVonMises FitVonMises.hpp FitVonMises.cpp  output-* temp.*;
            fit    = stan('file','FitVonMises.stan','data',data,'verbose',true,'iter',1000,'init',init);
            cd ~/Desktop/fancycarp;
        end
        
        function ModelRatings(self,run,funtype)
            %create a tuning object and fits FUNTYPE to it.
            self.tunings.rate{run} = Tuning(self.Ratings(run));%create a tuning object for the RUN for ratings.
            self.tunings.rate{run}.SingleSubjectFit(funtype);%call fit method from the tuning object
        end
        
        function ModelSCR(self,run,funtype)
            %create a tuning object and fits FUNTYPE to it.
            self.tunings.scr = Tuning(self.getSCRbars(run));%create a tuning object for the RUN for SCRS.
            %             self.tunings.scr.SingleSubjectFit(funtype);%call fit method from the tuning object
        end
        
        function getSCRtunings(self,run,funtype)
            self.ModelSCR(run,funtype);
        end
        
        function getSI(self,funtype)
            %fits FUNTYPE to behavioral ratings and computes Sharpening
            %Index.
            self.ModelRatings(3,funtype);
            self.ModelRatings(4,funtype);
            self.sigma_cond = [];
            self.sigma_test = [];
            for s = 1:length(self.subject)
                if funtype==3
                    self.SI         = [self.SI; self.tunings.rate{3}.singlesubject{s}.Est(:,2) - self.tunings.rate{4}.singlesubject{s}.Est(:,2)];%take the diff of sigma parameters.
                elseif funtype == 8
                    self.SI         = [self.SI; self.tunings.rate{4}.singlesubject{s}.Est(:,2) - self.tunings.rate{3}.singlesubject{s}.Est(:,2)];%take the diff of sigma parameters.
                else
                    self.SI = [];
                end
                self.sigma_cond = [self.sigma_cond; self.tunings.rate{3}.singlesubject{s}.Est(:,2)];
                self.sigma_test = [self.sigma_test; self.tunings.rate{4}.singlesubject{s}.Est(:,2)];
            end
        end
        function [ratings,ratings_sd] = getRatings(self,phases)
            for ph = phases(:)';
                xc=0;
                for x  = unique(self.Ratings(ph).x(:)')
                    xc=xc+1;
                    i             = self.Ratings(ph).x(1,:) == x;
                    ratings(:,xc,ph)      = mean(self.Ratings(ph).y(:,i),2);
                    ratings_sd(:,xc,ph)   = std(self.Ratings(ph).y(:,i),0,2);
                end
            end
            %clean the zeros
            ratings    = ratings(:,:,phases);
            ratings_sd = ratings_sd(:,:,phases);
        end
        %%
        function [out] = getSCR(self,varargin)
            valid = [];
            data = NaN(9*3,length(self.ids));
            for sc = 1:length(self.ids)
                try
                    if ~isempty(varargin)
                        data(:,sc) = self.subject{sc}.scr.ledalab_summary(varargin{:});
                    else
                        data(:,sc) = self.subject{sc}.scr.ledalab_summary; %takes default timewindow defined in SCR object
                    end
                    valid = [valid sc];
                catch
                    cprintf([1 0 0],'SCR failed for subject %03d...\n',sc);
                end
                
            end
            out.y   = data'; % comes already nanzscored from SCR object
            out.x   = repmat([-135:45:180 NaN]',3,size(out.y,1))';
            out.ids = self.ids;
        end
        
        function [out, tags] = parameterMat(self)
            tags = {'kappa_cond' ...
                'kappa_test' ...
                'kapp_SI'...
                'mu_cond'...
                'mu_test'...
                'fwhm_cond'...
                'fwhm_test'...
                'fwhm_SI'...
                };
            out = [];
            for s = 1:length(self.ids);
                a = self.subject{s}.get_fit('rating',3).params;
                b = self.subject{s}.get_fit('rating',4).params;
                out = [out; a(2) b(2) b(2)-a(2) a(3) b(3) vM2FWHM(a(2)) vM2FWHM(b(2))  vM2FWHM(a(2))-vM2FWHM(b(2))];
            end
        end
        function PlotRatingFit(self,subject)%this is a single subject option...
            if ~isempty(self.tunings.rate)
                i    =  find(self.ids == subject);
                ave  = mean(reshape(self.tunings.rate{3}.y(i,:),2,8));
                x    = mean(reshape(self.tunings.rate{3}.x(i,:),2,8));
                x_HD = linspace(min(x),max(x),1000);
                h    = figure(100);clf
                subplot(1,2,1)
                title(sprintf('Sub: %i, Likelihood: %03g (p = %5.5g)',subject,self.tunings.rate{3}.singlesubject{i}.Likelihood,self.tunings.rate{3}.singlesubject{i}.pval));
                hold on;
                plot(x_HD,self.tunings.rate{3}.singlesubject{i}.fitfun(x_HD,self.tunings.rate{3}.singlesubject{i}.Est),'ro','linewidth',3);
                plot(x,ave, 'b','linewidth', 3);
                ylabel('Cond')
                drawnow;
                grid on;
                
                subplot(1,2,2)
                ave  = mean(reshape(self.tunings.rate{4}.y(i,:),2,8));
                title(sprintf('CSP: %i, Likelihood: %03g (p = %5.5g)',self.subject{i}.csp,self.tunings.rate{4}.singlesubject{i}.Likelihood,self.tunings.rate{4}.singlesubject{i}.pval));
                hold on;
                plot(x_HD,self.tunings.rate{4}.singlesubject{i}.fitfun(x_HD,self.tunings.rate{4}.singlesubject{i}.Est),'ro','linewidth',3);
                
                plot(x,ave, 'b','linewidth', 3);
                ylabel('Test')
                EqualizeSubPlotYlim(h);
                drawnow;
                grid on;
                pause
            else
                fprintf('No tuning object found here yet...\n');
                
            end
        end
        %%
        function [rating] = PlotRatings(self,runs) % uses imagesc, shows single and group ratings
            hvfigure;
            trun = length(runs);
            crun = 0;
            for run = runs(:)'%for each run make a subplot column
                crun    = crun + 1;
                %
                subplot(2,trun,crun);
                rating  = self.Ratings(run);%collect group ratings
                imagesc(rating.y,[0 10]);%thincolorbar('vert');%single subject data
                set(gca,'xticklabel',{'CS+' 'CS-'},'xtick',[4 8],'fontsize',20,'yticklabel',{''});
                colormap hot
                %
                subplot(2,trun,crun+trun);
                [y x] = hist(rating.y);
                y     = y./repmat(sum(y),size(y,1),1)*100;%make it a percentage
                imagesc(rating.x(1,:),x,y,[0 75]);axis xy;
                thincolorbar('vert');
                hold on
                h     = errorbar(mean(rating.x),mean(rating.y),std(rating.y),'g-');
                axis xy;
                set(gca,'xticklabel',{'CS+' 'CS-'},'xtick',[0 180]);
                hold off;
            end
        end
        function fit_results = PlotRatingResults(self)%plots conditioning and test, in the usual bar colors. With GroupFit Gauss/Mises Curve visible
            %%
            ratings = self.getRatings(2:4);
            f=figure;
            for n = 1:3
                subplot(1,3,n);
                h = bar(-135:45:180,mean(ratings(:,:,n),1));
                SetFearGenBarColors(h);
                hold on;
                errorbar(-135:45:180,mean(ratings(:,:,n),1),std(ratings(:,:,n),0,1)./sqrt(size(ratings,1)),'k.','LineWidth',2);
                xlim([-160 200]);
                box off
                set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'});
                x = linspace(-135,180,100);
                fprintf('Estimating Groupfit, phase %g.\n',n)
                out = self.Ratings(n+1);
                t = Tuning(out);
                self.tunings.rate{n+1} = t;
                t.GroupFit(8);
                fit_results(n) = t.groupfit;
                figure(1);
                plot(x,t.groupfit.fitfun(x,t.groupfit.Est(:,1:end-1)) ,'k','linewidth',2);
                axis square
            end
            subplot(1,3,1);title('Base');
            subplot(1,3,2);title('Cond');
            subplot(1,3,3);title('Test');
            s = supertitle(sprintf('Ratings results from n = %g subjects.',length(self.ids)));
            set(s,'FontSize',14);
            EqualizeSubPlotYlim(gcf);
        end
        function PlotSCRResults(self)
            out = self.getSCR;
            M = mean(out.y);
            S = std(out.y);
            f=figure;
            %
            for ph = 1:3
                ind = ph-1;
                subplot(1,3,ph);
                h = bar(-135:45:180,M([1:8]+8*ind));SetFearGenBarColors(h);
                hold on;
                errorbar(-135:45:180,M([1:8]+8*ind),S([1:8]+8*ind)./sqrt(length(self.ids)),'k.','LineWidth',2);
                xlim([-180 225]);
                box off
                set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'});
                hold off
                axis square
            end
            subplot(1,3,1)
            t=title('Base');set(t,'FontSize',14);
            subplot(1,3,2)
            t=title('Cond');set(t,'FontSize',14);
            subplot(1,3,3)
            t=title('Test');set(t,'FontSize',14);
        end
        function fit_results = PlotResults(self) %plots both ratings and scr
            ratings = self.getRatings(2:4);
            scr = self.getSCR;
            f=figure(1);
            for n = 1:3
                subplot(2,3,n);
                h = bar(-135:45:180,mean(ratings(:,:,n)));
                SetFearGenBarColors(h);
                hold on;
                errorbar(-135:45:180,mean(ratings(:,:,n)),std(ratings(:,:,n))./sqrt(size(ratings,1)),'k.','LineWidth',2);
                xlim([-160 200]);
                box off
                set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'});
                x = linspace(-135,180,100);
                fprintf('Estimating Groupfit, phase %g.\n',n)
                out = self.Ratings(n+1);
                t = Tuning(out);
                self.tunings.rate{n+1} = t;
                t.GroupFit(8);
                fit_results(n) = t.groupfit;
                figure(1);
                plot(x,t.groupfit.fitfun(x,t.groupfit.Est(:,1:end-1)) ,'k','linewidth',2);
                ylim([0 10]);
                axis square
                %SCR for this phase
                %correct for nans resulting from invalid scr
                scr.y = scr.y(~isnan(scr.y(:,1)),:);
                scr.x = scr.x(~isnan(scr.y(:,1)),:);
                scr.ids = scr.ids(~isnan(scr.y(:,1)));
                M = mean(scr.y);
                S = std(scr.y);
                figure(1);
                ind = n-1;
                subplot(2,3,n+3);
                h = bar(-135:45:180,M([1:8]+8*ind));SetFearGenBarColors(h);
                hold on;
                errorbar(-135:45:180,M([1:8]+8*ind),S([1:8]+8*ind)./sqrt(length(self.ids)),'k.','LineWidth',2);
                xlim([-180 225]);
                box off
                set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'});
                if n~=2 %cond can't be fitted
                    dummy.y = scr.y(:,[1:8]+8*ind);
                    dummy.x = scr.x(:,[1:8]+8*ind);
                    dummy.ids = self.ids;
                    t = Tuning(dummy);
                    t.GroupFit(8);
                    fit_results(n+3) = t.groupfit;
                    figure(1);
                    plot(x,t.groupfit.fitfun(x,t.groupfit.Est(:,1:end-1)) ,'k','linewidth',2);
                end
                hold off
                axis square
            end
            subplot(2,3,1);title('Base');
            subplot(2,3,2);title('Cond');
            subplot(2,3,3);title('Test');
            s = supertitle(sprintf('Results from n = %g / %g (Rating/SCR) subjects.',length(self.ids),sum(~isnan(scr.y(:,1)))));
            set(s,'FontSize',14);
        end
        
        %%
        function [scr] = getSCRbars(self,run)
            %will collect the SCR tunings from single subjects
            scr.y = [];
            scr.x = [];
            scr.ids = [];
            for s = 1:length(self.subject)
                if ~isempty(self.subject{s})
                    dummy = self.subject{s}.GetSubSCRbars(run);
                    if ~isempty(dummy)
                        scr.y   = [scr.y; dummy.y];
                        scr.x   = [scr.x; dummy.x];
                        scr.ids = [scr.ids; self.ids(s)];
                    end
                end
            end
        end
        function [scr] = getSCRgraphs(self,run,cond)
            %will collect the SCR response over time from single subjects
            if nargin < 3
                cond = 1:8;
            end
            scr.y = [];
            scr.conds = [];
            scr.phase = [];
            scr.ids = [];
            for s = 1:length(self.subject)
                if ~isempty(self.subject{s})
                    dummy = self.subject{s}.GetSubSCRgraphs(run,cond);
                    if ~isempty(dummy)
                        scr.y   = cat(3,scr.y,dummy.y(1:800,:));
                        scr.conds   = [scr.conds; dummy.conds];
                        scr.phase   = [scr.phase; dummy.phase];
                        scr.ids = [scr.ids; self.ids(s)];
                    end
                end
            end
        end
        function [rating] = Ratings(self,run)
            %will collect the ratings from single subjects
            rating.y  = [];
            rating.x  = [];
            rating.ids = [];
            c = 0;
            for s = 1:length(self.subject)
                if ~isempty(self.subject{s})
                    dummy = self.subject{s}.get_rating(run);
                    if ~isempty(dummy)
                        c = c+1;
                        if self.mean_correction
                            dummy.y_mean = dummy.y_mean-mean(dummy.y_mean);
                            dummy.y      = dummy.y - mean(dummy.y(:));
                        end
                        rating.y   = [rating.y ; dummy.y(:)'];
                        rating.x   = [rating.x ; dummy.x(:)'];
                        rating.ids  = [rating.ids; self.ids(s)];
                    end
                end
            end
        end
        function [info]  = getInfo(self,varargin)
            % please enter which info you need as input string.
            % gender is 1 = male, 2 = female,
            % bdnf   is 1 = G, 2 = A (where A is the mutation)
            % CTQ    is 1 = yes, trauma, 0 = no trauma
            % LTE    is 1 = yes, 0 = no.
            
            if isempty(varargin)
                prompt = 'Which groupinfo do you want to collect? gender, bdnf, CTQ or LTE? enter as string! \n';
                x = input(prompt);
            else
                x = varargin{:};
                ind = strcmp(self.subject{1}.groupinfo.tags,x);
                info = [];
                for ns = 1:length(self.subject)
                    info = [info self.subject{ns}.groupinfo.groups(ind)];
                end
            end
        end
       
    end
end