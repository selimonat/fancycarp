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
        end
        %%
        function getPMF(self)
            c = 0;
            for s = 1:length(self.subject)
                c = c + 1;
                self.pmf.csp_before_alpha(c,1) = self.subject{s}.pmf.params1(1,1);
                self.pmf.csp_after_alpha(c,1)  = self.subject{s}.pmf.params1(3,1);
                self.pmf.csp_before_beta(c,1)  = self.subject{s}.pmf.params1(1,2);
                self.pmf.csp_after_beta(c,1)   = self.subject{s}.pmf.params1(3,2);
                %
                self.pmf.csn_before_alpha(c,1) = self.subject{s}.pmf.params1(2,1);
                self.pmf.csn_after_alpha(c,1)  = self.subject{s}.pmf.params1(4,1);
                self.pmf.csn_before_beta(c,1)  = self.subject{s}.pmf.params1(2,2);
                self.pmf.csn_after_beta(c,1)   = self.subject{s}.pmf.params1(4,2);
            end
        end
        function out = getSCR(self,varargin)
            data = NaN(8*3,length(self.ids));
            for sc = 1:length(self.ids)
                if ~isempty(varargin)
                    data(:,sc) = self.subject{sc}.getSCR(varargin{:});
                else
                    data(:,sc) = self.subject{sc}.getSCR; %take default timewindow defined in SCR object
                end
            end
            out = nanzscore(data);
        end
        
        function plotPMFbars(self)
            means     = reshape(mean(self.pmf.params1(:,1,:),3),2,2);%compute the mean
            stds      = reshape(std(self.pmf.params1(:,1,:),0,3),2,2);
            sem       = stds/sqrt(length(self.ids));
            
            fig=figure;
            [h,e] = barwitherr(sem,means);
            set(gca,'XTickLabel',{'CS+','CS-'})
            set(e,'LineWidth',1.5)
            set(h(1), 'FaceColor','r')
            set(h(2), 'FaceColor',[143/255 0 0 ])
            ylim([20 80])
            ylabel('threshold \alpha (degrees)')
            legend('before','after','orientation','horizontal','location','southoutside')
        end
        function [out, tags] = parameterMat(self)
            tags = {'csp_before_alpha' 'csp_after_alpha' 'csn_before_alpha' 'csn_after_alpha' ...
                'csp_before_beta' 'csp_after_beta' 'csn_before_beta' 'csn_after_beta' ...
                'csp_improvmt' 'csn_improvmnt' ...
                'csp_imprvmtn_cted' ...
                'kappa_cond' ...
                'kappa_test' ...
                'kapp_SI'...
                'mu_cond'...
                'mu_test'...
                'initial_alpha'...
                'fwhm_cond'...
                'fwhm_test'...
                'fwhm_SI'...
                };
            out = [];
            for s = 1:length(self.ids);
                out = [out;self.subject{s}.parameterMat];
            end
        end
        function PlotRatingFit(self,subject)
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
        function [rating] = PlotRatings(self,runs)
            hvfigure;
            trun = length(runs);
            crun = 0;
            for run = runs(:)'%for each run make a subplot column
                crun    = crun + 1;
                %
                subplot(2,trun,crun);
                rating  = self.Ratings(run);%collect group ratings
                imagesc(rating.y,[0 10]);thincolorbar('vert');%single subject data
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
        function PlotRatingResults(self)%plots conditioning and test, in the usual bar colors. With GroupFit Gauss/Mises Curve visible
            %%
            f=figure;
            subplot(1,2,1);
            h = bar(unique(self.tunings.rate{3}.x(1,:)),self.tunings.rate{3}.y_mean);SetFearGenBarColors(h);
            hold on;
            errorbar(unique(self.tunings.rate{3}.x(1,:)),self.tunings.rate{3}.y_mean,self.tunings.rate{3}.y_std./sqrt(length(self.ids)),'k.','LineWidth',2);
            xlim([-160 200]);
            box off
            set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'});
            x = linspace(self.tunings.rate{3}.x(1,1),self.tunings.rate{3}.x(1,end),100);
            plot(x ,  self.tunings.rate{3}.groupfit.fitfun( x,self.tunings.rate{3}.groupfit.Est(:,1:end-1)) ,'k--','linewidth',2);
            %             plot(x ,  self.tunings.rate{3}.singlesubject{1}.fitfun( x,mean(self.tunings.rate{3}.params(:,1:end-1))) ,'k--','linewidth',1);
            hold off
            %             set(gca,'fontsize',14);
            axis square
            t=title('Conditioning');set(t,'FontSize',14);
            %
            subplot(1,2,2);
            h = bar(unique(self.tunings.rate{4}.x(1,:)),self.tunings.rate{4}.y_mean);SetFearGenBarColors(h);hold on;
            errorbar(unique(self.tunings.rate{4}.x(1,:)),self.tunings.rate{4}.y_mean,self.tunings.rate{4}.y_std./sqrt(length(self.ids)),'k.','LineWidth',2);
            EqualizeSubPlotYlim(gcf);
            box off
            xlim([-160 200]);
            set(gca,'xtick',[0 180],'xticklabel',{'CS+' 'CS-'});
            x = linspace(self.tunings.rate{4}.x(1,1),self.tunings.rate{4}.x(1,end),100);
            %             plot(x ,  self.tunings.rate{4}.singlesubject{1}.fitfun( x,mean(self.tunings.rate{4}.params(:,1:end-1))) ,'k','linewidth',1);
            plot(x ,  self.tunings.rate{4}.groupfit.fitfun( x,self.tunings.rate{4}.groupfit.Est(:,1:end-1)) ,'k','linewidth',2);
            x = linspace(self.tunings.rate{3}.x(1,1),self.tunings.rate{3}.x(1,end),100);
            % % %            plot(x ,  self.tunings.rate{3}.singlesubject{1}.fitfun( x,mean(self.tunings.rate{3}.params(:,1:end-1))) ,'k--','linewidth',1);
            %             plot(x , self.tunings.rate{3}.groupfit.fitfun( x,self.tunings.rate{3}.groupfit.Est(:,1:end-1)) ,'k--','linewidth',2);
            %             set(gca,'fontsize',14);
            axis square
            t=title('Test');set(t,'FontSize',14);
            annotation(f,'textbox',[0.78 0.65 0.1 0.1],'String',['SI = ' num2str(nanmean(self.SI))],'FitBoxToText','off','LineStyle','none');
            hold off
        end
        function PlotSCR(self)
            out = self.getSCR;
            M = mean(out,2);
            S = std(out,0,2);
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
            subplot(1,3,1)
            t=title('Cond');set(t,'FontSize',14);
            subplot(1,3,1)
            t=title('Test');set(t,'FontSize',14);
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
                    dummy = self.subject{s}.GetRating(run);
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
        
    end
end
