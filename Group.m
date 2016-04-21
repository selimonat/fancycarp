classdef Group < Project
    properties (Hidden,Constant)
        mean_correction = 0;%decides if mean correction should be applied
        align_tunings   = 1;%should ratings be aligned to CS+ face
    end
    properties
        subject
        ids
        csps
        table_pmf
        ratings
        total_subjects
        tunings
        SI        
        sigma_cond
        sigma_test
        SCR_ampl
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
            end
        end
        %%
        function out = get.ratings(self)                        
            %will collect group ratings as a matrix in out(run).y,
            %out(run).x etc.
            for s = 1:self.total_subjects
                for fields = fieldnames(self.subject{s}.ratings)'
                    out.(fields{1})(s,:) = self.subject{s}.ratings.(fields{1})(:);
                end
            end
        end
        %%        
        function out = get.table_pmf(self)            
            %returns pmf parameters as a table
            out = [];
            for s = self.subject                
                out = [out ;s{1}.pmf_parameters];
            end            
        end
        %%        
        function csps = get.csps(self)   
            %returns the csp face for all the group... In the future one
            %could make a get_subject_field function
            for s = 1:length(self.subject)
                csps(s) = self.subject{s}.csp;
            end
        end   
        %% 
        function out = get.total_subjects(self)
            out = length(self.ids);
        end
        %%
        function plot_ratings(self)
            %%
            [y x]  =  GetSubplotNumber(self.total_subjects);
            for ns = 1:self.total_subjects
                subplot(y,x,ns);
                self.feargen_plot(self.ratings.y_mean(1,:,ns));
                box off;
                title(sprintf('s: %d, cs+: %d',self.ids(ns),self.subject{ns}.csp),self.font_style{:}); 
            end
            EqualizeSubPlotYlim(gcf);
            supertitle(self.path_project,1);
        end        
        %%
        function model_ratings(self,run,funtype)
            %create a tuning object and fits FUNTYPE to it.
            self.tunings.rate{run} = Tuning(self.ratings(run));%create a tuning object for the RUN for ratings.
            self.tunings.rate{run}.SingleSubjectFit(funtype);%call fit method from the tuning object
        end               
         %%
        function feargen_plot(self,data)
            %elementary function to make feargen plots
            h = bar(data,1);            
            self.set_feargen_colors(h,2:9);
            set(gca,'xtick',[4 8],'xticklabel',{'cs+' 'cs-'},'xgrid','on',self.font_style{:});
            axis tight;
        end
         
        %%
        function ModelSCR(self,run,funtype)
            %create a tuning object and fits FUNTYPE to it.
            self.tunings.scr = Tuning(self.getSCRs(run));%create a tuning object for the RUN for SCRS.
%             self.tunings.scr.SingleSubjectFit(funtype);%call fit method from the tuning object
        end
        function getSCRtunings(self,run,funtype)
            self.ModelSCR(run,funtype);
        end
        
        function [out] = getSCRmeans(self,phase)
            for n = 1:length(self.ids)
                    ind = self.subject{n}.scr.findphase(phase);
                    self.subject{n}.scr.cut(ind);
                    self.subject{n}.scr.run_ledalab;
                    out(n,:) = mean(self.subject{n}.scr.ledalab.mean(1:800,:));
            end
        end
        function [out] = loadmises(self)
            a = load(sprintf('%smidlevel%smisesmat.mat',self.path_project,filesep));
            out = a.misesmat(self.ids,:);
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
        function [out labels] = parameterMat(self)
            labels = {'csp_before_alpha' 'csp_after_alpha' 'csn_before_alpha' 'csn_after_alpha' ...
                      'csp_before_beta' 'csp_after_beta' 'csn_before_beta' 'csn_after_beta' ...                     
                      'csp_improvmt' 'csn_improvmnt' ...
                      'csp_imprvmtn_cted' ...
                      'rating_cond' ... 
                      'rating_test' ... 
                      'SI'...
                      'SCR ampl'};
            out = [self.pmf.csp_before_alpha,...
                   self.pmf.csp_after_alpha,...              
                   self.pmf.csn_before_alpha,...
                   self.pmf.csn_after_alpha,...   
                   self.pmf.csp_before_beta,...
                   self.pmf.csp_after_beta,...              
                   self.pmf.csn_before_beta,...
                   self.pmf.csn_after_beta,... 
                   self.pmf.csp_before_alpha - self.pmf.csp_after_alpha,...
                   self.pmf.csn_before_alpha - self.pmf.csn_after_alpha,...                   
                   (self.pmf.csp_before_alpha-self.pmf.csp_after_alpha)-(self.pmf.csn_before_alpha-self.pmf.csn_after_alpha),...
                   self.sigma_cond,...
                   self.sigma_test,...
                   self.SI,...
                   self.SCR_ampl];
            
               try
               if strcmp(self.tunings.rate{3}.singlesubject{1}.funname,'vonmisses_mobile')
                   
                   for s = 1:size(out,1)
                       out(s,15) = self.tunings.rate{3}.singlesubject{s}.Est(3);
                       out(s,16) = self.tunings.rate{4}.singlesubject{s}.Est(3);
                   end
                   out(:,14) = out(:,13) - out(:,12);
                   labels = [labels(1:14) { 'mu_cond' 'mu_test'}];
               end
               end
        end
        function [out labels] = misesMat(self)
            labels = {'rating_cond' ...
                'rating_test' ...
                'SI'...
                'Mu_cond'...
                'Mu_test'};
            for s = 1:length(self.ids)
                mu_cond(s) = self.tunings.rate{3}.singlesubject{s}.Est(3);
                mu_test(s) = self.tunings.rate{4}.singlesubject{s}.Est(3);
            end
            out = [self.sigma_cond,...
                self.sigma_test,...
                self.SI,...
                mu_cond',...
                mu_test'];
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
        %%
        function [scr]    = getSCRs(self,run)
            %will collect the ratings from single subjects 
            scr.y = [];
            scr.x = [];
            scr.ids = [];
            for s = 1:length(self.subject)
                if ~isempty(self.subject{s})
                  dummy = self.subject{s}.GetSubSCR(run);
                    if ~isempty(dummy)
                        scr.y   = [scr.y; dummy.y];
                        scr.x   = [scr.x; dummy.x];
                        scr.ids = [scr.ids; self.ids(s)];
                    end
                end
            end
        end
        
        
    end
end
