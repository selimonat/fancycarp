classdef Group < Project
    properties (Hidden,Constant)
        mean_correction = 0;%decides if mean correction should be applied
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
            for s = subjects
                fprintf('subject: %03d\n',s)
                c = c+1;
                dummy            = Subject(s);
                group.subject{c} = dummy;
                group.ids = subjects;
                group.getPMF;
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
            self.tunings{run} = Tuning(self.RatingsDemeaned(run));%create a tuning object for the RUN for ratings.
            self.tunings{run}.SingleSubjectFit(funtype);%call fit method from the tuning object
        end

        function getSI(self,funtype)
            self.ModelRatings(3,funtype);
            self.ModelRatings(4,funtype);
            self.sigma_cond = [];
            self.sigma_test = [];
            for s = 1:length(self.subject)
                self.SI = [self.SI; self.tunings{3}.singlesubject{s}.Est(:,2) - self.tunings{4}.singlesubject{s}.Est(:,2)];%take the diff of sigma parameters.
                self.sigma_cond = [self.sigma_cond; self.tunings{3}.singlesubject{s}.Est(:,2)];
                self.sigma_test = [self.sigma_test; self.tunings{4}.singlesubject{s}.Est(:,2)];
            end
        end
        %%
        function getPMF(self)
            self.pmf.params1   = [];
            self.pmf.subject_alpha = [];
            self.pmf.subject_beta  = [];
            for s = 1:length(self.subject)
                self.pmf.params1 = cat(3,self.pmf.params1,self.subject{s}.pmf.params1);%concatenate all subjects (third dim)
                self.pmf.subject_alpha =  [self.pmf.subject_alpha; self.subject{s}.pmf.subject_alpha];%mean alpha for CS+/CS- before exp
                self.pmf.subject_beta  =  [self.pmf.subject_beta; self.subject{s}.pmf.subject_beta];%mean beta for CS+/CS- before exp
            end
        end
        
        function out = parameterMat(self)
            %1/subject alpha
            %subject beta
            %
            %3/CS+ alpha pre
            %CS- alpha pre
            %CS+ alpha post
            %CS- alpha post
            %
            %7/CS+ beta pre
            %CS- beta pre
            %CS+ beta post
            %CS- beta post
            %
            %11/CS+ alpha improvement
            %CS+ beta improvement
            %CS- alpha improvement
            %CS- beta improvement
            %
            %15/sigma cond
            %sigma test
            %17/SI
            %SInorm
            
            out = [self.pmf.subject_alpha,...
                self.pmf.subject_beta,...
                squeeze(self.pmf.params1(1,1,:)),...
                squeeze(self.pmf.params1(2,1,:)),...
                squeeze(self.pmf.params1(3,1,:)),...
                squeeze(self.pmf.params1(4,1,:)),...
                squeeze(self.pmf.params1(1,2,:)),...
                squeeze(self.pmf.params1(2,2,:)),...
                squeeze(self.pmf.params1(3,2,:)),...
                squeeze(self.pmf.params1(4,2,:)),...
                (squeeze(self.pmf.params1(1,1,:))-squeeze(self.pmf.params1(3,1,:))),...
                (squeeze(self.pmf.params1(3,2,:))-squeeze(self.pmf.params1(1,2,:))),...
                (squeeze(self.pmf.params1(2,1,:))-squeeze(self.pmf.params1(4,1,:))),...
                (squeeze(self.pmf.params1(4,2,:))-squeeze(self.pmf.params1(2,2,:))),...
                self.sigma_cond,...
                self.sigma_test,...
                self.SI];
        end
        
        function PlotRatingFit(self,subject)
            i    =  find(self.ids == subject);
            x_HD = linspace(min(self.tunings{3}.x(1,:)),max(self.tunings{3}.x(1,:)),100);
            h    = figure(100);clf
            
            subplot(1,2,1)
            title(sprintf('Likelihood: %03g (p = %5.5g)',self.tunings{3}.singlesubject{i}.Likelihood,self.tunings{3}.singlesubject{i}.pval));
            plot(x_HD,self.tunings{3}.singlesubject{i}.fitfun(x_HD,self.tunings{3}.singlesubject{i}.Est),'ro','linewidth',3);
            hold on;
            plot(self.tunings{3}.x(i,:),self.tunings{3}.y(i,:), 'b','linewidth', 3);
            ylabel('Cond')
            drawnow;
            grid on;
            
            subplot(1,2,2)
            title(sprintf('Likelihood: %03g (p = %5.5g)',self.tunings{4}.singlesubject{i}.Likelihood,self.tunings{4}.singlesubject{i}.pval));
            plot(x_HD,self.tunings{4}.singlesubject{i}.fitfun(x_HD,self.tunings{4}.singlesubject{i}.Est),'ro','linewidth',3);
            hold on;
            plot(self.tunings{3}.x(i,:),self.tunings{4}.y(i,:), 'b','linewidth', 3);
            ylabel('Test')
            EqualizeSubPlotYlim(h);
            s = supertitle(sprintf('Rating Fits Subject %03d',subject),1);
            set(s,'FontSize',14);
            drawnow;
            grid on;
            
        end
                
        
        function [rating] = PlotRatingsDemeaned(self,runs,varargin)
            hvfigure;
            trun = length(runs);
            crun = 0;
            for run = runs(:)'%for each run make a subplot column
                crun    = crun + 1;
                %
                subplot(2,trun,crun);                
                rating  = self.RatingsDemeaned(run,varargin{:});%collect group ratings                
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
        %%
       
        function [rating] = PlotRatings(self,runs,varargin)
            hvfigure;
            trun = length(runs);
            crun = 0;
            for run = runs(:)'%for each run make a subplot column
                crun    = crun + 1;
                %
                subplot(2,trun,crun);                
                rating  = self.RatingsRaw(run,varargin{:});%collect group ratings                
                imagesc(rating.y,[0 10]);thincolorbar('vert');%single subject data
                set(gca,'xticklabel',{'CS+' 'CS-'},'xtick',[4 8],'fontsize',20,'yticklabel',{''});
                colormap hot
                %
                subplot(2,trun,crun+trun);
                [y x] = hist(rating.y,1:10);
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
        %%
        function [rating] = Ratings(self,run,varargin)
            %will collect the ratings from single subjects 
            rating.y = [];
            rating.x = [];
            c = 0;
            for s = 1:length(self.subject)
                if ~isempty(self.subject{s})
                    dummy = self.subject{s}.GetRating(run,varargin{:});
                    if ~isempty(dummy)
                        c = c+1;
                        if self.mean_correction
                            dummy.y_mean = dummy.y_mean-mean(dummy.y_mean);
                        end
                        rating.y   = [rating.y ; dummy.y_mean];
                        rating.x   = [rating.x ; mean(dummy.x)];
                    end
                end
            end
        end
        
    end
end
