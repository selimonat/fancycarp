classdef Group < Project
    properties
        subject
    end
    
    methods
        %%
        function group = Group(subjects)
            c = 0;
            for s = subjects
                c = c+1;
                dummy            = Subject(s);
                group.subject{c} = dummy;
            end
        end
        function csps = getcsp(self)
            csps = [];
            for s = 1:length(self.subject)
                csps = [csps self.subject{s}.csp];
            end
        end
        %%
        function ModelRatings(self,run)            
            T = Tuning(self.Ratings(run));
            keyboard
        end
        %%
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
        function [rating] = RatingsRaw(self,run,varargin)
            %will collect the ratings from single subjects 
            rating.y = [];
            rating.x = [];
            c = 0;
            for s = 1:length(self.subject)
                if ~isempty(self.subject{s})
                    dummy = self.subject{s}.GetRating(run,varargin{:});
                    if ~isempty(dummy)
                        c = c+1;
                        rating.y   = [rating.y ; dummy.y_mean];
                        rating.x   = [rating.x ; mean(dummy.x)];
                    end
                end
            end
        end
        %%
        function [rating] = RatingsDemeaned(self,run,varargin)
            %will collect the ratings from single subjects 
            rating.y = [];
            rating.x = [];
            c = 0;
            for s = 1:length(self.subject)
                if ~isempty(self.subject{s})
                    dummy = self.subject{s}.GetRating(run,varargin{:});
                    if ~isempty(dummy)
                        c = c+1;
                        rating.y   = [rating.y ; dummy.y_mean-mean(dummy.y_mean)];
                        rating.x   = [rating.x ; mean(dummy.x)];
                    end
                end
            end
        end
    end
end
