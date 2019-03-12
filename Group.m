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
        fit_results
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
        function relief = get_relief(self,varargin)
            
            force = 1;
            phases2collect = [1 2 5];
            filename = fullfile(Project.path_project,'midlevel','ratings',sprintf('relief_N%02d_phases_%s.mat',self.total_subjects,sprintf('%d',phases2collect)));
            fprintf('Filename for called relief ratings is: \n %s\n',filename)
            
            
            if ~exist(filename) || force == 1
                fprintf('Ratings not found for this group, collecting them!\n');
                relief = nan(self.total_subjects,10,length(phases2collect));
                phc = 0;
                for ph = phases2collect
                    phc = phc + 1;
                    for ns = 1:self.total_subjects;
                        relief(ns,:,phc) = self.subject{ns}.get_reliefmeans(ph,self.allconds,varargin{1});
                    end
                end
                fprintf('Saving...');
                save(filename,'relief');
                fprintf('done\n');
            else
                fprintf('Ratings for this group found in midlevel, loading them.\n');
                load(filename);
            end
        end
        function pain = get_pain(self,varargin)
            for ns = 1:self.total_subjects
                pain(ns,:,:) = self.subject{ns}.get_pain;
            end
            if nargin > 1
                runs = varargin{1};
                %mean of both test runs
                if ismember(5,runs)
                    pain  = cat(3,pain,nanmean(pain(:,:,3:4),3));
                end
                % we can kick the nan bc only short runs selected
                if all(runs<3) % 3 ratings
                    pain = pain(:,1:3,runs);
                else
                    pain = pain(:,:,runs);
                end
            end
            
        end
        function thresh = get_threshold(self)
            thresh = nan(1,self.total_subjects);
            for ns = 1:self.total_subjects
                thresh(ns) = self.subject{ns}.get_threshold;
            end
        end
        
        function plot_grouprelief_bar(self)
            plottype = 'bar';
            
            phases2plot = [1 2 5];
            dofit = 1;
            fitmethod = self.selected_fitfun;
            titles = {'Base','Cond','Test1','Test2','Test'};
            f=figure;
            f.Position = [73 410 1763 476];
            clf;
            spc = 0;
            for ph = phases2plot(:)'
                spc = spc + 1;
                
                for ns = 1:self.total_subjects;
                    relief(ns,:) = self.subject{ns}.get_reliefmeans(ph,self.allconds,'raw');
                end
                M = nanmean(relief);
                S = nanstd(relief)./sqrt(sum(~isnan(relief(:,end)))); %this way we get correct SEM even if missing data from one subject.
                subplot(1,length(phases2plot),spc)
                self.plot_bar(self.plotconds,M,S)
                hold on
                axis square
                set(gca,'XTickLabels',{'' '' '' 'CS+' '' '' '' 'CS-' 'UCS' 't0'},'XTickLabelRotation',45,'FontSize',12);
                xlim([min(self.plotconds)-40 max(self.plotconds+40)])
                title(titles{ph},'FontSize',14)
                if dofit == 1
                    if ph ~= 2
                        data.x = -135:45:180;
                        data.y = M(1:8);
                        data.ids = self.ids;
                        t = Tuning(data);
                        t.SingleSubjectFit(fitmethod);
                        self.fit_results{ph} = t.fit_results;
                        self.tunings{ph} = t;
                        
                        
                        if 10^-self.fit_results{ph}.pval < .05
                            plot(self.fit_results{ph}.x_HD,self.fit_results{ph}.fit_HD,'k-','LineWidth',3)
                        else
                            txt= text(min(xlim)+20,M(1)+S(1).*1.2,sprintf('p = %4.3f',self.fit_results{ph}.pval));
                            set(txt, 'rotation', 90)
                            l=line([-135,180],repmat(mean(self.tunings{ph}.y),1,2));
                            set(l,'LineWidth',3,'Color','k')
                        end
                    end
                end
                if spc == 1
                    ylabel('Relief rating M +/- SEM')
                    
                end
            end
            st = supertitle(sprintf('N = %02d subs',self.total_subjects));
            set(st,'FontSize',16);
            EqualizeSubPlotYlim(gcf);
            
            %             savefig(gcf,savefigf)
            %                 export_fig(gcf,savebmp)
            %                 export_fig(gcf,savepng,'-transparent')
        end
        function [relief tb tt] = plot_grouprelief_pirate(self)
            relief = self.get_relief('zscore');
            
            %             a=reshape(relief,size(relief,1),30);
            %             reliefZ = nanzscore(a')';
            %             relief = reshape(reliefZ,size(relief,1),10,3);
            
            fs = 14;
            vio = 0;
            baryn = 1;
            erb = 0;
            ml = 0;
            dotyn = 0;
            
            figure;
            clf;
            cols10 = self.GetFearGenColors;
            
            xcenters = -135:45:225;
            
            subplot(1,3,1)
            D = relief(:,1:9,1);
            COL = cols10;
            pirateplot(xcenters,D,'color',repmat([.3 .3 .3],9,1),'violin',vio,'bar',baryn,'errorbar',erb,'meanline',ml,'dots',dotyn);
            hold on;
            errorbar(xcenters,nanmean(D),nanstd(D)./sqrt(35),'.','Color',[.5 .5 .5],'LineWIdth',2,'MarkerFaceColor',[.5 .5 .5])
            title('Baseline','FontSize',fs)
            set(gca,'XTick',xcenters([4 8]),'XTickLabel',{'CS+','CS-'},'XTickLabelRotation',45,'FontSize',fs,'Ydir','reverse');
            ylabel('Relief ratings [zscore]')
            
            %             Publication_NiceTicks(gca,5);
            
            subplot(1,3,2)
            %             D = relief(:,:,2);
            %             COL = cols10;
            D = relief(:,[9 8],2);
            COL = cols10([9 8],:);
            xc = xcenters([8 9]);
            pirateplot(xc,D,'color',COL,'violin',vio,'bar',baryn,'errorbar',erb,'meanline',ml,'dots',dotyn);
            hold on;
            for n = 1:2
                errorbar(xc(n),nanmean(D(:,n)),nanstd(D(:,n))./sqrt(35),'.','Color',COL(n,:),'LineWidth',2,'MarkerFaceColor',COL(n,:))
            end
            title('Conditioning','FontSize',fs)
            xlim([135 270])
            set(gca,'XTick',xc,'XTickLabel',{'UCS','CS-'},'XTickLabelRotation',45,'FontSize',fs,'Ydir','reverse');
            
            subplot(1,3,3);
            D = relief(:,1:8,3);
            COL = cols10;
            pirateplot(xcenters,D,'color',COL,'violin',vio,'bar',baryn,'errorbar',erb,'meanline',ml,'dots',dotyn);
            hold on;
            for n = 1:8
                errorbar(xcenters(n),nanmean(D(:,n)),nanstd(D(:,n))./sqrt(35),'.','Color',COL(n,:),'LineWidth',2,'MarkerFaceColor',COL(n,:))
            end
            title('Test','FontSize',fs)
            set(gca,'XTick',xcenters([4 8 9]),'XTickLabel',{'CS+','CS-','UCS'},'XTickLabelRotation',45,'FontSize',fs,'Ydir','reverse');
            
            base.x = repmat(-135:45:180,self.total_subjects,1);
            base.y = relief(:,1:8,1);
            base.ids = self.ids;
            test.x = repmat(-135:45:180,self.total_subjects,1);
            test.y = relief(:,1:8,3);
            test.ids = self.ids;
            
            tb = Tuning(base);
            tt = Tuning(test);
            tb.GroupFit(8);
            tt.GroupFit(8);
            linecol = [.3 .3 .3];
            
            subplot(1,3,1);
            hold on;
            if 10.^-tb.groupfit.pval < .001
                plot(tb.groupfit.x_HD,tb.groupfit.fit_HD,'k','LineWidth',3,'Color',linecol,'LineStyle',':');
            else
                plot([-135 180],repmat(mean(tb.y_mean),1,2),'k','LineWidth',3,'Color',linecol);
                
            end
            
            subplot(1,3,3);
            hold on;
            if 10.^-tt.groupfit.pval < .001
                plot(tt.groupfit.x_HD,tt.groupfit.fit_HD,'k','LineWidth',3,'Color',linecol,'LineStyle','-');
            else
                plot([-180 225],repmat(mean(tt.y_mean),1,2),'k','LineWidth',3,'Color',linecol);
                
            end
            for n=2:3;
                subplot(1,3,n);
                hold on;
                %                 ylim([-2.5 2.5]);
                %                 set(gca,'YTick',-2:2);
                Publication_RemoveYaxis(gca)
            end
            
            set(gcf,'Color','w')
            EqualizeSubPlotYlim(gcf);
        end
        function [relief, tb, tt]= plot_grouprelief_pirate_BT(self,varargin)
            
            lwt = 4;%linewidth tuning
            lweb = 3; %linewidth errorbar;
            fs = 15;
            vio = 0;
            baryn = 1;
            erb = 0;
            ml = 0;
            dotyn = 0;
            
            COL = Project.GetFearGenColors;
            xlev = -135:45:180;
            relief = self.get_relief('zscore');
            %             yticki = [-1.5 0 1.5];
            %
            
            if nargin > 1
                fighand = varargin{1};
                set(fighand,'Color','w')
                spn = varargin{2};
                sp1 = varargin{3};
            else
                fighand = figure;
                set(fighand,'Color','w')
                spn = 2;
                sp1 = 1;
            end
            
            subplot(1,spn,sp1);
            D = relief(:,1:8,1);
            pirateplot(xlev,D,'color',repmat([.3 .3 .3],8,1),'violin',vio,'bar',baryn,'errorbar',erb,'meanline',ml,'dots',dotyn);
            hold on;
            for n = 1:8
                errorbar(xlev(n),nanmean(D(:,n)),nanstd(D(:,n))./sqrt(35),'.','Color',[.3 .3 .3],'LineWidth',lweb,'MarkerFaceColor',COL(n,:))
            end
            %             set(gca,'XTick',[0 180],'XTickLabel',{'CS+','CS-'},'FontSize',fs,'YTick',yticki,'Ydir','reverse');
            set(gca,'XTick',[0 180],'XTickLabel',{'CS+','CS-'},'FontSize',fs,'Ydir','reverse');
            ylabel('pain relief [VAS, z-score]','FontSize',fs+1)
            %             ylim([-2 2])
            hold on;
            axis square
            title('Baseline','FontSize',fs+1)
            
            
            
            subplot(1,spn,sp1+1);
            %             yticki = [-1.5 0 1.5];
            D = relief(:,1:8,3);
            pirateplot(xlev,D,'violin',vio,'bar',baryn,'errorbar',erb,'meanline',ml,'dots',dotyn);
            %             set(gca,'XTick',[0 180],'XTickLabel',{'CS+','CS-'},'FontSize',fs,'YTick',yticki,'Ydir','reverse');
            set(gca,'XTick',[0 180],'XTickLabel',{'CS+','CS-'},'FontSize',fs,'Ydir','reverse');
            hold on;
            for n = 1:8
                errorbar(xlev(n),nanmean(D(:,n)),nanstd(D(:,n))./sqrt(35),'.','Color',COL(n,:),'LineWidth',lweb,'MarkerFaceColor',COL(n,:))
            end
            
            hold on;
            axis square
            title('Test','FontSize',fs+1)
            
            %%
            base.x = repmat(-135:45:180,self.total_subjects,1);
            base.y = relief(:,1:8,1);
            base.ids = self.ids;
            test.x = repmat(-135:45:180,self.total_subjects,1);
            test.y = relief(:,1:8,3);
            test.ids = self.ids;
            
            tb = Tuning(base);
            tt = Tuning(test);
            tb.GroupFit(8);
            tt.GroupFit(8);
            linecol = [.3 .3 .3];
            
            subplot(1,spn,sp1);
            hold on;
            if 10.^-tb.groupfit.pval < .001
                plot(tb.groupfit.x_HD,tb.groupfit.fit_HD,'k','LineWidth',lwt,'Color',linecol,'LineStyle',':');
            else
                plot([-135 180],repmat(mean(tb.y_mean),1,2),'k','LineWidth',lwt,'Color',linecol);
                
            end
            
            subplot(1,spn,sp1+1);
            hold on;
            if 10.^-tt.groupfit.pval < .001
                plot(tt.groupfit.x_HD,tt.groupfit.fit_HD,'k','LineWidth',lwt,'Color',linecol,'LineStyle','-');
            else
                plot([-180 225],repmat(mean(tt.y_mean),1,2),'k','LineWidth',lwt,'Color',linecol);
                
            end
            %             EqualizeSubPlotYlim(gcf);
            
            yticki{1} = -.2:.1:.2;
            yticki{2} = [];
            spnn=0;
            for n = sp1:(sp1+1)
                spnn=spnn+1;
                subplot(1,spn,n);
                ylim([-.301 .243])
                set(gca,'YTick',yticki{spnn});
            end
        end
        function [relief] = plot_relief_cond_placebo(self,varargin)
            
            lweb = 3; %linewidth errorbar;
            fs = 15;
            vio = 0;
            baryn = 1;
            erb = 0;
            ml = 0;
            dotyn = 0;
            
            COL = Project.GetFearGenColors;
            xlev = [-135:45:180 225];
            relief = self.get_relief('raw');
            %
            if nargin > 1
                fighand = varargin{1};
                set(fighand,'Color','w')
                spn = varargin{2};
                sp1 = varargin{3};
            else
                fighand = figure;
                set(fighand,'Color','w')
                spn = 2;
                sp1 = 1;
            end
            set(fighand,'Color','w')
            subplot(1,spn,sp1);
            D = relief(:,1:9,2);
            pirateplot(xlev(:,[8 9]),D(:,[8 9]),'color',COL([8 9],:),'violin',vio,'bar',baryn,'errorbar',erb,'meanline',ml,'dots',dotyn);
            hold on;
            for n = 1:9
                errorbar(xlev(n),nanmean(D(:,n)),nanstd(D(:,n))./sqrt(35),'.','Color',COL(n,:),'LineWidth',lweb,'MarkerFaceColor',COL(n,:))
            end
            set(gca,'XTick',[0 180 225],'XTickLabel',{'','CS-','UCS'},'XTickLabelRotation',45,'FontSize',fs,'Ydir','reverse');
            ylabel('pain relief [VAS]','FontSize',fs)
            
            hold on;
            %             axis square
            title('Conditioning','FontSize',fs+1)
            
            
            
            subplot(1,spn,sp1+1);
            D = relief(:,1:9,3);
            pirateplot(xlev([7 8 9]),D(:,[4 8 9]),'color',COL([4 8 9],:),'violin',vio,'bar',baryn,'errorbar',erb,'meanline',ml,'dots',dotyn);
            set(gca,'XTick',xlev([7 8 9]),'XTickLabel',{'CS+','CS-','UCS'},'XTickLabelRotation',45,'FontSize',fs,'Ydir','reverse');
            hold on;
            for n = [4 8 9]
                xlev(4) = xlev(7);
                errorbar(xlev(n),nanmean(D(:,n)),nanstd(D(:,n))./sqrt(35),'.','Color',COL(n,:),'LineWidth',lweb,'MarkerFaceColor',COL(n,:))
            end
            
            hold on;
            %             axis square
            title('Test','FontSize',fs)
            
            %             EqualizeSubPlotYlim(gcf);
            ylimmi = [0 90];
            yticki{1} = 0:20:80;
            yticki{2} = [];
            xlimmi{1} = [90+22.5  270+22.5];
            xlimmi{2} = [90 270];
            spp = 0;
            for n = sp1:(sp1+1)
                spp = spp+1;
                sh(n)=subplot(1,spn,n);
                ylim(ylimmi);
                set(gca,'YTick',yticki{spp});
                xlim(xlimmi{spp});
            end
            [mean(relief(:,10,2)) mean(relief(:,9,2));sqrt([std(relief(:,10,2)) std(relief(:,9,2))])]
            [mean(relief(:,4,3)) mean(relief(:,8,3)) mean(relief(:,9,3));sqrt([std(relief(:,4,3)) std(relief(:,8,3)) std(relief(:,9,3))])]
        end
        
        function relief = plot_painNrelief_pirate(self)
            relief = self.get_relief('raw');
            pain   = self.get_pain([1 2 5]);
            pain   = squeeze(nanmean(pain,2));
            vio = 0;
            barf = 0;
            erbar = 0;
            fs = 14;
            colorpain = [184,70,64]./255;%zeros(3,3);
            colorrelief = [76,73,127]./255;%repmat([0 0 1],3,1);
            xpain = [1 5 9];
            xrelief = xpain+1;
            
            R = squeeze(nanmean(relief,2));
            
            
            plot(10,120,'o','MarkerSize',10,'MarkerFaceColor',colorpain,'MarkerEdgeColor',colorpain)
            hold on
            plot(10,120,'o','MarkerSize',10,'MarkerFaceColor',colorrelief,'MarkerEdgeColor',colorrelief)
            lg = legend('Pain','Relief');
            legend boxoff
            set(lg,'FontSize',14)
            
            pirateplot(xpain,pain,'color',repmat(colorpain,3,1),'violin',vio,'bar',barf,'errorbar',erbar);
            pirateplot(xrelief,-R,'color',repmat(colorrelief,3,1),'violin',vio,'bar',barf,'errorbar',erbar);
            
            hold on;
            %             b=bar(xpain,mean(pain));set(b,'FaceColor','none','EdgeColor',colorpain(1,:),'LineWidth',2,'Barwidth',.6);
            b=bar(xpain,mean(pain));set(b,'FaceColor',colorpain(1,:),'FaceAlpha',.6,'EdgeColor','none','LineWidth',2,'Barwidth',.6);
            errorbar(xpain,mean(pain),std(pain)./sqrt(self.total_subjects),'.','Color',colorpain,'LineWidth',2);
            set(gca,'XTickLabel',{'Base','Cond','Test'},'FontSize',fs);
            ylabel('VAS Rating');
            hold on
            
            hold on
            %             b=bar(xrelief,-mean(R));set(b,'FaceColor','none','EdgeColor',colorrelief,'LineWidth',2,'Barwidth',.6);
            b=bar(xrelief,-mean(R));set(b,'FaceColor',colorrelief,'FaceAlpha',.6,'EdgeColor','none','LineWidth',2,'Barwidth',.6);
            errorbar(xrelief,-mean(R),std(R)./sqrt(self.total_subjects),'.','Color',colorrelief,'LineWidth',2);
            set(gca,'XTickLabel',{'Base','Cond','Test'},'FontSize',fs);
            ylim([-100 100])
            set(gca,'YTick',-100:20:100,'YTickLabels',{'100','80','60','40','20','0','20','40','60','80','100'});
            l=line(xlim,[0 0]);set(l,'Color','k')
            set(gcf,'Color','w')
        end
        
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
        
        function plot_pain_vs_relief(self,ph)
            
            cmap = jet(self.total_subjects);
            %% plot interpl pain vs relief single trials, single subs as subplot
            forcexlim = 0;
            figure(100);clf;
            for ns = 1:self.total_subjects;
                [~, pain, relief] = self.subject{ns}.get_pmod_PR(ph);
                rss(ns) = corr(pain',relief);
                subplot(6,6,ns);
                plot(pain,relief','o','Color',cmap(ns,:),'MarkerFaceColor',cmap(ns,:))
                ls = lsline;set(ls,'Color','k','LineWidth',2)
                axis square
                box off
                if forcexlim
                    xlim([0 100]);
                end
                title(sprintf('s%02d, r=%04.2f',self.ids(ns),rss(ns)))
            end
            st = supertitle(sprintf('pain vs relief (phase %d)',ph));set(st,'FontSize',14,'Position',[.05 .3 0]);
            % ks density
            figure(101);clf;
            for ns = 1:self.total_subjects;
                [~, pain, relief] = self.subject{ns}.get_pmod_PR(ph);
                subplot(6,6,ns);
                % ksdensity
                %                    if forcexlim
                %                 xlim([0 100]);
                %                 end
                %                 [f,xi] = ksdensity(relief);
                %                 plot(xi,f,'Color',cmap(ns,:))
                %                 UCS = self.subject{ns}.paradigm{ph}.presentation.ucs;
                %                 MxTemp(1) = mean(relief(UCS));
                %                 MxTemp(2) = mean(relief(~UCS));
                %                 l1 = line(repmat(MxTemp(1),1,2),ylim);set(l1,'Color','r','LineWidth',2,'LineStyle',':')
                %                 l2 = line(repmat(MxTemp(2),1,2),ylim);set(l2,'Color','c','LineWidth',2)
                % UCS and other trials as lines, single ratings as dots
                
                UCS = self.subject{ns}.paradigm{ph}.presentation.ucs;
                MxTemp(1) = mean(relief(UCS));
                MxTemp(2) = mean(relief(~UCS));
                hold on;
                l1 = line([1 length(relief)],repmat(MxTemp(1),1,2),ylim);set(l1,'Color','r','LineWidth',2,'LineStyle',':');
                l2 = line([1 length(relief)],repmat(MxTemp(2),1,2),ylim);set(l2,'Color','k','LineWidth',2);
                scatter(1:length(relief),relief,20,[0 0 0],'filled')
                ylim([0 100]);
                xlim([0 length(UCS)])%length trials
                axis square
                box off
                title(sprintf('s%02d, r=%04.2f',self.ids(ns),rss(ns)))
            end
            st = supertitle(sprintf('Relief ratings per sub (phase %d)',ph));set(st,'FontSize',14,'Position',[.05 .3 0]);
            %% plot group, 4 phases, one dot per sub (Mpain vs Mrelief)
            figure(102);clf;
            tstr = {'Base','Cond','Test1','Test2'};
            for ph = 1:4
                clear pain
                clear relief
                for ns = 1:self.total_subjects
                    [~, pain(ns,:), relief(ns,:)] = self.subject{ns}.get_pmod_PR(ph);
                    rho(ns,ph) = corr(squeeze(pain(ns,:))',squeeze(relief(ns,:))');
                end
                subplot(1,5,ph);
                scatter(mean(pain(:,:),2),mean(relief(:,:),2),30,cmap,'filled')
                ls=lsline;set(ls,'LineWidth',2,'Color','k')
                axis square
                box off
                xlabel('Mpain')
                ylabel('Mrelief')
                xlim([0 100])
                ylim([0 100])
                title(tstr{ph});
                avecorr(ph) = fisherz_inverse(mean(fisherz(rho(:,ph))));
                text(50,90,sprintf('mean corr\n r = %04.3f',avecorr(ph)))
                [rr(ph),pp(ph)] = corr(mean(pain,2),mean(relief,2));
                text(10,10,sprintf('r_{mean} = %04.3f',rr(ph)))
            end
            subplot(1,5,5);
            for ns = 1:self.total_subjects
                for ph = 1:4
                    plot(ph,rho(ns,ph),'o','Color',cmap(ns,:),'MarkerFaceColor',cmap(ns,:))
                    hold on
                end
            end
            xlim([0 5])
            hold on;
            for n = 1:4
                l=line([n-.3 n+.3],repmat(avecorr(n),1,2));set(l,'LineWidth',2','Color','k')
            end
            axis square
            box off
            set(gca,'XTick',1:4,'XTicklabel',tstr,'XTickLabelRotation',45)
            xlabel('phase')
            ylabel('correlation r')
            title('subjects'' corr per phase')
            
            rhoz=fisherz(rho);
            for n=1:4;[rrz(n) ppz(n)] = ttest(rhoz(:,n),0);end
        end
        function plot_temp_pain_relief(self)
            vio   = 0;
            barf  = 1;
            erbar = 1;
            
            fs = 16;
            colorpain   = [184,70,64]./255;%[100 100 100]./255;%[184,70,64]./255;%zeros(3,3);
            colorrelief = [76,73,127]./255;%[170 170 170 ]./255;%[76,73,127]./255;%repmat([0 0 1],3,1);
            colortemp   = [30 30 30]./255;
            
            xpain = [1 5 9];
            xrelief = xpain+1;
            
            for ns = 1:self.total_subjects
                for ph = 1:4
                    tempr(ns,ph) = self.subject{ns}.paradigm{ph}.presentation.pain.tonic(1);
                end
            end
            
            tempr(:,3) = mean(tempr(:,3:4),2);
            tempr(:,4) = [];
            
            
            figure(1005);
            clf
            subplot(1,2,1)
            pirateplot(1:3,tempr,'color',repmat(colortemp,size(tempr,2),1),'violin',vio,'bar',barf,'errorbar',erbar);
            hold on;
            %% resulting pain and relief
            relief = self.get_relief('raw');
            pain   = self.get_pain([1 2 5]);
            pain   = squeeze(nanmean(pain,2));
            
            R = squeeze(nanmean(relief,2));
            subplot(1,2,2)
            plot(10,120,'o','MarkerSize',10,'MarkerFaceColor',colorpain,'MarkerEdgeColor',colorpain)
            hold on
            plot(10,120,'o','MarkerSize',10,'MarkerFaceColor',colorrelief,'MarkerEdgeColor',colorrelief)
            lg = legend('Pain','Relief');
            
            set(lg,'Box','Off','Location','northeast','FontSize',14)
            
            pirateplot(xpain,pain,'color',repmat(colorpain,3,1),'violin',vio,'bar',barf,'errorbar',erbar);
            pirateplot(xrelief,-R,'color',repmat(colorrelief,3,1),'violin',vio,'bar',barf,'errorbar',erbar);
            
            hold on;
            %             b=bar(xpain,mean(pain));set(b,'FaceColor','none','EdgeColor',colorpain(1,:),'LineWidth',2,'Barwidth',.6);
            %             errorbar(xpain,mean(pain),std(pain)./sqrt(self.total_subjects),'.','Color',colorpain,'LineWidth',2);
            set(gca,'XTickLabel',{'Base','Cond','Test'},'FontSize',fs);
            ylabel('VAS Rating','FontSize',fs+2);
            hold on
            
            hold on
            %             b=bar(xrelief,-mean(R));set(b,'FaceColor','none','EdgeColor',colorrelief,'LineWidth',2,'Barwidth',.6);
            %             errorbar(xrelief,-mean(R),std(R)./sqrt(self.total_subjects),'.','Color',colorrelief,'LineWidth',2);
            for n = 1:2
                subplot(1,2,n)
                set(gca,'XTickLabel',{'Base','Cond','Test'},'FontSize',fs);
                if n == 1
                    xlim([0 4])
                    ylabel(sprintf('Temperature [%cC]',char(176)),'FontSize',fs+2)
                    ylim([20 40])
                    set(gca,'YTick',20:5:50);
                elseif n == 2
                    ylim([-100 100])
                    %                     set(gca,'YTick',-100:20:100,'YTickLabels',{'100','80','60','40','20','0','20','40','60','80','100'});
                    
                    set(gca,'YTick',-100:50:100,'YTickLabels',{'100','50','0','50','100'});
                    l=line(xlim,[0 0]);set(l,'Color','k')
                end
            end
            set(gcf,'Color','w')
            
            
            
            
            
            
            
        end
        %%
        function selectedface = get_selectedface(self)
            selectedface = nan(self.total_subjects,1);
            for ns = 1:self.total_subjects
                selectedface(ns) = self.subject{ns}.selectedface;
            end
            figure;
            pirateplot(1:8,histc(selectedface,-135:45:180)','color',self.GetFearGenColors,'violin',0,'errorbar',0,'CI',0,'dots',0)
            set(gca,'XTick',[4 8],'XTickLabel',{'CS+' 'CS-'},'YTick',0:2:10,'FontSize',14)
            set(gcf,'Color','w')
            ylabel('N','FontSize',14);
            title('face selected as most effective','FontSize',14)
        end
        %%
        function model_ratings(self,run,funtype)
            %will fit to ratings from RUN the FUNTYPE and cache the result
            %in the midlevel folder.
            T = [];%future tuning object
            filename               = sprintf('%smidlevel/Tunings_Run_%03d_FunType_%03d_N%s.mat',self.path_project,run,funtype,sprintf('%s\b\b',sprintf('%ito',self.ids([1 end]))));
            if exist(filename) == 0
                %create a tuning object and fits FUNTYPE to it.
                T  = Tuning(self.ratings(run));%create a tuning object for the RUN for ratings.
                T.SingleSubjectFit(funtype);%call fit method from the tuning object
                save(filename,'T');
            else
                fprintf('Will load the tuning parameters from the cache:\n%s\n',filename);
                load(filename);
            end
            %get the relevant data from the tuning object
            self.fit_results = T.fit_results;
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
                'csp_improvmt' 'csn_improvmnt' ...2
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
    methods %(mri))
        function Run1stlevel(self,modelnum,varargin)
            if nargin > 2
                phases = varargin{1};
            end
            %% wrapper Creating the Onsets, Fitting HRF model to model, and computing contrasts
            fprintf('Preparing CondFiles\n')
            for ns = 1:self.total_subjects
                for ph = phases(:)'
                    if self.ids(ns)==15 && ph == 4
                        continue
                    else
                        self.subject{ns}.PrepareCondFile(ph,modelnum);
                    end
                end
            end
            for ns = 1:self.total_subjects
                for ph = phases(:)'
                    if self.ids(ns)~=15 && ph ==3
                        ph = [3 4];
                    end
                    fprintf('Fitting Subject %02d, run %d.\n',self.ids(ns),ph)
                    self.subject{ns}.FitHRF(ph,modelnum);
                    fprintf('ConImages for Subject %02d, run %d.\n',self.ids(ns),ph)
                    self.subject{ns}.Con1stLevel(ph(1),modelnum);
                end
            end
        end
        function Fit2ndlevel(self,nrun,modelnum,namestring,varargin)
            %% 2ndlevel ANOVA
            dependencies = 1;
            unequalvar   = 1;
            
            clear2ndlevel = 1;
            versiontag = 0;
            prefix = 's6_wCAT_';
            versiontag = 0;
            foldersuffix = sprintf('_N%02d',self.total_subjects);
            
            if nargin == 5 %one varargin is given
                foldersuffix = varargin{1};
            elseif nargin == 6
                foldersuffix = varargin{1};
                versiontag = varargin{2};
            end
            %covariate business
            if modelnum == 12 && ~isempty(strfind(foldersuffix,'_painCOVAR'))
                [cov,covstrct] = self.get_2ndlevel_covariate(nrun,modelnum);
            else
                [cov,covstrct] = self.get_2ndlevel_covariate(nrun,0);
            end
            
            if strcmp(namestring,'CSdiff')
                cons2collect   = 1;
            elseif strcmp(namestring,'8conds')
                %go through all conditions, i.e. 8 (B/T) or 2 (C)
                if ismember(nrun,[1 3])
                    cons2collect = 2:9;
                elseif nrun == 2
                    cons2collect = 2:3;
                end
                
            elseif strcmp(namestring,'8conds_rate')
                switch nrun
                    case 1
                        cons2collect = 1:8;%WIP
                    case 2
                        cons2collect = 2:4;
                    case 3
                        cons2collect = 2:4;
                end
            elseif strcmp(namestring,'VMdVM')
                cons2collect = 1:2; %main effect was not included in 1stlevel contrasts. So 1=VM, 2=dVM
            elseif strcmp(namestring,'VMdVM_BT')
                cons2collect = 1:2;
                nrun = [1 3];
            elseif strcmp(namestring,'PMOD')
                cons2collect = 1; %on 1stlevel, only pmod was computed as con. so we just take this.
            elseif strcmp(namestring,'SBS_down_PMOD')
                cons2collect = 2;
            elseif strcmp(namestring,'SBS_up_PMOD')
                cons2collect = 5;
            elseif strcmp(namestring,'SBS_down_ME')
                cons2collect = 1;
            elseif strcmp(namestring,'SBS_up_ME')
                cons2collect = 4;
            elseif strcmp(namestring,'SBS_Plateau')
                cons2collect = 3;
            end
            
            start = tic;
            fprintf('Starting 2nd Level for model %02d, run %02d, named %s, versiontag %d with foldersuffix ''%s''...\n',modelnum,nrun(1),namestring,versiontag,foldersuffix);
            
            path2ndlevel = fullfile(self.path_second_level,sprintf('model_%02d_chrf_%01d_%s_%s%s',modelnum,versiontag,namestring,self.nrun2phase{nrun(1)},foldersuffix));
            
            
            if exist(path2ndlevel) && (clear2ndlevel==1);
                system(sprintf('rm -fr %s*',strrep(path2ndlevel,'//','/')));
            end%this is AG style.
            if ~exist(path2ndlevel)
                mkdir(path2ndlevel);
            end
            clear matlabbatch
            
            load(self.subject{1}.path_spmmat(nrun,modelnum));
            % collect all subjects' con images for every cond
            c = 0;
            for runloop = nrun(:)'
                for ncon = cons2collect(:)'
                    c = c+1;
                    fprintf('Run %d: Getting con_%04d, %s from sub ',runloop,ncon, SPM.xCon(ncon).name)
                    clear files
                    for ns = 1:numel(self.ids)
                        files(ns,:) = strrep(self.subject{ns}.path_con(runloop,modelnum,prefix,ncon),'sub004',sprintf('sub%03d',self.ids(ns)));
                        fprintf('%d..',self.ids(ns))
                    end
                    fprintf('.done. (N=%02d).\n',ns)
                    matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(c).scans = cellstr(files); %one cond at a time, but all subs
                end
            end
            
            % specify rest of the model
            matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(path2ndlevel);
            matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = dependencies;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = unequalvar;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            if cov == 1
                matlabbatch{1}.spm.stats.factorial_design.cov = covstrct;
            else
                matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            end
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = -Inf;
            
            groupmask = [self.path_groupmeans sprintf('/thr20_ave_wCAT_s3_ss_data_N%02d.nii',self.total_subjects)];
            if exist(groupmask)
                matlabbatch{1}.spm.stats.factorial_design.masking.em = {groupmask};
            else
                matlabbatch{1}.spm.stats.factorial_design.masking.em = {[self.path_groupmeans '/ave_wCAT_s3_ss_data.nii']};
            end
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {[path2ndlevel filesep 'SPM.mat']};
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
            done = toc(start);
            fprintf('\n\nDone estimating 2nd Level for HRF model %d, called %s %s, version %d, N = %02d subs in %05.2f mins. (Simple one-way anova)\n',modelnum,namestring,foldersuffix,versiontag,length(self.ids),done./60)
            fprintf('Output folder is: %s\n',path2ndlevel)
            
        end
        function Con2ndlevel(self,nrun,modelnum,namestring,varargin)
            
            foldersuffix = sprintf('_N%02d',self.total_subjects);
            versiontag = 0;
            
            if nargin == 5 %one varargin is given
                foldersuffix = varargin{1};
            elseif nargin == 6
                foldersuffix = varargin{1};
                versiontag = num2str(varargin{2});
            end
            
            nF = 0;
            nT = 0;
            n  = 0;
            
            path_spmmat = fullfile(self.path_second_level,sprintf('model_%02d_chrf_%01d_%s_%s%s',modelnum,versiontag,namestring,self.nrun2phase{nrun},foldersuffix),'SPM.mat');
            
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(path_spmmat);
            
            if strcmp(namestring,'CSdiff')
                n = n + 1; nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'CSP>CSN';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                matlabbatch{1}.spm.stats.con.delete = 0;
                
            elseif strcmp(namestring,'8conds')
                if ismember(nrun,[1 3])
                    nconds = 8;
                else
                    nconds = 2;
                end
                
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_F_8conds';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'main_allconds';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                if ismember(nrun,[1 3])
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 0 0 1 0 0 0 -1];
                else
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [1 -1];
                end
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'CSP>CSN';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                if ismember(nrun,[1 3])
                    [VM, dVM] = self.compute_VM(-135:45:180,1,1,.001);
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = VM;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'VMtuning';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = dVM;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'dVMtuning';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [-repmat(1/7,1,3) 1 -repmat(1/7,1,3)];
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'CSP>rest';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                end
                
            elseif strcmp(namestring,'VMdVM')
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'pp';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(2);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                n  = n + 1;
                nF = nF + 1;
                vec = eye(2); vec(logical(eye(2))) = [1 -1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'pn';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = vec;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'VM';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [1 0];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'dVM';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
            elseif strcmp(namestring,'VMdVM_BT')
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_eye(4)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(4);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n  = n + 1;
                nF = nF + 1;
                vec = eye(4); vec(logical(eye(4))) = [1 -1 1 -1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'main_VMvsdVM';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = vec;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                n  = n + 1;
                nF = nF + 1;
                vec = eye(4); vec(logical(eye(4))) = -[1 -1 1 -1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'main_dVMvsVM';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = vec;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                n  = n + 1;
                nF = nF + 1;
                vec = eye(4); vec(logical(eye(4))) = [-1 -1 1 1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'test>base';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = vec;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
            elseif strcmp(namestring,'PMOD') && isempty(strfind(foldersuffix,'COVAR'))
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F_any';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'pmod';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            elseif strcmp(namestring,'PMOD') && ~isempty(strfind(foldersuffix,'COVAR'))
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F_pmod';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [0 1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'pmod';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            elseif any(strfind(namestring,'SBS'))
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 't';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            end
            
            
            if nT > 0
                for tc = 1:nT
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.name =      [matlabbatch{1}.spm.stats.con.consess{nF+tc}.tcon.name '_neg'];
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.weights =   -matlabbatch{1}.spm.stats.con.consess{nF+tc}.tcon.weights;
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.sessrep =   'none';
                    nT = nT + 1;
                end
            end
            
            matlabbatch{1}.spm.stats.con.delete = 1;
            
            spm_jobman('run',matlabbatch);
            
            ntotal = nT + nF;
            fprintf('Done creating %d 2ndlevel contrasts (%d F, %d T) for model %d, run %s, modelname %s %s.\n',ntotal,nF,nT,modelnum,self.nrun2phase{nrun},namestring,foldersuffix)
            if nF > 0
                for nnF = 1:nF
                    disp(['(F) ' matlabbatch{1}.spm.stats.con.consess{nnF}.fcon.name])
                end
            end
            for nnT = 1:nT
                disp(['(T) ' matlabbatch{1}.spm.stats.con.consess{nF+nnT}.tcon.name])
            end
        end
        function Fit2ndlevel_FIR(self,nrun,modelnum,order,namestring,varargin)
            start = tic;
            %% 2ndlevel 8conds ANOVA
            clear2ndlevel = 1;
            
            dependencies = 0;
            unequalvar   = 0;
            
            versiontag = 0;
            prefix = 's6_wCAT_';
            namestring1stlevel = '10conds';
            foldersuffix = sprintf('_N%02d',self.total_subjects);
            bins2take = 1:self.orderfir;
            
            
            if nargin == 6 %one varargin is given
                foldersuffix = varargin{1}; %here you can pass whatever extension you want, like test, or N39 or so
            elseif nargin == 7
                foldersuffix =  varargin{1}; %here you can pass whatever extension you want, like test, or N39 or so
                versiontag = varargin{2};   %probably never needed, better name then with foldersuffix, easier to remember and document
            elseif nargin > 7
                fprintf('Too many inputs. Please debug.')
                keyboard;
            end
            
            fprintf('Starting 2nd Level for FIR model %02d, run %02d, named %s, version %d with foldersuffix ''%s''...\n',modelnum,nrun,namestring,versiontag,foldersuffix);
            
            path2ndlevel = fullfile(self.path_second_level,'FIR',sprintf('model_%02d_FIR_%02d_%s_b%02dto%02d_%s%s',modelnum,versiontag,namestring,bins2take(1),bins2take(end),self.nrun2phase{nrun},foldersuffix));
            
            if exist(path2ndlevel) && (clear2ndlevel==1)
                system(sprintf('rm -fr %s*',strrep(path2ndlevel,'//','/')));
            end%this is AG style.
            if ~exist(path2ndlevel)
                mkdir(path2ndlevel);
            end
            
            
            % information, which contrasts to collect.
            %
            % for models without PMOD:
            % on firstlevel, we have bin 1-14 cond 1, bin 1-14 cond 2, etc,
            % then bin 1-14 CSdiff
            %
            % For VMdVM, con images on 1stlevel were skipped for main
            % effect. So just go for contrast 1:Npmod
            % For PMOD from ratings, same. so just 1:Npmod(which is =1)
            %
            % Here: defining chunks of contrasts to go for (+ N bins is done by loop later).
            
            if strcmp(namestring,'8conds')
                switch nrun
                    case 1
                        conds2collect = 1:8;
                    case 2
                        conds2collect = 1:2;
                    case 3
                        conds2collect = 1:8;
                end
            elseif strcmp(namestring,'9conds')
                switch nrun
                    case 1
                        conds2collect = 1:8;
                    case 2
                        conds2collect = 1:2;
                    case 3
                        conds2collect = [1:8 10];
                end
            elseif strcmp(namestring,'CSPCSN')
                switch nrun
                    case 1
                        conds2collect = [4 8];
                    case 2
                        conds2collect = [1 2];
                    case 3
                        conds2collect = [4 8];
                end
            elseif strcmp(namestring,'CSdiff')
                switch nrun
                    case 1
                        conds2collect = 9;
                    case 2
                        conds2collect = 3;
                    case 3
                        conds2collect = 9;
                end
            elseif strcmp(namestring,'VMdVM')
                conds2collect = [1 2];
            elseif strcmp(namestring,'VMdVM_BT')
                conds2collect = [1 2];
                nrun = [1 3];
            elseif strcmp(namestring,'PMOD_rate')
                conds2collect = 1;
            elseif strcmp(namestring,'PMOD_PRind_both')
                conds2collect = 1:2;
            elseif strcmp(namestring,'PMOD_PRind_ME')
                conds2collect = 1;
            elseif strcmp(namestring,'PMOD_PRind_pmod')
                conds2collect = 2;
            elseif any(strfind(namestring,'win')) && any(strfind(namestring,'CSdiff'))
                load(fullfile(self.subject{1}.path_FIR(nrun,modelnum,self.orderfir,'10conds'),'SPM.mat'))
                if any(strfind(namestring,'bin4'))
                    if nrun == 2
                        cons2collect = 45;
                    elseif nrun == 3
                        cons2collect = 150;%bin4win4
                    end
                elseif any(strfind(namestring,'bin5'))
                    if nrun == 2
                        cons2collect = 48;
                    elseif nrun == 3
                        cons2collect = 160;%bin5win4
                    end
                end
                fprintf('You chose contrasts named:\n')
                SPM.xCon(cons2collect).name
                fprintf('.................................................\n')
            elseif any(strfind(namestring,'bin4win4_CSPCSN_BCT'))
                load(fullfile(self.subject{1}.path_FIR(1,modelnum,self.orderfir,'10conds'),'SPM.mat'))
                nrun = [1 2 3];
                cons2collectstruct = {[130 134], 58:59 ,[144 148]};%bin4win4 %before:  {[130 134], 43:44,[144 148]};%but con 43:44 somehow wrong.
                fprintf('You chose contrasts named:\n')
                SPM.xCon(cons2collectstruct{1}).name
                fprintf('.................................................\n')
            elseif any(strfind(namestring,'win')) && any(strfind(namestring,'UCSCSN'))
                load(fullfile(self.subject{1}.path_FIR(nrun,modelnum,self.orderfir,'10conds'),'SPM.mat'))
                if any(strfind(namestring,'bin4'))
                    if nrun == 2
                        cons2collect = 58:59;%43:44;
                    elseif nrun == 3
                        cons2collect = [149 148];%bin4win4
                    end
                elseif any(strfind(namestring,'bin5'))
                    if nrun == 2
                        cons2collect = 46:47;
                    elseif nrun == 3
                        cons2collect = [159 158];%bin5win4
                    end
                end
                fprintf('You chose contrasts named:\n')
                SPM.xCon(cons2collect).name
                fprintf('.................................................\n')
            elseif any(strfind(namestring,'win')) && any(strfind(namestring,'CSPCSN'))
                load(fullfile(self.subject{1}.path_FIR(nrun,modelnum,self.orderfir,'10conds'),'SPM.mat'))
                if any(strfind(namestring,'bin4'))
                    if nrun == 2
                        cons2collect = 43:44;
                    elseif nrun == 3
                        cons2collect = [144 148];%bin4win4
                    end
                elseif any(strfind(namestring,'bin5'))
                    if nrun == 2
                        cons2collect = 43:44;
                    elseif nrun == 3
                        cons2collect = 154:158;%bin5win4
                    end
                end
                fprintf('You chose contrasts named:\n')
                SPM.xCon(cons2collect).name
                fprintf('.................................................\n')
            elseif any(strfind(namestring,'win')) && any(strfind(namestring,'allconds'))
                load(fullfile(self.subject{1}.path_FIR(nrun,modelnum,self.orderfir,'10conds'),'SPM.mat'))
                if any(strfind(namestring,'bin4'))
                    if nrun == 2
                        cons2collect = 43:44;
                    elseif nrun == 3
                        cons2collect = 141:149;%bin4win4
                    elseif nrun == 1
                        cons2collect = 127:134;%bin4win4
                    end
                elseif any(strfind(namestring,'bin5'))
                    if nrun == 2
                        cons2collect = 46:47;
                    elseif nrun == 3
                        cons2collect = 151:159;%bin4win4
                    end
                end
                fprintf('You chose contrasts named:\n')
                SPM.xCon(cons2collect).name
                fprintf('.................................................\n')
            elseif any(strfind(namestring,'win')) && any(strfind(namestring,'8conds_BT'))
                load(fullfile(self.subject{1}.path_FIR(nrun,modelnum,self.orderfir,'10conds'),'SPM.mat'))
                if any(strfind(namestring,'bin4win4'))
                    nrun = [1 3];
                    cons2collectstruct = {127:134, [],141:148};%bin4win4
                elseif any(strfind(namestring,'bin2win4'))
                    nrun = [1 3];
                    cons2collectstruct = {164:171, [],202:209};%bin2win4
                    %...
                elseif any(strfind(namestring,'bin4win6'))
                    nrun = [1 3];
                    cons2collectstruct = {173:180, [],212:219};%bin2win4
                elseif any(strfind(namestring,'bin5win6'))
                    nrun = [1 3];
                    cons2collectstruct = {182:189, [],222:229};%bin2win4
                end
                fprintf('You chose contrasts named:\n')
                SPM.xCon(cons2collectstruct{1}).name
                fprintf('.................................................\n')
            elseif any(strfind(namestring,'win2')) && any(strfind(namestring,'_BT'))
                load(fullfile(self.subject{1}.path_FIR(nrun,modelnum,self.orderfir,'10conds'),'SPM.mat'))
                if any(strfind(namestring,'bin4'))
                    nrun = [1 3];
                    cons2collectstruct = {137:144, [],161:168};%bin4win2
                elseif any(strfind(namestring,'bin5'))
                    nrun = [1 3];
                    cons2collectstruct = {146:153, [],171:178};%bin5win2
                elseif any(strfind(namestring,'bin6'))
                    nrun = [1 3];
                    cons2collectstruct = {155:162, [],181:188};%bin6win2
                end
                fprintf('You chose contrasts named:\n')
                SPM.xCon(cons2collectstruct{1}).name
                fprintf('.................................................\n')
            elseif strfind(namestring,'bin4win4_9conds_BCT')
                load(fullfile(self.subject{1}.path_FIR(nrun,modelnum,self.orderfir,'10conds'),'SPM.mat'))
                if any(strfind(namestring,'bin4'))
                    nrun = [1 2 3];
                    cons2collectstruct = {127:134, 58:59, 141:149};%bin4win4 %before: {127:134, 43:44,141:149};, but 43:44 was wrong somehow
                elseif any(strfind(namestring,'bin5'))
                    %...
                end
                fprintf('You chose contrasts named (excerpt):\n')
                SPM.xCon(cons2collectstruct{3}).name
                fprintf('.................................................\n')
            elseif any(strfind(namestring,'bin4win4')) && any(strfind(namestring,'Gauss_BT'))
                load(fullfile(self.subject{1}.path_FIR(nrun,modelnum,self.orderfir,'10conds'),'SPM.mat'))
                
                nrun = [1 3];
                cons2collectstruct = {136, [],191};%bin4win2
                
                fprintf('You chose contrasts named:\n')
                SPM.xCon(cons2collectstruct{1}).name
                fprintf('.................................................\n')
                
            end
            clear matlabbatch
            
            %load(self.subject{1}.path_spmmat(nrun,modelnum)); %allows lookup of things
            % collect all subjects' con images for every cond
            if any(strfind(namestring,'win')) %here we don\t take single bins to 2ndlevel, so we loop through the CONS, not conds identified above
                bc = 0;
                for runloop = nrun(:)'
                    if numel(nrun) > 1
                        cons2collect = cons2collectstruct{runloop};
                    end
                    for con = cons2collect(:)'
                        bc = bc+1;
                        fprintf('\nCollecting con images of con %04d:\n',con)
                        clear files
                        fprintf('Factor %02d - Looping through subs: ',bc);
                        for ns = 1:self.total_subjects
                            fprintf('%02d - ',ns)
                            files(ns,:) = cellstr([self.subject{ns}.path_FIR(runloop,modelnum,order,namestring1stlevel), sprintf('%scon_%04d.nii',prefix,con)]);
                        end
                        fprintf('completo.')
                        matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(bc).scans = cellstr(files); %one bin at a time, but all subs
                        
                    end
                end
            else
                bc = 0;
                for runloop = nrun(:)'
                    for cond = conds2collect(:)'
                        fprintf('\nCollecting con images for run %01d, cond %04d:\n',runloop,cond)
                        for bin = bins2take(:)' %loop through bins, then subjects. sub01_bin1 sub02_bin1 ... sub01_bin2 sub02_bin etc..
                            ind = self.findcon_FIR(order,cond,bin);
                            fprintf('\nbin %02d, i.e. con %03d...',bin,ind)
                            bc = bc + 1;
                            clear files
                            fprintf('Looping through subs: ');
                            for ns = 1:self.total_subjects
                                fprintf('%02d - ',ns)
                                files(ns,:) = cellstr([self.subject{ns}.path_FIR(runloop,modelnum,order,namestring1stlevel), sprintf('%scon_%04d.nii',prefix,ind)]);
                            end
                            fprintf('completo.')
                            matlabbatch{1}.spm.stats.factorial_design.des.anova.icell(bc).scans = cellstr(files); %one bin at a time, but all subs
                        end
                    end
                end
            end
            
            % specify rest of the model
            matlabbatch{1}.spm.stats.factorial_design.dir = cellstr(path2ndlevel);
            matlabbatch{1}.spm.stats.factorial_design.des.anova.dept = dependencies;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.variance = unequalvar;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.gmsca = 0;
            matlabbatch{1}.spm.stats.factorial_design.des.anova.ancova = 0;
            matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = -Inf;
            groupmask = [self.path_groupmeans sprintf('/thr20_ave_wCAT_s3_ss_data_N%02d.nii',self.total_subjects)];
            if exist(groupmask)
                matlabbatch{1}.spm.stats.factorial_design.masking.em = {groupmask};
            else
                matlabbatch{1}.spm.stats.factorial_design.masking.em = {[self.path_groupmeans '/ave_wCAT_s3_ss_data.nii']};
            end
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
            
            
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {[path2ndlevel filesep 'SPM.mat']};
            %             matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
            
            spm_jobman('run',matlabbatch);
            done = toc(start);
            fprintf('\n\nDone estimating 2nd Level for FIR model %d, called %s %s, version %d, N = %02d subs in %05.2f mins. (Simple one-way anova)\n',modelnum,namestring,foldersuffix,versiontag,length(self.ids),done./60)
            fprintf('Output folder is: %s\n',path2ndlevel)
        end
        
        
        function Con2ndlevel_FIR(self,nrun,modelnum,namestring,varargin)
            deletecons = 1;
            
            versiontag = 0;
            foldersuffix = sprintf('_N%02d',self.total_subjects);
            bins2take = 1:self.orderfir;
            
            
            if nargin == 5 %one varargin is given
                foldersuffix = varargin{1}; %here you can pass whatever extension you want, like test, or N39 or so
            elseif nargin == 6
                foldersuffix =  varargin{1}; %here you can pass whatever extension you want, like test, or N39 or so
                versiontag   = varargin{2};   %probably never needed, better name then with foldersuffix, easier to remember and document
            elseif nargin > 6
                fprintf('Too many inputs. Please debug.')
                keyboard;
            end
            
            fprintf('Starting 2nd Level Contrasts for FIR model %02d, run %02d, named %s, with foldersuffix ''%s''...\n',modelnum,nrun,namestring,foldersuffix);
            
            path2ndlevel = fullfile(self.path_second_level,'FIR',sprintf('model_%02d_FIR_%02d_%s_b%02dto%02d_%s%s',modelnum,versiontag,namestring,bins2take(1),bins2take(end),self.nrun2phase{nrun},foldersuffix));
            
            if ~exist(path2ndlevel)
                fprintf('Folder not found, please debug, or run 2ndlevel estimation first.\n')
            end
            
            nF = 0;
            nT = 0;
            n  = 0;
            
            path_spmmat = fullfile(path2ndlevel,'SPM.mat');
            
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(path_spmmat);
            
            if ismember(nrun,[1 3])
                nconds = 8;
            else
                nconds = 2;
            end
            %% handy contrast tiles
            square0 = zeros(self.orderfir,self.orderfir);
            square1 = ones(self.orderfir,self.orderfir);
            eye1    = eye(self.orderfir);
            vec0    =  zeros(1,self.orderfir);
            vec1    =  ones(1,self.orderfir);
            
            %% Tunings
            gauss = spm_Npdf(1:8,4)-mean(spm_Npdf(1:8,4));
            
            
            conds = -135:45:180;
            
            amp   = 1;
            kappa = 1;
            delta = .01;
            VMlookup = zscore(Tuning.VonMises(conds,amp,kappa,0,0));
            dVMlookup = -zscore((Tuning.VonMises(conds,amp,kappa+delta,0,0)-Tuning.VonMises(conds,amp,kappa-delta,0,0))./(2*delta)); %central difference formula
            
            %% hrf contrast vecs - CAVE: only valid if onsets are on RampDown, not on Face (thus offset_Face is negative)
            defaultparam = [6 16 1 1 6 0 32];%seconds!
            prebins   = 2;
            offset_face = -1.7;%secs
            offset_rate =  6.0;%secs
            param2change = 6; %p(1) - delay of response (relative to onset),p(6) - onset {seconds}
            
            param_ramp = defaultparam; param_ramp(param2change) = param_ramp(param2change) + prebins;
            hrf_Ramp = spm_hrf(self.TR,param_ramp);
            hrf_Ramp = hrf_Ramp(1:max(bins2take))' - mean(hrf_Ramp(1:max(bins2take))');
            
            param_face = defaultparam; param_face(param2change) = param_face(param2change) + prebins + offset_face;
            hrf_Face = spm_hrf(self.TR,param_face);
            hrf_Face = hrf_Face(1:max(bins2take))' - mean(hrf_Face(1:max(bins2take))');
            
            param_rate = defaultparam; param_rate(param2change) = param_rate(param2change) + prebins + offset_rate;
            hrf_Rating = spm_hrf(self.TR,param_rate);
            hrf_Rating = hrf_Rating(1:max(bins2take))' - mean(hrf_Rating(1:max(bins2take))');
            %%
            
            if strcmp(namestring,'CSdiff')
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = sprintf('F_eye(%02d)',max(bins2take));
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(max(bins2take));
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n = n + 1; nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 't_CSP>CSN';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,max(bins2take))./max(bins2take);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                
            elseif strcmp(namestring,'CSPCSN')
                nconds = 2;
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = sprintf('F_eye(%02dx%02d)',nconds,max(bins2take));
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(nconds*max(bins2take));
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = sprintf('F_%01dxeye(%02d)',nconds,max(bins2take));
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = repmat(eye(max(bins2take)),1,nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F_CSP>CSN';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [eye(max(bins2take)) -eye(max(bins2take))];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 't_main_bothCSPCSN';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,nconds*max(bins2take))./(nconds*max(bins2take));
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 't_CSP>CSN';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [ones(1,max(bins2take)) -ones(1,max(bins2take))]./max(bins2take);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'HRF_CSPCSN_Ramp';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = repmat(hrf_Ramp,1,nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'HRF_CSPCSN_Face';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = repmat(hrf_Face,1,nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'HRF_CSPCSN_Rating';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = repmat(hrf_Rating,1,nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'HRF_CSP>CSN_Ramp';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [hrf_Ramp -hrf_Ramp];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'HRF_CSP>CSN_Face';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [hrf_Face -hrf_Face];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'HRF_CSP>CSN_Rating';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [hrf_Rating -hrf_Rating];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                
            elseif strcmp(namestring,'8conds')
                
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = sprintf('F_eye(%02dx%02d)',nconds,max(bins2take));
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(nconds*max(bins2take));
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                n  = n + 1;
                nF = nF + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = sprintf('F_%01dxeye(%02d)',nconds,max(bins2take));
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = repmat(eye(max(bins2take)),1,nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                n  = n + 1;
                nT = nT + 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'main_allconds';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,nconds);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                if ismember(nrun,[1 3])
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 0 0 1 0 0 0 -1];
                else
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [1 -1];
                end
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'CSP>CSN';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                if ismember(nrun,[1 3])
                    [VM, dVM] = self.compute_VM(-135:45:180,1,1,.001);
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = VM-mean(VM);
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'VMtuning';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = dVM-mean(dVM);
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'dVMtuning';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [-repmat(1/7,1,3) 1 -repmat(1/7,1,4)];
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'CSP>rest';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                end
            elseif strcmp(namestring,'VMdVM_BT')
                
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(4*self.orderfir);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(4*allbins)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %%both pmods test
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [zeros(2*self.orderfir,2*self.orderfir),eye(2*self.orderfir)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_T_(2*allbins)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %%both pmods test vs base
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [-eye(2*self.orderfir),eye(2*self.orderfir)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_T_vs_B';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                %%both pmods test vs base
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [-eye1 -eye1 eye1 eye1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'VM_dVM_T_vs_B';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                % only VM T vs B
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [-eye1,square0,eye1,square0];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'VM_T_vs_B';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                % only dVM T vs B
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [square0,-eye1,square0,eye1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'dVM_T_vs_B';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                %% T Tests
                %only VM in test
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [vec0 vec0 vec1 vec0];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'VM_T';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                %only dVM in test
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [vec0 vec0 vec0 vec1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'dVM_T';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
            elseif strcmp(namestring,'PMOD_rate')
                
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(self.orderfir);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(allbins)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %% T Tests
                % all bins positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = vec1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_allbins';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                % HRF ramp
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = hrf_Ramp;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_HRF_RampOn';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = hrf_Face;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_HRF_FaceOn';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = hrf_Rating;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_HRF_RateOn';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            elseif strcmp(namestring,'PMOD_PRind_both')
                
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(self.orderfir*2);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(allbins)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %eoi ME
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [eye(self.orderfir) square0];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(ME)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                %eoi PMOD
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [square0 eye(self.orderfir)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(PMOD)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                
                %% T Tests
                % all bins positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [vec1 vec1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_allbins';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                % only main effect pos
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [vec1 vec0];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_ME';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                % only PMOD pos
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [vec0 vec1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_PRind';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
            elseif strcmp(namestring,'PMOD_PRind_ME')|| strcmp(namestring,'PMOD_PRind_pmod')
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(self.orderfir);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(allbins)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %% T Tests
                % all bins positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = vec1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_allbins';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            elseif any(strfind(namestring,'win')) && any(strfind(namestring,'CSdiff'))
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = 1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %% T Tests
                % all bins positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = 1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            elseif (any(strfind(namestring,'win')) && any(strfind(namestring,'UCSCSN'))) && ~any(strfind(namestring,'BCT'))|| (any(strfind(namestring,'win')) && any(strfind(namestring,'CSPCSN'))) && ~any(strfind(namestring,'BCT'))
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(2);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(2)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [1 -1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F_UCSvsCSN';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %% T Tests
                % all bins positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [1 1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_both';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [1 0];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_UCS';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_CSN';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [1 -1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_CSdiff';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            elseif any(strfind(namestring,'win')) && any(strfind(namestring,'allconds'))
                if nrun == 3
                    %% F TESTS
                    %eoi
                    n = n + 1;
                    nF = nF+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(9);
                    matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(allconds)';
                    matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                    n = n + 1;
                    nF = nF+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [0 0 0 1 0 0 0 -1];
                    matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F_CSPvsCSN';
                    matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                    
                    %% T Tests
                    % all bins positive
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,9);
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_all';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 0 0 1];
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_CSP';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 0 0 0 0 0 0 1];
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_CSN';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [-1/7 -1/7 -1/7 1 -1/7 -1/7 -1/7 -1/7 0];
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_CSP>rest';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    
                    n = n + 1;
                    nT = nT+1;
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [repmat(-1/8,1,8) 1];
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_UCS>rest';
                    matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                end
            elseif any(strfind(namestring,'win')) && any(strfind(namestring,'8conds_BT'))
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(16);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(allconds)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                %testphase
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [eye(8) zeros(8)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Base)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %baseline
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [zeros(8) eye(8)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Test)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %baseline vs test
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [eye(8) -eye(8)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Test)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                
                %% T Tests
                % all bins positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,16);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_all';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [ones(1,8) zeros(1,8)];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_Base';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights =[zeros(1,8) ones(1,8)];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_Test';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                %% VonMises
                [VM, dVM] = self.compute_VM(-135:45:180,1,1,.001);
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [VM VM];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_VM_both';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [VM zeros(1,8)];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_VM_base';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [zeros(1,8) VM];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_VM_test';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [-VM VM];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_VM_test>base';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                %% Gauss
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [gauss gauss];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_gauss_both';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [gauss zeros(1,8)];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_gauss_base';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [zeros(1,8) gauss];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_gauss_test';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [-gauss gauss];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_Gauss_test>base';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [zeros(1,8) dVMlookup];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'dVM_test';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
            elseif any(strfind(namestring,'bin4win4_Gauss_BT'))
                %%
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(2);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_BT';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                %testphase
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [1 0];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Base)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %baseline
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [0 1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Test)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %baseline vs test
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [1 -1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Test)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                
                %% T Tests
                % all conds positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,2);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_BT';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [1 0];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_Base';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [0 1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_Test';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = [-1 1];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name = 'T_Test>Base';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
            elseif strfind(namestring,'bin4win4_9conds_BCT')
                %%
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(19);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_19_BCT';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                %base
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = ones(1,8);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Base)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %cond
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [zeros(1,8) ones(1,2)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Cond)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %test
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [zeros(1,10) ones(1,9)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Test)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %baseline vs test
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [ones(1,8) zeros(1,2) -ones(1,9)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F_Base_vs_Test';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %UCS vs noUCS
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [-ones(1,8) 1 -1 -ones(1,8) 1];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F_UCS_vs_rest';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                %% T Tests
                % all conds positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,19)./19;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_main_relief';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                % all conds positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights =  [zeros(1,8) 1 0 zeros(1,8) 1]./2;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_UCS_main';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                % all conds positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights =  [ones(1,8) 0 1 ones(1,8) 0]./sum([ones(1,8) 0 1 ones(1,8) 0]);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_noUCS_main';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                % all conds positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights =  [-ones(1,8)./17 1/2 -1./17  -ones(1,8)./17  1/2];
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_UCS>rest';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
            elseif strfind(namestring,'bin4win4_CSPCSN_BCT')
                %%
                %% F TESTS
                %eoi
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = eye(6);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_6_BCT';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                %base
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = ones(1,2);
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Base)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %cond
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [zeros(1,2) ones(1,2)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Cond)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %test
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [zeros(1,4) ones(1,2)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'eoi_(Test)';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %baseline vs test
                n = n + 1;
                nF = nF+1;
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights = [ones(1,2) zeros(1,2) -ones(1,2)];
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.name = 'F_Base_vs_Test';
                matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                
                %% T Tests
                % all conds positive
                n = n + 1;
                nT = nT+1;
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights = ones(1,6);
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    = 'T_main_relief';
                matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                
            else
                fprintf('No 2ndlevel contrasts defined here, please check code.\n')
                keyboard
                
            end
            
            nn = n; %final counter
            if nT > 0
                for tc = 1:nT
                    nn = nn + 1;
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.name =      [matlabbatch{1}.spm.stats.con.consess{nF+tc}.tcon.name '_neg'];
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.weights =   -matlabbatch{1}.spm.stats.con.consess{nF+tc}.tcon.weights;
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.sessrep =   'none';
                    nT = nT + 1;
                end
            end
            
            figure(100);
            for nFc = 1:nF
                subplot(ceil(sqrt(nF)),ceil(sqrt(nF)),nFc);
                imagesc(matlabbatch{1}.spm.stats.con.consess{nFc}.fcon.weights,[-1 1]);
                title(matlabbatch{1}.spm.stats.con.consess{nFc}.fcon.name)
            end
            figure(101);
            cc = 0;
            for nTc = (nF+1):nn
                cc = cc + 1;
                subplot(ceil(sqrt(nT)),ceil(sqrt(nT)),cc);
                imagesc(matlabbatch{1}.spm.stats.con.consess{nTc}.tcon.weights,[-1 1]);
                title(matlabbatch{1}.spm.stats.con.consess{nTc}.tcon.name)
            end
            
            
            matlabbatch{1}.spm.stats.con.delete = deletecons;
            
            spm_jobman('run',matlabbatch);
            
            ntotal = nT + nF;
            fprintf('Done creating %d 2ndlevel contrasts (%d F, %d T) for model %d, run %s, modelname %s %s.\n',ntotal,nF,nT,modelnum,self.nrun2phase{nrun},namestring,foldersuffix)
            if nF > 0
                for nnF = 1:nF
                    disp(['(F) ' matlabbatch{1}.spm.stats.con.consess{nnF}.fcon.name])
                end
            end
            for nnT = 1:nT
                disp(['(T) ' matlabbatch{1}.spm.stats.con.consess{nF+nnT}.tcon.name])
            end
        end
        
        
        function Con2ndlevel_FIR_rollwin(self,nrun,modelnum,namestring,binwin,varargin)
            deletecons = 0;
            t_con = 0;
            F_con = 1;
            
            versiontag = 0;
            foldersuffix = sprintf('_N%02d',self.total_subjects);
            bins2take = 1:14;
            
            
            if nargin == 6 %one varargin is given
                foldersuffix = varargin{1}; %here you can pass whatever extension you want, like test, or N39 or so
            elseif nargin == 7
                foldersuffix =  varargin{1}; %here you can pass whatever extension you want, like test, or N39 or so
                versiontag   = varargin{2};   %probably never needed, better name then with foldersuffix, easier to remember and document
            elseif nargin > 7
                fprintf('Too many inputs. Please debug.')
                keyboard;
            end
            
            fprintf('Starting 2nd Level Contrasts for FIR model %02d, run %02d, named %s, with foldersuffix ''%s''...\n',modelnum,nrun,namestring,foldersuffix);
            if modelnum > 3
                path2ndlevel = fullfile(self.path_second_level,'FIR',sprintf('model_%02d_FIR_%02d_%s_b%02dto%02d_%s%s',modelnum,versiontag,namestring,bins2take(1),bins2take(end),self.nrun2phase{nrun},foldersuffix));
            else %former FIR models were with 39 subs, so no suffix
                path2ndlevel = fullfile(self.path_second_level,sprintf('model_%02d_FIR_%02d_%s_b%02dto%02d_%s',modelnum,versiontag,namestring,bins2take(1),bins2take(end),self.nrun2phase{nrun}));
            end
            
            
            if ~exist(path2ndlevel)
                fprintf('Folder not found, please debug, or run 2ndlevel estimation first.\n')
            end
            
            nF = 0;
            nT = 0;
            n  = 0;
            
            path_spmmat = fullfile(path2ndlevel,'SPM.mat');
            
            matlabbatch{1}.spm.stats.con.spmmat = cellstr(path_spmmat);
            
            
            if strcmp(namestring,'CSdiff')
                ncond = 1;
            elseif strcmp(namestring,'CSPCSN')
                ncond =2;
            elseif strcmp(namestring,'8conds')
                ncond = 8;
                if nrun == 2
                    ncond = 2;
                end
            end
            
            
            for bin = bins2take(:)'
                ind = bin:(bin+binwin-1);
                if max(ind) <= (max(bins2take)) %if it falls out of number of bins
                    
                    vec = zeros(1,max(bins2take)); vec(ind) = 1./binwin;
                    vec = repmat(vec,1,ncond);
                    vec = vec./ncond;
                    
                    if t_con == 1
                        n = n + 1; nT = nT + 1;
                        matlabbatch{1}.spm.stats.con.consess{n}.tcon.name    =  sprintf('bin%02d_binwin%d',bin,binwin);
                        matlabbatch{1}.spm.stats.con.consess{n}.tcon.weights =  vec;
                        matlabbatch{1}.spm.stats.con.consess{n}.tcon.sessrep = 'none';
                    end
                    if F_con == 1
                        n = n + 1; nF = nF + 1;
                        matlabbatch{1}.spm.stats.con.consess{n}.fcon.name    =  sprintf('bin%02d_win%d',bin,binwin);
                        matlabbatch{1}.spm.stats.con.consess{n}.fcon.weights =  vec;
                        matlabbatch{1}.spm.stats.con.consess{n}.fcon.sessrep = 'none';
                    end
                end
            end
            
            
            if nT > 0
                for tc = 1:nT
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.name =      [matlabbatch{1}.spm.stats.con.consess{nF+tc}.tcon.name '_neg'];
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.weights =   -matlabbatch{1}.spm.stats.con.consess{nF+tc}.tcon.weights;
                    matlabbatch{1}.spm.stats.con.consess{n+tc}.tcon.sessrep =   'none';
                    nT = nT + 1;
                end
            end
            
            matlabbatch{1}.spm.stats.con.delete = deletecons;
            
            spm_jobman('run',matlabbatch);
            
            ntotal = nT + nF;
            fprintf('Done creating %d 2ndlevel contrasts (%d F, %d T) for model %d, run %s, modelname %s %s.\n',ntotal,nF,nT,modelnum,self.nrun2phase{nrun},namestring,foldersuffix)
            if nF > 0
                for nnF = 1:nF
                    disp(['(F) ' matlabbatch{1}.spm.stats.con.consess{nnF}.fcon.name])
                end
            end
            for nnT = 1:nT
                disp(['(T) ' matlabbatch{1}.spm.stats.con.consess{nF+nnT}.tcon.name])
            end
        end
        function [covyesno, cov] = get_2ndlevel_covariate(self,run,model_num)
            covyesno = 0;
            cov = struct('c', {}, 'cname',{}, 'iCFI', {}, 'iCC', {});
            switch model_num
                case 12
                    covyesno = 1;
                    Mpain = nan(1,self.total_subjects);
                    for sc = 1:self.total_subjects;
                        Mpain(sc) = nanmean(self.subject{sc}.get_pain(run));
                    end
                    cov = struct('c', {Mpain}, 'cname','Mpain', 'iCFI', {1}, 'iCC', {1});
            end
        end
        %% ROI methods
        function [out, base, test,tb,tt,Y] = get_betas_ROI(self,coords,namestr,varargin)
            path2ndlevel = fullfile(self.path_project,'spm/FIR/model_04_FIR_00_bin4win4_8conds_BT_b01to14_B_N35/');
            cd(path2ndlevel);
            
            vis = 1;
            dofit = 1;
            fitmethod = 3;
            plot_timeband = 1;
            
            if plot_timeband == 1
                sp = 3;
            else
                sp = 2;
            end
            
            if isempty(coords) %chose atlas manually
                if ~exist(fullfile(path2ndlevel,['VOI_',namestr,'.mat']))
                    choice = spm_select(1,'nii','Select Atlas or VOI',[]);
                    matlabbatch{1}.spm.util.voi.spmmat                   = {[path2ndlevel 'SPM.mat']};
                    matlabbatch{1}.spm.util.voi.adjust                   = 0;
                    matlabbatch{1}.spm.util.voi.session                  = 1;
                    matlabbatch{1}.spm.util.voi.name                     = namestr;
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.image        = {choice};
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold    = 0.5;
                    matlabbatch{1}.spm.util.voi.expression               = 'i1';
                    spm_jobman('run',matlabbatch);
                    if nargout > 0
                        load(['VOI_',namestr,'.mat'])
                    end
                else
                    cd(path2ndlevel)
                    load(['VOI_',namestr,'.mat'])
                end
                
            else
                r_sphere = 4;
                if nargin >3
                    r_sphere = varargin{1};
                end
                if ~exist(fullfile(path2ndlevel,['VOI_',namestr,'.mat']))
                    
                    matlabbatch{1}.spm.util.voi.spmmat                   = {[path2ndlevel 'SPM.mat']};
                    matlabbatch{1}.spm.util.voi.adjust                   = 0;
                    matlabbatch{1}.spm.util.voi.session                  = 1;
                    matlabbatch{1}.spm.util.voi.name                     = namestr;
                    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre     = coords;
                    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius     = r_sphere;
                    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
                    matlabbatch{1}.spm.util.voi.expression               = 'i1';
                    spm_jobman('run',matlabbatch);
                    
                    if nargout > 0
                        load(['VOI_',namestr,'.mat'])
                    end
                else
                    cd(path2ndlevel)
                    load(['VOI_',namestr,'.mat'])
                end
            end
            
            out =  reshape(Y,self.total_subjects,16);
            base.x = repmat(-135:45:180,self.total_subjects,1);
            base.y = out(:,1:8);
            base.ids = self.ids;
            test.x =  repmat(-135:45:180,self.total_subjects,1);
            test.y = out(:,9:16);
            test.ids = self.ids;
            
            if vis
                X = -135:45:180;
                bar_Y = nanmean(out);
                SEM = nanstd(out)./sqrt(self.total_subjects);
                
                graycol = [.3 .3 .3];
                figure(1000);
                clf
                subplot(1,sp,1);
                bar(X,bar_Y(:,1:8),'facecolor',graycol,'edgecolor','none','facealpha',.8);
                hold on
                errorbar(X,bar_Y(1:8),SEM(1:8),'.','Color',graycol,'LineWidth',2)
                set(gca,'xtick',[0 180],'xticklabel',{'CS+','CS-'},'TickLength',[0.03 0.025],'FontSize',14,'LineWidth',2);
                box off;
                set(gca,'color','none');
                drawnow;
                axis tight;box off;axis square;drawnow;alpha(.5);
                xlim([min(X)-mean(diff(X)) max(X)+mean(diff(X))])
                
                ylabel('Contrast Estimate [a.u.]')
                titl= title('Baseline');set(titl,'FontWeight','normal');
                if dofit==1
                    tb = Tuning(base);
                    tb.GroupFit(fitmethod);
                    if (10^-tb.groupfit.pval)<.05
                        plot(tb.groupfit.x_HD,tb.groupfit.fit_HD,'Color',graycol,'LineWidth',3)
                    else
                        plot(tb.groupfit.x_HD,ones(1,numel(tb.groupfit.fit_HD))*mean(tb.y_mean),'Color',graycol,'LineWidth',3)
                    end
                end
                figure(1000);
                subplot(1,sp,2);
                cmap  = Project.GetFearGenColors;
                cc= 0;
                for i = 9:16
                    cc = cc+1;
                    try
                        h(cc)    = bar(X(cc),bar_Y(i),40,'facecolor',cmap(cc,:),'edgecolor','none','facealpha',.8);
                    catch
                        h(cc)    = bar(X(cc),bar_Y(i),40,'facecolor',cmap(cc,:),'edgecolor','none');
                    end
                    hold on
                    errorbar(X(cc),bar_Y(i),SEM(i),'.','Color',cmap(cc,:),'LineWidth',2);
                end
                hold on
                if dofit==1
                    tt = Tuning(test);
                    tt.GroupFit(fitmethod);
                    if (10^-tt.groupfit.pval)<.05
                        plot(tt.groupfit.x_HD,tt.groupfit.fit_HD,'Color',graycol,'LineWidth',3)
                    else
                        plot(tt.groupfit.x_HD,ones(1,numel(tt.groupfit.fit_HD))*mean(tt.y_mean),'Color',graycol,'LineWidth',3)
                    end
                end
                
                %
                box off;
                set(gca,'color','none','xtick',[0 180],'xticklabel',{'CS+','CS-'},'TickLength',[0.03 0.025],'FontSize',14,'LineWidth',2);
                drawnow;
                axis tight;box off;axis square;drawnow;alpha(.5);
                xlim([min(X)-mean(diff(X)) max(X)+mean(diff(X))])
                title('Test','FontWeight','normal');
                
                EqualizeSubPlotYlim(gcf);
                
                if plot_timeband == 1
                    out = self.get_betas_ROI_fullFIR(coords,namestr,r_sphere);
                end
                
                if ~isempty(coords)
                    stit = supertitle(sprintf('[%3.1f %3.1f %3.1f] %s , p_{test} = %05.3f',coords(1),coords(2),coords(3),namestr,10.^-tt.groupfit.pval));
                else
                    stit = supertitle(sprintf('%s , p_{test} = %05.3f',namestr,10.^-tt.groupfit.pval));
                end
                cd(path2ndlevel)
                set(stit,'FontSize',16)
                set(gcf,'color','white')
                print([namestr '_r600.png'],'-dpng','-r600')
                
                if plot_timeband == 1
                    print([namestr '_svg.svg'],'-dsvg')
                end
            end
            
        end
        function [out] = get_betas_ROI_fullFIR(self,coords,namestr,varargin)
            
            %              path_SPM_fullFIR{1} = [Project.path_second_level,'FIR/model_04_FIR_01_8conds_b01to14_B_N35/'];
            path_SPM_fullFIR = [Project.path_second_level,'FIR/model_04_FIR_00_8conds_b01to14_T_N35/'];
            
            vis = 1;
            
            
            path2ndlevel = path_SPM_fullFIR;
            cd(path2ndlevel)
            if isempty(coords) %chose atlas manually
                if ~exist(fullfile(path2ndlevel,['VOI_',namestr,'.mat']))
                    choice = spm_select(1,'nii','Select Atlas or VOI',[]);
                    matlabbatch{1}.spm.util.voi.spmmat                   = {[path2ndlevel 'SPM.mat']};
                    matlabbatch{1}.spm.util.voi.adjust                   = 0;
                    matlabbatch{1}.spm.util.voi.session                  = 1;
                    matlabbatch{1}.spm.util.voi.name                     = namestr;
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.image        = {choice};
                    matlabbatch{1}.spm.util.voi.roi{1}.mask.threshold    = 0.5;
                    matlabbatch{1}.spm.util.voi.expression               = 'i1';
                    spm_jobman('run',matlabbatch);
                    if nargout > 0
                        load(['VOI_',namestr,'.mat'])
                    end
                else
                    cd(path2ndlevel)
                    load(['VOI_',namestr,'.mat'])
                end
                
            else
                r_sphere = 4;
                if nargin>3
                    r_sphere = varargin{1};
                end
                if ~exist(fullfile(path2ndlevel,['VOI_',namestr,'.mat']))
                   
                    
                    matlabbatch{1}.spm.util.voi.spmmat                   = {[path2ndlevel 'SPM.mat']};
                    matlabbatch{1}.spm.util.voi.adjust                   = 0;
                    matlabbatch{1}.spm.util.voi.session                  = 1;
                    matlabbatch{1}.spm.util.voi.name                     = namestr;
                    matlabbatch{1}.spm.util.voi.roi{1}.sphere.centre     = coords;
                    matlabbatch{1}.spm.util.voi.roi{1}.sphere.radius     = r_sphere;
                    matlabbatch{1}.spm.util.voi.roi{1}.sphere.move.fixed = 1;
                    matlabbatch{1}.spm.util.voi.expression               = 'i1';
                    spm_jobman('run',matlabbatch);
                    
                   
                    load(['VOI_',namestr,'.mat'])
                    
                else
                    cd(path2ndlevel)
                    load(['VOI_',namestr,'.mat'])
                end
            end
            
            out =  reshape(Y,self.total_subjects,14,8);
            
            
            if vis
                nrun = 3;
                X = 1:size(out,2);
                Y = squeeze(nanmean(out));
                SEM = squeeze(nanstd(out)./sqrt(self.total_subjects));
                CSP = [4 1 4];
                CSN = [8 2 8];

                graycol = [.3 .3 .3];
                
                cmap = self.GetFearGenColors;
                lwe = 3; %line width errorsGetFear
                lwa = 1.5; %line width axes
                
                figure(1000);
%                 clf
                
                subplot(1,3,3);
                ylims = 1.2*[min(min(Y(:,[4 8])))-max(max((SEM(:,[4 8])))) max(max(Y(:,[4 8])))+max(max((SEM(:,[4 8]))))];
                yticki = unique(sort([ceil(ylims(1)*100)./100 floor(ylims(2)*100)./100,0]));
                
                binbar = [[2:5 5.9;4:7 7.9;9:12 12.9]-3]...
                    *Project.TR; %8 is only to have the timewindow from bins.
                bins2sec = [[1:14]-3]*self.TR;
                %% timecourse CSP CSN
                
                ylabel('MRI signal [a.u.]');
                %         xlabel('secs since Ramp Onset');
                  title(Project.plottitles_BCT{nrun},'FontWeight','normal');
                ylim(ylims);
                xlim([min(bins2sec)-1 max(bins2sec)+1]);
                set(gca,'XTick',[0 11],'XTickLabel',{'TreatOn','TreatOff'},'YTick',yticki,'TickLength',[0.03 0.025],'XTickLabelRotation',0,'LineWidth',lwa);
                hold on;
                colrec = repmat([.6 .1 .6]',1,3);
                for tt =2
                    %             rectangle('Position',[binbar(tt,1),min(ylim)+range(ylim)*.05,binbar(tt,end)-binbar(tt,1),range(ylim)*.8],'EdgeColor',colrec(tt,:),'LineStyle',':','LineWidth',2);
                    pp= patch('vertices',...
                        [binbar(tt,1), min(Y(:))-max(SEM(:))-max(SEM(:)).*.5;...
                        binbar(tt,1) ,yticki(end);...
                        binbar(tt,end),yticki(end);...
                        binbar(tt,end),min(Y(:))-max(SEM(:))-max(SEM(:)).*.5],...
                        'faces',[1,2,3,4],'EdgeColor','none', 'FaceColor',graycol*.5,'FaceAlpha',.15);
                end
                %add shaded Errorbar for CSP and CSN
                if nrun > 1
                    shadedErrorBar(bins2sec,Y(:,CSP(nrun)),SEM(:,CSP(nrun)),'lineprops',{'Color',cmap(4,:),'LineWidth',lwe});hold on;
                    shadedErrorBar(bins2sec,Y(:,CSN(nrun)),SEM(:,CSN(nrun)),'lineprops',{'Color',cmap(8,:),'LineWidth',lwe});
                end
                box off;axis square;
                ylabel('MRI signal [a.u.]');
                set(gca,'FontSize',14,'LineWidth',2);
                if ~isempty(coords)
                    stit = supertitle(sprintf('[%3.1f %3.1f %3.1f] %s ',coords(1),coords(2),coords(3),namestr));
                else
                    stit = supertitle(sprintf('%s',namestr));
                end
                set(stit,'FontSize',16)
                set(gcf,'color','white')
                print([namestr '.png'],'-dpng')
            end
        end
        function overlay_roi(self,varargin)
            clear matlabbatch
            path2ROI = '/projects/crunchie/treatgen/data/spm/FIR/model_04_FIR_00_bin4win4_8conds_BT_b01to14_B_N35/';
            
            skullstrip_mean = '/projects/crunchie/treatgen/data/spm/groupmeans/ave_wCAT_ss_data.nii';
            if nargin > 1
                namestr = varargin{1};
                if strcmp(namestr,'bothHippo')
                    matlabbatch{1}.spm.util.imcalc.input = {
                        fullfile(path2ROI,['VOI_' 'rHippoKahntWimmerLissekOnat' '_mask.nii']); ...
                        fullfile(path2ROI,['VOI_' 'lHippoKahntWimmerLissek' '_mask.nii']); ...
                        skullstrip_mean
                        };
                    matlabbatch{1}.spm.util.imcalc.expression = 'i3.*(i1+i2+.5)';
                    matlabbatch{1}.spm.util.imcalc.output = ['overlay_mean_ss_ROI_' namestr ];
                else
                    matlabbatch{1}.spm.util.imcalc.input = {
                        fullfile(path2ROI,['VOI_' namestr '_mask.nii']); ...
                        skullstrip_mean
                        };
                    matlabbatch{1}.spm.util.imcalc.expression = 'i2.*(i1+.5)';
                    matlabbatch{1}.spm.util.imcalc.output = ['overlay_mean_ss_ROI_' namestr ];
                    
                end
            else
                
                matlabbatch{1}.spm.util.imcalc.input = {
                    fullfile(path2ROI,['VOI_' 'ACCGeuterBingel10_x0' '_mask.nii']); ...
                    fullfile(path2ROI,['VOI_' 'rHippoKahntWimmerLissekOnat' '_mask.nii']); ...
                    fullfile(path2ROI,['VOI_' 'lHippoKahntWimmerLissek' '_mask.nii']); ...
                    fullfile(path2ROI,['VOI_' 'lAmyDunsmoorBingel' '_mask.nii']); ...
                    fullfile(path2ROI,['VOI_' 'rAmyBingelOnat' '_mask.nii']); ...
                    skullstrip_mean
                    };
                matlabbatch{1}.spm.util.imcalc.expression = 'i6.*(i1+i2+i3+i4+i5+.5)';
                
                matlabbatch{1}.spm.util.imcalc.output = ['overlay_mean_ss_5ROIs'];
                
            end
            matlabbatch{1}.spm.util.imcalc.outdir = {path2ROI};
            
            matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
            matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
            matlabbatch{1}.spm.util.imcalc.options.mask = 0;
            matlabbatch{1}.spm.util.imcalc.options.interp = 1;
            matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
            spm_jobman('run',matlabbatch);
        end
        function [pvals, niftimat,idx,params] = ROI_singlevoxelfit(self,namestr)
            
            path2ndlevel = fullfile(self.path_project,'spm/FIR/model_04_FIR_00_bin4win4_8conds_BT_b01to14_B_N35/');
            filename =    fullfile(path2ndlevel,['VOI_',namestr,'.mat']);        

            fitmethod = 3;
            vis       = 1;
            if exist(filename)
                load(filename)
            else
                fprintf('No VOI.mat found, check your namestring or run it first! \n');
                keyboard
            end
            
            n_vox         = size(xY.y,2);
            coords_of_vox = xY.XYZmm;
            data_per_vox  = reshape(xY.y,self.total_subjects,16,n_vox);
            
            pvals    = nan(n_vox,1);
            LL       = nan(n_vox,1);
            LLnull   = nan(n_vox,1);
            params   = nan(n_vox,2);
            niftimat = nan(xY.spec.dim);
            niftimat_neg = nan(xY.spec.dim);
            niftimat_Ampl = nan(xY.spec.dim);            
            niftimat_Ampl_neg = nan(xY.spec.dim);
            niftimat_Sigma = nan(xY.spec.dim);
            
            fprintf('\nStarting single voxel fit...0...');
            start = tic;
            for n = 1:n_vox
               if mod(n,10)==0
                   fprintf('%d...',n)
               end
               data.x   = repmat(-135:45:180,self.total_subjects,1);
               data.y   = squeeze(data_per_vox(:,9:16,n));
               data.ids = 1:35;
               t = Tuning(data);
               t.GroupFit(fitmethod);
               pvals(n) = 10.^-t.groupfit.pval;
               params(n,:) = t.groupfit.Est(1:2);
               LL(n) = t.groupfit.Likelihood;
               LLnull(n) = t.groupfit.null_Likelihood;
               
            end
            done = toc(start);
            dur = done./60;
            fprintf('\ndone after %04.2f mins.\n',dur)
            
            % load a con image to find out which coords in mm correspond to
            % which coords in vox
            
            [~,xyz] = spm_read_vols(spm_vol(fullfile(path2ndlevel,'con_0010.nii'))); %random con number
                  
            for n = 1:n_vox
                idx(n) = find(xyz(1,:)==coords_of_vox(1,n) & xyz(2,:)==coords_of_vox(2,n) & xyz(3,:)==coords_of_vox(3,n));
            end
            niftimat(idx) = pvals;
            niftimat_neg(idx) = 1-pvals;
            niftimat_Ampl(idx) = params(:,1);
            
            niftimat_Ampl_neg(idx) = -params(:,1);
            niftimat_Sigma(idx) = params(:,2);
            
            
            
            dummynii = spm_vol(fullfile(path2ndlevel,'con_0010.nii')); %random con number
            dummynii = rmfield(dummynii,'pinfo');
            dummynii.fname = strrep(dummynii.fname,'con_0010.nii',sprintf('nifti_pvals_method%d_ROI_%s.nii',fitmethod, namestr));
            V =spm_write_vol(dummynii,niftimat);
            clear dummynii
            dummynii = spm_vol(fullfile(path2ndlevel,'con_0010.nii')); %random con number
            dummynii = rmfield(dummynii,'pinfo');
            dummynii.fname = strrep(dummynii.fname,'con_0010.nii',sprintf('nifti_neg_pvals_method%d_ROI_%s.nii',fitmethod, namestr));
            Vneg =spm_write_vol(dummynii,niftimat_neg);
            dummynii = spm_vol(fullfile(path2ndlevel,'con_0010.nii')); %random con number
            dummynii = rmfield(dummynii,'pinfo');
            dummynii.fname = strrep(dummynii.fname,'con_0010.nii',sprintf('nifti_Ampl_method%d_ROI_%s.nii',fitmethod, namestr));
            VAmpl =spm_write_vol(dummynii,niftimat_Ampl);
            
            dummynii = spm_vol(fullfile(path2ndlevel,'con_0010.nii')); %random con number
            dummynii = rmfield(dummynii,'pinfo');
            dummynii.fname = strrep(dummynii.fname,'con_0010.nii',sprintf('nifti_Ampl_neg_method%d_ROI_%s.nii',fitmethod, namestr));
            VAmplneg =spm_write_vol(dummynii,niftimat_Ampl_neg);
            
            dummynii = spm_vol(fullfile(path2ndlevel,'con_0010.nii')); %random con number
            dummynii = rmfield(dummynii,'pinfo');
            dummynii.fname = strrep(dummynii.fname,'con_0010.nii',sprintf('nifti_Sigma_method%d_ROI_%s.nii',fitmethod, namestr));
            VSigma =spm_write_vol(dummynii,niftimat_Sigma);
            
            save(fullfile(self.path_project,'midlevel',sprintf('pvals_method%d_ROI_%s.mat',fitmethod,namestr)),'pvals','niftimat','niftimat_neg','niftimat_Ampl','niftimat_Sigma','idx','params')
          
            
        end
        function plot_ROI_peakvox(self,namestr,varargin)
            
%             if nargin > 2
%                 r_sphere = varargin{1};
%             else
%                 r_sphere = 4;
%             end
            
           
            
            plot_timeband = 1;
            sp = 2;
            dofit     = 1;
            fitmethod = 3;
            criterion = 'min_p';
            
            path2ndlevel = fullfile(self.path_project,'spm/FIR/model_04_FIR_00_bin4win4_8conds_BT_b01to14_B_N35/');
            filename =    fullfile(path2ndlevel,['VOI_',namestr,'.mat']);
            
            fitmethod = 3;
            vis       = 1;
            
            if exist(filename)
                load(filename)
            else
                fprintf('No VOI.mat found, check your namestring or run it first! \n');
                keyboard
            end
            
            if exist(fullfile(self.path_project,'midlevel',sprintf('pvals_method%d_ROI_%s.mat',fitmethod,namestr)))
                load(fullfile(self.path_project,'midlevel',sprintf('pvals_method%d_ROI_%s.mat',fitmethod,namestr)));
            else
                [pvals, niftimat,idx,params] = self.ROI_singlevoxelfit(namestr);
            end
            
            n_vox              = length(pvals);
            data_per_vox       = reshape(xY.y,self.total_subjects,16,n_vox);
            if strcmp(criterion,'min_p')
                [pval_min,ind_peak] = min(pvals);
            elseif strcmp(criterion,'max_Ampl')
                [max_ampl,ind_peak] = max(abs(params(:,1).*(pvals<.05)));
            end
            coords_min_mm      = xY.XYZmm(:,ind_peak);
                        
            data = squeeze(data_per_vox(:,:,ind_peak));
            
            
            
            if plot_timeband == 1
                path_SPM_fullFIR = [Project.path_second_level,'FIR/model_04_FIR_00_8conds_b01to14_T_N35/'];
                a =  load(fullfile(path_SPM_fullFIR,['VOI_',namestr,'_FIR.mat']),'xY','Y');
                xY_FIR = a.xY;
                Y_FIR  = a.Y;
                data_mean_ROI_FIR =  reshape(Y_FIR,self.total_subjects,14,8);
                data_per_vox_FIR  = reshape(xY_FIR.y,self.total_subjects,14,8,n_vox);
                
            end
                        
            data_FIR = squeeze(data_per_vox_FIR(:,:,:,ind_peak));
            
            %% plotting bin4win4 Tuning from peak voxel
            
            
            base.x = repmat(-135:45:180,self.total_subjects,1);
            base.y = data(:,1:8);
            base.ids = self.ids;
            test.x =  repmat(-135:45:180,self.total_subjects,1);
            test.y = data(:,9:16);
            test.ids = self.ids;
            
            if vis
                if plot_timeband == 1
                    sp = 3;
                end
                % plot tuning
                % baseline
                X = -135:45:180;
                bar_Y = nanmean(data);
                SEM = nanstd(data)./sqrt(self.total_subjects);
                
                graycol = [.3 .3 .3];
                figure(1000);
                clf
                subplot(1,sp,1);
                  hold on
                bar(X,bar_Y(:,1:8),'facecolor',graycol,'edgecolor','none','facealpha',.8);
                hold on
                errorbar(X,bar_Y(1:8),SEM(1:8),'.','Color',graycol,'LineWidth',2)
                set(gca,'xtick',[0 180],'xticklabel',{'CS+','CS-'},'TickLength',[0.03 0.025],'FontSize',14,'LineWidth',2);
                box off;
                set(gca,'color','none');
                drawnow;
                axis tight;box off;axis square;drawnow;alpha(.5);
                xlim([min(X)-mean(diff(X)) max(X)+mean(diff(X))])
%                 yyy(1,:) = ylim;
                
                ylabel('Contrast Estimate [a.u.]')
                titl= title('Baseline');set(titl,'FontWeight','normal');
                if dofit==1
                    tb = Tuning(base);
                    tb.GroupFit(fitmethod);
                    if (10^-tb.groupfit.pval)<.05
                        plot(tb.groupfit.x_HD,tb.groupfit.fit_HD,'Color',graycol,'LineWidth',3)
                    else
                        plot(tb.groupfit.x_HD,ones(1,numel(tb.groupfit.fit_HD))*mean(tb.y_mean),'Color',graycol,'LineWidth',3)
                    end
                end
                % test 
                figure(1000);
                subplot(1,sp,2);
                hold on
                cmap  = Project.GetFearGenColors;
                cc= 0;
                for i = 9:16
                    cc = cc+1;
                    try
                        h(cc)    = bar(X(cc),bar_Y(i),40,'facecolor',cmap(cc,:),'edgecolor','none','facealpha',.8);
                    catch
                        h(cc)    = bar(X(cc),bar_Y(i),40,'facecolor',cmap(cc,:),'edgecolor','none');
                    end
                    hold on
                    errorbar(X(cc),bar_Y(i),SEM(i),'.','Color',cmap(cc,:),'LineWidth',2);
                end
                hold on
                if dofit==1
                    tt = Tuning(test);
                    tt.GroupFit(fitmethod);
                    if (10^-tt.groupfit.pval)<.05
                        plot(tt.groupfit.x_HD,tt.groupfit.fit_HD,'Color',graycol,'LineWidth',3)
                    else
                        plot(tt.groupfit.x_HD,ones(1,numel(tt.groupfit.fit_HD))*mean(tt.y_mean),'Color',graycol,'LineWidth',3)
                    end
                end
                
                %
                box off;
                set(gca,'color','none','xtick',[0 180],'xticklabel',{'CS+','CS-'},'TickLength',[0.03 0.025],'FontSize',14,'LineWidth',2);
                drawnow;
                axis tight;box off;axis square;drawnow;alpha(.5);
                xlim([min(X)-mean(diff(X)) max(X)+mean(diff(X))])
%                 yyy(2,:) = ylim;
                title('Test','FontWeight','normal');
                
                EqualizeSubPlotYlim(gcf);
                
                 %% plot timeband
                if plot_timeband == 1
                    nrun = 3;
                    X = 1:size(data_per_vox_FIR,2);
                    Y = squeeze(nanmean(data_FIR));
                    SEM = squeeze(nanstd(data_FIR)./sqrt(self.total_subjects));
                    CSP = [4 1 4];
                    CSN = [8 2 8];
                    
                    graycol = [.3 .3 .3];
                    
                    cmap = self.GetFearGenColors;
                    lwe = 3; %line width errorsGetFear
                    lwa = 1.5; %line width axes
                    
                    figure(1000);
                    
                    subplot(1,sp,3);
                    ylims = 1.2*[min(min(Y(:,[4 8])))-max(max((SEM(:,[4 8])))) max(max(Y(:,[4 8])))+max(max((SEM(:,[4 8]))))];
                    yticki = unique(sort([ceil(ylims(1)*100)./100 floor(ylims(2)*100)./100,0]));
                    
                    binbar = [[2:5 5.9;4:7 7.9;9:12 12.9]-3]...
                        *Project.TR; %8 is only to have the timewindow from bins.
                    tw = 2;
                    bins2sec = [[1:14]-3]*self.TR;
                    %% timecourse CSP CSN
                    ylabel('MRI signal [a.u.]');
                    title(Project.plottitles_BCT{nrun},'FontWeight','normal');
                    ylim(ylims);
                    xlim([min(bins2sec)-1 max(bins2sec)+1]);
                    set(gca,'XTick',[0 11],'XTickLabel',{'TreatOn','TreatOff'},'YTick',yticki,'TickLength',[0.03 0.025],'XTickLabelRotation',0,'LineWidth',lwa);
                    hold on;
                    for tw =2
                        %             rectangle('Position',[binbar(tt,1),min(ylim)+range(ylim)*.05,binbar(tt,end)-binbar(tt,1),range(ylim)*.8],'EdgeColor',colrec(tt,:),'LineStyle',':','LineWidth',2);
                        pp= patch('vertices',...
                            [binbar(tw,1), min(Y(:))-max(SEM(:))-max(SEM(:)).*.5;...
                            binbar(tw,1) ,yticki(end);...
                            binbar(tw,end),yticki(end);...
                            binbar(tw,end),min(Y(:))-max(SEM(:))-max(SEM(:)).*.5],...
                            'faces',[1,2,3,4],'EdgeColor','none', 'FaceColor',graycol*.5,'FaceAlpha',.15);
                    end
                    %add shaded Errorbar for CSP and CSN
                    if nrun > 1
                        shadedErrorBar(bins2sec,Y(:,CSP(nrun)),SEM(:,CSP(nrun)),'lineprops',{'Color',cmap(4,:),'LineWidth',lwe});hold on;
                        shadedErrorBar(bins2sec,Y(:,CSN(nrun)),SEM(:,CSN(nrun)),'lineprops',{'Color',cmap(8,:),'LineWidth',lwe});
                    end
                    box off;axis square;
                    ylabel('MRI signal [a.u.]');
                    set(gca,'FontSize',14,'LineWidth',2);
                end
                %% format figure
                stit = supertitle(sprintf('Peak voxel (%s) at [%3.1f %3.1f %3.1f] mm %s , p_{test} = %05.3f',criterion,coords_min_mm(1),coords_min_mm(2),coords_min_mm(3),strrep(namestr,'_',' '),10.^-tt.groupfit.pval));
               set(stit,'FontSize',16);
                for spp = 1:2
                    subplot(1,sp,spp)
                    yyy = ylim;
                    set(gca,'YTick',[round(min(yyy)*100)./100 0 floor(max(yyy)*100)./100])
                end
                cd(path2ndlevel)
                set(stit,'FontSize',16)
                set(gcf,'color','white')
                print(sprintf('%s_%s_peakVOX_r600.png',namestr,criterion),'-dpng','-r600')
                
                if plot_timeband == 1
                    print(sprintf('%s_%s_peakVOX_svg.svg',namestr,criterion),'-dsvg')
                end
            end
        end
    end
end