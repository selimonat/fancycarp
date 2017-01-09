classdef Fixmat < Project
    properties (Hidden,Constant)
        window       = 250;%half size of the fixation maps
    end
    
    properties (Hidden,SetAccess = public)
        %related to map computation
        bc                  = 0;%cocktail blank correction
        unitize             = 1;%sums to 1 or not
        maptype             = 'conv';%conv or bin
        kernel_fwhm         = Fixmat.PixelPerDegree*.8;
        binsize             = 25;
        linkage_method      = 'average';
        linkage_metric      = 'correlation';
        similarity_metric   = 'correlation';
        dendro_tree         = [];%results of linkage analysis will be stored here.
        dendro_D            = [];
        dendro_leafOrder    = [];
        maps%current maps;
        mds%mds results
        cmap_limits         = 7;
    end
    
    properties (Hidden,SetAccess = private)
        %internal
        query
        all_queries
        map_titles
        selection
        realcond%all conditions that are not ucs,odd,or nulltrial
    end
    
    properties (SetAccess = private,Dependent,Hidden)
        %rect: [y x y_size x_size]
        rect       %the aperture
        binedges   %bin edges for binned fixation maps
        stimulus   %contains average stimulus
    end
    
    properties (SetAccess = private)
        %experimental observations
        subject = [];
        phase   = [];
        start   = [];
        stop    = [];
        x       = [];
        y       = [];
        eye     = [];
        deltacsp= [];
        file    = [];
        oddball = [];
        trialid = [];
        ucs     = [];
        fix     = [];
        chain   = [];
        isref   = [];
    end
    
    events
        UpdateGraph
    end
    
    methods
        function obj = Fixmat(subjects,runs)%constructor
            %%
            %initialize
            for run = runs(:)'
                for subject = subjects(:)'
                    %                     if exist(regexprep(obj.path2data(subject,run,'eye'),'.mat','.edf')) ~= 0 %does the subject has an edf file?
                    
                    %% edf2mat conversion
                    if ~exist(obj.path2data(subject,run,'eye'))%is the edf file converted?
                        fprintf('-edfread not yet ran (subject:%03d, run:%03d)\n-Roger that, will do it ASAP, sir!\n',subject,run);
                        [trials info] = edfread(regexprep(obj.path2data(subject,run,'eye'),'.mat','.edf'),'TRIALID');
                        save(obj.path2data(subject,run,'eye'),'trials','info');
                        %if the conversion fails the code will fail, it is recommended to rename the data.edf file to something else then.
                    else
                        fprintf('edf is converted (subject:%03d, run:%03d)\n',subject,run);
                    end
                    
                    dummy = load(obj.path2data(subject,run,'eye'));
                    %this is necessary to expand the PTB message to
                    %something that is understandable by the fixations
                    %method.
                    for nt = 1:length(dummy.trials)
                        %extract relevant info first from the text message
                        T = dummy.trials(nt).TRIALID.msg(1,:);
                        %add them as msg so that .getfixmat can create the fields
                        pairs = regexp(T,'[\w-]*','match');
                        for np = 1:2:length(pairs)-1
                            dummy.trials(nt).(lower(pairs{np})).msg = pairs{np+1};
                        end
                    end
                    %
                    dummy.trials = rmfield(dummy.trials,{'TRIALID'});%make .getfixmat happy
                    dummy        = obj.getfixmat(dummy.trials,dummy.info);%get the fixmat
                    dummy.subject= repmat(uint32(subject),1,length(dummy.x));
                    %and append it to the previous fixmat
                    for fns = properties(obj)'
                        %% take fixations which are only coming from the required phase.
                        %(e.g. ratings in baseline are coded as 5)
                        valid_fix = dummy.phase == run;
                        if isfield(dummy,fns{1})%if it is not a property dont even consider
                            obj.(fns{1}) = [obj.(fns{1}) dummy.(fns{1})(valid_fix)];%append it to the previous
                        else
                            obj.(fns{1}) = [obj.(fns{1}) zeros(1,sum(valid_fix))];%append it to the previous
                        end
                    end
                end
                %                 end
            end
            %% remove fixations outside of the image border
            obj.selection = ~(obj.x < obj.rect(2) | obj.x > (obj.rect(2)+obj.rect(4)-1) | obj.y < obj.rect(1) | obj.y > (obj.rect(1)+obj.rect(3)-1) );
            obj.ApplySelection;
            %% round coordinates to pixels
            obj.x = round(obj.x);
            obj.y = round(obj.y);
            obj.realcond = unique(obj.deltacsp(~ismember(obj.deltacsp,[500 1000 3000]))); %500 UCS, 1000 Odd, 3000 Null
        end
        function UpdateSelection(obj,varargin)
            %takes VARARGIN pairs to update the selection vektor
            obj.selection = logical(ones(1,length(obj.subject)));
            for n = 1:2:length(varargin)
                obj.selection = obj.selection .* ismember(obj.(varargin{n}),varargin{n+1});
            end
            %             fprintf('Selection vector updated...\n');
            obj.selection = logical(obj.selection);
            obj.query     = varargin;
        end
        function ApplySelection(obj)
            %removes fixations based on the fixation vecktor, only used in
            %the constructor.
            %removes fixations which are FALSE in selection
            for p = properties(obj)'
                if ~strcmp(p{1},'rect') && ~strcmp(p{1},'maps') && ~strcmp(p{1},'selection');%all the properties unrelated to fixation data.
                    if ~isempty(obj.(p{1}))%sometimes the property is not filled in
                        obj.(p{1})(~obj.selection) = [];
                    end
                end
            end
            obj.query = [];
            fprintf('Selection (%04d fixations) removed from the object...\n',sum(~obj.selection));
        end
        function [cormat] = condcorrmat(obj,ph,plot)
            tsub = length(unique(obj.subject));
            cormat = nan(8,8,tsub);
            subc            = 0;
            for sub = unique(obj.subject);
                subc = subc + 1;
                %create the query cell
                v = [];
                c = 0;
                for cond = -135:45:180
                    c    = c+1;
                    v{c} = {'phase', ph, 'deltacsp' cond 'subject' sub};
                end
                % plot and save fixation maps
                obj.getmaps(v{:});
                maps = obj.maps;
                obj.maps   = maps(:,:,1:8) - repmat(mean(maps(:,:,1:8),3),[1 1 8]);%take out the average
                cormat(:,:,subc) = obj.corr;
            end
            if plot ==1
                figure;
                imagesc(median(cormat,3),[-.1 .1])
                axis square;colorbar
                set(gca,'fontsize',15)
                axis off
            end
        end
        function out = entropy(obj)
            
            obj.maps = obj.maps + eps;
            out = nan(size(obj.maps,3),1);
            maps = obj.vectorize_maps;
            out = -diag(maps'*log2(maps));
            %normalize by maximal entropy possible
            map0 = repmat(1/length(maps),[length(maps) 1]);
            entr0  = -diag(map0'*log2(map0));
            out = out./entr0;
            
            
        end
        
        function mds_run(obj,subjects,phase)
            %will run an MDS analysis on subjects
            mds_dim = 2;
            if mds_dim == 2
                initial_positions = [cosd(double(obj.realcond));sind(double(obj.realcond))]';
                initial_positions = reshape(zscore(initial_positions(:)),size(initial_positions));%shifts and scale all coordinates to have same center and scale
            else
                initial_positions = [];
            end
            mds = [];
            sc  = 0;
            for ns   = subjects(:)'
                sc          = sc +1;
                %
                c = 0;
                for cond = unique(obj.realcond)
                    c    = c+1;
                    v{c} = {'subject' ns 'deltacsp' cond 'phase' phase};
                end
                obj.getmaps(v{:});
                vecmaps     = obj.vectorize_maps;
                %
                vecmaps     = vecmaps - repmat(mean(vecmaps,2),1,size(vecmaps,2));
                %                 S           = squareform(pdist(vecmaps',obj.similarity_metric));
                S           = squareform(pdist(vecmaps','correlation'));
                coor        = mdscale(S,mds_dim,'criterion','metricstress');%,'start',initial_positions);
                %                 coor        = cmdscale(S,mds_dim);
                %                 coor        = zscore(coor(:));%shifts and scale all coordinates to have same center and scale
                mds(:,sc)   = coor(:); %set to (im,2) for 2-dimensional scaling.
            end
            %now we procrustes-align
            E_mean      = mean(mds,2);
            E_mean      = reshape(mean(mds,2),size(E_mean,1)./mds_dim,mds_dim);
            mds_aligned = [];
            for ns = 1:size(mds,2)
                coor              = mds(:,ns);
                coor              = reshape(coor,size(coor,1)./mds_dim,mds_dim)
                [d z transform]   = procrustes(E_mean,coor,'Reflection',false);
                mds_aligned(:,ns) = z(:);
            end
            obj.mds = mds_aligned;
        end
        function mds_plot(obj)
            
            x = obj.mds(1:8,:);
            y = obj.mds(9:end,:);
            plot(mean(x,2),mean(y,2),'r--');
            hold on;
            plot(mean(x,2),mean(y,2),'ro');
            
            %             for node = 1:size(x,1)
            %                 error_ellipse([x(node,:) ;y(node,:)]','color',Project.colors(node+1,:),'linewidth',1.5);
            %             end
            hold off
            
            
            
        end
        
        function getcondmaps(obj,varargin)
            %will return maps for different conditions
            if nargin > 1
                subjects = varargin{1};
            else
                subjects = unique(obj.subject);
                
            end
            %will get single subject maps
            v = [];
            c = 0;
            for cond = unique(obj.realcond)
                c    = c+1;
                v{c} = {'subject' subjects 'deltacsp' cond};
            end
            obj.getmaps(v{:});
        end
        
        function getsubmaps(obj,varargin)
            %will return all subjects' maps. VARARGIN cd contain specific
            %subjects, otherwise all are returned.
            if nargin > 1
                subjects = varargin{1};
            else
                subjects = unique(obj.subject);
            end
            %will get single subject maps
            v = [];
            c = 0;
            for sub = subjects
                c    = c+1;
                v{c} = {'subject' sub 'deltacsp' obj.realcond};
            end
            obj.getmaps(v{:});
        end
        function linkage(obj)
            fprintf('Conducting linkage analysis with linkage method: _%s_ and linkage metric: _%s_\n',char(obj.linkage_method),obj.linkage_metric);
            %provides data for a dendrogram analysis.
            vecmaps     = obj.vectorize_maps;
            obj.dendro_tree        = linkage(vecmaps',obj.linkage_method,obj.linkage_metric);
            obj.dendro_D           = pdist(vecmaps',obj.linkage_metric);
            obj.dendro_leafOrder   = optimalleaforder(obj.dendro_tree,obj.dendro_D);
        end
        
        function [branch_id,order0,pmat]=dendrogram(obj,k,varargin)
            %will plot fixmaps as a dendrogram. VARARGIN limits the number
            %of leafs, use 0 to have as many leafs as number of fixmaps
            
            if isempty(obj.dendro_tree)
                obj.linkage;
            end
            
            % plotting business
            %plot the dendro
            figure;
            [H,T,order]    = dendrogram(obj.dendro_tree,k);
            [H2,~,order0]   = dendrogram(obj.dendro_tree,0,'colorthreshold',obj.dendro_tree(end-k+2,3),'reorder',obj.dendro_leafOrder);
            title('optimal leaforder')
            set(H2,'LineWidth',2);
            axis off
            branch_id = T;
            % plot the single subject fixation maps
            figure;
            maps = obj.maps(:,:,order0);
            %
            imagesc(reshape(maps,[size(maps,1) size(maps,2)*size(maps,3)]));
            axis image;
            set(gca,'xtick',0:size(maps,1):max(xlim)-1,'ytick',[],'xticklabel',[],'xcolor','r');
            grid on
            subplotChangeSize(gca,.15,.25);
            if ~isempty(varargin)
                data = varargin;
                data = data{:};
                % plot cluster averages
                figure;
                tcluster = max(T(:));
                c = 0;
                for ncluster = unique(T(order0),'stable')'
                    c = c + 1;
                    subplot(3,tcluster,c)
                    imagesc(mean(obj.maps(:,:,T == ncluster),3))
                    axis image
                    axis off;
                    title(sprintf('(N  = %03d)\n',sum(T==ncluster)));
                end
                %
                %plot the data vector for each cluster also
                subplot(3,tcluster,k+1:k*2)
                clusterord = unique(T(order0),'stable')';
                c=0;
                for ncluster = clusterord
                    c=c+1;
                    N(ncluster) = sum(~isnan(data(T == ncluster)));
                    m(ncluster) = nanmean(data(T == ncluster));
                    s(ncluster) = nanstd(data(T == ncluster))./sqrt(N(ncluster));
                    text(c-0.1,m(ncluster)+s(ncluster)*1.1, sprintf('(N = %03d)\n',N(ncluster)));
                    hold on;
                end
                errorbar(1:tcluster,m(clusterord),s(clusterord),'LineWidth',2);
                %
                subplot(3,tcluster,(k+1:k*2)+k)
                for ncluster1 = unique(T(order0),'stable')'
                    for ncluster2 = unique(T(order0),'stable')'
                        [h p stats] = ttest2(data(T == ncluster1),data(T == ncluster2));
                        tmat(ncluster1,ncluster2) = -log10(p);
                        pmat(ncluster1,ncluster2) = p;
                    end
                end
                imagesc(tmat);
                set(gca,'XTick',1:tcluster,'XTickLabel',clusterord,'YTick',1:tcluster,'YTickLabel',clusterord)
                colorbar
                axis square
                fprintf('t-test p = %g. \n',10.^-max(tmat))
            end
        end
        
        function plotband(obj,varargin)%varargin can reorder the subjmaps
            if nargin > 1
                order = varargin{1};
            else
                order = unique(obj.subject);
            end
            N = length(order);
            figure;
            c=0;
            for sub = order(:)';
                c=c+1;
                obj.getmaps({'subject' sub 'deltacsp' obj.realcond})
                h = subplot(1,N,c);imagesc(obj.maps);
                %                 title(num2str(sub))
                axis square
                axis off
                subplotChangeSize(h,.01,.01);
                colorbar off
            end
        end
        function [out] = clustertest(obj,cluster,criterion)
            k = size(cluster,2);
            N = length([cluster(1:k).subs]);
            fprintf('Basing analysis on %d clusters.\n',k)
            clf
            group =nan(N,1);
            for n = 1:k
                group(cluster(n).subs) = n;
                m(n) = mean(criterion(cluster(n).subs));
                s(n) = std(criterion(cluster(n).subs));
                bar(n,m(n))
                hold on
                errorbar(n,m(n),s(n)./sqrt(length(cluster(n).subs)),'k.')
                xlim([0 k+1])
            end
            [out.anova.p,out.anova.tab,out.anova.stats] = anova1(criterion,group); %compares k groups
            imax = find(max(m));
            imin = find(min(m));
            [~,out.ttest.p,out.ttest.ci] = ttest2(criterion(cluster(imax).subs),criterion(cluster(imin).subs));%compares higgest to smallest bar
            
        end
        function plot(obj,varargin)
            
            M = obj.maps;
            if nargin > 1
                cmap = varargin{1};
            else
                cmap = 'linear';
            end
            %get colormap limits.
            if strcmp(cmap,'linear')
                [d u] = GetColorMapLimits(M,obj.cmap_limits);
                if sum(obj.maps(:) < 0) == 0%if there are no negative values
                    d = 0;
                end
            elseif strcmp(cmap,'log')
                u = max(log10(M(:)));
                d = u - 1;
                M = log10(M);
            end
            %
            %             ffigure;
            clf
            if nargin > 2
                nsp = varargin{2};
            else
                nsp     = obj.subplot_number;
            end
            for nc = 1:size(M,3)
                h   = subplot(nsp(1),nsp(2),nc);
                
                %plot the image;
                imagesc(obj.bincenters_x(500),obj.bincenters_y(500),obj.stimulus);
                hold on;
                h     = imagesc(obj.bincenters_x,obj.bincenters_y,M(:,:,nc),[d u]);
                set(h,'alphaData',Scale(abs((obj.maps(:,:,nc))))*.5+.5);
                axis image;
                axis off;
                            
%                 try
%                     t     = sprintf('%s %d ',obj.map_titles{nc}{1:2});
%                     title(t,'interpreter','none');
%                 end
                drawnow
            end
%                         thincolorbar('vert');
colorbar
        end
        function getmaps_split(obj,varargin)
            %will duplicate the number of varargin with odd and even trials
            c= 0;cc=0;
            v_new{1} = {};
            for v = varargin
                c = c+1;
                obj.UpdateSelection(v{1}{:});
                trials          = sort(unique(obj.trialid(obj.selection)));
                ttrial          = length(trials);
                if ttrial >=2;%there sd be at least two trials
                    cc              = cc +1;
                    mid             = round(ttrial/2);
                    v_new{1}{c}     = [v{1} 'trialid' trials(1:mid)];                
                    v_new2{1}{c}    = [v{1} 'trialid' trials(mid+1:end)];                
                end
            end 
            if cc == c
                v_new{1} = [v_new{:} v_new2{:}];            
                obj.getmaps(v_new{1}{:});       
            else
                obj.maps = [];
                cprintf([1 0 0],'not enough trials for oddevenings...\n');
            end
        end    
        function getmaps(obj,varargin)
            %will populate the maps property based on the filter in
            %VARARGIN, example:
            %fs.getmaps({'phase',4,'deltacsp',-135});
            
            obj.maps = [];%clear whatever is there
            c  = 0;
            obj.all_queries = varargin;
            for v = varargin
                c                 = c+1;
                obj.UpdateSelection(v{1}{:});
                if strcmp(obj.maptype,'conv')
                    %accum and conv
                    FixMap            = accumarray([obj.current_y-obj.rect(1)+1 obj.current_x-obj.rect(2)+1],1,[obj.rect(3) obj.rect(4)]);
                    FixMap            = conv2(sum(obj.kernel),sum(obj.kernel,2),FixMap,'same');
                elseif strcmp(obj.maptype,'bin')
                    %divide by a factor (i.e. binning) and accum, but no conv
                    [FixMap]          = hist3([obj.current_y obj.current_x],'edges',obj.binedges);
                    FixMap            = FixMap/obj.current_ttrial;%divide by the number of trials
                    %remove the last column and row
                    FixMap(end,:)     = [];
                    FixMap(:,end)     = [];
                end
                %
                if obj.unitize
                    FixMap            = FixMap./sum(FixMap(:));
                end
                obj.maps(:,:,c)       = FixMap;
                %
                obj.map_titles{c}     = obj.query(:);
            end
            %correct for baseline if wanted.
            if obj.bc
                obj.maps = obj.maps - repmat(nanmean(obj.maps,3),[1 1 size(obj.maps,3)]);
                nnan = sum(isnan(obj.maps(1,1,:)));%total number of maps with nans;
                if nnan ~= 0
                    cprintf('*[1 0 0]','ATTENTION: %i maps was just full of NaNs\n',nnan);
                end
            end
        end
        function bild = get.stimulus(obj)
            
            bild    = imread(obj.find_stim);
            bild    = bild( obj.rect(1):obj.rect(1)+obj.rect(3)-1,  obj.rect(2):obj.rect(2)+obj.rect(4)-1);
            bild    = repmat(bild,[1 1 3]);
        end
        function out = cropmaps(obj,map)
            if strcmp(obj.maptype,'conv')
                out      = map(obj.rect(2)/2-obj.mapsize : obj.rect(2)/2+obj.mapsize, obj.rect(4)/2-obj.mapsize: obj.rect(4)/2+obj.mapsize);
            elseif strcmp(obj.maptype,'bin')
                new_size = obj.mapsize./obj.binsize;
                
                out      = map(obj.rect(2)/obj.binsize/2-new_size : obj.rect(2)/obj.binsize/2+new_size, obj.rect(4)/obj.binsize/2-new_size: obj.rect(4)/obj.binsize/2+new_size);
            end
        end
        function out = kernel(obj)
            dummy       = make_gaussian2D(obj.kernel_fwhm*2.65,obj.kernel_fwhm*2.65,obj.kernel_fwhm,obj.kernel_fwhm,obj.kernel_fwhm*2.65/2,obj.kernel_fwhm*2.65/2);
            out         = dummy./sum(dummy(:));
        end
        function maps   = vectorize_maps(obj)
            if ~isempty(obj.maps)
                maps = reshape(obj.maps,size(obj.maps,1)*size(obj.maps,2),size(obj.maps,3));
            else
                fprintf('no maps stored here\n');
            end
        end
        function [out] = subplot_number(self);
            % how many subplots
            [row col] = GetSubplotNumber(size(self.maps,3));
            out       = sort([row col]);
        end
        function out = cov(obj)
            %             figure(3);clf
            %             set(gcf,'position',[1104         152         337         299]);
            out = cov(obj.vectorize_maps);
            %             imagesc(out);
            %             axis image;
            %             colorbar
        end
        function [out pval] = corr(obj)
            %plots the corrmat between maps
            %             figure(2);clf
            %             set(gcf,'position',[1104         152         337         299])
            [out pval] = corr(obj.vectorize_maps);
            %             imagesc(out);
            %             axis image;
            %             colorbar
        end
        function [x] = current_x(obj)
            %returns the current x y coordinates
            x = obj.x(obj.selection)';
        end
        function [y] = current_y(obj)
            %returns the current x y coordinates
            y = obj.y(obj.selection)';
        end
        function ttrial = current_ttrial(obj)
            %computes number of trials included in the current selection
            ttrial = length(unique([obj.trialid(obj.selection) ;obj.subject(obj.selection)]','rows'));
        end
        function out = get.rect(obj)
            out = [obj.screen_resolution(1)/2-obj.window obj.screen_resolution(2)/2-obj.window [obj.window obj.window]*2];
        end
        function out = get.binedges(obj)
            out = {obj.rect(1):obj.binsize:(obj.rect(1)+obj.rect(3)) ...
                obj.rect(2):obj.binsize:(obj.rect(2)+obj.rect(4))};
        end
        function [out] = bincenters_y(obj,varargin)
            %
            if nargin < 2
                resolution = length(obj.binedges{1});
            else
                resolution = varargin{1};
            end
            %transforms edges to centers with a required resolution.
            out = linspace( min(obj.binedges{1}), max(obj.binedges{1}),resolution);
            out = out + mean(unique(diff(out)))/2;
            out(end) =[];
        end
        function [out] = bincenters_x(obj,varargin)
            if nargin < 2
                resolution = length(obj.binedges{2});
            else
                resolution = varargin{1};
            end
            %transforms edges to centers with a required resolution.
            out = linspace( min(obj.binedges{2}), max(obj.binedges{2}),resolution);
            out = out + mean(unique(diff(out)))/2;
            out(end) =[];
        end
        function [out] = mapsize(obj,varargin)
            if ~isempty(obj.maps)
                %returns size of fixmap along VARARGIN dimension
                out = size(obj.maps,varargin{1});
            else
                out =[];
            end
        end
        function [count index] = histogram(obj)
            %will generate an histogram of fixation numbers
            subject    = unique(obj.subject);
            phase      = unique(obj.phase);
            condition  = unique(obj.deltacsp);
            count      = nan(length(subject),length(condition));
            sub = 0;
            for ns = subject
                sub = sub + 1;
                index.sub(sub) = ns;
                ph  = 0;
                for np = phase
                    ph = ph + 1;
                    cond = 0;
                    for nc = condition
                        cond = cond + 1;
                        index.cond(cond) = double(nc);%column identities
                        obj.UpdateSelection('subject',ns,'phase',np,'deltacsp',nc);
                        %number of trials
                        repet              = obj.current_ttrial;
                        %average number of fixation
                        count(sub,cond,ph) = sum(obj.selection)./repet;
                    end
                end
            end
            figure;
            set(gcf,'position',[440   393   920   405]);
            for phases = 1:ph
                subplot(1,ph,phases)
                imagesc(count(:,:,phases),[0 7]);
                colorbar
                ylabel('subjects')
                xlabel('condition');
                axis image;
                box off;
                title(sprintf('Phase: %03d',phases))
            end
            supertitle('Fixation Counts',1);
        end
        function fixmat  = getfixmat(obj,edfdata, edfmeta)
            % fixmat = fixations(edfdata, edfmeta)
            %   where [edfdata, edfmeta] = edfread('file.edf')
            %   returns a flat fixation matrix
            %   returns user defined meta data per trial and per experiment
            %   removes pre-stimulus-onset fixations
            
            % parameter with_pupil: return pupil data
            with_pupil = 0;
            % number of trials
            trials = length(edfdata);
            
            % get field names of custom per-trial meta data
            function my_fields = trial_metadata()
                fields = fieldnames(edfdata);
                my_fields = setdiff(fields, {'left', 'right', 'button'});
            end
            
            % get field names of custom per-experiment meta data
            function my_fields = experiment_metadata()
                fields = fieldnames(edfmeta);
                my_fields = setdiff(fields, {'calib', 'header'});
            end
            
            % get empty trials for each eye
            function empty = empty_trials()
                empty = zeros(2, trials);
                for i=(1:trials)
                    empty(1, i) = isstruct(edfdata(i).left);
                    empty(2, i) = isstruct(edfdata(i).right);
                end
            end
            
            % get eye with best calibration
            function better = eye_calib()
                % init result:
                better = zeros(2, trials);
                % number of calibrations realized
                ncalib = size(edfmeta.calib, 1);
                % running from the back, in which trial we calibrated the last time
                lastc  = size(edfdata, 2) + 1;
                % run over calibrations, from the last one to the first one
                for c=(ncalib:-1:1)
                    % get current calibration
                    cal = edfmeta.calib(c);
                    % is the calibration more recent than the last one?
                    if cal.trial ~= lastc
                        % we have a fresh calibration
                        if isstruct(cal.left) && ((~isstruct(cal.right) || cal.left.err_avg <  cal.right.err_avg) || isnan(cal.right.err_avg)) && ( ~ isnan(cal.left.err_avg) )
                            % left eye has better calibration
                            better(1, max(1, cal.trial):lastc-1) = 0.5;
                            lastc = cal.trial;
                        elseif ~ isnan(cal.right.err_avg)
                            % right eye rocks
                            better(2, max(1, cal.trial):lastc-1) = 0.5;
                            lastc = cal.trial;
                        end % no suitable calibration. we will skip this
                    end
                end
                if lastc~=0
                    display(['Error: no calibration before trial ', int2str(lastc), ' - using left eye']);
                end
            end
            
            % initialize the output structure
            function fixmat = init_struct(fields)
                fixmat.start = zeros(0, 0, 'int32');
                fixmat.stop  = zeros(0, 0, 'int32');
                fixmat.x     = zeros(0, 0, 'single');
                fixmat.y     = zeros(0, 0, 'single');
                fixmat.fix   = zeros(0, 0, 'int32');
                if with_pupil
                    fixmat.pupil = zeros(0, 0, 'single');
                end
                fixmat.eye   = zeros(0, 0, 'uint8');
                for i=(1:length(fields))
                    fixmat.(fields{i}) = zeros(0, 0, 'int32');
                end
            end
            
            % determine which eye to use
            eye = eye_calib() + empty_trials();
            want_eye = 1+(eye(1,:)<eye(2,:));
            eyes = {'left'; 'right'};
            
            % read user defined meta data
            trial_fields = trial_metadata();
            trial_fields_length = length(trial_fields);
            exp_fields = experiment_metadata();
            exp_fields_length = length(exp_fields);
            
            % prepare output structure
            fixmat = init_struct(union(trial_fields, exp_fields));
            
            % main loop: iterate over single trials
            for trial=(1:trials)
                % get fixations of 'better' eye
                fix = edfdata(trial).(eyes{want_eye(trial)}).fixation;
                tstart = fix.start;
                % skip pre-stimulus fixations
                idx = find(tstart>0);
                fix_in_trial = length(idx);
                % output: eye (1 - left, 2 - right)
                fixmat.eye = [fixmat.eye, ones(1, fix_in_trial) * want_eye(trial)];
                % add user defined meta data
                for f=(1:trial_fields_length)
                    % check for empty fields, mark them with -maxint
                    if isempty(edfdata(trial).(trial_fields{f}))
                        edfdata(trial).(trial_fields{f}).msg = '-2147483648';
                    end
                    field_val = str2double(edfdata(trial).(trial_fields{f}).msg);
                    fixmat.(trial_fields{f}) = [fixmat.(trial_fields{f}), ones(1, fix_in_trial) * field_val];
                    
                end
                % add default data fields
                fixmat.start = [fixmat.start, tstart(idx)];
                fixmat.stop  = [fixmat.stop, fix.end(idx)];
                fixmat.x     = [fixmat.x, fix.x(idx)];
                fixmat.y     = [fixmat.y, fix.y(idx)];
                fixmat.fix   = [fixmat.fix,1:fix_in_trial];
                if with_pupil
                    fixmat.pupil = [fixmat.pupil, fix.pupil(idx)];
                end
            end;
            % add user defined meta data
            for f=(1:exp_fields_length)
                fixmat.(exp_fields{f}) = ones(1, length(fixmat.start), 'int32') * str2double(edfmeta.(exp_fields{f}));
            end
            
        end
        function [C, Cr] = get_corrFixbyFix(obj,phases,conditions,fixations)
            %Will run across all subjects and collect fixation maps for
            %each FIXATION, CONDITION and PHASE. This will result in a
            %similarity matrix of
            %[FIXATIONxCONDITIONxPHASE FIXATIONxCONDITIONxPHASE N).
            %Fixations maps are always cocktail blank corrected.
            
            %set blank correction to 1, but will revert it at the end to
            %its initial value.
            old_value = obj.bc;
            obj.bc    = 1;
            c_sub     = 0;
            %
            for ns = unique(obj.subject);
                c_sub   = c_sub +1;
                fprintf('Processing subject %d...\n',ns);
                M = [];
                for phase = phases
                    v       = [];
                    counter = 0;
                    for nfix = fixations
                        for ncond = conditions
                            counter      = counter + 1;
                            v{counter}   = {'phase' phase 'deltacsp' ncond 'subject' ns 'fix' nfix};
                        end
                    end
                    obj.getmaps(v{:});
                    M = cat(3,M,obj.maps);
                end
                obj.maps          = M;
                C(:,:,c_sub)      = obj.cov;
                Cr(:,:,c_sub)     = obj.corr;
            end
            obj.bc = old_value;
        end
         function [count]=EyeNoseMouth(obj,map,normalize)                       
           %% count number of fixations in each roi.
           roi = obj.GetFaceROIs;
           for n = 1:size(roi,3)               
            count(n) = nanmean(map(:).*Vectorize(roi(:,:,n)));
           end
           if nargin>2
               if normalize;
                   fprintf('Will Normalize...\n');
                   count  = count./sum(count);
               end           
           end

         end
         function roi = GetFaceROIs(obj)
             %to visualize
             %roi=fixmat.GetFaceROIs;h=imagesc(fixmat.stimulus);set(h,'alphadata',roi(:,:,1)+roi(:,:,2));
             
             [x y] = meshgrid(1:size(obj.stimulus,2),1:size(obj.stimulus,1));
             %% build rois.
             %coordinates of ROI centers.
             coor = [[140 172 22 14];[360 172 22 14];[255 269 16 25]; [257 425 30 15]];%x and y coordinates for left eye (from my perspective), right eye, nose and mouth.
             for n = 1:size(coor,1)
                 roi(:,:,n) = sqrt(((x-coor(n,1))./coor(n,3)).^2 + ((y-coor(n,2))./coor(n,4)).^2)<4;
             end
             roi(:,:,n+1) = sum(roi,3) == 0;
             %            figure(2);imagesc(obj.stimulus);alpha(sum(roi,3));
             %            roi(:,:,5) = sum(roi(:,:,1:4),3);
             %            for n = 1:size(roi,3);figure(n);h=imagesc(obj.stimulus);set(h,'alphaData',roi(:,:,n));drawnow;end
         end
    end
end

