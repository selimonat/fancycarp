classdef Tuning < handle
    properties (Hidden)
        visualization = 0;%visualization of fit results
        gridsize      = 20;%resolution per parameter for initial estimation.
        options       = optimset('Display','none','maxfunevals',10000,'tolX',10^-12,'tolfun',10^-12,'MaxIter',10000,'Algorithm','interior-point');
    end
    %tuning object, can contain any kind of fear-tuning SCR, rating etc.
    properties
        x =[];
        y =[];
        y_mean = [];
        y_std  = [];
        groupfit
        singlesubject
        params;
        pval
    end
    
    methods
        function tuning = Tuning(data,varargin)
            %data is anything that has a x and y fields.
            tuning.x = data.x;
            tuning.y = data.y;
            for x  = unique(tuning.x(:)')
                i             = tuning.x == x;
                tuning.y_mean = [tuning.y_mean mean(tuning.y(i))];
                tuning.y_std  = [tuning.y_std  std(tuning.y(i))];
            end
        end
        
        function SingleSubjectFit(self,funtype)
            %fit FUNTYPE to each individual subject
            for ns = 1:size(self.x,1)
                fprintf('Fitting subject %03d\n',ns)
                self.singlesubject{ns} = self.Fit(self.x(ns,:),self.y(ns,:),funtype);
            end
            self.FitGetParam;
        end
        
        function FitGetParam(self)
            self.params = NaN(length(self.singlesubject),size(self.singlesubject{1}.Est,2));
            for unit = 1:length(self.singlesubject)
                self.params(unit,:) = self.singlesubject{unit}.Est;
                self.pval(unit)     = self.singlesubject{unit}.pval;
            end
        end
        
        function GroupFit(self,funtype)
            %pools different subjects and fit FUNTYPE
            self.groupfit = self.Fit(self.x(:),self.y(:),funtype);
            %
        end
        
        function result = Fit(self,x,y,funtype)
            %Fits FUNTYPE to a tuning defined in x and y
            %for later plotting, check if Group or SingleSubject Fit:
            if size(x,1)>length(unique(x))
                isgroup = 1;
            else
                isgroup = 0;
            end
            %% set the function to be fitted
            x        = x(:);
            y        = y(:);%make it sure to have columns
            y        = y + rand(length(y),1)*eps;
            CONSTANT = 0;%will be added to all the data points
            if funtype == 1
                result.fitfun = @(x,p) repmat(p(1),length(x),1);
                L           = [ min(y)  0];%mean sigma
                U           = [ max(y)  2*std(y)];
                result.dof    = 1;
                result.funname= 'null';
            elseif funtype == 2
                result.fitfun = @(x,p) make_gaussian_fmri(x,p(1),p(2),p(3));%2 amp, std, offset
                L           = [-range(y)*2    0       mean(y)-range(y)*2          .01    ];
                U           = [range(y)*2     180      mean(y)+range(y)*2    std(y(:)+rand(length(y),1).*eps)*2 ];%                
                result.dof    = 3;
                result.funname= 'gaussian';
            elseif funtype == 3
                result.fitfun = @(x,p) self.make_gaussian_fmri_zeromean(x,p(1),p(2));%2 amp fwhm
                L           = [ -range(y)*2  5          .01    ];
                U           = [  range(y)*2  180       std(y(:)+rand(length(y),1).*eps)*2 ];
                result.dof    = 3;
                result.funname= 'gaussian_ZeroMean';
                %detect the mean, store it and subtract it
                CONSTANT    = mean(y);
                y           = y-CONSTANT;%we are not interested in the baseline, just remove it so we don't need to estimated it.
            elseif funtype == 4
                result.fitfun = @(x,p) make_gaussian_fmri_tau(x,p(1),p(2),p(3));%amp, tau, offset
                L           = [-range(y)*2    0    mean(y)-range(y)*2          .01    ];
                U           = [range(y)*2     180  mean(y)+range(y)*2    std(y(:)+rand(length(y),1).*eps)*2 ];
                result.dof    = 3;
                result.funname= 'gaussian_tau';
            elseif funtype == 5
                result.fitfun = @(x,p) VonMises_fmri(x, p(1), p(2), p(3));
                L           = [ -range(y)*2  10.^-16  mean(y)-range(y)*2   0];%amp kappa offset
                U           = [  range(y)*2  180       mean(y)+range(y)*2   std(y)];
                result.dof    = 3;
                result.funname= 'vonmises';
            elseif funtype == 6
                result.fitfun = @(x,p) make_gabor1d_ZeroMean(x,p(1),p(2),p(3),p(4));%amp std freq baseline
                L           = [ -range(y)*2  0    1  -range(y)*2    .01    ];
                U           = [  range(y)*2  180   4   range(y)*2  std(y(:)+rand(length(y),1).*eps)*2 ];
                result.dof    = 6;
                result.funname= 'gabor';
            elseif funtype == 7
                result.fitfun = @(x,p) p(1)*cos(x*p(2)) + p(3);%amp std freq baseline
                L           = [ -range(y)*2  1  -range(y)*2    .01    ];
                U           = [  range(y)*2  4   range(y)*2 std(y(:)+rand(length(y),1).*eps)*2 ];
                result.dof    = 3;
                result.funname= 'cosine';
            elseif funtype == 8
                result.fitfun = @(x,p) self.VonMises(x,p(1),p(2),p(3),p(4));%amp,kappa,centerX,offset
                L             = [ min(y(:))-std(y)     0.1        min(x)   min(y(:))-std(y)   eps ];
                U             = [ max(y(:))+std(y)  range(x)*1.5  max(x)   max(y(:))+std(y)   std(y(:)+rand(length(y),1).*eps)*2 ];
                %                 L      = [ eps                   0.1   eps     -pi   eps ];
                %                 U      = [ min(10,range(y)*1.1)  20   2*pi   pi   10];
                result.dof    = 3;
                result.funname= 'vonmisses_mobile';
            end
            %% set the objective function
            result.likelihoodfun  = @(params) sum(-log( normpdf( y - result.fitfun( x,params(1:end-1)) , 0,params(end)) ));
            
            %% Initial estimation of the parameters
            %if gabor or gaussian, make a grid-estimatation
            if funtype > 1
                Init = self.RoughEstimator(x,y,result.fitfun,L(1:end-1),U(1:end-1));%[7.1053 15.5556 1.0382];
            else %null model
                Init = mean(y);
            end
            
            %% estimate initial values for sigma_noise
            %based on the likelihood of sigma given the data points assuming a Gaussian normal distribution.
            tsample      = 1000;
            sigmas       = linspace(0.001,std(y)*3,tsample);
            PsigmaGiveny = sigmas.^-length(y) .* exp( - (1./(2.*sigmas.^2) ) .* sum((y - result.fitfun(x,Init)).^2) );
            [m i]        = max(PsigmaGiveny);
            Init         = [Init sigmas(i)];%
            %% Optimize!
            try
                [result.Est, result.Likelihood, result.ExitFlag]  = fmincon(result.likelihoodfun, Init, [],[],[],[],L,U,[],self.options);
                result.Likelihood = result.likelihoodfun(result.Est);
            catch
                result.Est        = Init;
                result.Likelihood = result.likelihoodfun(result.Est);
                result.ExitFlag   = 1;
            end
            %% get the null fit
            null.fitfun     = @(x,p) repmat(p(1),length(x),1);
            null.L          = [ min(y)  0];%mean sigma
            null.U          = [ max(y)  2*std(y)];
            null.Init       = Init;
            null.funny      = @(params) sum(-log( normpdf( y - null.fitfun( x, params(1:end-1)) , 0,params(end)) ));
            y               = y + rand(size(y))*eps;
            [null.Est, result.null_Likelihood, null.ExitFlag]  = fmincon( null.funny, [mean(y) std(y)], [],[],[],[], null.L , null.U , [], self.options);
            %% Compute a residual
            %This is the RMS differences between the model and data points...
            result.residuals      = y-result.fitfun(x,result.Est);
            result.ss_residuals   = sqrt(mean((y-result.fitfun(x,result.Est)).^2));
            %% Get the likelihood
            if ~isempty(result.null_Likelihood)
                result.dof   = result.dof - 1;%the DOF of the null hypothesis is 0.
                result.pval  = -log10(1-chi2cdf(-2*(result.Likelihood - result.null_Likelihood),result.dof) + eps);
            end
            result.x       = unique(x);
            x_HD           = linspace(min(result.x),max(result.x),100);
            result.fit     = result.fitfun(result.x,result.Est);
            
            %% show fit if wanted
            if self.visualization
                
                c       = 0;
                tsub       = size(self.x,1);%to compute SEM we need the number of subjects.
                xs = sort(unique(self.x));
                for nx = xs(:)'%for each unique X variable compute a the average/std/sem values.
                    c               = c+1;
                    Y_ave(c)        = mean(self.y(self.x==nx));
                    Y_ave_trim(c)   = trimmean(self.y(self.x==nx),5);
                    Y_std(c,1)      = std(self.y(self.x==nx));%will be used to plot errorbars
                    Y_sem(c,1)      = Y_std(c,1)./sqrt(tsub);
                end
                
                figure(100);clf
                plot(x_HD,result.fitfun(x_HD,result.Est),'ro','linewidth',3);
                hold on
                plot(x_HD, result.fitfun(x_HD,Init)  ,'color',[.3 .3 .3] ,'linewidth',3);
                
                if isgroup
                    errorbar(result.x, Y_ave+CONSTANT, Y_sem   , 'b'   ,'linewidth', 3);
                else
                    plot(result.x, y+CONSTANT, 'b'   ,'linewidth', 3);
                end
                
                hold off
                if funtype > 1
                    title(sprintf('Likelihood: %03g (p = %5.5g)',result.Likelihood,result.pval));
                end
                xlim([min(x(:)) max(x(:))]);
                drawnow;
                grid on;
            end
        end
        
        function [params]=RoughEstimator(self,x,y,fun,L,U)
            %will roughly estimate free parameter values bounded between L
            %and U for FUN and x values. Error is computed according to Y.        
            %% get a grid            
            tparam   = length(L);
            for n = 1:tparam
                grid_in{n} = linspace(L(n),U(n),self.gridsize);
            end
            G      = cell(1,tparam);
            [G{:}] = ndgrid(grid_in{:});%generate a grid
            G      = cat(tparam+1,G{:});%cell to matrix
            G      = permute(G,[tparam+1 1:tparam]);%change dimension orders so that
            G      = reshape(G,[tparam,self.gridsize^tparam])';%we can make Nx3 matrix with reshape
            %% compute the error for each parameter combination
            error  = zeros(1,self.gridsize^tparam);
            for npoint = 1:self.gridsize^tparam;
                error(npoint) = sum((y - fun(x,G(npoint,:))).^2);%residual error
            end
            [m i]  = min(error);
            params = double(G(i,:));            
        end
    end
    methods (Static)
        function  [out] = VonMises(X,amp,kappa,centerX,offset)
            %[out] = VonMises(X,amp,centerX,kappa,offset)
            %
            %This is a variation of the von Mises distribution explained here:
            %http://en.wikipedia.org/wiki/Von_Mises_distribution.
            out    =  (exp(kappa*cos(deg2rad(X-centerX)))-exp(-kappa))./(exp(kappa)-exp(-kappa));%put it btw [0 and 1]
            out    =  amp*out + offset;%and now scale it
        end
        function [out] = make_gaussian_fmri_zeromean(x,amp,sd)            
            %
            %	Generates a Gaussian with AMP, FWHM and offset parameters. The center location is
            %	fixed to zero. 
                                                
            d = mean(diff(x));                
            %for speed length(X) can be precomputed
            XT = length(x);
            %as well as d
            % d  = 0.02;
            %out      = amp.*exp(-tau*(x.^2)/2) - amp*sqrt(2*pi/tau)./d./XT;
            
            out      = amp.*exp(-(x./sd).^2/2) - amp*sqrt(2*pi*sd.^2)./d./XT;
        end
    end
end