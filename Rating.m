classdef Tuning < handle
    properties (Hidden)
        visualization =0;
        options = optimset('Display','none','maxfunevals',10000,'tolX',10^-12,'tolfun',10^-12,'MaxIter',10000,'Algorithm','interior-point');
    end
    %tuning object, can contain any kind of fear-tuning SCR, rating etc.
    properties
        x =[];
        y =[];
        p =[];
    end
    
    methods
        function tuning = Tuning(x,y,varargin)
            tuning.x = x;
            tuning.y = y;
        end
        function fit(tuning,model)
            
        end
        function plot(self)
            plot(
        end
        
        function Fit(self,funtype)
            %Fits FUNTYPE to a tuning 
            %% make sure we have columns
            X        = self.x(:);
            Y        = self.y(:);
            CONSTANT = 0;%will be add to all the data points
            
            %% x-axis, original and upscaled
            xs      = sort(unique(X));            
            xsori   = linspace(xs(1),xs(end),8);
            %% set the function to be fitted
            if funtype == 1
                fitfun = @(X,p) repmat(p(1),length(X),1);
                L      = [ min(Y)  0];%mean sigma
                U      = [ max(Y)  2*std(Y)];
                dof    = 1;
                funname= 'null';
            elseif funtype == 2                
                fitfun = @(X,p) make_gaussian_fmri(X,p(1),p(2),p(3));%2 amp, std, offset
                L      = [-range(Y)*2    0       mean(Y)-range(Y)*2          .01    ];
                U      = [range(Y)*2     20      mean(Y)+range(Y)*2    std(Y(:)+rand(length(Y),1).*eps)*2 ];%
                dof    = 3;
                funname= 'gaussian';                
            elseif funtype == 3                
                fitfun = @(X,p) make_gaussian_fmri_zeromean(X,p(1),p(2));%2 amp fwhm
                L      = [ -range(Y)*2  0     .01    ];
                U      = [  range(Y)*2  20       std(Y(:)+rand(length(Y),1).*eps)*2 ];
                dof    = 3;
                funname = 'gaussian_ZeroMean';
                %detect the mean, store it and subtract it
                CONSTANT = mean(Y);
                Y        = Y-CONSTANT;%we are not interested in the baseline, just remove it so we don't need to estimated it.                
            elseif funtype == 4                
                fitfun = @(X,p) make_gaussian_fmri_tau(X,p(1),p(2),p(3));%amp, tau, offset
                L      = [-range(Y)*2    0    mean(Y)-range(Y)*2          .01    ];
                U      = [range(Y)*2     20  mean(Y)+range(Y)*2    std(Y(:)+rand(length(Y),1).*eps)*2 ];
                dof    = 3;
                funname= 'gaussian_tau';                
            elseif funtype == 5
                fitfun = @(X,p) VonMises_fmri(X, p(1), p(2), p(3));
                L      = [ -range(Y)*2  10.^-16  mean(Y)-range(Y)*2   0];%amp kappa offset
                U      = [  range(Y)*2  20       mean(Y)+range(Y)*2   std(Y)];
                dof    = 3;
                funname= 'vonmises';                
            elseif funtype == 6
                fitfun = @(X,p) make_gabor1d_ZeroMean(X,p(1),p(2),p(3),p(4));%amp std freq baseline
                L      = [ -range(Y)*2  0    1  -range(Y)*2    .01    ];
                U      = [  range(Y)*2  20   4   range(Y)*2  std(Y(:)+rand(length(Y),1).*eps)*2 ];
                dof    = 6;
                funname= 'gabor';
            elseif funtype == 7
                fitfun = @(X,p) p(1)*cos(X*p(2)) + p(3);%amp std freq baseline
                L      = [ -range(Y)*2  1  -range(Y)*2    .01    ];
                U      = [  range(Y)*2  4   range(Y)*2 std(Y(:)+rand(length(Y),1).*eps)*2 ];
                dof    = 3;
                funname= 'cosine';                
            end
            %% set the objective function            
            self.likelihoodfun  = @(params) sum(-log( normpdf( Y - fitfun( X, params(1:end-1)) , 0,params(end)) ));
            
            %% Initial estimation of the parameters            
            %if gabor or gaussian, make a grid-estimatation
            if funtype > 1
                Init = RoughEstimator(X,Y,funname,L,U)';
                %     Init = GaussianEstimate(X,Y);
            else %null model
                Init = mean(Y);
            end
            %% once we have the rough parameter values, create a sigma grid, and
            % estimate the sigma as well.
            tsample = 10000;
            sigmas  = linspace(0.001,std(Y)*3,tsample);            
            %% alternative method: based on the likelihood of sigma given the data points assuming a Gaussian normal distribution.
            PsigmaGiveny = sigmas.^-length(Y) .* exp( - (1./(2.*sigmas.^2) ) .* sum((Y - fitfun(X,Init)).^2) );
            [m i]        = max(PsigmaGiveny);
            %%
            Init   = [Init sigmas(i)];
            U(end) = Init(end)*10;
            %% Optimize!            
            % set optim options            
            [self.Est, self.Likelihood, self.ExitFlag]  = fmincon(self.likelihoodfun, Init, [],[],[],[],L,U,[],self.options);
            %% get the null fit
            null.fitfun     = @(X,p) repmat(p(1),length(X),1);
            null.L          = [ min(Y)  0];%mean sigma
            null.U          = [ max(Y)  2*std(Y)];
            null.Init       = Init;
            null.funny      = @(params) sum(-log( normpdf( Y - null.fitfun( X, params(1:end-1)) , 0,params(end)) ));
            Y               = Y + rand(size(Y))*eps;
            self.null_Likelihood = [];
            [null.Est, self.null_Likelihood, null.ExitFlag]  = fmincon(null.funny, [mean(Y) std(Y)], [],[],[],[], null.L , null.U , [], self.options);
            %% Compute a residual
            %This is the RMS differences between the model and data points...            
            residuals     = Y-self.fitfun(X,Est);
            ss_Residual   = sqrt(mean((Y-self.fitfun(X,Est)).^2));            
            %
            %% produce the fit with the estimated parameters
            fit               = fitfun( xsori  , Est) + CONSTANT;%add the mean we subtracted                        
            %% Get the likelihood
            pval = NaN;
            if ~isempty(null.Likelihood)
                dof   = dof - 1;%the DOF of the null hypothesis is 0.
                pval  = -log10(1-chi2cdf(-2*(self.Likelihood - self.null_Likelihood),dof) + eps);
            end
        end
        function plot_fit
            %% Visualization
            if self.visualization
                c       = 0;
                tsub    = length(X)./8;%to compute SEM we need the number of subjects.
                for nx = xs(:)'%for each unique X variable compute a the average/std/sem values.
                    c               = c+1;
                    Y_ave(c)        = mean(Y(X==nx));
                    Y_ave_trim(c)   = trimmean(Y(X==nx),5);
                    Y_std(c,1)      = std(Y(X==nx));%will be used to plot errorbars
                    Y_sem(c,1)      = Y_std(c,1)./sqrt(tsub);
                end
                figure(100);clf
                plot(xsup,fitup,'ro','linewidth',3);
                hold on
                plot(xsup,fit_init,'color',[.3 .3 .3],'linewidth',3);
                errorbar(xs, Y_ave+CONSTANT, Y_sem, 'b'  ,'linewidth',3);
                hold off
                if funtype > 1
                    title(sprintf('Likelihood: %03g',Likelihood));
                end
                drawnow;
                grid on;
            end
        end
    end
end