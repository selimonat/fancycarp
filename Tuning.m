classdef Tuning < handle
    properties (Hidden)
        visualization =1;
        options = optimset('Display','none','maxfunevals',10000,'tolX',10^-12,'tolfun',10^-12,'MaxIter',10000,'Algorithm','interior-point');
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
        
        %         Est
        %         Likelihood
        %         ExitFlag
        %         pval = NaN;
        %         dof
        %
        %         ss_residuals
        %         residuals
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
                result.fitfun = @(x,p) make_gaussian_fmri_zeromean(x,p(1),p(2));%2 amp fwhm
                L           = [ -range(y)*2  0     .01    ];
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
                Init = self.RoughEstimator(x,y,result.funname,L,U)';
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
                
                %
                %                 tsub       = size(self.x,1);
                %                 Y_ave      = mean(self.y);
                %                 Y_std      = std(self.y);
                %                 Y_sem      = Y_std./sqrt(tsub);
                
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
                drawnow;
                grid on;
            end
        end
        
        function [params]=RoughEstimator(self,x,y,funname,L,U)
            %%[params]=RoughEstimator(x,y,funtype,L,U)
            %
            %   Will roughly estimate the parameters of a Gaussian or Gabor function
            %   using a grid approach. Generally these initial estimates are a bit
            %   slower compared to more heuristic approaches however they are nearly as
            %   good as the result of the gradient descent which makes the optimization
            %   use less iteration. X and Y are the data points and FUNTYPE follows the
            %   convention in FitGauss.m and it works for either Gabor or Gaussian
            %   functions. The estimation uses all the data points, not average at
            %   given angles. L and U are the lower/upper bounds of the grid. Pay
            %   attention that the argument orders corresponds to the function selected
            %   by the funtype.
            %
            %   Selim Onat
            
            grid = 10;
            %%
            
            
            %create a reportoire of Gaussians or Gabors.
            if strcmp(funname,'gaussian_ZeroMean')
                %%
                g    = single(zeros(length(x),50*25));
                lut  = single(zeros(length(L)-1,size(g,2)));
                %
                amps      = linspace(L(1),U(1),grid*2);
                widths     = linspace(L(2),U(2),grid);
                %
                c       = 0;
                for amp = amps
                    for width = widths
                        c        = c+1;
                        lut(:,c) = [amp;width];
                        dummy    = make_gaussian_fmri_zeromean(x,amp,width);
                        g(:,c)   = dummy;
                    end
                end
                
            elseif strcmp(funname,'gaussian')
                %%
                g    = single(zeros(length(x),50*25));
                lut  = single(zeros(length(L)-1,size(g,2)));
                %
                amps      = linspace(L(1),U(1),grid*2);
                widths     = linspace(L(2),U(2),grid);
                baselines = linspace(L(3),U(3),grid);
                %
                c       = 0;
                for amp = amps
                    for width = widths
                        for baseline = baselines
                            c        = c+1;
                            lut(:,c) = [amp;width;baseline];
                            dummy    = make_gaussian_fmri(x,amp,width,baseline);
                            g(:,c)   = dummy;
                        end
                    end
                end
                
            elseif strcmp(funname,'gaussian_tau')
                %%
                g    = single(zeros(length(x),50*25));
                lut  = single(zeros(length(L)-1,size(g,2)));
                %
                amps       = linspace(L(1),U(1),grid*2);
                widths     = linspace(L(2),U(2),grid);
                baselines  = linspace(L(3),U(3),grid);
                %
                c       = 0;
                for amp = amps
                    for width = widths
                        for baseline = baselines
                            c        = c+1;
                            lut(:,c) = [amp;width;baseline];
                            dummy    = make_gaussian_fmri_tau(x,amp,width,baseline);
                            g(:,c)   = dummy;
                        end
                    end
                end
                
                
            elseif strcmp(funname,'vonmises')
                %%
                g    = single(zeros(length(x),50*25));
                lut  = single(zeros(length(L)-1,size(g,2)));
                %
                amps      = linspace(L(1),U(1),grid*2);
                widths     = logspace(L(2),U(2),grid);
                baselines = linspace(L(3),U(3),grid);
                %
                c       = 0;
                for amp = amps
                    for width = widths
                        for baseline = baselines
                            c        = c+1;
                            lut(:,c) = [amp;width;baseline];
                            dummy    = make_gaussian_fmri_tau(x,amp,width,baseline);
                            g(:,c)   = dummy;
                        end
                    end
                end
                
                
                
                
            elseif strcmp(funname,'gabor')
                %%
                grid = 15;
                g    = single(zeros(length(x),grid*grid*10*1*grid));
                lut  = single(zeros(length(L)-1,size(g,2)));
                %
                amps      = linspace(L(1),U(1),grid);
                widths    = linspace(L(2),U(2),grid);
                freqs     = linspace(L(3),U(3),grid);
                baselines = linspace(L(4),U(4),grid);
                %
                c       = 0;
                for amp = amps
                    for width = widths
                        for freq = freqs;
                            for baseline = baselines;
                                c        = c+1;
                                lut(:,c) = [amp;width;freq;baseline];
                                dummy    = make_gabor1d_ZeroMean(x,amp,width,freq,baseline);
                                g(:,c)   = dummy;
                            end
                        end
                    end
                end
                
            elseif strcmp(funname,'cosine')
                %%
                grid = 15;
                g    = single(zeros(length(x),grid*grid*10*1*grid));
                lut  = single(zeros(length(L)-1,size(g,2)));
                %
                amps       = linspace(L(1),U(1),grid);
                freqs      = linspace(L(2),U(2),grid);
                baselines  = linspace(L(3),U(3),grid);
                %
                c       = 0;
                for amp = amps
                    for freq = freqs
                        for baseline = baselines
                            c        = c+1;
                            lut(:,c) = [amp;freq;baseline];
                            dummy    = amp*cos(x*freq)+baseline;
                            g(:,c)   = dummy;
                        end
                    end
                end
                
            elseif strcmp(funname,'vonmisses_mobile')
                %%
                g    = single(zeros(length(x),grid*grid*10*1*grid));
                lut  = single(zeros(length(L)-1,size(g,2)));
                %
                amps      = linspace(L(1),U(1),grid);
                widths    = linspace(L(2),U(2),grid);
                freqs     = linspace(L(3),U(3),grid);
                baselines = linspace(L(4),U(4),grid);
                %
                c       = 0;
                for amp = amps
                    for width = widths
                        for freq = freqs;
                            for baseline = baselines;
                                c        = c+1;
                                lut(:,c) = [amp;width;freq;baseline];                                
                                dummy    = self.VonMises(x,amp,width,freq,baseline);
                                g(:,c)   = dummy;
                            end
                        end
                    end
                end
                
            end
            
            
            %%
            res   = (sum( (g - repmat((y(:)),1,size(g,2))).^2 ));
            [m i] = min(res);
            
            params = double(lut(:,i));
            
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
    end
end