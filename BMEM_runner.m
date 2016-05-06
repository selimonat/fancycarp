function BMEM_runner(inputs)


for input = inputs
    if input == 1
        clear all;g = Group(62:100);
        magnet = BMEM(g);magnet.prior_function = 'female';
        magnet.fit_ss(1:3);magnet.fit_constantprior;
        r = magnet.fit_quality;
        save ~/Desktop/magnet_Gmale_Pfemale.mat r
        %%
    elseif input == 2
        clear all;g = Group(62:100);
        magnet = BMEM(g);magnet.prior_function = 'male';
        magnet.fit_ss(1:3);magnet.fit_constantprior;
        r = magnet.fit_quality;
        save ~/Desktop/magnet_Gmale_Pmale.mat r
        %%
    elseif input == 3
        clear all;g = Group(101:143);
        magnet = BMEM(g);magnet.prior_function = 'female';
        magnet.fit_ss(1:3);magnet.fit_constantprior;
        r = magnet.fit_quality;
        save ~/Desktop/magnet_Gfemale_Pfemale.mat r
        %%
    elseif input == 4
        clear all;g = Group(101:143);
        magnet = BMEM(g);magnet.prior_function = 'male';
        magnet.fit_ss(1:3);magnet.fit_constantprior;
        r = magnet.fit_quality;
        save ~/Desktop/magnet_Gfemale_Pmale.mat r
        %%
    elseif input == 5
        clear all;g = Group(62:203);
        magnet = BMEM(g);
        magnet.prior_function = 'male';
        magnet.fit_ss(1:3);magnet.fit_constantprior;
        r = magnet.fit_quality;
        save ~/Desktop/magnet_Pmale.mat r
        %%
    elseif input == 6
        clear all;g = Group(62:203);
        magnet = BMEM(g);
        magnet.prior_function = 'female';
        magnet.fit_ss(1:3);magnet.fit_constantprior;
        r = magnet.fit_quality;
        save ~/Desktop/magnet_Pfemale.mat r
        %%
    elseif input == 7
        clear all;g           = Group(62:203);
        magnet                = BMEM(g);
        magnet.prior_function = 'double_category';
        magnet.fit_ss(1:3);magnet.fit_constantprior;
        r = magnet.fit_quality;
        save ~/Desktop/magnet_Pdouble_category.mat r
    end
end
%%
rdouble = load('~/Desktop/magnet_Pdouble_category.mat');
rfemale = load('~/Desktop/magnet_Pfemale.mat');
rmale = load('~/Desktop/magnet_Pmale.mat');
%%
phase = 2;
bar([squeeze(nanmedian(rmale.r.r(:,phase,:))) squeeze(nanmedian(rfemale.r.r(:,phase,:))) squeeze(nanmedian(rdouble.r.r(:,phase,:)))]);
bar([squeeze(nanmedian(rmale.r.NonExpVar(:,phase,:))) squeeze(nanmedian(rfemale.r.NonExpVar(:,phase,:))) squeeze(nanmedian(rdouble.r.NonExpVar(:,phase,:)))])
%%
%%
clear all;
g      = Group(62);
magnet = BMEM(g);
p      = [];
k1c    = 0;
magnet.prior_function  = 'male';
magnet.param_kappa_csp = 3;
magnet.csp             = 0;
c = 0;
params = [];
for k1 = linspace(0.001,25,10)
    magnet.param_kappa_csp = k1;
    k1c = k1c + 1;
    k2c = 0;
    for k2 = linspace(0.001,25,10)
        magnet.param_kappa_perceptual = k2;
        k2c                  = k2c + 1;
        data                 = magnet.p_f_given_csp1;
        magnet.current_data  = data;
        magnet.current_model = 2;
        magnet.SetBoundary;
%         magnet.find_params([1 1 1]);
        c = c +1;
%         subplot(5,5,c)        
        p = FitGauss(magnet.x(:),magnet.p_f_given_csp1,9);
        params(k1c,k2c,:) = p.Est(1:3);        
        drawnow;
    end
end
%%
figure(101);subplot(1,2,1);
imagesc(linspace(0.001,25,25),linspace(0.001,25,25),params(:,:,2));
drawnow;
ylabel('ActiveGeneralization (Kappa)');xlabel('Perceptual Accuracy (Kappa)')
title('Width of Behavioral Generalization')
axis square
axis xy;colorbar







