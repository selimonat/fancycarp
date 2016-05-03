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