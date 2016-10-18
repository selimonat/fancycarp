classdef ROI < Project
    
    properties
        hash
        base_atlas
        index
        name
        model
        pattern_evoked %evoked spatial pattern of activity in this roi
        pattern %spatial pattern of activity for single subjects.
        similarity
        sk   = 6;        
        beta = [];
        selected_subjects
    end
    
    
    
    methods
        function roi = ROI(roi_index,model,betas,subjects)
            %collects the data for the roi and construct the ROI object.
            
            %create variables            
            roi.index                = roi_index;
            roi.model                = model;
            roi.beta                 = betas;%p.model2validbetas;
            roi.selected_subjects    = subjects;
            roi.base_atlas           = roi.path_atlas(roi.index);%path to the atlas
            roi.name                 = roi.get_atlasROIname(roi.index);%the name of the roi
            roi.atlas2mask_threshold = 50;
            %
            
            roi.hash               = DataHash([roi.index roi.atlas2mask_threshold roi.model roi.beta roi.sk roi.selected_subjects.list(:)']);            
            filename               = sprintf('%smidlevel/ROI_%s.mat',roi.path_project,roi.hash);
            if exist(filename) == 0;
                sc                     =  0;
                for ns = roi.selected_subjects.list
                    sc                 = sc + 1;
                    %
                    s                  = Subject(ns);
                    betas              = s.path_beta(1,roi.model,'w_',roi.beta);%(run, model)
                    vol                = spm_vol(betas);
                    vol                = spm_smoothto16bit(vol,roi.sk);%smooth on the fly, slow but disk efficient.
                    
                    XYZmm              = roi.get_XYZmmNormalized(roi.index);
                    XYZvox             = roi.get_mm2vox(XYZmm,vol(1));%in EPI voxel space.
                    %get the BOLD at these coordinates.
                    roi.pattern(:,:,sc)= spm_get_data(vol,XYZvox)';%[voxel condition];                    
                end
                save(filename,'roi');
            else
                load(filename)                
            end
        end        
    end
    methods
        function out = get.pattern_evoked(self)
            %this only makes sense if the data is normalized.
            %so check here if the beta images have a w_ prefix, otherwise
            %return an empty matrix.
            out = mean(self.pattern,3);%compute a mean across subjects.
        end
        
        function out     = get.similarity(self)
            Patterns        = self.pattern_evoked;
            out.euclidian   = squareform(pdist(Patterns','euclidean'));
            out.seuclidian  = squareform(pdist(Patterns','seuclidean'));
%             out.maha        = squareform(pdist(Patterns','mahalanobis',cov));
%             out.cosine      = squareform(pdist(Patterns','cosine'));
            out.rank = squareform(pdist(Patterns','spearman'));
            out.cov         = cov(Patterns);
            out.corr        = corr(Patterns);
            out.mean        = mean(Patterns);
            out.median      = median(Patterns);
            out.var         = var(Patterns);
        end
    end
end