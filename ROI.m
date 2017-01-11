classdef ROI < Project
    %ROI object. Can collect betas values for different atlas based or
    %voxel based specified areas.
    properties
        hash
        index
        name
        model
        evoked_pattern%evoked spatial pattern of activity in this roi
        pattern %spatial pattern of activity for single subjects.
        sk   = [];
        beta = [];
        selected_subjects
        evoked_similarity
        induced_similarity
        XYZmm
    end
    
    
    
    methods
        function roi = ROI(voxel_selection,ROI_OR_Voxels,default_model,model_number,betas,subjects,sk)
            %roi = ROI(voxel_selection,ROI_OR_Voxels,default_model,model_number,betas,subjects)
            %
            %input_type: roi_based or voxel_based as a string.
            %
            %collects the data for the roi and construct the ROI object.
            %input
            %
            %
            %EXAMPLE:
            %r = ROI('voxel_based',[0 0 0 1]','chrf_0_0',0,1:8,5)
            tic;
            force = 1;
            %store some variables
            roi.model                = model_number;
            roi.sk                   = sk;
            roi.beta                 = betas;%p.model2validbetas;
            roi.selected_subjects    = subjects(:)';
            %get a cache name
            roi.hash                 = DataHash(getByteStreamFromArray([uint8(voxel_selection),ROI_OR_Voxels(:)',uint8(default_model),roi.atlas2mask_threshold roi.model roi.beta(:)' roi.sk]));
            %
            %compute and save
            if GetVoxelsXYZmm;%get coordinates and continue if success
                roi.XYZmm              = XYZmm;
                sc                     =  0;
                for ns = roi.selected_subjects%run across subjects
                    sc                   = sc + 1;
                    s                    = Subject(ns);
                    s.default_model_name = default_model;
                    filename             = sprintf('%sROI/ROI_%s.mat',s.path_midlevel(1),roi.hash);
                    if exist(fileparts(filename)) == 0
                        mkdir(fileparts(filename));
                    end
                    cprintf([0 1 0],'File name:\n%s\n',filename);
                    if exist(filename) == 0 | force;
                        cprintf([1 0 0],'Not yet cached, will compute...\n');
                        betas                = s.path_beta(1,roi.model,'w_',roi.beta)%(run, model)
                        vol                  = spm_vol(betas);
                        vol                  = spm_smoothto16bit(vol,roi.sk);%smooth on the fly, slow but disk efficient.                        
                        %get the BOLD at these coordinates.
                        for NNN = 1:size(XYZmm,2)
                            XYZvox               = roi.get_mm2vox(XYZmm(:,NNN),vol(1));%in EPI voxel space.
                            data(NNN,:)          = spm_get_data(vol,XYZvox)';
                        end
                        save(filename,'data');
                        cprintf([0 1 0],'Saved under %s...\n',roi.hash);
                    else
                        cprintf([0 1 0],'Loading from cache...\n',roi.hash);
                        dummy = load(filename);
                        data  = dummy.data;
                    end
                    roi.pattern(:,:,sc)= data;%[voxels conditions (i.e beta)];                                        
                end
            else
                cprintf([1 0 0],'Not doing anything here\n');
            end
            cprintf([0 0 1],'Finished in %1.2g seconds.\n',toc);
                                  
            function ok = GetVoxelsXYZmm
                %returns voxels in MM space.
                ok = 1;
                if strcmp(voxel_selection,'roi_based');
                    cprintf([0 0 1],'ROI (%03d) based voxel selection.\n',ROI_OR_Voxels);
                    %roi given, find all voxels which are valid based on
                    %atlas2mask threshold, and fill in ROI info.
                    XYZmm                    = roi.get_XYZmmNormalized(ROI_OR_Voxels);
                    roi.index                = ROI_OR_Voxels;
                    %roi.base_atlas           = roi.path_atlas(roi.index);%path to the atlas
                    roi.name                 = roi.get_atlasROIname(roi.index);%the name of the roi
                    %
                elseif strcmp(voxel_selection,'voxel_based');
                    cprintf([0 0 1],'Direct voxel selection.\n');
                    %voxel is directly given, no need to feedle around with
                    %the atlas. Added possibility to give more than one
                    %voxel
                    XYZmm              = ROI_OR_Voxels;
                    roi.name = '';
                    for coor = XYZmm
                        roi.name           = strvcat(roi.name,roi.get_atlaslabel(coor(1:3)).name{1});%the name of the roi
                    end
                    if ~all(size(XYZmm,1) == 4);
                        cprintf([1 0 0],'WRONG size of voxel inputs\n');
                        XYZmm = [];
                        ok = 0;
                    end
                else
                    fprintf('Unknown option\n')
                    return
                end
            end
        end
    end
    methods
        
        function out = get.evoked_pattern(self)
            %this only makes sense if the data is normalized.
            %so check here if the beta images have a w_ prefix, otherwise
            %return an empty matrix.
            out = mean(self.pattern,3);%compute a mean across subjects.
        end
        
        function out     = get.evoked_similarity(self)
            %computes an evoked similarity matrix
            
            out = self.similarity(self.evoked_pattern);
        end
        
        function out     = get.induced_similarity(self)
            %will do it separately for each subject first and then average
            %the results.
            for ns = 1:size(self.pattern,3)
                out2(ns) = self.similarity(self.pattern(:,:,ns));
            end
            out = [];
            for fn = fieldnames(out2)'
                out.(fn{1}) = mean(cat(3,out2(:).(fn{1})),3);
            end
        end
        
    end
    
    methods (Static)
        function out = similarity(Patterns)
            if size(Patterns,1) > 1
                out.rank        = squareform(pdist(Patterns','spearman'));
                out.cov         = cov(Patterns);
                out.corr        = squareform(pdist(Patterns','correlation'));
            end
            out.mean        = mean(Patterns,1);
            out.var         = var(Patterns,1,1);
            %             out.median      = median(Patterns);
            %             out.euclidian   = squareform(pdist(Patterns','euclidean'));
            %             out.seuclidian  = squareform(pdist(Patterns','seuclidean'));
            %             out.maha        = squareform(pdist(Patterns','mahalanobis',cov));
            %             out.cosine      = squareform(pdist(Patterns','cosine'));
        end
    end
end