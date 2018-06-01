# fancycarp
A repository for streamlining fMRI preprocessing and analysis based on SPM.

## Fancycrap

![Logo](https://lh3.googleusercontent.com/JiBFuSBsphKFLgHn3DkIR2YQCpP-B8Spvoo3wrW2Rm3HiyC4yHNlSWmxEoLTLkB8Zw=w300)

Repository for streamlining fMRI preprocessing and analysis based on SPM.

## What can I do with this repo?

Create a folder hierarchy for a project.
Download the fMRI data from the internal dicom server (only for internals of IFSN).
Preprocess fmri data with SPM.
Conduct first- and second-level analyses.

It can be integrated together with SCR, Pupil objects for the analysis of autonomic responses.

## How do I proceed?

### 1/ Get the repo
First clone this repository and switch to mrt/main branch.
```shell
onat@neocortex:/tmp$ git clone https://github.com/selimonat/fancycarp.git
Cloning into 'fancycarp'...
remote: Counting objects: 2154, done.
remote: Compressing objects: 100% (4/4), done.
remote: Total 2154 (delta 0), reused 0 (delta 0), pack-reused 2150
Receiving objects: 100% (2154/2154), 986.53 KiB | 801.00 KiB/s, done.
Resolving deltas: 100% (1414/1414), done.
Checking connectivity... done.

onat@neocortex:/tmp/fancycarp$ git checkout mrt/main
Branch mrt/main set up to track remote branch mrt/main from origin.
Switched to a new branch 'mrt/main'
```

If you would like to version-track your repository (which you should) create a new branch following /mrt/XXXX, where XXXX is the name of your project.

### 2/ Add paths
Fire up Matlab, add Fancycarp to Matlab's path. Don't forget to add SPM to your repository.
```
addpath('/tmp/fancycarp')
```

### 3/ Project Object

Fancycarp is OOP. This is convenient as it helps one to define a Project with a set of properties and methods. 
Functions that are specific to a project are declared in the ProjectObject, whereas other functions, which are for example specific to the participants are declared accordingly at their respective objects. 
This removes clutter and helps defining clear functional separation of code. 
ProjectObject contains all sort of functions used for low-level house keeping (e.g. finding the path to a subject at a certain run you can run 
```matlab
>> Project().pathfinder(2,2)
ans =
/tmp/mynextsciencepaper/data//sub002/run002/
```
The role of the ProjectObject is to provide functionality that spans multiple participants at once. For example, at the bottom second level analysis is defined, this is because a project-wide operation that uses all the participants.
Similarly, path finding, dicom transformation, sanity checks, mysql requests are all defined here at the ProjectObject.
<img src="https://github.com/selimonat/fancycarp/images/ProjectObject.png" height="400">


Now, when you start a new project the first thing is to change the parameters that are defined in the Project object, namely the file called Project.m.
```matlab
properties (Hidden, Constant)%adapt these properties for your project
        %All these properties MUST BE CORRECT and adapted to one owns project
        ====> Path to your project folder.
        path_project          = '/tmp/mynextsciencepaper/data/';
        ====> Path to your SPM installation.
        path_spm              = '/home/onat/Documents/Code/Matlab/spm12-6685/';
        ====> Here enter your participants number.
        trio_sessions         = {'PRISMA_19873' 'PRISMA_19875'};
        ====> For every participant enter the number of acquisition runs of TRIO/PRISMA scanner to be downloaded.
        dicom_serie_selector  = {[8 19 20 21 6 7 17 18 ] [8 19 20 21 6 7 17 18 ] };
        ====> Runs [8, 19, 20...] will be distributed to folders [1, 2, 3, ...]
        dicom2run             = repmat({[1:8]},1,length(Project.dicom_serie_selector));
        ====> For field map corrections, enter their runs as stored in the dicom server.
        runs_fieldmap         = [{5 6} {7 8}];%The order is important: first one is the magnitude and the second one is the phase.
        ====> Related to VDM correction. 
        apply_vdm             = [{1}   {2 3 4}];        
        data_folders          = {'midlevel' 'mrt' 'design'};%if you need another folder, do it here.
        ====> Your TR, High-pass filtering values, wheterh to exclude surface extraction (requires CAT12)
        TR                    = 0.99;              
        HParam                = 128;%parameter for high-pass filtering
        surface_wanted        = 0;%do you want CAT12 toolbox to generate surfaces during segmentation (0/1)                
        smoothing_factor      = 4;%how many mm images should be smoothened when calling the SmoothVolume method        
    end
```

### 4/ Data Folders 
Standardization of data folders in a project is the first step to facilitate code sharing and reproduction.
To create a folder hiararchy, run

```matlab
>> Project().CreateFolderHierarchy()
```

This will use Properties of the ProjectObject and loop over all  participants and runs to create a folder structure to download anatomical and functional data. Check it out with the shell command tree.
```matlab
>> !tree                      
.
└── data
    ├── sub001
    │   ├── run000
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run001
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run002
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run003
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run004
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run005
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run006
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run007
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   └── run008
    │       ├── design
    │       ├── midlevel
    │       └── mrt
    ├── sub002
    │   ├── run000
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run001
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run002
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run003
    │   │   ├── design
    │   │   ├── midlevel
    │   │   └── mrt
    │   ├── run004
    │   │   ├── design
```

### 5/ Downloading Anatomical Data

Now we have the folder hierarcy, we can proceed with downloading the anatomical and functional data from the internal server. 
To download data of a given participant, we need to create a Subject instance. 
Thanks to OOP inheritence, SubjectObject receives all the methods defined in the ProjectObject in addition to methods defined in the SubjectObject.
Running the Subject().get_hr method will download the subject's anatomical scans, convert that to NifTi and store them at the 
```subXXX/run000/mrt``` folder.

```matlab
>> s = Subject(1)
Subject Constructor for id:1 is called:
s = 
  Subject with properties:
              id: 1
            path: '/tmp/mynextsciencepaper/data//sub001/'
    trio_session: 'PRISMA_19873'
       total_run: 8
>> s.get_hr                                 
get_hr:
Will now dump the latest HR (16:05:01)
Dicom Server returns:
=====
Database: prisma
Patient: XXXXXX
  * Examination: PRISMA_XXXXX (XXXX, F)                          [ 1|  23|  5192]
    + Study:  1 (2017-11-08 11:15:00)                            [  |  23|  5192] 
      - Series:   1 {localizer (expectb-a, blank)              } [  |    |    11] 
      - Series:   2 {AAHead_Scout_64ch-head-coil               } [  |    |   128] 
      - Series:   3 {AAHead_Scout_64ch-head-coil_MPR_sag       } [  |    |     5] 
      - Series:   4 {AAHead_Scout_64ch-head-coil_MPR_cor       } [  |    |     3] 
      - Series:   5 {AAHead_Scout_64ch-head-coil_MPR_tra       } [  |    |     3] 
      - Series:   6 {gre_field_map, 2mm, filter M              } [  |    |    90] 
      - Series:   7 {gre_field_map, 2mm, filter M              } [  |    |    45] 
      - Series:   8 {ep2d_bold, mb3, loc                       } [  |    |  1193] 
      - Series:   9 {ep2d_diff, mb3, ref                       } [  |    |    45] 
      - Series:  10 {ep2d_diff, mb3, blip inv                  } [  |    |    45] 
      - Series:  11 {mprage, HR64                              } [  |    |   240] 
      - Series:  12 {localizer (expectb-b, blank)              } [  |    |    11] 
      - Series:  13 {AAHead_Scout_64ch-head-coil               } [  |    |   128] 
      - Series:  14 {AAHead_Scout_64ch-head-coil_MPR_sag       } [  |    |     5] 
      - Series:  15 {AAHead_Scout_64ch-head-coil_MPR_cor       } [  |    |     3] 
      - Series:  16 {AAHead_Scout_64ch-head-coil_MPR_tra       } [  |    |     3] 
      - Series:  17 {gre_field_map, 2mm, filter M              } [  |    |    90] 
      - Series:  18 {gre_field_map, 2mm, filter M              } [  |    |    45] 
      - Series:  19 {ep2d_bold, mb3, part 1                    } [  |    |  1061] 
      - Series:  20 {ep2d_bold, mb3, part 2                    } [  |    |   887] 
      - Series:  21 {ep2d_bold, mb3, part 3                    } [  |    |  1061] 
      - Series:  22 {ep2d_diff, mb3, ref                       } [  |    |    45] 
      - Series:  23 {ep2d_diff, mb3, blip inv                  } [  |    |    45] 
=====
The latest recorded HR data:
Series:  11 {mprage, HR64                              } [  |    |   240] 
DicomDownload:
Calling system's COPY function to dump the data...16:05:04
source:/common/mrt32/prisma/images/XXXX
destination:mrt/
COPY finished successully 16:05:04
ConvertDicom:
Found 240 files...
Dicom conversion s#1... (16:05:10)
Will call spm_jobman...
Running SPM jobman 1...
------------------------------------------------------------------------
Running job #1
------------------------------------------------------------------------
Running 'DICOM Import'
   Changing directory to: mrt/
   Changing back to directory: /tmp/fancycarp
Done    'DICOM Import'
Done
Finished... (16:05:15)
Deleting DICOM images in (16:05:15)
mrt/
Finished... (16:05:15)
```
### 6/ Downloading Functional Data
Similarly, ```get_epi``` method downloads the EPIs based on ProjectObject properties.
```matlab
>> s.get_epi
Making a dicom query, sometimes this might take long (so be patient)...(16:05:14)
This is what I found for you:
You told me to download the following series: 8,19,20,21,6,7,17,18,
Double check if everything is fine.
Making a dicom query, sometimes this might take long (so be patient)...(16:05:15)
Will now dump series (16:05:17)
DicomDownload:
Calling system's COPY function to dump the data...16:05:17
source:/common/mrt32/prisma/images/1.3.12.2.1107.5.2.43.167015.30000017110806004710100000013/1.3.12.2.1107.5.2.43.167015.201711081124596897912160.0.0.0
destination:/tmp/mynextsciencepaper/data//sub001/run001/mrt/
```
### 7/ Preprocessing

```s.preprocess_pipeline``` takes care of the preprocessing steps.
It tries to make a fieldmap correction if required files are present. And continues with 
a/ Surface Segmentation of the anatomical data,
b/ Normalization to MNI space,
c/ Realignment and Co-registration,
d/ Gray to White matter segmentation with newSegment,
e/ Strips away skull voxels

```
 function preprocess_pipeline(self,runs)
            %meta method to run all the required steps for hr
            %preprocessing. RUNS specifies the functional runs, make it a
            %vector if needed. RUNS will be used with Re_Coreg.
            if nargin > 1
                try
					fprintf('Will attempt field corrections\n');
					self.ComputeVDM;
                	self.ApplyVDM;
					self.epi_prefix = 'u_';%if we came that far, change the default EPI prefix.
				catch
					fprintf('Failed... Will work on non-field corrected EPIs\n');
				end
                %
	    		self.SegmentSurface_HR;%cat12 segmentation
            	self.SkullStrip;%removes non-neural voxels
            	self.MNI2Native;%brings the atlas (if present) to native space
            	self.Re_Coreg(runs);%realignment and coregistration
            	self.Segment_meanEPI;%segments mean EPI with new segment
				self.SkullStrip_meanEPI;%creates a native mask
	    	else
				fprintf('One input argument is required!\n');
		    end
        end
```        
8/ Analysis of functional data.
So far, Fancycarp covered all the basic preprocessing steps, which should be fairly common to anybody, independent of their projects.
However, first-level analyses in fMRI has certainly project specific flavors, that cannot be fully automated with additional data.
To this end, the CreateFolderHierarchy method has created the ```design``` folder, which is supposed to contain parameters necessary for the building of a design matrix.
One may create as many design models, and use their number as argument to the Subject instance's first-level analysis method. 

```shell
>> pwd
ans =
/mnt/data/project_fearamy/data/sub005/run001/design
>> ls -l
total 136
drwxr-xr-x 2 onat onat 4096 Feb 24  2017 model01
drwxr-xr-x 2 onat onat 4096 Feb 24  2017 model02
drwxrwxr-x 2 onat onat 4096 Feb 13 11:29 model03
drwxr-xr-x 2 onat onat 4096 Feb 24  2017 model04
....
```
In each of these folders there should be a ```data.mat``` file that contains SPM compatible timeing information for different stimuli and conditions.
For example the model below contains 10 condition as a Matlab structure which could be directly fed to the SPM gui.
```
>> cond
cond = 
1x10 struct array with fields:
    name
    onset
    duration
    tmod
    pmod
```

Calling ```s.analysis_firstlevel(1,1)``` should fit a first-level GLM using the model 1 to run 1.

