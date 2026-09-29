All scripts have been architeced and edited by Sean Ryan and generated using Claude. 


# Instructions for replication

## Setup

1. On your local machine and in your profile on remission remote desktop, clone the repository. The script is designed so that the repository operates from F:/ since this maximizes compatibility with various software and minimizes the amount that data needs to be moved between remote drives. 

`git clone <>` 

## 1. Export from Ethos

1. Navigate to Ethos treatment manager
2. Navigate to the patient
3. Right click the square indicating the treatment session
4. Click the three dots > Export session data. Export to ETHOS_Simulations/EthosExports
5. Do not let the program time out until the export is complete. This shouldn't be a problem on campus wifi but remote is much slower. 

Parallelization: You can open multiple instances of Ethos and export at once. I think the bottleneck is the Ethos interface and not the user download speed, so having multiple instances seems to be faster than just one. 

## 2. Process Ethos data 

1. Navigate to ./pipeline_setup.m on your local computer. 
2. Edit the list of patients and list of sessions to match the ones you want to process. 
*Make sure the id and session match the directory strings exactly*
3. Run pipeline_setup. Mind any errors. 
4. Verify output: 
a. The directory ETHOS_Simulations/Raystation_input/[id]/[session] should now contain a bunch of CT files, 10-20 RTPLAN files, 2 RTSTRUCT files, and registrations. 
b. The directory ETHOS_Simulations/RayStationFiles/[id]/[session] should now contain an RTPlan file. 

Parallelization: I recommend doing this in batches of patients of sessions. 


## 3. Upload to RS 

### IF it is the first patient: 
1. Import new patient
2. Import all files from ./RayStation_input
3. Rename the patient in RS to add (research) after the last name


### Otherwise

1. Navigate to existing patient 
2. In patient data management, navigate to import
3. In import directory, navigate to ./RayStation_input/[id]/[session]
4. Select all files in the directory. 
5. At the bottom, select "Create new case for existing patient"
6. Import

A lot of errors will probably appear. This should not affect the dose calculation. If there are problems with the dose calculation just make it on the tracking form and I will try to fix it. If you want to debug you can try deleting all case data, reimporting, and saving the import error.

Parallelization: RayStation will let you upload data for multiple patients at once, but not multiple cases for an existing patient. However, you can open many instances of RayStation and have a lot uploading at once. 


## 4. Edit RS parameters for dose calculation

The goal is now to try to do a dose calculation with any plan. Errors will pop up that need to be fixed. Known errors / solutions are: 

1. Double click each CT In patient data management. Set imaging system to "generic ct"
2. Navigate to patient modeling. Select new ROI geometry > Create external ROI. By default it should create an ROI for the Body contour. The default threshold should work. Make sure all CT's are selected when creating ROI.
3. Navigate to plan design. Select dosegrid. Set the default grid (should be 2.5mm resolution.) 
4. Double click the couch contours. Change the material override to the default one supported by raystation i.e. For some reason the ethos imported contours give an undefined material override of the same name, so you need to manually change material Water (1) to Water etc. If the dose calculation keeps giving errors for material overrides, just delete whatever contours you need to until it stops

## 5. Dose calculation

1. Navigate to the sidebar: scripting > Script creation. 
2. Search for ./calc_beam_plan_doses.py. Run the script. 
3. Upon script completion, loosely verify the output. The directory ./RayStationFiles/[id]/[session] should now contain a bunch of CT files, 165 (0 to 164) RTDOSE files for each beam, 2 RTSTRUCT files. 

TO see execution details, navigate to the execution details tab. The script should now run. If execution stops, the error message will be at the very top of the error log and you may need to scroll up to find it. If it is not trivial to resolve, flag it and let me know. 

Parallelization: Same as the last step. Most of the time is spent on the file export working with the RayStation interface. Usually doing parallel dose calculations you would worry about using up clinic resources, but since most of the script time is just the file export and there are only 17 dose calculations over about 2 hours per run of the script, you are safe to open multiple instances and run in parallel. 

## 6. Upload all RayStationFiles files to Remission

I find the best way is to use git bash and scp the session directory to the patient directory on Remission. YOu can also upload via OOD but this is sensitive to interruption when dealing with a lot of files. There may be other better solutions. 

## 7. Process doses on remission. 

1. Navigate to ./ETHOS_Simulations in matlab. 
2. Run pipeline_compress.m
3. Verify there are no output red flags. You can run `verify_pipeline_compress_output.m` as a heuristic. Common ones include zero dose or zero body voxels beacuse it didn't load the RTSTRUCts properly. If there are errors, they will become apparent in the next step 

This should be all that is necessary for this step. Matlab will do some post processing of RS dose files and turn them into .mat arrays in RayStationFiles/[patient]/[session]/processed and keep a log.

## 8. Run kwave simulation

1. Run pipeline_simulate.m

A simulation costs about 40 gpu hours. I find that checking out 4 gpu's is a good balance between job time, accessibility, and not hogging resources. You can divide up these jobs however you see fit. 

