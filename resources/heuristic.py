import os
import glob

# Heudiconv cannot extract correct information from the dcm, it switches dim4 and dim3
# that's why you have to put dim3 > 30, which is reasonable for functional data (we can go even to 40, just to be sure)
# the 30 here is the number of slices => does not work
# I found out for 4D images, the image_type has to be ('ORIGINAL', 'PRIMARY', 'NON_PARALLEL', 'NONE')
# and the 3d volumes => ('ORIGINAL', 'PRIMARY', 'VOLUME', 'NONE')


# give the DICOM directory and it will navigate to where the files are
# use --files flag with the folder you get from unzipping
# do not use -d flag
# try:
#     # keep going down the tree as long as you are facing dirs
#     while os.path.isdir(os.getcwd()):
#         curr_dir = glob.glob("*")
#         os.chdir(curr_dir[-1])
#
#     # once you face dicoms, go up one level (that's where all data is)
# except NotADirectoryError:
#     os.chdir(os.path.dirname(os.getcwd()))


# ======================================================================================================================
def create_key(template, outtype=('nii.gz',), annotation_classes=None):
    if template is None or not template:
        raise ValueError('Template must be a valid format string')

    return template, outtype, annotation_classes
# ======================================================================================================================
def filter_files(fl):
    if fl.endswith(".dcm"):
        # better to add the SatckId and the rest of info separately, otherwise it messes up the header
        pass
    return fl

# ======================================================================================================================
def infotodict(seqinfo):
    """Heuristic evaluator for determining which runs belong where
    allowed template fields - follow python string module:
    item: index within category
    subject: participant id
    seqitem: run number during scanning
    subindex: sub index within group
    """
    # anatomical
    t2w_tse = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-TSE_run-{item:01d}_T2w')
    t2w_gre = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-GRE_run-{item:01d}_T2w')
    t1w_pregad = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-pregad_run-{item:01d}_T1w')
    t1w_postgad = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-postgad_run-{item:01d}_T1w')
    t2w_rare_orig = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-RARE_rec-orig_run-{item:01d}_T2w')
    t2w_rare_den = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-RARE_rec-denoised_run-{item:01d}_T2w')

    t2starw_mag = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_part-mag_run-{item:01d}_T2starw')
    t2starw_swi = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_rec-SWI_run-{item:01d}_T2starw')
    t2starw_phase = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_part-phase_run-{item:01d}_T2starw')

    tof = create_key('sub-{subject}/{session}/anat/sub-{subject}_{session}_acq-3dtof_run-{item:01d}_angio')
    # ==================================================================================================================
    # resting-state
    func_resting_R = create_key(
            'sub-{subject}/{session}/func/sub-{subject}_task-rest_dir-PA_{session}_run-{item:01d}_bold')
    func_resting_RV = create_key(
            'sub-{subject}/{session}/func/sub-{subject}_task-rest_dir-AP_{session}_run-{item:01d}_bold')

    # ==================================================================================================================
    # diffusion
    dwi = create_key('sub-{subject}/{session}/dwi/sub-{subject}_{session}_run-{item:01d}_dwi')


    # ==================================================================================================================

    info = {t2w_tse: [],
            t2w_gre: [],
            t2w_rare_orig: [],
            t2w_rare_den: [],
            t2starw_mag: [],
            t2starw_swi: [],
            t2starw_phase: [],
            tof: [],

            t1w_pregad:[],
            t1w_postgad:[],
            func_resting_R: [],
            func_resting_RV: [],

            dwi: [],

            }



    # ---------------------------------------------------------------------------------------
    # T2starw SWI have 3 series (Mag, SWI, Phase, then complex data). Complex data (NON_PARALLEL datatype) 
    #  seems to be misaligned spatially so leaveing it out, but others are found by sequential order

    # First: collect candidate T2star / SWI series but exclude NON_PARALLEL types
    t2starw_candidates = []
    for s in seqinfo:
        desc = s.series_description.lower()
        if ('swi' in desc or 't2star' in desc):
            # Exclude those where s.image_type contains 'NON_PARALLEL'
            if not any('NON_PARALLEL' in t for t in s.image_type):
                t2starw_candidates.append(s)

    # Sort by series_id to ensure acquisition order
    t2starw_candidates = sorted(t2starw_candidates, key=lambda x: x.series_id)

    # Assign triplets sequentially: mag, swi, phase
    for i, s in enumerate(t2starw_candidates):
        pos = i % 3
        if pos == 0:
            info[t2starw_mag].append(s.series_id)
        elif pos == 1:
            info[t2starw_swi].append(s.series_id)
        elif pos == 2:
            info[t2starw_phase].append(s.series_id)


    # T2w RARE is also similar, but we have the original, then the denoised
    # First: collect candidate T2w RARE 
    t2w_rare_candidates = []
    for s in seqinfo:
        desc = s.series_description.lower()
        if ('RARE' in desc):
           t2w_rare_candidates.append(s)

    # Sort by series_id to ensure acquisition order
    t2w_rare_candidates = sorted(t2w_rare_candidates, key=lambda x: x.series_id)

    # Assign triplets sequentially: mag, swi, phase
    for i, s in enumerate(t2w_rare_candidates):
        pos = i % 2
        if pos == 0:
            info[t2w_rare_orig].append(s.series_id)
        elif pos == 1:
            info[t2w_rare_den].append(s.series_id)



    # ---------------------------------------------------------------------------------------


    # extract the digits of the name to separate
    # you can even add the videos
    # magnitude data has name: 170001, phase: 170002
    # reverse phase and normal phase
    # we need the json files as well
    # TODO: no of volumes
    # trying to get the denopised version which is usually x0002
    for idx, s in enumerate(seqinfo):
        if 'tse2d' in s.series_description:
            info[t2w_tse].append(s.series_id)

        if 'Gre3D' in s.series_description:
            info[t2w_gre].append(s.series_id)



        elif 'PreGad' in s.series_description:
            info[t1w_pregad].append(s.series_id)

        elif 'PostGad' in s.series_description:
            info[t1w_postgad].append(s.series_id)



        elif 'TOF3D' in s.series_description:
            info[tof].append(s.series_id)
        # ==================================================rest========================================================
        # if the name does not contain "_RV_" then it is a normal phase
        elif  'rsfMRI' in s.series_description:
            if 'RPE' in s.series_description:
                info[func_resting_RV].append(s.series_id)
            else:
                info[func_resting_R].append(s.series_id)

        # ==================================================diffusion========================================================
        if 'DWI' in s.series_description:
            info[dwi].append(s.series_id)


    return info
