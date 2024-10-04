function [urls, target_directory] = aj_download_default()

% url_list: cell array of URLs to download (ex: {'url1', 'url2'})
% target_dir: the base BIDS directory where files will be saved

% Liste des URL des fichiers à télécharger
urls = {
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_acq-mprage_T1w.nii.gz?versionId=L7NqLanecid9geIgqDFn113Q7HNL5ykW',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-1_echo-1_FLASH.nii.gz?versionId=Rhx.tPSEmbycRWxA15iA_PswLOz0fuMQ',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-1_echo-2_FLASH.nii.gz?versionId=RFfAPsGJnYGKdoDNcX.ReTy9eHZuV2yn',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-1_echo-3_FLASH.nii.gz?versionId=f76hIrwanLoox2bj6jaUnxqPE0S4neDO',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-1_echo-4_FLASH.nii.gz?versionId=IxU1Uk5jOZ8015ydo8VCjF4K5rTpX_gf',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-1_echo-5_FLASH.nii.gz?versionId=n54jU1WaforbKqr7vKXgmjWwElyuutGr',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-1_echo-6_FLASH.nii.gz?versionId=lJXqy6m40pRtFBoY6rluNexzt5xdQOip',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-1_echo-7_FLASH.nii.gz?versionId=SqIuvVZ3HgXEGaBmzXXx6L9dvIKMMgqy',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-2_echo-1_FLASH.nii.gz?versionId=meM6HPFIjW_skukjOTB.XuNUkptFbxiY',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-2_echo-2_FLASH.nii.gz?versionId=kQwUAnls8mGn9sj9Vm5AEFAtymzFH_gZ',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-2_echo-3_FLASH.nii.gz?versionId=AjvBxtUuA4RRRDTmeQ9kuaLvzVBbqsu.',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-2_echo-4_FLASH.nii.gz?versionId=oQAU0YdkP8SPDs7x26IsvOFlWxEGDQQa',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-2_echo-5_FLASH.nii.gz?versionId=ZcM4MH7pmeU4R_Wlikqby.uJfSHtYUV3',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-2_echo-6_FLASH.nii.gz?versionId=51XytQAOlpfjmIOs1aJ_tzPpVo8XBrC3',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/anat/sub-01_ses-mri_run-2_echo-7_FLASH.nii.gz?versionId=UhvpBuu8Hxc6vCx09QHelwAcDv7VxFwV',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/dwi/sub-01_ses-mri_dwi.nii.gz?versionId=xjT2BuEYylKlqP_T8mrAuomci3tZBK0i',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/fmap/sub-01_ses-mri_magnitude1.nii?versionId=BszE62t5q2ewvX4WB1xEL022IXuoTEXD',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/fmap/sub-01_ses-mri_magnitude2.nii?versionId=ljF9pcc6.iV8x7.aI_GY42hk5JaaY2EI',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/fmap/sub-01_ses-mri_phasediff.nii?versionId=PElStVzxl_c.6QdXDSaz4sxSGfqSj9OR',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/func/sub-01_ses-mri_task-facerecognition_run-01_bold.nii.gz?versionId=0aMsjQSUxgCQedX5FXBysEzTEuQL12fu',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/func/sub-01_ses-mri_task-facerecognition_run-02_bold.nii.gz?versionId=fehHFhVYkgVQ5eWqkXPDJIqUg1j0mC5B',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/func/sub-01_ses-mri_task-facerecognition_run-03_bold.nii.gz?versionId=xVRPWwpqNrbW8gqatHITBhZRxq30npUo',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/func/sub-01_ses-mri_task-facerecognition_run-04_bold.nii.gz?versionId=LK8F8pKUJsl2MQE5NIwnigCYIDaUWYch',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/func/sub-01_ses-mri_task-facerecognition_run-05_bold.nii.gz?versionId=7a2oo8rPmvpzDy.yoy8OVJ91OQKTVvFr',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/func/sub-01_ses-mri_task-facerecognition_run-06_bold.nii.gz?versionId=3LdGTAfKt5dHvWot7mSEWcW1.KE_mW.i',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/func/sub-01_ses-mri_task-facerecognition_run-07_bold.nii.gz?versionId=8zbJmhnQeYMxSTw8gBa.oskhwI44YFaG',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/func/sub-01_ses-mri_task-facerecognition_run-08_bold.nii.gz?versionId=ufc0R73tPc9MX474yxWbhlQNL1R2KDlD',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-01/ses-mri/func/sub-01_ses-mri_task-facerecognition_run-09_bold.nii.gz?versionId=FM.Yr3P7SoJNCbSFVWFT_MYvtOiK0EH5'
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_acq-mprage_T1w.nii.gz?versionId=iPlD.IqXiHfcM.i10d0rk3UXNawqUUoH',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-1_echo-1_FLASH.nii.gz?versionId=NQ_YcyXtpsG15cyeMKw2Ov0FozJTP_hk',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-1_echo-2_FLASH.nii.gz?versionId=HkEQEnScgq1K.uR.LOTcALmOnLXZi8RD',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-1_echo-3_FLASH.nii.gz?versionId=Cf6rzx4LThetKzdm0LMtHtDYkkm_nuBH',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-1_echo-4_FLASH.nii.gz?versionId=17CmaZa1GWyqcecs0di0Ca18klQ5qqvu',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-1_echo-5_FLASH.nii.gz?versionId=czISbo1nZa9a54eTCDj.DYNf0YoLsQYO',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-1_echo-6_FLASH.nii.gz?versionId=eje6yER_WMYJ6q5swHTZ2.BSpvUkO5lj',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-1_echo-7_FLASH.nii.gz?versionId=qciPe3B2X.JZfZMKG5NUyjz19BYU7J7S',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-2_echo-1_FLASH.nii.gz?versionId=bZ_DtjpSUJEKukzySKJBApVK8qyOonc.',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-2_echo-2_FLASH.nii.gz?versionId=HVyyw6oUq9PiCTDK_gRbjYhdpAxYtmP7',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-2_echo-3_FLASH.nii.gz?versionId=Pzn1b01h1PwTEpAcnNz.anmO2Ey_MvwM',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-2_echo-4_FLASH.nii.gz?versionId=qfpOlPGhE.mphoCTwSCYsexSjBB1I5fQ',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-2_echo-5_FLASH.nii.gz?versionId=ZLikisZ7sSRw9OjTvwnDaQ2afVMR4.K5',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-2_echo-6_FLASH.nii.gz?versionId=ufaztB1JaiM775DflJTPisHZh4.iV4YS',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/anat/sub-02_ses-mri_run-2_echo-7_FLASH.nii.gz?versionId=ZYK53Bya9SrKq2W2V3fB_Kz_1gA_LxnW',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/dwi/sub-02_ses-mri_dwi.nii.gz?versionId=iIRbyFkwtw19pWEkWk8x2k87000kAJap',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/fmap/sub-02_ses-mri_magnitude1.nii?versionId=bc5xpsAnB24ac0v1u6PlH6TjR5kwzM1I',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/fmap/sub-02_ses-mri_magnitude2.nii?versionId=mzJ9q0WGw19apvAa5Kx.wAANuMInUF9h',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/fmap/sub-02_ses-mri_phasediff.nii?versionId=vIL6t5kOJU2VfIHHTftBVE0x8dBW4GYl',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/func/sub-02_ses-mri_task-facerecognition_run-01_bold.nii.gz?versionId=UfEVSq6ajLRiy_hXa__wp8VgqX2YX2eC',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/func/sub-02_ses-mri_task-facerecognition_run-02_bold.nii.gz?versionId=ZVKJ903UBQ3XlFoOIpJBWUIdxqMLBiK6',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/func/sub-02_ses-mri_task-facerecognition_run-03_bold.nii.gz?versionId=ZOShdS.yKc0Bxik10hfYgrQMlLmdBBvG',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/func/sub-02_ses-mri_task-facerecognition_run-04_bold.nii.gz?versionId=i2qxeWVNZOnRV1cdj_Z_CO6VmGszpG3E',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/func/sub-02_ses-mri_task-facerecognition_run-05_bold.nii.gz?versionId=63yCzzURKi1D75kKXwlsJ5AarQOH7AHY',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/func/sub-02_ses-mri_task-facerecognition_run-06_bold.nii.gz?versionId=pBj1GWDT5o8LGNNfZXrf91rELIaHOulZ',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/func/sub-02_ses-mri_task-facerecognition_run-07_bold.nii.gz?versionId=cvEO_JaieTPyU02HBvjRhXe0swpvIq4q',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/func/sub-02_ses-mri_task-facerecognition_run-08_bold.nii.gz?versionId=w20s.PnWZBHGV4yKBKC6ax6ND_eBsY6I',
    'https://s3.amazonaws.com/openneuro.org/ds000117/sub-02/ses-mri/func/sub-02_ses-mri_task-facerecognition_run-09_bold.nii.gz?versionId=OMrxe.Id8C7eiBFlKeKrQWEQhGFQ1Bn6'
    };

% Répertoire cible
target_directory = 'C:\Users\antoi\Documents\master_thesis\Data\test2';

end
