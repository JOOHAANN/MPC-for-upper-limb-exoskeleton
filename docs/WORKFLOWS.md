# Workflow catalog

Each configuration preserves the original local function versions and input filenames.
Use `prepare_workflow(id)` before running an entry script. Do not add `src` recursively to the MATLAB path.

## `demo_calibration__geometry`

Original: `Good demo for training BSE model/Step 1 from Torque to determinate the ODE model`

Scripts: `compare_fd_activation_only_demo`, `demo_theta_ls`, `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_Calibrate_Geo_model`

## `demo_calibration__muscle_tendon`

Original: `Good demo for training BSE model/Step 2 from ODE model to muscle tendon unit`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `demo_joint_damping__geometry`

Original: `Good demo for training BSE model determining jointdamping/Step 1 from Torque to determinate the ODE model`

Scripts: `compare_fd_activation_only_demo`, `demo_theta_ls`, `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_Calibrate_Geo_model`

## `demo_joint_damping__muscle_tendon`

Original: `Good demo for training BSE model determining jointdamping/Step 2 from ODE model to muscle tendon unit`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `map_predictor`

Original: `MAP predictor`

Scripts: `convert_const`, `data_combination`, `EA_convertor_test`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `Lumped_model_check_with_coeffs_thesis_combination`, `MAP_output_model_fusion`, `NNMF_component`, `offline_simulator_NMPC_quasimodo_incremental_MPC_with_MAP`, `test_MAP`, `untitled`

## `control_experimental`

Original: `Simulate the controller with the experimental activation`

Scripts: `Lumped_model_check_with_coeffs_thesis`, `NNMF_component`, `offline_simulator_NMPC_quasimodo_incremental_MPC`, `simplify_pulley_model_fusion_MTU_test_coeffs`

## `control_experimental__period_correct_model_controller`

Original: `Simulate the controller with the experimental activation/Period_correct_model_controller`

Scripts: `offline_simulator_NMPC_quasimodo_period_model`

## `control_map`

Original: `Simulate the controller with the MAP activation`

Scripts: `convert_const`, `data_combination`, `EA_convertor_test`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `Lumped_model_check_with_coeffs_thesis_combination`, `MAP_output_model_fusion`, `NNMF_component`, `offline_simulator_NMPC_quasimodo_incremental_MPC_with_MAP`, `test_MAP`, `untitled`

## `training_bse`

Original: `Train the lumped muscle BSE model`

Scripts: `data_combination`, `data_combination_Qualisys`, `Lumped_model_check_with_coeffs_thesis_combination`, `Lumped_System_Dynamics`, `parameters`, `Sewing_model_with_coeffs`, `Train_Lumped_muscle_model`

## `training_bse__from_theta_to_torque`

Original: `Train the lumped muscle BSE model/From theta to Torque`

Scripts: `Lumped_model_check_with_coeffs`, `Lumped_model_check_with_coeffs_ICNR`, `Lumped_model_check_with_coeffs_thesis`, `Lumped_model_check_with_coeffs_thesis_test`, `Lumped_model_training_from5to15`, `Lumped_System_Dynamics`, `parameters`, `Sewing_model_with_coeffs`

## `training_bse__geometry`

Original: `Train the lumped muscle BSE model/Step 1 From LMT to MTF Calibrate Geometric model`

Scripts: `compare_fd_activation_only_demo`, `demo_theta_ls`, `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_Calibrate_Geo_model`

## `training_bse__muscle_tendon`

Original: `Train the lumped muscle BSE model/Step 2 From EMG and LMT to MTF Calibrate MTU`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_bse__closed_loop`

Original: `Train the lumped muscle BSE model/Step 3 Closed-loop musculoskeletal model`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model`

## `training_qualisys`

Original: `Train the lumped muscle BSE model Qualisys`

Scripts: `data_combination`, `Lumped_model_check_with_coeffs_thesis_combination`, `Lumped_System_Dynamics`, `parameters`, `Sewing_model_with_coeffs`, `Train_Lumped_muscle_model`

## `training_qualisys__from_lmt_to_mtf`

Original: `Train the lumped muscle BSE model Qualisys/From LMT to MTF`

Scripts: `Lumped_model_check_with_coeffs`, `Lumped_model_training_from5to15`, `Lumped_System_Dynamics`

## `training_qualisys__from_mtf_to_lmt`

Original: `Train the lumped muscle BSE model Qualisys/From MTF to LMT`

Scripts: `Lumped_model_check_with_coeffs`, `Lumped_model_training_from5to15`, `Lumped_System_Dynamics`

## `training_qualisys__from_theta_to_torque`

Original: `Train the lumped muscle BSE model Qualisys/From theta to Torque`

Scripts: `Lumped_model_check_with_coeffs`, `Lumped_model_check_with_coeffs_ICNR`, `Lumped_model_check_with_coeffs_thesis`, `Lumped_model_check_with_coeffs_thesis_test`, `Lumped_model_training_from5to15`, `Lumped_System_Dynamics`, `parameters`, `Sewing_model_with_coeffs`

## `training_map__old__geometry`

Original: `Train the lumped muscle BSE model with MAP/old/Step 1 From LMT to MTF Calibrate Geometric model`

Scripts: `compare_fd_activation_only_demo`, `demo_theta_ls`, `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_Calibrate_Geo_model`

## `training_map__old__muscle_tendon`

Original: `Train the lumped muscle BSE model with MAP/old/Step 2 From EMG and LMT to MTF Calibrate MTU`

Scripts: `data_combination`, `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__geometry`

Original: `Train the lumped muscle BSE model with MAP/Step 1`

Scripts: `compare_fd_activation_only_demo`, `demo_theta_ls`, `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_Calibrate_Geo_model`

## `training_map__muscle_tendon__5_35_10_55_15_45_done`

Original: `Train the lumped muscle BSE model with MAP/Step 2/5_35_10_55_15_45_done`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__muscle_tendon__5_35_10_55_15_45_done__build_data`

Original: `Train the lumped muscle BSE model with MAP/Step 2/5_35_10_55_15_45_done/Build data combination`

Scripts: `data_combination`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `MAP_output_model_fusion`

## `training_map__muscle_tendon__5_45_10_35_15_55_done`

Original: `Train the lumped muscle BSE model with MAP/Step 2/5_45_10_35_15_55_done`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__muscle_tendon__5_45_10_35_15_55_done__build_data`

Original: `Train the lumped muscle BSE model with MAP/Step 2/5_45_10_35_15_55_done/Build data combination`

Scripts: `data_combination`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `MAP_output_model_fusion`

## `training_map__muscle_tendon__5_55_10_45_15_35_done`

Original: `Train the lumped muscle BSE model with MAP/Step 2/5_55_10_45_15_35_done`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__muscle_tendon__5_55_10_45_15_35_done__build_data`

Original: `Train the lumped muscle BSE model with MAP/Step 2/5_55_10_45_15_35_done/Build data combination`

Scripts: `data_combination`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `MAP_output_model_fusion`

## `training_map__muscle_tendon__general_model`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General Model`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__muscle_tendon__general_model__build_data`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General Model/Build data combination`

Scripts: `data_combination`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `MAP_output_model_fusion`

## `training_map__muscle_tendon__general_model_op_damp`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_damp`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `MAP_output_model_fusion`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__muscle_tendon__general_model_op_damp__build_data`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_damp/Build data combination`

Scripts: `data_combination`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `MAP_output_model_fusion`

## `training_map__muscle_tendon__general_model_op_scale`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_scale`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `MAP_output_model_fusion`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__muscle_tendon__general_model_op_scale__5_35_10_55_15_45_done`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_scale/5_35_10_55_15_45_done`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__muscle_tendon__general_model_op_scale__5_35_10_55_15_45_done__build_data`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_scale/5_35_10_55_15_45_done/Build data combination`

Scripts: `data_combination`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `MAP_output_model_fusion`

## `training_map__muscle_tendon__general_model_op_scale__5_45_10_35_15_55_done`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_scale/5_45_10_35_15_55_done`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__muscle_tendon__general_model_op_scale__5_45_10_35_15_55_done__build_data`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_scale/5_45_10_35_15_55_done/Build data combination`

Scripts: `data_combination`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `MAP_output_model_fusion`

## `training_map__muscle_tendon__general_model_op_scale__5_55_10_45_15_35_done`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_scale/5_55_10_45_15_35_done`

Scripts: `Lumped_model_check_with_coeffs_thesis_combination`, `Train_Lumped_muscle_model_EMG_and_LMT_to_MTF`

## `training_map__muscle_tendon__general_model_op_scale__5_55_10_45_15_35_done__build_data`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_scale/5_55_10_45_15_35_done/Build data combination`

Scripts: `data_combination`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `MAP_output_model_fusion`

## `training_map__muscle_tendon__general_model_op_scale__build_data`

Original: `Train the lumped muscle BSE model with MAP/Step 2/General_Model_op_scale/Build data combination`

Scripts: `data_combination`, `Fitting_Double_gauss_model`, `Fitting_exp`, `Fitting_Poly`, `MAP_output_model_fusion`, `untitled`
