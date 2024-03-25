

       % pt fin data_chb14_2_cb01_cb03_age_rrmarch2021_mht560_final_age_model_2023_03_13
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')
data_chb     = readtable('data_chb14_2_cb01_cb03_age_rrmarch2021_mht560.txt');
data_chb     = table2array(data_chb);
%data_chb(:,1)= data_chb(:,1) ./ 1000;
data_chb_string_cb0103_rr_560_rr = [
    "data_chb14_2_cb01_cb03_age_rrmarch2021_mht560"
];

% Interpolating data to evenly spaced time axis.
data_chb_age_cb0103_rr_560_rr(:,1) = agemodelmin : agemodelres : agemodelmax;


       data_chb_age_cb0103_rr_560_rr(:,2) = interp1(data_chb(:,1),...
        data_chb(:,2),data_chb_age_cb0103_rr_560_rr(:,1),inttype);
               
clear  data_chb 
