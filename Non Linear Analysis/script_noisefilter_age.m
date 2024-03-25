%% Script to call function for nonlinear noise reduction
%
% Trauth, M.H., Asrat, A., Cohen, A., Duesing, W., Foerster, V.,
% Kaboth-Bahr, S., Kraemer, H.,  Lamb, H., Marwan, N., Maslin, M.,
% Schaebitz, F. (2021) Recurring types of variability and transitions in
% the ~620 kyr record of climate change from the Chew Bahir basin, southern
% Ethiopia, Quaternary Science Reviews.
%
% https://doi.org/10.1016/j.quascirev.2020.106777
%
% Required functions and external scripts:
% - nonlin_noise_reduction.m
%
% Expected variables:
% - x (vector with proxy values)
%
% Delivers:
% - x (vector with filtered proxy values)

x = nonlin_noise_reduction(x,11,30);
