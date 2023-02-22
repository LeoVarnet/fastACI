function [templ_tar,templ_ref, c] = Normalise_two_templates(templ_tar,templ_ref,subfs,templ_num)
% function [templ_tar,templ_ref, c] = Normalise_two_templates(templ_tar,templ_ref,subfs,templ_num)
%
% Used by model_template.m, segmentationSIMU_template.m
%
% Programmed by Alejandro Osses, ENS, 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
templ_tar = templ_tar / templ_num; 
templ_ref = templ_ref / templ_num; 

% The following is a method to simulateously scale both templates to unit energy
[~,c] = Normalise_signal(1/sqrt(2)*[templ_tar; templ_ref],subfs);
templ_tar = c*templ_tar;
templ_ref = c*templ_ref;
