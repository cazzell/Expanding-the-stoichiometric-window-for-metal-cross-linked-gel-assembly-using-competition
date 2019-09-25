%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the master driver that calls other functions to generate
% theoretical predictions for metal-coordinated hydrogels of linear,
% telechelic nitrocatechol coordinated with iron or aluminum and histidine
% coordinated with nickel. Predicts gel plateau moduli as a function of
% metal concentration and pH.

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% This code calls functions modified or taken from the work of Donald Scott Smith,
% which is covered by the BSD 2 Clause License. Those functions include his BSD 2 license.

% Functions called in this code are directly or modified from Donald Scott Smith
% see his work
% Smith, D. Scott, "Solution of Simultaneous Chemical Equilibria in Heterogeneous Systems: Implementation in Matlab" (2019).
% Chemistry Faculty Publications. 14. https://scholars.wlu.ca/chem_faculty/14

% see reference 25 in "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% cazzell.lbi@gmail.com

%% This is the master program. Here the variable to modify is the model_number.
% Pick the model_number that corresponds to which ligand-metal combination you want to examine

% model_number descriptions

% model_number = 1
% iron and nitrocatechol without hydroxide interactions

% model_number = 2
% iron and nitrocatechol with hydroxide interactions

% model_number = 3
% aluminum and nitrocatechol without hydroxide interactions

% model_number = 4
% aluminum and nitrocatechol with hydroxide interactions

% model_number = 5
% nickel and histidine without hydroxide interactions

% model_number = 6
% nickel and histidine with hydroxide interactions

%%
close all
clear
clc

for model_number = 1:6
	
	warning('off','all')
	
	% For Production Run
 	%metal_ratios = linspace(0,5/3,251);
 	%pH = linspace(0,14,281);
	% For Quick Run
 	metal_ratios = 0:  (1/3)/10  :5/3;
 	pH = 0:0.1:14;
	
	
	cd Models
	
	model_str = num2str(model_number);
	cd_str = [num2str(0), model_str];
	cd(cd_str)
	
	% Generates model used to describe the desired system for the speciation program
	[species_input_model] = generate_models(model_number)
	
	% Runs speciation program to predict species concentrations vs. pH and metal concentration
	[speciation] = calc_speciation(species_input_model,pH,metal_ratios,model_number)
	
	% Runs mechanical prediction program to predict plateau moduli from speciation prediction
	[speciation] = calc_Gp_linear_Miller_Macosko(speciation,species_input_model,metal_ratios,pH)
	
	% Function to organize and generate plots
	[spec_gp_surf] = plot_spec_gp_surf(metal_ratios,model_number,species_input_model,speciation,pH)
	
	% Code to save generated data
	save('structures.mat','species_input_model','speciation','spec_gp_surf')
	
	close all
	
	cd ..
	cd ..
	
end