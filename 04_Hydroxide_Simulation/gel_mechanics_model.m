%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the master driver that calls other functions to generate
% theoretical predictions for metal-coordinated hydrogels of simulated linear,
% telechelic ligand/metal combinations. Predicts gel plateau moduli as a function of
% metal concentration and pH as a function of MOH competitor strength.

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

%% This program is used to study the role of a single hydroxide competitor's equilibrium constant and stoichiometry
% on the gel's plateau moduli as a function of metal concentration and pH

% This is the master program. Here the variable to modify is the model_number, the equilibrium constant for the competitor,
% and the stoichiometry of the competitor

close all
clear
clc

warning('off','all')

% model_number descriptions

%K_competitor_list = [1:60]
%stoichiometry_list = [1:50]

% K_competitor_list = [0:50]
% stoichiometry_list = [1:3]

K_competitor_list = 0:50;
stoichiometry_list = 1:3;
pH = 0:0.1:14;


for model_number = [1,3,20]
    
    metal_ratios = [0:  (1/3)/10  :5/3];
        
    cd Models
        
	model_str = num2str(model_number);
	
	if model_number < 10
		cd_str_mod = [num2str(0), model_str];
	else
		cd_str_mod = [model_str];
	end
		
	cd(cd_str_mod)
        
    for stoichiometry = stoichiometry_list
        
        stoich_str = num2str(stoichiometry);
        cd_str_stoich = [num2str(0), stoich_str];
        cd(cd_str_stoich)
        
        for K_competitor = K_competitor_list
            
            K_c_str = num2str(ceil(K_competitor));
            if K_competitor < 10
                cd_str_K_c = [num2str(0), K_c_str];
            else
                cd_str_K_c = [K_c_str];
            end
            cd(cd_str_K_c)
            
            % Generates model used to describe the desired system for the speciation program
            [species_input_model] = generate_models(model_number,K_competitor,stoichiometry)
    
            % Runs speciation program to predict species concentrations vs. pH and metal concentration
			[speciation] = calc_speciation(species_input_model,pH,metal_ratios,model_number)
						
            % Runs mechanical prediction program to predict plateau moduli from speciation prediction
            [speciation] = calc_Gp_linear_Miller_Macosko(speciation,species_input_model,metal_ratios,pH)
    
            % Function to organize and generate plots
            [spec_gp_surf] = plot_spec_gp_surf(metal_ratios,model_number,species_input_model,speciation,pH,stoichiometry,K_competitor)
            
            % Code to save generated data
			save('structures.mat','spec_gp_surf')
            
            close all
            cd ..
            end
        cd ..
    end
    cd ..
    cd ..
end