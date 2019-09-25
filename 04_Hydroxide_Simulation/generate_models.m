%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% This function is used to generate a matrix indicating component and species stoichiometry,
% and a vector giving the stability constant of formation for that species or component

% This takes randomly generated constants in matrix "Constants" for
% beta_mono, beta_bis and beta_tris, and rasters log[beta_M(OH)z] from 1-50 as the
% simulation iterates.

% Also runs an intial simulation where there is no competition.

function [species_input_model] = generate_models(model_number,K_competitor,stoichiometry)

Constants = [   15.2178   22.7129   30.1492
   14.0381   20.5300   24.8822
   18.0369   25.9098   29.2689
   11.8528   21.8804   29.7877
   12.7511   16.9667   20.3143
    7.9247   14.8964   19.5875
    6.4594   12.5088   17.0104
   18.8370   23.9409   27.8380
   17.0292   21.8637   25.7777
   14.6326   25.1159   30.5402
   17.5458   26.8506   33.3400
   22.4224   34.0285   37.1844
   24.5027   34.6175   39.5980
   22.9864   29.5456   33.8071
   15.2178   22.7129   30.1492
   23.2579   30.7612   37.1325
   16.6812   22.5992   25.6125
   15.8608   22.2668   27.7013
    8.2270   15.5309   22.0711
   23.0958   31.7863   36.8439
   11.8203   21.1139   27.2071
   16.4637   24.1055   30.0208
   12.3098   20.1555   25.8724
   14.6876   22.1555   27.2026
   17.5615   21.7225   25.0102
   26.8575   33.4868   38.8600
   14.9539   23.2633   29.0766
   15.8516   22.3424   28.2031
   10.2734   14.7637   19.1390
   24.9969   29.5424   32.7600
   24.7240   33.4035   37.2873
   24.0792   34.2867   39.5311
   22.9189   33.9889   38.2538
   19.3615   25.5612   28.8137
   26.1265   36.0539   39.8153
   22.0456   28.0340   32.2162
   17.6636   27.9653   31.8503
   13.0555   20.6508   26.0212
   20.0680   26.0325   31.9348
   24.2577   30.4174   34.2139
    8.7598   17.3369   20.8893
   11.1932   17.8803   23.5029
   18.7164   29.8474   34.4331
    9.7669   17.3374   23.5347
   19.0142   27.2360   30.8883
    9.5951   19.1762   23.9062
   19.8258   26.0547   29.7482
   21.1898   29.0655   33.1632
   22.6364   31.3241   37.2646
   12.6428   20.3981   24.5961];

if K_competitor ~= 0 

    log_beta_mono = Constants(model_number,1);
    log_beta_bis = Constants(model_number,2);
    log_beta_tris = Constants(model_number,3);
    log_beta_competitor = K_competitor - (13.997 .*stoichiometry);
        
    spec_names =		{'H'   'M'   'L'      'OH'          'LM'        'L2M'              'L3M'               'MOH'		};
    Model      =        [ 1      0      0       -1           0            0                  0             -1 .* stoichiometry; ...  % H
                          0      1      0       0            1            1                  1                       1;  ...  % M
                          0      0      1      	0            1            2                  3                       0]; ...  % L
    log_beta   =        [ 0     0      0       -14      log_beta_mono    log_beta_bis     log_beta_tris       log_beta_competitor];
    beta         = 10.^log_beta;

	description = {'Hypothetical Hydroxide Competition'};
	num_solids = 0;

end

if K_competitor == 0 % This is the no competition case. 
    
    log_beta_mono = Constants(model_number,1);
    log_beta_bis = Constants(model_number,2);
    log_beta_tris = Constants(model_number,3);
        
    spec_names =		{   'H'    'M'   'L'        'OH'        'LM'                   'L2M'					 'L3M'};
    Model      =        [	1		0      0		-1           0                         0                    0; ...  % H
                            0		1      0		 0           1                         1                    1;  ...  % M
							0		0      1		 0           1                         2                    3]; ...  % L
    log_beta   =         [	0		0      0	   -14    log_beta_mono				 log_beta_bis     log_beta_tris];
    beta         = 10.^log_beta;

	description = {'Hypothetical Hydroxide Competition'};
	num_solids = 0;
	
end
                
% Get the indices for unbound, mono, bis and tris ligand.
unbound_inc = 1;
mono_inc = 1;
bis_inc = 1;
tris_inc = 1;

[~,species_num] = size(Model);
for current_species = 1:species_num
    L_value = Model(3,current_species);
    M_value = Model(2,current_species);
    if M_value == 0
        if L_value == 1
            unbound_L_indices(unbound_inc) = current_species;
            unbound_inc = unbound_inc + 1;
        end
    end
    if M_value == 1
        if L_value == 1
            mono_L_indices(mono_inc) = current_species;
            mono_inc = mono_inc + 1;
        end
        if L_value == 2
            bis_L_indices(bis_inc) = current_species;
            bis_inc = bis_inc + 1;
        end
        if L_value == 3
            tris_L_indices(tris_inc) = current_species;
            tris_inc =tris_inc + 1;
        end
    end
end

species_input_model.description = description;
species_input_model.spec_names = spec_names;
species_input_model.Model = Model;
species_input_model.log_beta = log_beta;
species_input_model.beta = beta;
species_input_model.unbound_L_indices = unbound_L_indices;
species_input_model.mono_L_indices = mono_L_indices;
species_input_model.bis_L_indices = bis_L_indices;
species_input_model.tris_L_indices = tris_L_indices;
species_input_model.num_solids = num_solids;

Model = Model';
log_beta = log_beta';

species_input_model.calc_spec.KSOLUTION = log_beta(1:end-num_solids,:);
species_input_model.calc_spec.ASOLUTION = Model(1:end-num_solids,:);
species_input_model.calc_spec.KSOLID = log_beta(end-num_solids+1:end,:);
species_input_model.calc_spec.ASOLID = Model(end-num_solids+1:end,:);
species_input_model.calc_spec.SOLUTIONNAMES = char(spec_names(1:end-num_solids));
species_input_model.calc_spec.SOLIDNAMES = char(spec_names(end-num_solids+1:end));

end
