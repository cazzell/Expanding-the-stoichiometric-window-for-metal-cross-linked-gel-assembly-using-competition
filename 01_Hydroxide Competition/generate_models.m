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

function [species_input_model] = generate_models(model_number)

% Stability Constants are taken from Smith, R. M. & Martell, A. E. "Critical Stability Constants," volume indicated

% Protonation
% Nitrocatechol
% HL/H.L
n.pKa2 = 10.85; % Volume 6
% H2L/HL.H
n.pKa1 = 6.69; % Volume 6
n.H1 = n.pKa2;
n.H2 = n.pKa1 + n.pKa2;

% Histamine
% HL/H.L
h.pKa2 = 9.83; % Volume 2
% H2L/HL.H
h.pKa1 = 6.07; % Volume 2
h.H1 = h.pKa2;
h.H2 = h.pKa1 + h.pKa2;

% Metals-Ligand
% Nitrocatechol
% Iron
n_Fe.K1 = 17.1; % Volume 5
n_Fe.K2 = 30.5; % Volume 5
n_Fe.K3 = 40;    % Volume 5
% Aluminum
n_Al.K1 = 13.74; % Volume 6
n_Al.K2 = 25.4;   % Volume 6
n_Al.K3 = 34.3;   % Volume 6

% Histamine
% Nickel
h_Ni.K1 = 6.78;   % Volume 2
h_Ni.K2 = 11.78; % Volume 2
h_Ni.K3 = 14.9;   % Volume 2

% Hydroxide-Metal
% Iron
Fe.MOH = -2.187;     % Volume 4
Fe.MOH2 = -5.694;   % Volume 4
Fe.MOH3 = -3.191;   % Volume 4
Fe.MOH4 = -21.588; % Volume 4
Fe.M2OH2 = -2.894; % Volume 4
Fe.M3OH4 = -6.288; % Volume 4


% Aluminum
Al.MOH = -4.987;        % Volume 4
Al.MOH2 = -9.294;      % Volume 4
Al.MOH3 =-14.991;     % Volume 4
Al.MOH4 = -22.998;    % Volume 4
Al.M2OH2 = -7.694;    % Volume 4
Al.M3OH4 = -13.888;  % Volume 4
Al.M1OH3 = -8.491;    % Volume 4

% Nickel
Ni.MOH = -9.897;      % Volume 4
Ni.MOH2 = -19.994;   % Volume 4
Ni.MOH3 = -30.991;   % Volume 4
Ni.M2OH = -10.697;   % Volume 4
Ni.M4OH4 = -27.688; % Volume 4
Ni.M1OH2 = -12.794; % Volume 4


if model_number == 1
    
    spec_names = {'H'  'M'    'L'       'OH'    'LH2'    'LH1'     'LM'       'L2M'     'L3M'};
    Model      = [1     0      0        -1         2         1      0           0         0; ...  % H
                  0     1      0        0          0         0      1           1         1;  ...  % M
                  0     0      1      	0          1         1      1           2         3]; ...  % L
    log_beta   = [0     0      0       -14       n.H2     n.H1   n_Fe.K1      n_Fe.K2     n_Fe.K3];
    beta         = 10.^log_beta;
    
    description = {'Fe-nCat without hydroxide interaction'};
    num_solids = 0;
end


if model_number == 2
    
    spec_names =	{'H'     'M'      'L'     'OH'       'LH2'    'LH1'     'LM'       'L2M'      'L3M'   'MOH'    'MOH2'   'MOH4'  'M2OH2'  'M3OH4'  'MOH3'};
    Model      =	[1        0        0       -1         2          1        0         0          0        -1     -2        -4      -2         -4		-3; ...  % H
					 0        1        0       0          0          0        1         1          1         1      1         1       2         3		 1;  ...  % M
					 0        0        1       0          1          1        1         2          3         0      0         0       0         0		 0]; ...  % L
    log_beta   =    [0        0        0       -14        n.H2     n.H1    n_Fe.K1   n_Fe.K2   n_Fe.K3   Fe.MOH  Fe.MOH2    Fe.MOH4  Fe.M2OH2  Fe.M3OH4 Fe.MOH3];
    beta         = 10.^log_beta;
    
    description = {'Fe-nCat with hydroxide interaction'};
    num_solids = 1;
end


if model_number == 3
    
    spec_names = {'H'   'M'    'L'   'OH'   'LH2'    'LH1'     'LM'     'L2M'     'L3M'};
    Model      = [1     0      0       -1       2    1           0        0         0; ...  % H
				  0     1      0        0       0    0           1        1         1;  ...  % M
				  0     0      1        0       1    1           1        2         3]; ...  % L
    log_beta =   [0     0      0      -14    n.H2   n.H1    n_Al.K1    n_Al.K2     n_Al.K3];
    beta         = 10.^log_beta;
    
    description = {'Al-nCat without hydroxide interaction'};
    num_solids = 0;
end


if model_number == 4
    
    spec_names = { 'H'   'M'    'L'    'OH'  'LH2'    'LH1'      'LM'        'L2M'      'L3M'       'MOH'    'MOH2'  'MOH3'  'MOH4'    'M2OH2'   'M3OH4'   'M1OH3'};
    Model      = [  1     0      0      -1       2     1           0           0         0           -1       -2      -3      -4         -2        -4         -3; ...  % H
					0     1      0       0       0     0           1           1         1            1        1       1       1          2         3          1;  ...  % M
					0     0      1       0       1     1           1           2         3            0        0       0       0          0         0          0]; ...  % L
    log_beta   = [  0     0      0     -14    n.H2     n.H1    n_Al.K1      n_Al.K2     n_Al.K3   Al.MOH  Al.MOH2  Al.MOH3  Al.MOH4  Al.M2OH2  Al.M3OH4  Al.M1OH3];
    beta         = 10.^log_beta;
    
    description = {'Al-nCat with hydroxide interaction'};
    num_solids = 1;
end


if model_number == 5
    
    spec_names = {'H'    'M'    'L'    'OH'   'LH2'       'LH1'     'LM'       'L2M'        'L3M'  };
    Model      = [  1     0      0      -1       2          1         0           0            0; ...  % H
					0     1      0       0       0          0         1           1            1;  ...  % M
					0     0      1       0       1          1         1           2            3]; ...  % L
    log_beta   = [  0     0      0     -14     h.H2       h.H1     h_Ni.K1      h_Ni.K2     h_Ni.K3   ];
    beta         = 10.^log_beta;
    
    description = {'Ni-His without hydroxide interaction'};
    num_solids = 0;
end


if model_number == 6
    
    spec_names = {   'H'   'M'   'L'   'OH'    'LH2'         'LH1'     'LM'         'L2M'      'L3M'    'MOH'     'MOH2'   'MOH3'    'M2OH'    'M4OH4'  'M1OH2'};
    Model      =  [ 1      0      0       -1       2          1          0            0          0        -1       -2        -3       -1        -4         -2; ...  % H
					0      1      0        0       0          0          1            1          1         1        1         1        2         4          1 ;  ...  % M
					0      0      1        0       1          1          1            2          3         0        0         0        0         0          0]; ...  % L
    log_beta   =  [ 0      0      0       -14    h.H2       h.H1     h_Ni.K1      h_Ni.K2     h_Ni.K3   Ni.MOH  Ni.MOH2  Ni.MOH3   Ni.M2OH   Ni.M4OH4  Ni.M1OH2 ];
    beta         = 10.^log_beta;
    
    description = {'Ni-His with hydroxide interaction'};
    num_solids = 1;
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

% Clean up for calc_speciation
Model = Model';
log_beta = log_beta';

species_input_model.calc_spec.KSOLUTION = log_beta(1:end-num_solids,:);
species_input_model.calc_spec.ASOLUTION = Model(1:end-num_solids,:);
species_input_model.calc_spec.KSOLID = log_beta(end-num_solids+1:end,:);
species_input_model.calc_spec.ASOLID = Model(end-num_solids+1:end,:);
species_input_model.calc_spec.SOLUTIONNAMES = char(spec_names(1:end-num_solids));
species_input_model.calc_spec.SOLIDNAMES = char(spec_names(end-num_solids+1:end));

end
