%% Modification of this code by Seth Allen Cazzell from original code by Scott Grindy.
% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.
 
%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
%% Copyright 2017 Scott Grindy All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Function to convert speciation data to plateau moduli predictions
% Based off of the methods developed by Miller and Macosko
% D. R. Miller, C. W. Macosko, A new derivation of post gel properties of network polymers. Macromolecules 9, 206-211 (1976).
% and Scott Grindy
% S. C. Grindy, "Complex mechanical design of bio-inspired model transient network hydrogels," PhD thesis, Massachusetts Institute of Technology (2017).

% Remember in MM that the higher functional molecule is A and the linear
% molecule is B. Here we flip that so all letter references in MM need to
% be reversed

function [speciation] = calc_Gp_linear_Miller_Macosko(speciation,C_ratios)

for C_list_number = 1:length(C_ratios)
    
    C_list_number
    
    'MM'
	
	p = speciation.n_C(C_list_number).pB;
	r = speciation.B_ratios;
	
	rpsquare = r .* p .^ 2;
	
	% MM eqn 14
	P_B_out = (1-rpsquare) ./ rpsquare;
	
	molar_ligand = speciation.molar_A;
	molar_B = speciation.molar_B;

    %P_B_in = P_B_out.^(3-1);
    %p_B_3 = (1/3) .* fraction_tris .* (1 - P_A_in).^3;
	p_B_3 = (1 - P_B_out).^3;
    
    % molar concentration of 3-functional crosslinks
    %molar_X3 = p_B_3 .* molar_ligand;
	molar_X3 = p_B_3 .* molar_B'  .* (1/3);
    % convert to per m^3 from per L
    c_X3 = molar_X3 .* 1000;
    % c represents concentration in mol/m^3
    c_crosslinks = c_X3;
    % concentration of elastically active chains
    c_chains = 1.5*c_X3;
    
    R = 8.314;
    % Absolute Temperature
    T = 25+273.15;
    RT = R * T;
    gp_phantom = (c_chains - c_crosslinks) .* RT;
    
    % Save Data
    speciation.n_C(C_list_number).gp_phantom = gp_phantom;
	speciation.n_C(C_list_number).P_A = gp_phantom;
    
end
end
