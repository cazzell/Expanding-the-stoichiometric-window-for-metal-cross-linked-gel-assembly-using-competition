%% Modification of this code by Seth Allen Cazzell from original code by Donald Scott Smith.
% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Copyright 2019 Donald Scott Smith All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% This code is modified from Donald Scott Smith
% see his work
% Smith, D. Scott, "Solution of Simultaneous Chemical Equilibria in Heterogeneous Systems: Implementation in Matlab" (2019).
% Chemistry Faculty Publications. 14. https://scholars.wlu.ca/chem_faculty/14

% see reference 25 in "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% cazzell.lbi@gmail.com

%% Carrayrou et al AIChE journal 2002, 48, 894-904. % implement their method using their notation
% try HFO ppte as an example calc and at fixed pH
%function II=chem_equilib_Fe_fixedpH(pH,T)

%% This function takes the input stoichiometric/thermodynamic model and outputs the species concentrations as a function of A and B concentration.

% Remember in MM that the higher functional molecule is A and the linear
% molecule is B. Here we flip that so all letter references in MM need to
% be reversed

function [speciation] = calc_speciation(K_AB,K_BC,B_ratios,C_ratios)

% molar concentrations (mol/L)
n_A = 200 ./ 5000; %5kDa, 2-functional 100 mg/mL aka typical experimental hydrogel
n_B = n_A.*B_ratios;
n_C = n_A.*C_ratios;
master_index = 1;

for current_n_C = n_C
	
	current_n_C
	
	% Build thermodynamic models when competitor is present
	if K_BC > 0.5
		SpeciesName =	{'A'  'B'    'C'  'AB',   'BC'  };
		Model =			 [1     0      0    1        0    ;  %A
						  0     1      0    1        1    ;  %B
						  0     0      1    0        1    ]; %C
		log_beta =		 [0     0      0   K_AB      K_BC   ];
		beta  = 10.^log_beta;
		
		B_index = 2;
		AB_index = 4;
		BC_index = 5;
		
		num_solids = 0;
		
	end
	
	%Build thermodynamic model when ignoring competitor
	if K_BC < 0.5
		SpeciesName =	{'A'  'B' 'C'      'AB'  };
		Model =         [1     0	0		1     ;  %A
						 0     1	0		1     ;  %B
						 0     0	1		0    ];  %C
		log_beta =	    [0     0	1		K_AB  ];
		beta  = 10.^log_beta;
		
		B_index = 2;
		AB_index = 4;
		
		num_solids = 0;
		
	end
	
	Model = Model';
	log_beta = log_beta';
	
	c_species = zeros(length(B_ratios), length(beta));
	
	numpts=size(n_C,2);
	%Ncp is number of condensed phases
	Ncp=0;
	% maximum number of times to run newton raphson
	iterations=100;
	% Acceptance criteria
	criteria=1e-16;
	
	Ksolution = log_beta(1:end-num_solids,:);
	Asolution = Model(1:end-num_solids,:);
	Ksolid = log_beta(end-num_solids+1:end,:);
	Asolid = Model(end-num_solids+1:end,:);
	
	% Loop at each nB value
	for i=1:size(B_ratios,2)
		
		T = [n_A, n_B(i), current_n_C];
		
		if master_index == 1
			guess = [0.04, 1e-16, 1e-16];
		else
			guess = speciation.n_C(master_index-1).species(1,1:3);
		end
		
		Asolid_SI_check=Asolid;
		Ksolid_SI_check=Ksolid;
		
		% number of remaining components
		Nx=size(Asolution,2);
		% number of condensed phase species Ncp
		Ncp=size(Asolid,1);
		% number of solution species Ncp
		Nc=size(Asolution,1);
		
		% calculate species using NR
		solids=zeros(1,Ncp);
		
		% Runs NR, but seeds guess from previous results when moving on to new pH values
		if i==1
			[species,err,SI]=NR_method_solution(Asolution,Asolid,Ksolid,Ksolution,T',[guess(1:Nx)]',iterations,criteria);
			[species,err,SI]=NR_method_solution(Asolution,Asolid,Ksolid,Ksolution,T', species, iterations,criteria);
		end
		if i>1
			[species,err,SI]=NR_method_solution(Asolution,Asolid,Ksolid,Ksolution,T', species, iterations,criteria);
		end
		
		%% Runs calculation for condensed phase, makes corrections if saturation occurs
		for qq=1:Ncp
			
			[Y,I]=max(SI);
			
			if Y>1.000000001 Iindex(qq)=I;
				Asolidtemp(qq,:)=Asolid_SI_check(I,:);
				Ksolidtemp(qq,:)=Ksolid_SI_check(I,:);
				solidguess(qq)=T(I)*0.5;
				if i>1
					txt=['solidguess(qq)=',SOLIDNAMES(I,:),'(i-1);']; eval(txt);
				end
				guess=[species(2:Nx+1)' solidguess];
				[species,err,SItst,solids]=NR_method(Asolution,Asolidtemp',Ksolidtemp,Ksolution,T',guess',iterations,criteria);
				for q=1:size(solids,1)
					txt=[SOLIDNAMES(Iindex(q),:),'(i)=solids(q);']; eval(txt)
				end
			end
			
			Q=Asolid*log10(species(2:Nx+1)); SI=10.^(Q+Ksolid); Ifirst=I;
			
		end
		
		Q=Asolid*log10(species(2:Nx+1)); SI=10.^(Q+Ksolid);
		SI_summary(i,:)=SI;
		
		species_summary(i,:)=species;
		mass_err_summary(i,:)=(err(1));
		
		Asolidtemp=[]; Ksolidtemp=[];
		
	end
		
	% Save raw speciation data
	speciation.n_C(master_index).species = species_summary;
	
 	B_frac_AB = species_summary(:,AB_index) ./ n_B';
 	speciation.n_C(master_index).pB = B_frac_AB;
	
	% Record molar n_A and n_B
	speciation.molar_A = n_A;
	speciation.molar_B = n_B;
	speciation.molar_C = n_C;
	
	% Record ratios n_B/n_A and n_C/n_A
	speciation.B_ratios = B_ratios';
	speciation.C_ratios = C_ratios';
	
	master_index = master_index + 1;
	
end

end
