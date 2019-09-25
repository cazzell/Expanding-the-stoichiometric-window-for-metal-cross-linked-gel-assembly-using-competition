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

%% This function takes the input stoichiometric/thermodynamic model and outputs the species concentrations as a function of pH and metal concentration.

function [speciation] = calc_speciation(species_input_model,pH,metal_concentrations,model_number)

speciation.pH = pH

% hydrogels were made at 100 mg/mL
c_0_L = 2.*(100./5000);
speciation.molar_L = c_0_L;

for concentration_number=1:length(metal_concentrations)
	
	% Metal Ratio
	'calc species'
	metal_ratio = metal_concentrations(concentration_number)
	c_0_M = c_0_L .* metal_ratio;
	
	% Total component concentration = [metal_conc,L_conc] ;
	T = [c_0_M, c_0_L];
	
	% Extract Model Data from species_input_model
	KSOLUTION = species_input_model.calc_spec.KSOLUTION;
	KSOLID = species_input_model.calc_spec.KSOLID;
	ASOLUTION = species_input_model.calc_spec.ASOLUTION;
	ASOLID = species_input_model.calc_spec.ASOLID;
	SOLUTIONNAMES = species_input_model.calc_spec.SOLUTIONNAMES;
	SOLIDNAMES = species_input_model.calc_spec.SOLIDNAMES;
			
	% iteration 1 guess concentration of free components guess
	if concentration_number > 1
		guess = speciation.metal_concentration(concentration_number-1).species(1,2:3);
	else
		guess = [1e-10, 1e-10];
	end
		
	%numpt = number of points
	numpts=size(pH,2);
	%Ncp is number of condensed phases
	Ncp=size(ASOLID,1);
	% maximum number of times to run newton raphson
	iterations=100;
	% Acceptance criteria
	criteria=1e-16;
	
	for i=1:size(SOLIDNAMES,1)
		txt=[SOLIDNAMES(i,:),'=zeros(numpts,1);']; eval(txt)
	end
	
	% Loop at each pH value
	for i=1:size(pH,2)
		
		% Recast equilibrium equations with fixed pH
		[Ksolution,Ksolid,Asolution,Asolid]=get_equilib_fixed_pH(KSOLUTION,KSOLID,ASOLUTION,ASOLID,pH(i));
		
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
		end
		if i>1
			[species,err,SI]=NR_method_solution(Asolution,Asolid,Ksolid,Ksolution,T',[species(2:Nx+1)],iterations,criteria);
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
		
		combined_species = [species',solids];
		combined_species_L = (combined_species .* species_input_model.Model(3,:)) ./ c_0_L;
		combined_species_M = (combined_species .* species_input_model.Model(2,:)) ./ c_0_M;
		
		%% Save data
		speciation.metal_concentration(concentration_number).species(i,:) = combined_species;
		speciation.metal_concentration(concentration_number).fraction.ligand(i,:) = combined_species_L;
		speciation.metal_concentration(concentration_number).fraction.metal(i,:) = combined_species_M;
		
	
				
	end
	
	speciation.metal_concentration(concentration_number).metal_concentration = metal_ratio;
	
	for i=1:size(species_summary,2)
		txt=[SOLUTIONNAMES(i,:),'=species_summary(:,i);']; eval(txt)
	end
	
end

end
