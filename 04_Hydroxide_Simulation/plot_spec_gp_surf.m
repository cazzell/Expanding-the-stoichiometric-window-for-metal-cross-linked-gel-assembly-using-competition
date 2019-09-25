%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% This function plots data generated from mechanical prediction

function [spec_gp_surf] = plot_spec_gp_surf(metal_ratios,model_number,species_input_model,speciation,pH,stoichiometry,K_competitor)


x = metal_ratios;
y = pH;

[X,Y] = meshgrid(x,y);

for metal_number = 1:length(metal_ratios)
    
	Z(:,metal_number) = speciation.metal_concentration(metal_number).gp_phantom;
    
end

K_1_str = num2str(round(species_input_model.log_beta(5)));
K_2_str = num2str(round(species_input_model.log_beta(6)));
K_3_str = num2str(round(species_input_model.log_beta(7)));
K_C_str = num2str(round(K_competitor));
stoich_str = num2str(round(stoichiometry));

description_str = ['K1_', K_1_str, '_K_2_', K_2_str, '_K_3_', K_3_str, '_K_C_', K_C_str,'_s_', stoich_str]; 

spec_gp_surf.X = X;
spec_gp_surf.Y = Y;
spec_gp_surf.Z = Z;

% Plot Contour Map 
fontsize = 9.5;

fig_contour_print = figure;
hold all
[~,h] = contourf(X,Y,Z./1000,6);
h.LineWidth = .5;
h.LineStyle = '-';
set(gca,'FontName','Helvetica Neue','FontSize',fontsize)
set(gca,'Box','on')
xlabel('Metal Equivalency','fontsize',fontsize,'FontName', 'Helvetica Neue','FontWeight','Bold')
ylabel('pH','fontsize',fontsize,'FontName', 'Helvetica Neue','FontWeight','Bold')
xticks([0; 1/3; 2/3; 3/3; 4/3; 5/3]);
xticklabels({'0/3'; '1/3'; '2/3'; '3/3'; '4/3'; '5/3'});
yticks([0; 2; 4; 6; 8; 10; 12; 14]);
yticklabels({0; 2; 4; 6; 8; 10; 12; 14});
hA = gca;
hA.YAxis.MinorTickValues = [0:1:14];

cmax = 15;
cmin = 0;
caxis([cmin cmax])

h.LevelList
c = colorbar;
c.FontSize = fontsize;
c.FontName = 'Helvetica Neue';
c.Label.String = 'Plateau Modulus (kPa)';
c.Label.FontSize = fontsize;
c.Label.FontName = 'Helvetica Neue';
c.Label.FontWeight = 'Bold';
c.Ticks = [0  5  10  15];
c.TickDirection = 'both';
width = 2.75;
height = 2;
set(fig_contour_print, 'Position', [400,400,width   *80,height    *76.53])

axis square

set(gcf,'renderer','Painters')

saveas(fig_contour_print,[description_str, 'contour_print'])
saveas(fig_contour_print,[description_str, 'contour_print'],'epsc')

if max(h.LevelList) < cmax
	disp('Well Scaled')
else
	disp('Poorly Scaled')
end

score = sum(sum(Z,1),2);
[pixels_row, pixels_column] = size(Z);
num_pixels = pixels_row .* pixels_column;
average_modulus = score ./ num_pixels;
spec_gp_surf.score = score;
spec_gp_surf.average_modulus = average_modulus;

end


