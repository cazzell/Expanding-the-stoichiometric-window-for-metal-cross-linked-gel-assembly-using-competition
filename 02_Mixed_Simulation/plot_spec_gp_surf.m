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

function [spec_gp_surf,speciation] = plot_spec_gp_surf(speciation,K_AB,K_BC,B_ratios,C_ratios)

r_K_A = round(K_AB);
r_K_C = round(K_BC);

K_A_str = num2str(r_K_A);
K_C_str = num2str(r_K_C);

description_str = ['K_A_B-', K_A_str, '-K_B_C-', K_C_str]; 

x = B_ratios;
y = C_ratios;

[X,Y] = meshgrid(x,y);

for N_C_number = 1:length(C_ratios)
    
    % Makes Z matrix of predicted phantom network moduli
    Z(:,N_C_number) = speciation.n_C(N_C_number).gp_phantom;
    
end

Z = Z';

Z = Z.*(Z>0);

Z(isnan(Z))=0;

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
xlabel('N_B/N_A','fontsize',fontsize,'FontName', 'Helvetica Neue','FontWeight','Bold')
ylabel('N_C/N_A','fontsize',fontsize,'FontName', 'Helvetica Neue','FontWeight','Bold')
yticks([0; 1; 2; 3; 4; 5]);
yticklabels({0; 1; 2; 3; 4; 5});
set(gca,'TickDir','both')

[cmin, cmax] = caxis;
cmax = ceil(cmax);
caxis([cmin cmax])
caxis([0 15])

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

saveas(fig_contour_print,['contour_print_', description_str])
saveas(fig_contour_print,['contour_print_',description_str],'epsc')

if max(h.LevelList) < cmax
    disp('Well Scaled')
else
    disp('Poorly Scaled')
end

[a,b] = size(Z);
num_pixels = a .* b;

score = nansum(nansum(Z,1),2);
spec_gp_surf.score = score;
spec_gp_surf.average = score./num_pixels;

end