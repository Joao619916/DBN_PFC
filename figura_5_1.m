clear; clc; close all;
% 1. Mapa de erro relativo da DBN em Rayleigh
N_valores = 20:20:400;
SNR_dB = -10:0.5:10;
[H,W] = canal(0,1);
erro_integral = zeros(numel(N_valores),numel(SNR_dB));
erro_simplificada = erro_integral;

% 2. Comparar cada aproximacao com a BEP exata no mesmo ponto
for j = 1:numel(N_valores)
    N = N_valores(j);
    [E,peso] = energyNodes(N,160);
    for i = 1:numel(SNR_dB)
        s = 10^(SNR_dB(i)/10)*H;
        exata = W'*betainc(1./(2*(1+s)),N,N);
        sinal = E*s';
        integral_DBN = peso'*(0.5*erfc(sinal./sqrt(2*sinal+N)))*W;
        simplificada = W'*(0.5*erfc(sqrt(N)*s./sqrt(2*s+1)));
        erro_integral(j,i) = 100*abs(integral_DBN/exata-1);
        erro_simplificada(j,i) = 100*abs(simplificada/exata-1);
    end
end

% 3. Mostrar os dois mapas com a mesma escala de cores
figure('Color','w');
tiledlayout(1,2,'TileSpacing','compact');
mapas = {erro_integral,erro_simplificada};
titulos = {'Aproximação gaussiana','Aproximação por energia média'};
for k = 1:2
    nexttile;
    surf(SNR_dB,N_valores,mapas{k},'EdgeColor','none','FaceColor','interp');
    view(2); axis tight; caxis([0 2.6]);
    xlabel('SNR média por amostra (dB)'); ylabel('Amostras por símbolo, N');
    title(titulos{k});
end
cores = [255 255 240;255 237 160;254 178 76;240 59 32;128 0 38]/255;
colormap(interp1(0:4,cores,linspace(0,4,256)));
barra = colorbar; barra.Layout.Tile = 'east'; barra.Label.String = 'Erro relativo (%)';

% 4. Salvar os valores que deram origem ao mapa
[snr,n] = meshgrid(SNR_dB,N_valores);
dados = table(snr(:),n(:),erro_integral(:),erro_simplificada(:), ...
    'VariableNames',{'SNR_dB','N','Erro_integral_pct','Erro_simplificada_pct'});
pasta = fileparts(mfilename('fullpath'));
writetable(dados,fullfile(pasta,'figura_5_1.dat'),'FileType','text','Delimiter',' ');
savefig(gcf,fullfile(pasta,'figura_5_1.fig'));
exportgraphics(gcf,fullfile(pasta,'figura_5_1.pdf'),'ContentType','vector');
exportgraphics(gcf,fullfile(pasta,'figura_5_1.png'),'Resolution',300);

function [H,W] = canal(kappa,mu)
% Quadratura de Gauss-Legendre em log(H), para E[H] = 1.
ordem = 640;
j = (1:ordem-1)';
b = j./sqrt(4*j.^2-1);
[V,D] = eig(diag(b,1)+diag(b,-1));
[t,idx] = sort(diag(D));
logH = 23.5*t-18.5;
H = exp(logH);
densidade = zeros(size(H));
if kappa == 0 && mu == 1
    densidade = H.*exp(-H);
else
    taxa = mu*(1+kappa);
    lambda = mu*kappa;
    for k = 0:40
        a = mu+k;
        densidade = densidade + exp(-lambda+k*log(lambda)-gammaln(k+1) ...
            + a*log(taxa)+a*logH-taxa*H-gammaln(a));
    end
end
W = 47*V(1,idx)'.^2.*densidade;
end

function [x,w] = energyNodes(N,ordem)
% Media sobre uma energia Gama(N,1) por quadratura de Gauss-Laguerre.
% Guardar os pesos evita repetir a decomposicao para o mesmo N.
persistent chaves nos pesos
if isempty(chaves)
    chaves = zeros(0,2);
    nos = {};
    pesos = {};
end
idx = find(all(chaves == [N ordem],2),1);
if ~isempty(idx)
    x = nos{idx}; w = pesos{idx};
    return;
end
i = (1:ordem)';
j = (1:ordem-1)';
b = sqrt(j.*(j+N-1));
[V,D] = eig(diag(2*i+N-2)+diag(b,1)+diag(b,-1));
[x,idx] = sort(diag(D));
w = V(1,idx)'.^2;
w = w/sum(w);
chaves(end+1,:) = [N ordem];
nos{end+1} = x;
pesos{end+1} = w;
end
