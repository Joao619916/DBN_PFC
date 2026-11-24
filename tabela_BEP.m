clear; clc; close all;
% 1. Dimensionamento para a meta de BEP do TCC
meta = 1e-3;
N_max = 10000;
L = 100;                          % Intervalos por quadro DBN, incluindo a referencia
alpha = 10;
SNR_dB = [-10 -5 0 5 10];
canais = {'AWGN','LoS','nLoS'};
kappa = [0 2.66 1.05]; mu = [1 0.83 0.84];
nomes = {'NoiseMod','OODN','DBN'};
pasta = fileparts(mfilename('fullpath'));

for c = 1:numel(canais)
    if c == 1
        H = 1; W = 1;
    else
        [H,W] = canal(kappa(c),mu(c));
    end
    tabela = table(SNR_dB(:),'VariableNames',{'SNR_dB'});

    % 2. Buscar o menor N inteiro que atinge a meta em cada modulacao
    for k = 1:3
        N_min = strings(numel(SNR_dB),1);
        taxa = NaN(numel(SNR_dB),1);
        for i = 1:numel(SNR_dB)
            s = 10^(SNR_dB(i)/10)*H;
            if bep(N_max,s,W,k,alpha) > meta
                N_min(i) = sprintf('>%d',N_max);
                continue;
            end
            baixo = 1; alto = N_max;
            while baixo < alto
                meio = floor((baixo+alto)/2);
                if bep(meio,s,W,k,alpha) <= meta
                    alto = meio;
                else
                    baixo = meio+1;
                end
            end
            N_min(i) = string(baixo);
            taxa(i) = 1/baixo;
            if k == 3
                taxa(i) = (L-1)/(L*baixo);
            end
        end
        tabela.([nomes{k} '_N_min']) = N_min;
        tabela.([nomes{k} '_R_fs']) = taxa;
    end

    % 3. Mostrar e salvar a tabela do canal
    disp(canais{c}); disp(tabela);
    writetable(tabela,fullfile(pasta,['tabela_BEP_' canais{c} '.csv']));
end

function p = bep(N,s,W,k,alpha)
if k == 3
    [E,peso] = energyNodes(N,160);
    sinal = E*s';
    p = peso'*(0.5*erfc(sinal./sqrt(2*sinal+N)))*W;
else
    pot0 = [2/(1+alpha) 0]; pot1 = [2*alpha/(1+alpha) 2];
    lambda = limiar(s,pot0(k),pot1(k),1);
    p0 = gammainc(N*lambda./(1+pot0(k)*s),N,'upper');
    p1 = gammainc(N*lambda./(1+pot1(k)*s),N,'lower');
    p = W'*(p0+p1)/2;
end
end

function lambda = limiar(s,p0,p1,q)
% Limiar da energia media, com potencia de ruido estimada q.
v0 = q+p0*s;
v1 = q+p1*s;
d = (p1-p0)*s./v0;
r = ones(size(d));
idx = d>0;
r(idx) = log1p(d(idx))./d(idx);
lambda = v1.*r;
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
