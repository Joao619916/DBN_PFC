clear; clc; close all;
rng(20261114,'twister');

% 1. Parametros da Figura 5.3
N = 50:5:400;
pontos_MC = 50:50:400;
SNR_dB = [-10 -5];
N_bits = 1e8;
bloco = 2e5;
alpha = 10;
nome_canal = 'AWGN';
pasta = fileparts(mfilename('fullpath'));
cores = [228 117 0;32 147 46;19 112 185;174 37 65;114 48 129;204 77 25]/255;
nomes = {'NoiseMod','OODN','DBN'};
marcas = {'^','s','o'};
estilos = {'-','--'};
figure('Color','w');
H = 1; W = 1;
hold on; grid on; box on; set(gca,'YScale','log');
dados = table(N(:),'VariableNames',{'N'});
dados_MC = table(pontos_MC(:),'VariableNames',{'N'});
for j = 1:numel(SNR_dB)
    % 2. Calcular as curvas teoricas
    bep = zeros(numel(N),3);
    for i = 1:numel(N)
        N_atual = N(i);
        bep(i,:) = teoria(N_atual,10^(SNR_dB(j)/10),H,W,alpha,1);
    end
    % 3. Simular em blocos e contar os erros de cada detector
    erros = zeros(numel(pontos_MC),3);
    for i = 1:numel(pontos_MC)
        N_atual = pontos_MC(i);
        for inicio = 1:bloco:N_bits
            m = min(bloco,N_bits-inicio+1);
            bits = randi([0 1],m,1);
            ganho = 1;
            G = randg(N_atual,m,1);
            A = randg(N_atual,m,1);
            C = randg(N_atual,m,1);
            s = 10^(SNR_dB(j)/10)*ganho;
            q = 1;
            % Variancias de sinal para bit 0 e bit 1, em cada modulacao.
            pot0 = [2/(1+alpha) 0]; pot1 = [2*alpha/(1+alpha) 2];
            for k = 1:2
                v = 1+s.*(pot0(k)+(pot1(k)-pot0(k))*bits);
                T = v.*G/N_atual;
                erros(i,k) = erros(i,k)+nnz((T >= limiar(s,pot0(k),pot1(k),q)) ~= bits);
            end
            erro_DBN = (1+2*s).*A < C;
            erros(i,3) = erros(i,3)+nnz(erro_DBN);
        end
    end
    % 4. Guardar os dados e tracar as curvas
    mc = erros/N_bits;
    legenda = sprintf('%g dB',SNR_dB(j));
    sufixo = sprintf('m%d',-SNR_dB(j));
    for k = 1:3
        nome = [nomes{k} '_' sufixo];
        dados.(nome) = bep(:,k);
        dados_MC.([nome '_MC']) = mc(:,k);
        y = mc(:,k); y(y==0) = NaN;
        cor = cores(3*(j-1)+k,:);
        plot(N,bep(:,k),'k','LineStyle',estilos{j}, ...
            'Marker',marcas{k},'MarkerIndices',1:10:numel(N), ...
            'MarkerFaceColor','w','HandleVisibility','off');
        plot(pontos_MC,y,'Color',cor,'LineStyle',estilos{j}, ...
            'Marker',marcas{k},'MarkerFaceColor',cor, ...
            'DisplayName',[nomes{k} ', ' legenda]);
    end
    g = 10^(SNR_dB(j)/10);
    aprox = 0.5*erfc(sqrt(N(:)).*g./sqrt(2*g+1));
    dados.(['DBN_aprox_' sufixo]) = aprox;
    curva = plot(N,aprox,'r:','DisplayName','DBN (aprox.)');
    if j==1, curva.HandleVisibility = 'off'; end
end
xlabel('Amostras por símbolo, N'); ylabel('BEP'); title(nome_canal);
xlim([N(1) N(end)]); ylim([1e-6 0.5]);
legend('Location','best','FontSize',8);
prefixo = fullfile(pasta,['figura_5_3_' nome_canal]);
writetable(dados,[prefixo '_teorica.dat'],'FileType','text','Delimiter',' ');
writetable(dados_MC,[prefixo '_MC.dat'],'FileType','text','Delimiter',' ');
savefig(gcf,fullfile(pasta,'figura_5_3.fig'));
exportgraphics(gcf,fullfile(pasta,'figura_5_3.pdf'),'ContentType','vector');
exportgraphics(gcf,fullfile(pasta,'figura_5_3.png'),'Resolution',300);

function bep = teoria(N,g,H,W,alpha,q)
% NoiseMod e OODN: BEP condicional integrada sobre o canal.
s = g*H;
pot0 = [2/(1+alpha) 0];
pot1 = [2*alpha/(1+alpha) 2];
bep = zeros(1,3);
for k = 1:2
    lambda = limiar(s,pot0(k),pot1(k),q);
    p0 = gammainc(N*lambda./(1+pot0(k)*s),N,'upper');
    p1 = gammainc(N*lambda./(1+pot1(k)*s),N,'lower');
    bep(k) = W'*(p0+p1)/2;
end
% DBN: media da aproximacao gaussiana sobre a energia e o canal.
[E,peso] = energyNodes(N,160);
sinal = E*s';
bep(3) = peso'*(0.5*erfc(sinal./sqrt(2*sinal+N)))*W;
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
