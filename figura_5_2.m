clear; clc; close all;
rng(20261114,'twister');

% 1. Parametros da Figura 5.2
SNR_dB = -10:0.25:10;
pontos_MC = -10:1:10;
N_valores = [20 100];
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
dados = table(SNR_dB(:),'VariableNames',{'SNR_dB'});
dados_MC = table(pontos_MC(:),'VariableNames',{'SNR_dB'});
for j = 1:numel(N_valores)
    % 2. Calcular as curvas teoricas
    bep = zeros(numel(SNR_dB),3);
    N = N_valores(j);
    for i = 1:numel(SNR_dB)
        bep(i,:) = teoria(N,10^(SNR_dB(i)/10),H,W,alpha,1);
    end
    % 3. Simular em blocos e contar os erros de cada detector
    erros = zeros(numel(pontos_MC),3);
    for inicio = 1:bloco:N_bits
        m = min(bloco,N_bits-inicio+1);
        bits = randi([0 1],m,1);
        ganho = 1;
        G = randg(N,m,1);
        A = randg(N,m,1);
        C = randg(N,m,1);
        for i = 1:numel(pontos_MC)
            s = 10^(pontos_MC(i)/10)*ganho;
            q = 1;
            % Variancias de sinal para bit 0 e bit 1, em cada modulacao.
            pot0 = [2/(1+alpha) 0]; pot1 = [2*alpha/(1+alpha) 2];
            for k = 1:2
                v = 1+s.*(pot0(k)+(pot1(k)-pot0(k))*bits);
                T = v.*G/N;
                erros(i,k) = erros(i,k)+nnz((T >= limiar(s,pot0(k),pot1(k),q)) ~= bits);
            end
            erro_DBN = (1+2*s).*A < C;
            erros(i,3) = erros(i,3)+nnz(erro_DBN);
        end
    end
    % 4. Guardar os dados e tracar as curvas
    mc = erros/N_bits;
    legenda = sprintf('N=%d',N);
    sufixo = sprintf('N%d',N);
    for k = 1:3
        nome = [nomes{k} '_' sufixo];
        dados.(nome) = bep(:,k);
        dados_MC.([nome '_MC']) = mc(:,k);
        y = mc(:,k); y(y==0) = NaN;
        cor = cores(3*(j-1)+k,:);
        plot(SNR_dB,bep(:,k),'k','LineStyle',estilos{j}, ...
            'Marker',marcas{k},'MarkerIndices',1:4:numel(SNR_dB), ...
            'MarkerFaceColor','w','HandleVisibility','off');
        plot(pontos_MC,y,'Color',cor,'LineStyle',estilos{j}, ...
            'Marker',marcas{k},'MarkerFaceColor',cor, ...
            'DisplayName',[nomes{k} ', ' legenda]);
    end
    g = 10.^(SNR_dB(:)/10);
    aprox = 0.5*erfc(sqrt(N(:)).*g./sqrt(2*g+1));
    dados.(['DBN_aprox_' sufixo]) = aprox;
    curva = plot(SNR_dB,aprox,'r:','DisplayName','DBN (aprox.)');
    if j==1, curva.HandleVisibility = 'off'; end
end
xlabel('SNR média por amostra (dB)'); ylabel('BEP'); title(nome_canal);
xlim([SNR_dB(1) SNR_dB(end)]); ylim([1e-5 0.5]);
legend('Location','best','FontSize',8);
prefixo = fullfile(pasta,['figura_5_2_' nome_canal]);
writetable(dados,[prefixo '_teorica.dat'],'FileType','text','Delimiter',' ');
writetable(dados_MC,[prefixo '_MC.dat'],'FileType','text','Delimiter',' ');
savefig(gcf,fullfile(pasta,'figura_5_2.fig'));
exportgraphics(gcf,fullfile(pasta,'figura_5_2.pdf'),'ContentType','vector');
exportgraphics(gcf,fullfile(pasta,'figura_5_2.png'),'Resolution',300);

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
