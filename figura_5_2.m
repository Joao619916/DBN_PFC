clear; clc; close all;
rng(20261114,'twister');

% 1. Parametros da simulacao
N_valores = [20 100];
SNR_dB = -10:0.25:10;
SNR_linear = 10.^(SNR_dB/10);
SNR_MC = -10:1:10;
N_bits = 1e8;
bloco = 2e5;  % Numero de bits por bloco para economizar memoria
alpha = 10;

% Potencia do ruido do canal igual a 1
% Variancias recebidas para os bits 0 e 1 da NoiseMod
v0_NoiseMod = 1 + 2*SNR_linear/(1+alpha);
v1_NoiseMod = 1 + 2*alpha*SNR_linear/(1+alpha);

% Variancias recebidas para os bits 0 e 1 da OODN
v0_OODN = ones(size(SNR_linear));
v1_OODN = 1 + 2*SNR_linear;

% 2. Calcular os limiares otimos dos detectores de energia
lambda_NoiseMod = v0_NoiseMod.*v1_NoiseMod.* ...
    log(v1_NoiseMod./v0_NoiseMod)./(v1_NoiseMod-v0_NoiseMod);
lambda_OODN = v0_OODN.*v1_OODN.* ...
    log(v1_OODN./v0_OODN)./(v1_OODN-v0_OODN);

% Pontos da curva teorica que correspondem aos pontos de Monte Carlo
pontos = 1:4:length(SNR_dB);

% Cores e marcadores das curvas: NoiseMod, OODN e DBN
cores = [228 117 0; 32 147 46; 19 112 185; ...
         174 37 65; 114 48 129; 204 77 25]/255;
nomes = {'NoiseMod','OODN','DBN'};
marcadores = {'^','s','o'};
figure; hold on; grid on; box on;
set(gca,'YScale','log');
barra = waitbar(0,'Simulando...');

% Criar as tabelas para salvar os resultados
dados_teoricos = table(SNR_dB(:),'VariableNames', {'SNR_dB'});
dados_MC = table(SNR_MC(:),'VariableNames', {'SNR_dB'});

% Repetir os calculos para N = 20 e N = 100
for j = 1:length(N_valores)
    N = N_valores(j);

    % 3. Calcular a BEP teorica da NoiseMod e da OODN
    BEP_NoiseMod = 0.5*(gammainc(N*lambda_NoiseMod./v0_NoiseMod,N,'upper') + ...
                       gammainc(N*lambda_NoiseMod./v1_NoiseMod,N,'lower'));
    BEP_OODN = 0.5*(gammainc(N*lambda_OODN./v0_OODN,N,'upper') + ...
                    gammainc(N*lambda_OODN./v1_OODN,N,'lower'));

    % 4. Calcular a BEP aproximada da DBN pela integral (Eq. 3.10)
    % x = Eu/(N*sigma_u^2) possui distribuicao Gama(N,1/N)
    PDF = @(x) exp(N*log(N)-gammaln(N)+(N-1)*log(x)-N*x);
    BEP_DBN = zeros(size(SNR_dB));

    for i = 1:length(SNR_dB)
        E = @(x) N*SNR_linear(i).*x;
        BEP_DBN(i) = integral(@(x) 0.5*erfc(E(x)./ ...
            sqrt(2*E(x)+N)).*PDF(x),0,Inf,'AbsTol',1e-14,'RelTol',1e-8);
    end

    % Aproximacao simplificada da DBN (Eq. 3.11)
    BEP_DBN_aprox = 0.5*erfc(sqrt(N)*SNR_linear./sqrt(2*SNR_linear+1));

    % 5. Inicializar os contadores de erros de Monte Carlo
    erros_NoiseMod = zeros(size(SNR_MC));
    erros_OODN = zeros(size(SNR_MC));
    erros_DBN = zeros(size(SNR_MC));

    for inicio = 1:bloco:N_bits
        n_bits = min(bloco,N_bits-inicio+1);

        % Gerar bits e energias aleatorias Gama(N,1)
        % As mesmas realizacoes sao usadas nos diferentes pontos de SNR
        bits_tx = randi([0 1],n_bits,1);
        G = randg(N,n_bits,1);
        A = randg(N,n_bits,1);
        C = randg(N,n_bits,1);

        for i = 1:length(SNR_MC)
            p = pontos(i);
            gamma = SNR_linear(p);

            % 6. NoiseMod: gerar a energia recebida e detectar os bits
            v = v0_NoiseMod(p)*(1-bits_tx) + v1_NoiseMod(p)*bits_tx;
            T = v.*G/N;
            bits_rx = T >= lambda_NoiseMod(p);
            erros_NoiseMod(i) = erros_NoiseMod(i) + sum(bits_rx ~= bits_tx);

            % 7. OODN: gerar a energia recebida e detectar os bits
            v = v0_OODN(p)*(1-bits_tx) + v1_OODN(p)*bits_tx;
            T = v.*G/N;
            bits_rx = T >= lambda_OODN(p);
            erros_OODN(i) = erros_OODN(i) + sum(bits_rx ~= bits_tx);

            % 8. DBN: evento de erro exato do detector diferencial
            % A e C representam energias independentes de soma e diferenca
            % dos simbolos recebidos, normalizadas pela variancia.
            erros_DBN(i) = erros_DBN(i) + sum((1+2*gamma)*A < C);
        end

        % Atualizar a barra de progresso
        waitbar(((j-1)+(inicio+n_bits-1)/N_bits)/length(N_valores),barra);
    end

    % 9. Calcular a BEP simulada
    BEP_MC = [erros_NoiseMod; erros_OODN; erros_DBN]/N_bits;
    BEP_MC(BEP_MC == 0) = NaN;  % Omitir zeros no eixo logaritmico
    BEP_teorica = [BEP_NoiseMod; BEP_OODN; BEP_DBN];

    % 10. Guardar os resultados nas tabelas
    for k = 1:3
        coluna = sprintf('%s_N%d',nomes{k},N);
        dados_teoricos.(coluna) = BEP_teorica(k,:).';
        coluna = sprintf('%s_N%d_MC',nomes{k},N);
        dados_MC.(coluna) = BEP_MC(k,:).';
    end
    coluna = sprintf('DBN_N%d_aprox',N);
    dados_teoricos.(coluna) = BEP_DBN_aprox(:);

    % 11. Plotar as curvas teoricas e simuladas
    if N == 20
        estilo = '-';
    else
        estilo = '--';
    end

    for k = 1:3
        cor = cores(3*(j-1)+k,:);
        semilogy(SNR_dB,BEP_teorica(k,:),'k','LineStyle',estilo, ...
            'Marker',marcadores{k},'MarkerIndices',pontos, ...
            'MarkerFaceColor','w','HandleVisibility','off');
        semilogy(SNR_MC,BEP_MC(k,:),'Color',cor,'LineStyle',estilo, ...
            'Marker',marcadores{k},'MarkerFaceColor',cor, ...
            'DisplayName',sprintf('%s, N=%d',nomes{k},N));
    end

    curva = semilogy(SNR_dB,BEP_DBN_aprox,'r:', ...
        'DisplayName','DBN (aprox.)');
    if j == 1
        curva.HandleVisibility = 'off';
    end
end

close(barra);
xlabel('$\bar{\gamma}$ (dB)','Interpreter','latex');
ylabel('BEP');
xlim([-10 10]); ylim([1e-5 0.5]);
legend('Location','northeast');


% 12. Salvar os resultados e a figura
% Os arquivos ficam na pasta em que o script esta salvo
pasta = fileparts(mfilename('fullpath'));
if isempty(pasta)
    pasta = pwd;
end

writetable(dados_teoricos,fullfile(pasta,'fig1_Teorica.dat'), ...
    'FileType','text','Delimiter',' ');
writetable(dados_MC,fullfile(pasta,'fig1_Pontos_MC.dat'), ...
    'FileType','text','Delimiter',' ');

savefig(gcf,fullfile(pasta,'figura_5_2.fig'));
exportgraphics(gcf,fullfile(pasta,'figura_5_2.pdf'),'ContentType','vector');
exportgraphics(gcf,fullfile(pasta,'figura_5_2.png'),'Resolution',300);

disp('Arquivos salvos em:');
disp(pasta);
