# Modulação Binária Diferencial por Ruído — Figura 5.2

Código MATLAB para a comparação da probabilidade de erro de bit (BEP) das modulações DBN, OODN e NoiseMod em canal AWGN, referente à Figura 5.2 do Projeto Final de Curso.

## Execução

Abra `figura_5_2.m` no MATLAB e execute-o. Não são necessários arquivos externos. Recomenda-se MATLAB R2020a ou posterior.

- Amostras por símbolo: N = 20 e 100.
- SNR média por amostra: −10 a 10 dB.
- Passo analítico: 0,25 dB; passo de Monte Carlo: 1 dB.
- Decisões de Monte Carlo por configuração: 10^8.
- Razão entre variâncias dos sinais transmitidos na NoiseMod: 10.
- Variância do ruído complexo: 1.

Para testar a execução, utilize temporariamente `N_bits = 1e5`. Para a simulação completa, restaure `N_bits = 1e8`.

## Metodologia

O código calcula a BEP exata dos detectores de energia da NoiseMod e da OODN e as duas aproximações da DBN. As simulações de Monte Carlo utilizam variáveis aleatórias Gama para representar as estatísticas suficientes dos detectores, com processamento em blocos. As realizações aleatórias são compartilhadas entre os pontos de SNR. Resultados sem erros observados são exportados como NaN para representação em escala logarítmica.

## Arquivos gerados

- `fig1_Teorica.dat`: resultados analíticos (81 linhas de dados e 9 colunas).
- `fig1_Pontos_MC.dat`: resultados simulados (21 linhas de dados e 7 colunas).
- `figura_5_2.fig`: figura editável no MATLAB.
- `figura_5_2.pdf`: figura vetorial.
- `figura_5_2.png`: figura com resolução de 300 dpi.

Os arquivos são salvos na pasta do script. Os dados originais e o código TikZ não estão incluídos. Novas realizações de Monte Carlo podem diferir dos resultados originais. A compatibilidade dos nomes das colunas com o TikZ original e a execução integral desta versão ainda precisam ser verificadas.

## Referência

A. A. dos Anjos, J. V. F. Borges, H. S. Silva e R. D. Vieira, “Differential Binary Noise Modulation,” IEEE Wireless Communications Letters, vol. 15, 2026.
