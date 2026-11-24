# Modulação Binária Diferencial por Ruído — Resultados em MATLAB

Scripts MATLAB para reproduzir as Figuras 5.1 a 5.7 e as tabelas de dimensionamento do Projeto Final de Curso. O trabalho avalia a probabilidade de erro de bit (BEP) da DBN e compara seu desempenho com o da OODN e da NoiseMod.

## Scripts

| Arquivo | Resultado |
|---|---|
| `figura_5_1.m` | Erro relativo de duas aproximações da BEP da DBN em relação à expressão exata em canal Rayleigh. |
| `figura_5_2.m` | BEP em função da SNR em AWGN, para 20 e 100 amostras por símbolo. |
| `figura_5_3.m` | BEP em função do número de amostras por símbolo em AWGN, para SNRs de −10 e −5 dB. |
| `figura_5_4.m` | BEP em função da SNR em canais κ–μ, nos cenários LoS e nLoS, para 20 e 100 amostras por símbolo. |
| `figura_5_5.m` | BEP em função do número de amostras por símbolo nos mesmos cenários LoS e nLoS, para SNRs de −10 e −5 dB. |
| `figura_5_6.m` | Efeito do erro na estimativa da potência de ruído em canal Rayleigh, com SNR de 0 dB. |
| `figura_5_7.m` | Efeito do desvio residual de frequência da portadora (CFO) em nLoS, com SNR de −3 dB. |
| `tabela_BEP.m` | Número mínimo de amostras por símbolo e taxa normalizada R/fₛ para a meta de BEP de 10⁻³ em AWGN, LoS e nLoS. |

## Execução

Abra o script desejado no MATLAB e clique em **Run**. Cada arquivo contém os parâmetros e as funções necessários para gerar o respectivo resultado.

As Figuras 5.2 a 5.7 combinam curvas analíticas com simulações de Monte Carlo. Os scripts usam `N_bits = 1e8` decisões por configuração, processadas em blocos de `2e5`. Para um teste rápido, altere `N_bits` para `1e5`. Com menos decisões, aumenta a incerteza na estimativa de probabilidades de erro pequenas. A Figura 5.1 e as tabelas usam cálculo analítico.

## Arquivos gerados

Cada script salva os resultados na própria pasta.

- **Figuras 5.1 a 5.7** — cada figura é exportada em `.fig` editável, PDF vetorial e PNG a 300 dpi.
- **Figura 5.1** — `figura_5_1.dat` contém a SNR, o número de amostras por símbolo e os dois erros relativos percentuais.
- **Figuras 5.2 a 5.7** — cada canal gera um arquivo `_teorica.dat` com as curvas analíticas e um `_MC.dat` com as estimativas de Monte Carlo. Por exemplo, a Figura 5.2 gera `figura_5_2_AWGN_teorica.dat` e `figura_5_2_AWGN_MC.dat`. Quando uma simulação termina com zero erros observados, o arquivo registra zero e o gráfico omite o ponto no eixo logarítmico.
- **Tabelas** — `tabela_BEP_AWGN.csv`, `tabela_BEP_LoS.csv` e `tabela_BEP_nLoS.csv`. A indicação `>10000` significa que a modulação não atingiu a meta de BEP até o limite de 10 000 amostras por símbolo.

As Figuras 5.1 a 5.7 do PFC foram geradas em LaTeX com TikZ a partir dos arquivos `.dat` exportados pelo MATLAB. Para usar os dados com o TikZ original do TCC, ajuste os nomes dos arquivos e das colunas no LaTeX. Os scripts também geram PDFs que podem ser incluídos diretamente no documento.

## Referência da DBN

A. A. dos Anjos, J. V. F. Borges, H. S. Silva e R. D. Vieira, “Differential Binary Noise Modulation”, *IEEE Wireless Communications Letters*, vol. 15, pp. 3746–3750, 2026. [DOI 10.1109/LWC.2026.3704137](https://doi.org/10.1109/LWC.2026.3704137)
