# Base de dados da tese

Instâncias sintéticas usadas nos experimentos da tese. Cada instância é um par de genomas (origem e alvo) com um cromossomo linear, genes com sinal e regiões intergênicas. As regiões da origem têm tamanho fixo, e as do alvo podem ser flexíveis (intervalos).

## Organização

```
db/N<n>/<mod>/N<n>_<mod><op>_I<indels>_L<alfabeto>_F<flex>.txt
```

| parâmetro | significado | valores |
|---|---|---|
| `n` | número de genes da origem (`n + 1` regiões intergênicas) | 200, 500 |
| `mod` | tipo de operação conservativa usada para gerar o alvo | `SDCJ` (DCJ com sinal), `SREV` (reversões com sinal), `SREV_TRANS` (reversões com sinal e transposições), `TRANS` (transposições) |
| `op` | número de operações conservativas aplicadas | 10, 50, 100, 150, 200 |
| `indels` | número de indels aplicados: `indels/2` inserções e `indels/2` deleções de um único gene | 0, 10, 50, 70, 100 |
| `alfabeto` | tamanho do alfabeto usado para gerar a origem; alfabetos menores geram mais cópias de cada gene | 10, 50, 100, 150, 200 |
| `flex` | porcentagem usada para gerar as regiões flexíveis do alvo: cada região vira um intervalo de `flex`% abaixo a `flex`% acima do seu tamanho | 0, 20, 40, 60, 100 |

Há um arquivo para cada combinação de valores (625 por tamanho e operação, 5000 no total). Cada arquivo tem 100 instâncias.

## Formato

Cada instância ocupa 4 linhas:

1. genes da origem;
2. regiões intergênicas da origem (tamanhos);
3. genes do alvo;
4. regiões intergênicas do alvo, no formato `mínimo:máximo`.

As linhas de genes começam com `L` (cromossomo linear) e terminam com `;`. Os genomas não têm caps, então um cromossomo com `m` genes tem `m + 1` regiões intergênicas. Exemplo de início de instância:

```
L 1 6 -1 -7 8 -3 ... ;
18 92 22 52 96 57 ...
L 6 -1 8 10 4 -9 ... ;
32:74 14:30 19:43 37:85 ...
```

## Subconjunto usado nos experimentos de empacotamento e distância

Os scripts `scripts/run_dec_tese.sh` e `scripts/run_dist_tese.sh` usam, por padrão, só os arquivos com `n = 200`, `mod` em `SDCJ`, `SREV` e `SREV_TRANS`, `op` em 100 e 200, `indels` em 0 e 50, `alfabeto` em 10, 100 e 200 e `flex` em 0, 40 e 100 (36 arquivos por operação, 108 no total). Os scripts de partição usam a base inteira.
