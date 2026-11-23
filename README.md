# Códigos da tese

Códigos usados nos experimentos da tese de doutorado de Gabriel Siqueira sobre distâncias de rearranjo com regiões intergênicas: heurísticas de empacotamento do grafo de adjacências por ciclos, algoritmos de aproximação para distâncias de rearranjo e algoritmos de partição de genomas.

## Estrutura

| pasta | conteúdo |
|---|---|
| `cyc/` | heurísticas de empacotamento do grafo de adjacências por ciclos (`dec`), distâncias de rearranjo (`dist`), partição de genomas a partir de um empacotamento (`part`) e conversão para permutações (`toperm`); veja `cyc/README.md` |
| `partition/` | algoritmos de partição de genomas (`MCFISP-exe`): aproximação, guloso, correspondência de adjacências (SOAR/SOARComb) e algoritmo FPT; veja `partition/README.md` e `partition/doc/` |
| `db/` | base de dados sintética usada nos experimentos; veja `db/README.md` |
| `scripts/` | scripts que rodam os experimentos da tese e resumem os resultados |

## Compilação

```sh
cd cyc
make            # gera cyc/dec, cyc/dist, cyc/part e cyc/toperm
```

Requer `cmake`, `g++` com OpenMP e as dependências listadas em `cyc/shell.nix` (com Nix: `nix-shell cyc/shell.nix`).

```sh
cd partition
make            # gera partition/MCFISP-exe e partition/DB-exe
```

Requer o `stack` (Haskell), configurado para usar Nix (`partition/default.nix`). Os testes rodam com `make test`.

Os scripts compilam os executáveis sozinhos antes de rodar. Os scripts Python usam Python 3, e `table_part_tese.py` usa também o `pandas`.

## Base de dados

A pasta `db/` tem as instâncias sintéticas usadas nos experimentos: pares de genomas lineares com 200 ou 500 genes, gerados com quatro tipos de operação conservativa (`SDCJ`, `SREV`, `SREV_TRANS` e `TRANS`), variando o número de operações, o número de indels, o tamanho do alfabeto e a flexibilidade das regiões intergênicas do alvo. São 5000 arquivos com 100 instâncias cada. O formato e o significado dos parâmetros estão em `db/README.md`.

Os arquivos da base são armazenados com [Git LFS](https://git-lfs.com). Instale o `git-lfs` antes de clonar o repositório; num clone feito sem ele, os arquivos de `db/` são só ponteiros, e `git lfs pull` baixa o conteúdo.

Os scripts leem a base de `db/` por padrão (a variável `INPUT_FOLDER` troca a pasta):

- os de partição (`run_part_tese.sh` e `run_part_fix_tese.sh`) usam a base inteira (o `run_part_fix_tese.sh` só `N200` com `SDCJ`, `SREV` e `SREV_TRANS`);
- os de empacotamento e distância (`run_dec_tese.sh` e `run_dist_tese.sh`) usam um subconjunto: `N200` com `SDCJ`, `SREV` e `SREV_TRANS`, 100 e 200 operações, 0 e 50 indels, alfabetos 10, 100 e 200 e flexibilidade 0, 40 e 100 (108 arquivos). As variáveis `OP_VALS`, `I_VALS`, `L_VALS` e `F_VALS` trocam esses valores; uma variável vazia mantém todos os valores do parâmetro.

## Saídas

Todas as saídas padrão ficam em `TESE_OUT/`, na raiz do repositório, e não são versionadas:

| pasta | conteúdo |
|---|---|
| `TESE_OUT/ANS_TESE_PART/` | partições, custos, instâncias reduzidas e partições fixas |
| `TESE_OUT/ANS_TESE_DEC[_<ALGO>]/` | empacotamentos do `dec` |
| `TESE_OUT/ANS_TESE_DIST[_<ALGO>]/` | distâncias do `dist` |

Todos os scripts podem ser interrompidos e retomados: os arquivos de entrada cujos resultados já existem são pulados. As variáveis de ambiente aceitas (pastas de entrada e saída, número de iterações, número de instâncias por arquivo, valores dos parâmetros a usar) estão descritas no cabeçalho de cada script.

## Experimentos de empacotamento e distância

| script | função |
|---|---|
| `scripts/run_dec_tese.sh` | roda `dec` (heurísticas `ga` e `rand`, com a decomposição `bfs`; `map` e `1c2c3c` também estão disponíveis) |
| `scripts/run_dist_tese.sh map\|dec` | roda `dist` com mapeamentos aleatórios (`map`) ou com os empacotamentos de `run_dec_tese.sh` (`dec`) |

Exemplos:

```sh
# empacotamentos e distâncias no subconjunto padrão
scripts/run_dec_tese.sh
scripts/run_dist_tese.sh dec
# distâncias com mapeamentos aleatórios, só nos arquivos com 100 operações e sem indels
OP_VALS=100 I_VALS=0 scripts/run_dist_tese.sh map
# empacotamentos em todos os valores de flexibilidade
F_VALS= scripts/run_dec_tese.sh
```

## Partição de genomas

O `partition/` implementa os algoritmos do capítulo de partição de genomas (opção `-a` do `MCFISP-exe`):

| algoritmo | descrição |
|---|---|
| `Approx` | algoritmo de aproximação (classes ativas), com fator 2tk para relações de compatibilidade simétricas e transitivas; genomas desbalanceados passam antes por uma heurística de balanceamento |
| `Greedy`, `GreedySin` | algoritmo guloso (`GreedySin` começa pelos blocos com genes que ocorrem uma vez em cada genoma) |
| `SOAR`, `SOARComb` | algoritmo de correspondência de adjacências (`SOARComb` constrói a partição minimal a partir da correspondência) |
| `FPT` | algoritmo exato parametrizado, com limite de tempo de 10 minutos por instância |
| `ApproxLB` | escreve só o limitante inferior da prova da aproximação |

O formato de entrada e de saída e as demais opções (`-s` para genes com sinal, `--direct` para partição direta) estão em `partition/doc/Main.md`. Com `-a Approx`, o arquivo `<saída>_phase2` registra quantas vezes a segunda etapa da aproximação foi usada em cada instância.

| script | função |
|---|---|
| `scripts/run_part_tese.sh` | roda os algoritmos de partição e grava, em `TESE_OUT/ANS_TESE_PART/`, as partições (`PART/`), os custos (`COST/`), as instâncias reduzidas (`RED/`, via `scripts/reduce_partition.py`) e, para o `Approx`, os usos da segunda etapa (`PHASE2/`) |
| `scripts/run_part_fix_tese.sh` | escreve as partições em posições dos genomas originais (`FIX/`, via `scripts/partition_fixed.py`), para fixar as adjacências no `dec` (opções `-f` e `-p`) |
| `scripts/table_part_tese.py` | resume os custos de `COST/` numa tabela |

A lista de algoritmos de cada script pode ser trocada com a variável `ALGOS` (por padrão, `Greedy GreedySin SOAR SOARComb` no primeiro e `Greedy SOARComb Approx` no segundo).

## Partição combinada com as heurísticas

As partições podem ser usadas pelo `cyc` de duas formas, escolhidas com variáveis de `run_dec_tese.sh` e `run_dist_tese.sh`:

- `PART=<ALGO>`: cada bloco da partição vira um único gene, e as heurísticas rodam sobre as instâncias reduzidas (`RED/<ALGO>`);
- `FIX=<ALGO>`: as heurísticas rodam sobre as instâncias originais, com as adjacências que não são breakpoints da partição fixas (`FIX/<ALGO>`); com `FIX_PAIRS=1`, a correspondência entre os genes dos blocos também é fixa.

Por exemplo:

```sh
scripts/run_part_tese.sh
PART=SOARComb scripts/run_dec_tese.sh
PART=SOARComb scripts/run_dist_tese.sh dec
```

Os detalhes da redução, dos arquivos `FIX` e do tratamento das instâncias sem partição estão em `scripts/PARTICAO_CYC.md`.
