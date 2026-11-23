# Partições de `partition/` como entrada para o cyc

Este arquivo resume como as partições calculadas pelo `partition/` (MCFISP-exe) são passadas para as heurísticas do `cyc`, para testar a combinação de algoritmos de partição com heurísticas descrita na tese (`adj_graph.tex`, "Fixando Adjacências com uma Partição de Genomas").

## Scripts

| script | função |
|---|---|
| `run_part_tese.sh` | roda os algoritmos de partição na base da tese e gera custos e instâncias reduzidas |
| `reduce_partition.py` | transforma as saídas do MCFISP-exe de um arquivo de entrada numa instância reduzida e num arquivo de custos (chamado pelo `run_part_tese.sh`) |
| `run_part_fix_tese.sh` | grava as mesmas partições sem redução, em posições dos genomas originais (arquivos `FIX`), para o grafo de adjacências com adjacências fixas |
| `partition_fixed.py` | converte um arquivo `PART` em `FIX` e `COST` (chamado pelo `run_part_fix_tese.sh`) |
| `run_dec_tese.sh`, `run_dist_tese.sh` | rodam `dec`/`dist` do `cyc`; com `PART=<ALGO>` rodam sobre as instâncias reduzidas desse algoritmo, e com `FIX=<ALGO>` sobre as instâncias originais com as adjacências fixas da partição (variáveis abaixo) |

A base é `db/`, na raiz do repositório (ver `db/README.md`): `N<n>/<mod>/N<n>_<mod><op>_I<indels>_L<alfabeto>_F<flex>.txt`, com 100 instâncias por arquivo, 4 linhas por instância, genomas com sinal, sem caps, n genes e n + 1 regiões intergênicas.

## `run_part_tese.sh`

Para cada arquivo de entrada e algoritmo (`-a` do MCFISP-exe, sempre com `-s`), grava em `../TESE_OUT/ANS_TESE_PART/` (a pasta `TESE_OUT/` da raiz do repositório):

- `PART/<ALGO>/N<n>/<mod>/<arquivo>.txt`: as saídas do MCFISP-exe juntas. Por instância: partição da origem (genes com caps `0`, blocos separados por `|`), as regiões intergênicas (as de dentro de cada bloco intercaladas com a do breakpoint), o mesmo para o alvo, a lista de blocos compatíveis de cada bloco da origem, e o tempo (`# Time:`).
- `COST/<ALGO>/N<n>/<mod>/<arquivo>.txt`: uma linha por instância:
  ```
  breakpoints_origem breakpoints_alvo custo blocos_origem blocos_alvo exclusivos_origem exclusivos_alvo tempo
  ```
  O custo é a soma dos breakpoints das duas partições (definição da tese). Uma instância sem partição (erro ou limite de tempo) tem a linha `NA`.
- `RED/<ALGO>/N<n>/<mod>/<arquivo>.txt`: as instâncias reduzidas, com o mesmo formato e a mesma estrutura de pastas da base. Uma instância sem partição vira 4 linhas `NA`, para manter a ordem; `dec` e `dist` escrevem `NA` como resultado dela (`Dec: NA` nos arquivos `-best` do `dec`, `Dist: NA` nos `-dcj` do `dec` e nos `-all` do `dist`, e as linhas `NA` e `Dist: NA` nos `-best` do `dist`).

Algoritmos padrão: Greedy, GreedySin, SOAR e SOARComb. O Approx levou cerca de 8 minutos por instância de N200, e o FPT é exponencial; por isso ficaram de fora. O Greedy levou cerca de 1 s de CPU por instância de N200, e o SOARComb cerca de 0,05 s.

O script pode ser interrompido e retomado (arquivos já completos são pulados). Ele compila o MCFISP-exe com `stack install --local-bin-path .` em `../partition`.

## Redução (`reduce_partition.py`)

Cada bloco vira um único gene:
- os blocos da origem e do alvo são pareados por um emparelhamento máximo entre blocos compatíveis (a compatibilidade é recalculada a partir dos genes e das regiões intergênicas, porque a linha de blocos compatíveis do MCFISP-exe lista candidatos, não uma bijeção);
- cada par recebe um rótulo novo; no alvo o sinal é invertido quando os blocos são compatíveis no sentido reverso;
- um bloco sem par (bloco exclusivo, nas instâncias com indels) recebe um rótulo que só aparece no seu genoma, e o `cyc` o trata como indel;
- ficam as regiões intergênicas dos breakpoints; as de dentro de um bloco são descartadas, exceto a que liga um cap ao primeiro (último) gene do bloco que contém o cap, que vira a primeira (última) região.

Os genomas reduzidos não têm repetições, então um único mapeamento basta (`dist -k 1`). Só cromossomos lineares são aceitos (a base só tem lineares).

Ressalva: a tese garante essa redução para relações de compatibilidade transitivas (igualdade). Com a relação de intervalo (F > 0), a redução fixa um pareamento específico entre os blocos, em vez de agrupar os blocos compatíveis.

## Variáveis de ambiente

| variável | scripts | padrão |
|---|---|---|
| `INPUT_FOLDER` | todos | `db/`, na raiz do repositório |
| `OP_VALS`, `I_VALS`, `L_VALS`, `F_VALS` | `run_dec_tese.sh`, `run_dist_tese.sh` | o subconjunto da tese (`100 200`, `0 50`, `10 100 200` e `0 40 100`); vazio mantém todos os valores |
| `OUT_FOLDER` | todos | `TESE_OUT/ANS_TESE_DEC` (`run_dec_tese.sh`), `TESE_OUT/ANS_TESE_DIST` (`run_dist_tese.sh`) e `TESE_OUT/ANS_TESE_PART` (`run_part_tese.sh`, `run_part_fix_tese.sh`), na raiz do repositório |
| `DEC_FOLDER` | `run_dist_tese.sh` | `TESE_OUT/ANS_TESE_DEC` |
| `PART` | `run_dec_tese.sh`, `run_dist_tese.sh` | vazio (instâncias originais) |
| `PART_FOLDER` | `run_dec_tese.sh`, `run_dist_tese.sh` | `TESE_OUT/ANS_TESE_PART`; com `FIX`, a pasta com `FIX/` |
| `FIX`, `FIX_PAIRS` | `run_dec_tese.sh`, `run_dist_tese.sh` | vazios (sem adjacências fixas), ver "Leitura do `FIX` no `cyc`" |
| `K` | `run_dec_tese.sh` | 500 (iterações do `dec`) |
| `K_MAP` | `run_dist_tese.sh` | 1000 |
| `N_INST` | todos | 100 (instâncias por arquivo) |

Com `PART=<ALGO>` (por exemplo `SOARComb`), os dois scripts usam como entrada `PART_FOLDER/RED/<ALGO>` e gravam em `TESE_OUT/ANS_TESE_DEC_<ALGO>` e `TESE_OUT/ANS_TESE_DIST_<ALGO>`; o `run_dist_tese.sh dec` lê as decomposições de `TESE_OUT/ANS_TESE_DEC_<ALGO>`. Como os genomas reduzidos não têm repetições, o empacotamento de ciclos fica determinado: o `run_dec_tese.sh` usa só `rand` com uma iteração e o `run_dist_tese.sh map` usa um único mapeamento (`K_MAP=1`). `INPUT_FOLDER`, `OUT_FOLDER`, `DEC_FOLDER` e `K_MAP` continuam sobrepondo esses valores.

Os scripts pulam padrões de arquivo sem correspondência (tamanho ou operação ausente na pasta de entrada).

## Exemplo

```
# partições e instâncias reduzidas
./run_part_tese.sh

# decomposições e distâncias sobre as instâncias reduzidas do SOARComb
PART=SOARComb ./run_dec_tese.sh
PART=SOARComb ./run_dist_tese.sh map
PART=SOARComb ./run_dist_tese.sh dec
```

Teste com 10 instâncias por arquivo (SREV, `reversal_indel25`): distância média 114,5 (balanceadas) e 166,9 (com indels) sobre as instâncias reduzidas do SOARComb, contra 170,7 e 206,8 sobre as instâncias originais com 100 mapeamentos aleatórios. Teste dos scripts com `PART=SOARComb` (10 instâncias de 4 arquivos de SREV e SREV_TRANS, uma delas trocada por `NA`): `run_dec_tese.sh`, `run_dist_tese.sh map` e `run_dist_tese.sh dec` completaram e juntaram todos os arquivos, a instância `NA` ficou na sua posição em todas as saídas, e as distâncias de `map` e `dec` foram iguais em todas as instâncias.

## Saída sem redução (`run_part_fix_tese.sh`)

Mantém os genomas originais da base como entrada do `cyc` e grava, à parte, a partição em posições desses genomas, para que o `cyc` monte o grafo de adjacências com adjacências fixas, e opcionalmente com caracteres fixos (`foundation.tex:401-403`, `adj_graph.tex:44-57`).

O script reaproveita os arquivos `PART/<ALGO>/N<n>/<mod>/<arquivo>.txt` do `run_part_tese.sh` quando existem; senão roda o MCFISP-exe e grava o `PART`. Depois grava:
- `FIX/<ALGO>/N<n>/<mod>/<arquivo>.txt`: 3 linhas por instância, na mesma ordem das instâncias do arquivo de entrada (que continua sendo a entrada do `cyc`):
  1. **breakpoints da origem**: índices das regiões intergênicas cuja adjacência é breakpoint. Todas as outras adjacências estão dentro de um bloco, ou seja, são fixas;
  2. **breakpoints do alvo**, da mesma forma;
  3. **correspondência de genes**: para cada gene j da origem (j = 1..n), a posição do gene pareado no alvo, **negativa** se os blocos foram pareados no sentido reverso, e **0** se o gene está num bloco exclusivo (sem par). Os genes do alvo que não aparecem nessa linha estão em blocos exclusivos do alvo.
  Uma linha sem breakpoints é `-`. Uma instância sem partição tem as três linhas `NA`.
- `COST/...`: os mesmos custos do `run_part_tese.sh`.

Convenção de posições (a do arquivo de entrada): os genes de cada cromossomo são numerados 1..n, sem caps, e a região i é a que fica depois do gene i.
- Cromossomo linear: regiões 0..n; a região 0 fica antes do gene 1 (depois do cap) e a região n depois do gene n (antes do cap).
- Cromossomo circular: regiões 1..n; a região n fica entre o gene n e o gene 1 (é a última região da entrada). Um circular inteiro (sem breakpoint) tem `-`; um circular cortado tem tantos breakpoints quanto blocos.
- Com vários cromossomos, as regiões são numeradas por cromossomo e separadas por `;`, como na entrada; os genes da linha 3 são numerados em sequência por todos os cromossomos.

O `partition_fixed.py` recebe também o arquivo de entrada (`partition_fixed.py FIX COST PART ENTRADA`). Para um circular cortado, o MCFISP-exe escreve primeiro o bloco que atravessa a origem; a posição em que esse bloco começa é encontrada procurando a rotação do cromossomo de entrada em que genes e regiões coincidem com os dos blocos (se houver mais de uma, por simetria rotacional, qualquer uma dá blocos com o mesmo conteúdo). Dois circulares inteiros são pareados a menos de rotação e reflexão, e a correspondência de genes segue essa rotação (negativa na reflexão).

Verificação: nas 80 instâncias da base (10 de 2 arquivos SREV, 4 algoritmos), em 6000 instâncias aleatórias com cromossomos lineares e circulares (vários cromossomos; saídas do Greedy, SOAR e SOARComb), em 3000 com circulares repetidos rotacionados e em casos feitos à mão (bloco que atravessa a origem, circular cortado uma vez, circulares inteiros casados por rotação e por reflexão), todas consistentes: os genes pareados têm o mesmo rótulo (com sinal trocado nos pares reversos), genes vizinhos dentro de um bloco vão para posições vizinhas do alvo (dando a volta nos circulares), a região intergênica interna da origem cabe no intervalo da região correspondente do alvo, a adjacência correspondente do alvo também não é breakpoint, e a correspondência é injetiva. Nos lineares, os arquivos `FIX` ficaram idênticos aos da versão anterior.

### Leitura do `FIX` no `cyc`

O `dec` lê o `FIX` com as opções:

- `-f FIX_FILE` (**adjacências fixas**, linhas 1 e 2): as adjacências que não são breakpoints são fixadas no grafo (`fix_partition` em `cycle/cycles.cpp`). Uma adjacência fixa não pode virar breakpoint: se um gene de uma adjacência fixa é pareado, o gene vizinho tem de ser pareado com o vizinho correspondente no outro genoma (as duas arestas de adjacência e as duas arestas de rótulo formam um ciclo trivial balanceado, testado por intervalo); se um é indel, o vizinho também é (o bloco exclusivo inteiro é removido/inserido). As heurísticas continuam escolhendo a bijeção entre blocos compatíveis; ao fixar uma aresta de rótulo ou um indel, o fecho dessas implicações (o bloco inteiro) é fixado junto (`gray_closure`, `indel_closure`), em `bfs`, `map`, `1c2c3c` e nos operadores do `ga`. A partição induzida tem breakpoints contidos nos do `FIX` (pode ter menos).
- `-f FIX_FILE -p` (**adjacências e caracteres fixos**, as três linhas): fixa também a correspondência de genes da linha 3 (genes com `0` e genes do alvo fora da linha são indels; com um cromossomo em cada genoma, os caps são pareados). O empacotamento fica determinado.

Detalhes:
- A instância k do `FIX` corresponde à instância k da entrada; o número de partições tem de ser o de instâncias. Os cromossomos lineares têm de ser estendidos (`-e`); o `dec` confere o número de genes da linha 3 e os índices das regiões.
- Uma instância `NA` no `FIX` tem resultado `NA` (`Dec: NA` e `Dist: NA`), como nas instâncias reduzidas; o `dist -D` também escreve `NA` para uma decomposição `NA`.
- Com `-f`, as regras de redução (`suboptimal_rule_*`) não são aplicadas, porque podem fixar pares que contradizem a partição.
- Os pareamentos e indels de blocos inteiros não seguem a contagem de cópias de cada rótulo, então pode sobrar um gene sem par disponível; depois da busca, esses genes são removidos/inseridos com o seu bloco e a busca é refeita.
- **Circulares.** O `cyc` lê o `FIX` com cromossomos circulares (convenção acima). No grafo, um circular não tem caps: o gene j começa no vértice `2(j − 1)` a partir do início do cromossomo e a região i liga o gene i ao gene i + 1 (a região n liga o gene n ao gene 1). Um circular inteiro casado com outro é uma sequência de ciclos contínuos que se fecha, o que o fecho dos blocos já aceita. Com `-p`, os caps só são pareados quando cada genoma tem um único cromossomo, e linear; com vários cromossomos, os caps não aparecem na linha 3 e a heurística os pareia.
- **Caps cruzados.** Com um único cromossomo linear em cada genoma, o grafo usa caps distintos (o cap esquerdo de um genoma só é pareado com o cap esquerdo do outro), mas a partição do MCFISP-exe rotula os dois caps com `0` e pode parear, no sentido reverso, o bloco do cap esquerdo da origem com o bloco do cap direito do alvo (gene 1 da origem com `-n` do alvo na linha 3), ou o do cap direito da origem com o do cap esquerdo do alvo (gene n da origem com `-1`). Essas adjacências de cap não podem ser mantidas no grafo, então o `fix_partition` as trata como breakpoints (as demais adjacências do bloco continuam fixas). Antes disso o bloco ficava sem par e caía no beco sem saída abaixo, e com `-p` a correspondência era rejeitada (`the correspondence breaks the fixed adjacency of vertices 0 and 1`). Nas partições da base isso ocorre em 167 instâncias do Greedy e 288 do SOARComb (de 135975 cada; nenhuma mistura os dois sentidos).
- **Becos sem saída.** Fora o caso dos caps cruzados, existe sempre um empacotamento que respeita as adjacências fixas, mas a escolha gulosa pode não chegar a ele (por exemplo, um bloco que termina num cap fica sem par porque o único gene da outra ponta ao lado de um cap livre virou indel). Depois das passadas de indels, as adjacências fixas dos vértices que ainda sobrarem são liberadas e a busca é refeita, para que a decomposição fique completa; o `dec` avisa no stderr (`instance k: m fixed adjacencies released`). Nesses casos a partição induzida pode ter breakpoints que não estão no `FIX`.
- Sem `-f`, o comportamento é o anterior (com `rand bfs -k 1`, que é determinístico, as saídas são idênticas às do commit anterior).

Nos scripts, `FIX=<ALGO>` faz o `run_dec_tese.sh` rodar sobre as instâncias originais com `-f PART_FOLDER/FIX/<ALGO>/N<n>/<mod>/<arquivo>.txt` e gravar em `TESE_OUT/ANS_TESE_DEC_FIX_<ALGO>`; com `FIX_PAIRS=1` usa também `-p`, grava em `TESE_OUT/ANS_TESE_DEC_FIX_<ALGO>_PAIRS` e usa só `rand` com uma iteração. O `run_dist_tese.sh dec` com as mesmas variáveis lê essas decomposições e grava em `TESE_OUT/ANS_TESE_DIST_FIX_<ALGO>[_PAIRS]` (`map` não usa a partição e não é aceito com `FIX`).

```
./run_part_fix_tese.sh
FIX=SOARComb ./run_dec_tese.sh
FIX=SOARComb ./run_dist_tese.sh dec
FIX=SOARComb FIX_PAIRS=1 ./run_dec_tese.sh
FIX=SOARComb FIX_PAIRS=1 ./run_dist_tese.sh dec
```

Verificação (10 instâncias de 4 arquivos de SREV e SREV_TRANS, com e sem indels, F0/F40/F100; partições dos 4 algoritmos; `rand bfs`, `rand map` e `ga bfs`; `-f` e `-f -p`; build de debug): nas 960 decomposições, todo vértice está num ciclo, nenhuma adjacência fixa da origem ou do alvo virou breakpoint (os dois genes foram pareados com genes adjacentes, ou os dois são indels) e, com `-p`, o pareamento é o da linha 3. Com `-f -p` as distâncias foram iguais às das instâncias reduzidas do mesmo algoritmo. Distância média (SOARComb, `rand bfs -k 20`, `reversal_indel25` em SREV e `reversal_transposition_indel6` em SREV_TRANS):

| arquivo | sem partição | `-f` | `-f -p` |
|---|---|---|---|
| N200_SREV100_I0_L100_F0 | 136,0 | 103,3 | 116,6 |
| N200_SREV100_I50_L10_F40 | 175,4 | 136,2 | 156,3 |
| N200_SREV_TRANS100_I50_L100_F100 | 187,1 | 166,8 | 171,9 |
| N200_SREV_TRANS200_I0_L10_F0 | 211,4 | 193,2 | 205,7 |

O tempo com `-f` variou de 0,7 a 1,9 vez o tempo sem partição, por instância.
