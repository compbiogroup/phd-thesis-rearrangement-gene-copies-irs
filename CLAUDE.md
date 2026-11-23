# Contexto para o Claude (code_tese)

## Origem e estado do repositório

- `cyc/` e os scripts `run_dec_tese.sh` e `run_dist_tese.sh` foram copiados de `Gabriel-Siqueira/Doutorado` no commit `ecbf3a1`. `chim_flex/` e os scripts de partição foram trazidos por outra instância.
- `code_tese` e `GitHub/cyc` (repositório do doutorado) evoluem separados: mudanças num não chegam no outro. Vale decidir qual dos dois é a referência.
- `cyc/quickcheck` está copiado como arquivos comuns (QuickCheck++, GPLv3, usado só nos testes); não é submódulo.
- `cyc/main_dec.cpp` está na versão commitada, com as regras subótimas do MSOAR **ativas**. No `GitHub/cyc`, o usuário tem uma alteração local, não commitada, que as comenta.
- `cyc/BUGS.md` tem a lista de bugs em aberto do cyc, com as correções sugeridas.

## Padrões dos scripts (decididos pelo usuário)

- Todas as saídas ficam em `code_tese/TESE_OUT/` (caminhos relativos a `cyc/` nos scripts do cyc).
- A decomposição (`run_dec_tese.sh`, `run_dist_tese.sh`) lê `~/Documents/TESE_SMALL`, um subconjunto de arquivos inteiros do `TESE` (mesmos nomes e conteúdo, veja `copy_db.sh` no repositório do doutorado). Por isso os arquivos `FIX` da partição batem pelo nome com as entradas do `dec`.
- A partição lê `TESE`. Os dois scripts de partição ainda apontam para lugares diferentes: `~/Dropbox/Context/doutorado/Dados/TESE` (`run_part_tese.sh`) e `~/Documents/TESE` (`run_part_fix_tese.sh`). Como o segundo reaproveita os arquivos `PART` do primeiro, vale alinhar.
- Operações (`mods`): `SDCJ SREV SREV_TRANS`.

## 1c2c3c (`cyc/cycle/cycles.cpp`)

- `decompose_1c_2c_3c`, em ordem:
  1. `mark_indel_edges`;
  2. `add_small_cycles(s)` para s = 1, 2, 3 (ciclos com 2, 4 e 6 arestas pretas, menores primeiro);
  3. `decompose_with_bfs` para completar.
- Na etapa de ciclos pequenos só entram ciclos **limpos** (sem arestas de indel) e **balanceados**. O empate é decidido pela ordem aleatória das arestas pretas.
- `add_small_cycles` usa meet-in-the-middle: para cada aresta preta, guarda numa tabela hash os caminhos de s arestas cinzas que chegam ao início e junta com os caminhos de s arestas cinzas pela frente. Custo O(d^s) por aresta, com d = número de cópias de cada gene.
- A versão anterior enumerava e guardava todos os caminhos com até 6 arestas pretas (~n·d⁵) e era morta por falta de memória com alfabeto pequeno (L ≤ 10, N=200). A nova leva menos de 1 s nesses casos.
- Restos que o 1c2c3c não usa mais: `get_all_tuptok_cycles`, `get_tuptok_cycle`, `select_cycles`, `find_2_cycles`.
- O `--help` do `dec` ainda diz que o 1c2c3c é só para genomas com uma cópia de cada gene; está desatualizado.

## A tese não bate com o código

Em `tese/tese_tex/adj_graph.tex` (repositório do doutorado), seção "Seleção de Ciclos Pequenos":

- O texto descreve MIS para os 2-ciclos, 6-SP para os ciclos até 3 e a união dos dois ordenada por tamanho. A implementação agora é incremental por tamanho, só com ciclos limpos e balanceados, e completa com BFS (não com random map).
- O Algoritmo `alg:scp` tem dois erros de texto: o laço final diz "ciclo C em ℰ" em vez de ℰ₄, e a frase antes dele diz "pseudocódigo da heurística de busca em largura".

## Como testar

- A base `~/Documents/TESE*` não existe nesta máquina; gere instâncias sintéticas.
- Formato da entrada: 4 linhas por instância (genes da origem, regiões intergênicas da origem, genes do alvo, regiões do alvo). Os genomas lineares têm n+1 regiões e devem ser rodados com `-e`. Use `-d` quando os genomas não têm os mesmos genes e `-z` para pesos zero (2 linhas por instância).
- O `rand` roda as iterações em paralelo com OpenMP, então a memória se multiplica pelo número de threads (16 núcleos, 27 GB nesta máquina).
- Para comparar versões, compile o `dec` do HEAD separadamente e rode os dois nos mesmos arquivos.
- Para conferir propriedades dos ciclos adicionados, um log temporário no `stderr`, removido depois, funcionou bem.

## Preferências do usuário

- Conversa em português; código, comentários e mensagens de commit em inglês.
- Commit só quando o usuário pedir, terminando com a linha `Co-Authored-By` do Claude.
- Não incluir em commits alterações locais do usuário sem perguntar.
