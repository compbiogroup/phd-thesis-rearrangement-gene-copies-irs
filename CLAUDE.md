# Contexto para o Claude (code_tese)

## Origem e estado do repositório

- `cyc/` e os scripts `run_dec_tese.sh` e `run_dist_tese.sh` foram copiados de `Gabriel-Siqueira/Doutorado` no commit `ecbf3a1`.
- `partition/` (antes `chim_flex/`) foi copiado do submódulo `GitHub/chim_flex` (branch `unbalanced`, commit `e3f0553`): todos os arquivos rastreados, menos `datasets/` (dados antigos, sem uso). Os scripts de partição vieram de `GitHub/scripts` (commit `ecbf3a1`), com as saídas padrão movidas para `TESE_OUT/`. Os dois lados evoluem separados, como o `cyc`.
- `code_tese` e `GitHub/cyc` (repositório do doutorado) evoluem separados: mudanças num não chegam no outro. Vale decidir qual dos dois é a referência.
- `cyc/quickcheck` está copiado como arquivos comuns (QuickCheck++, GPLv3, usado só nos testes); não é submódulo.
- `cyc/main_dec.cpp` está na versão commitada, com as regras subótimas do MSOAR **ativas**. No `GitHub/cyc`, o usuário tem uma alteração local, não commitada, que as comenta.
- `cyc/BUGS.md` tem a lista de bugs em aberto do cyc, com as correções sugeridas.

## Padrões dos scripts (decididos pelo usuário)

- Todas as saídas ficam em `code_tese/TESE_OUT/` (caminhos relativos a `cyc/` nos scripts do cyc).
- A base da tese está em `db/` (5000 arquivos, 2,6 GB), e todos os scripts a leem por padrão. Os `.txt` de `db/` vão para o Git LFS (`.gitattributes`); o `db/README.md` não. Só dê `git add db` com o `git-lfs` instalado e configurado (`git lfs install`), senão os arquivos entram como texto comum.
- A partição usa a base inteira. A decomposição (`run_dec_tese.sh`, `run_dist_tese.sh`) usa só o subconjunto que antes ficava em `~/Documents/TESE_SMALL` (`copy_db.sh` no repositório do doutorado): N200, op 100 e 200, I 0 e 50, L 10, 100 e 200, F 0, 40 e 100. Os padrões de `OP_VALS`, `I_VALS`, `L_VALS` e `F_VALS` (`${VAR-padrão}`, vazio mantém todos) selecionam esse subconjunto, que é o mesmo das saídas antigas em `~/Dropbox/Context/doutorado/Dados/TESE/ANS/` (108 arquivos). Os nomes são os mesmos da base, então os arquivos `FIX` e `RED` da partição batem com as entradas do `dec`.
- Operações (`mods`): `SDCJ SREV SREV_TRANS`.
- Partição: `run_part_tese.sh` roda `Greedy GreedySin SOAR SOARComb` e `run_part_fix_tese.sh` roda `Greedy SOARComb Approx`. Os dois aceitam `ALGOS` e gravam em `TESE_OUT/ANS_TESE_PART` (`PART`, `COST`, `RED`, `FIX`, `PHASE2`).

## Partição de genomas (`partition/`)

### Código e tese

- O texto é o capítulo `tese/tese_tex/partition.tex` do repositório do doutorado. O módulo `partition/src/PApprox.hs` segue o Algoritmo 1 atual: conjuntos B_C por caso do Lema `tmin_comp`, classes ativas e o Lema `tmin_ativa`.
- O fator 2tk do Teorema `teo:aprox_part` só vale para relações simétricas e transitivas (trivial e igualdade). A segunda etapa, que corta um bloco D, só roda com relações não transitivas (intervalo, F > 0), e a tese não garante fator nesse caso. O número de usos dessa etapa vai para `<saída>_phase2`.
- O `ApproxLB` devolve q, o número de ativações de B_C, que é a cota inferior da prova. Para genomas desbalanceados, o `Balancer` faz com que ele não seja um limitante de verdade.
- O C escolhido em cada passo é a componente imperfeita com **menos adjacências**, não com menos genes. Um circular inteiro de n genes contém subcromossomos de n genes e não está em T^min; escolher por genes quebrava o Lema `tmin_ativa`.
- O `partition/CompatibilityWithText.md` compara o código com a tese. A seção "Em aberto", no começo, lista o que falta no texto (Lema `tmin_comp`, circulares inteiros no grafo, índices das relações, SOAR) e no código. `partition/BUGS.md` registra os bugs corrigidos. Atualize os dois, e `partition/doc/*.md`, quando mudar o comportamento.
- Notas de pesquisa sobre o fator tk (`ProvaTk.md`) ficaram só em `GitHub/chim_flex`, sem versionamento.
- O FPT (`PFpt.hs`) aceita genomas desbalanceados: Regra 2 da tese, com o conjunto D no campo `removed` do grafo. Os genes de rótulos ausentes no outro genoma vão para D já no grafo inicial. A poda soma a `2|A| − n_c` um limitante de blocos exclusivos (`exclusiveRuns`), calculado só quando algum rótulo tem ocorrências diferentes. A busca começa com a melhor partição do Guloso e do SOARComb (o custo continua o ótimo, mas a partição devolvida pode ser a da heurística quando ela já é ótima). O tempo limite do FPT é a opção `-t` do `MCFISP-exe` (padrão 600 s).
- Experimento do FPT: com 200 genes ele não termina (nem com a solução inicial). Calibração com 60 s por instância (`gen_db_fpt.sh`, F = 0): resolve 86% com N = 25, 33% com N = 50 e 1% com N = 100; o número de cópias (N/L) é o que mais pesa. `run_fpt_small.sh` (N = 25 e 50, 100 instâncias por arquivo, 600 s) foi estimado em cerca de 30 h. O travamento raro do `MCFISP-exe` está no #30 do `partition/BUGS.md`.
- A tese (`partition.tex:497`) diz que os blocos exclusivos de D só crescem ao longo de um ramo, mas duas sequências de D podem se juntar; o código conta só as que nunca se juntam (item 9a do `CompatibilityWithText.md`).

### Comportamentos que outros arquivos supõem

- O `Main` iguala o número de cromossomos lineares acrescentando lineares vazios (`equalizeLinear`), como em `foundation.tex:109`. Eles aparecem na saída e contam no custo. `reduce_partition.py` e `partition_fixed.py` não os tratam; na base da tese G e H sempre têm o mesmo número de lineares, então isso não acontece.
- O nome `<saída>_phase2` não pode casar com `<saída>_[0-9]*`, que os scripts usam para contar as saídas das instâncias.
- O `Approx` usa a compatibilidade calculada uma vez em G₀ e atualiza o grafo incrementalmente (`subgraphOf`, `removeAdjacency`). Isso depende de a compatibilidade entre dois subcromossomos não depender da partição.

### Como testar

- `cd partition && stack test`. Os scripts compilam com `stack install --local-bin-path .`, e o `stack.yaml` usa Nix. O `stack` não está instalado nesta máquina; use `nix-shell -p stack --run "stack test"`.
- O `PFptCheck` compara o FPT com uma força bruta (`bruteForceCost`) em instâncias pequenas, balanceadas e desbalanceadas, com a relação de igualdade.
- Ao otimizar, as saídas não devem mudar. Compile o executável antes e depois e compare as saídas, sem as linhas `# Time`, em muitas instâncias, com e sem sinal, com circulares, com intervalo (`-f`) e com números diferentes de lineares.
- O `DB-exe` gera instâncias. O `RandDB` exige `-d`. No `RepDB`, a opção `-r` serve para duas coisas e conflita. O `DB-exe` nunca gera números diferentes de lineares; para isso, escreva um gerador próprio.
- O FPT é exato. Para conferir, compare com uma força bruta que enumera os conjuntos de breakpoints de instâncias pequenas, incluindo circulares inteiros e circulares cortados uma vez. Foi assim que apareceu o bug #29.
- Para contar breakpoints na saída, use a linha de genes nos lineares (número de `|`). Nos circulares, use a linha das regiões intergênicas (metade do número de `|`, arredondada para cima), porque a linha de genes não mostra o corte de fechamento.
- Números de referência do `Approx` depois da otimização: 3 instâncias com 150 genes levam cerca de 2 s (antes, 63 s).

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

- A base está em `db/` (cópia de `~/Dropbox/Context/doutorado/Dados/TESE/`, que também tem as saídas antigas em `ANS/`); gere instâncias sintéticas também.
- Formato da entrada: 4 linhas por instância (genes da origem, regiões intergênicas da origem, genes do alvo, regiões do alvo). Os genomas lineares têm n+1 regiões e devem ser rodados com `-e`. Use `-d` quando os genomas não têm os mesmos genes e `-z` para pesos zero (2 linhas por instância).
- O `rand` roda as iterações em paralelo com OpenMP, então a memória se multiplica pelo número de threads (16 núcleos, 27 GB nesta máquina).
- Para comparar versões, compile o `dec` do HEAD separadamente e rode os dois nos mesmos arquivos.
- Para conferir propriedades dos ciclos adicionados, um log temporário no `stderr`, removido depois, funcionou bem.

## Preferências do usuário

- Conversa em português; código, comentários e mensagens de commit em inglês.
- Commit só quando o usuário pedir, terminando com a linha `Co-Authored-By` do Claude.
- Não incluir em commits alterações locais do usuário sem perguntar.
