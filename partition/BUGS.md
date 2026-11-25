# Relatório de bugs

Revisão do estado atual da árvore de trabalho (branch `unbalanced`, com as modificações ainda não commitadas). Li todos os arquivos de `src/`, `app/` e `test/`, compilei com `stack build --fast` e rodei os executáveis em instâncias pequenas (feitas à mão e geradas pelo `DB-exe`).

**Legenda**
- **[confirmado]**: reproduzido rodando o código.
- **[leitura]**: encontrado lendo o código, sem reprodução.

A maior parte dos problemas listados em `AddCircular.md` já foi corrigida (`nextIdx`/`prevIdx`/`adjacencies`, `subGenome` com volta na origem, `intergenicFullReversal`, `canonicPosition`, índices de `ChromList` etc.). Os que ainda restam estão incluídos abaixo.

---

## Resumo

| # | Sev. | Local | Problema |
|---|------|-------|----------|
| 1 | ✅ corrigido | `Partition.hs:204` | `getBlocksFit` quebra com `!!: index too large` (todo genoma desbalanceado falha) |
| 2 | ✅ corrigido | `Genomes.hs:573,585` | `mkRGenome`/`mkFGenome` quebram com cromossomo vazio (`maximum []`) |
| 3 | ✅ corrigido | `PFpt.hs:183` | FPT não põe breakpoint no fechamento de cromossomo circular → partição não comum → crash |
| 4 | ✅ corrigido | `PFpt.hs:335` | FPT (que deveria ser exato) devolve partições subótimas |
| 5 | ✅ corrigido | `PSOAR.hs:40` | SOAR devolve a partição dos genomas **modificados** (genes trocados na saída) |
| 6 | ✅ corrigido | `PApprox.hs:160` | `approxPartition` entra em loop infinito quando `findBreakpoint` não progride |
| 7 | ✅ corrigido | `PApprox.hs:354` | `updateEntry` mistura os lados esquerdo/direito quando a orientação é invertida |
| 8 | ✅ corrigido | `Genomes.hs:754` | `intergenicDeletion` circular com `j == size` descarta a IR de fechamento |
| 9 | ✅ corrigido (assert) | `Genomes.hs:685,724` | Reversão/transposição circular que cobre a origem escreve duas vezes na mesma IR |
| 10 | ✅ corrigido | `PFpt.hs:324` | FPT emparelha genes de orientação incompatível |
| 11 | ✅ corrigido | `app/DB.hs:302` | Reversões/transposições nunca atingem os 2 últimos genes |
| 12 | ✅ corrigido | `test/` | A suíte de testes não compila |
| 13 | ✅ corrigido | `app/Main.hs:87` | Linha em branco no arquivo de entrada derruba o programa |
| 14 | ✅ corrigido | `Genomes.hs:646` | `writeMultiC` escreve os cromossomos em ordem inversa |
| 15 | ✅ corrigido | — | Diversos (ver seção 15) |
| 16 | ✅ corrigido | `Genomes.hs:1030` | `randomGenomeWithReplicas` gera mais de `n` genes (entrada inválida) |
| 17 | ✅ corrigido | `Genomes.hs` (`shuffleGenome`) | `DB-exe -r -1` embaralha os caps junto com os genes |
| 18 | ✅ corrigido (entrada) | `Genomes.hs` (`readRGenome`) | Instância com sinal lida sem `-s` quebra com `logicError` |
| 19 | ✅ corrigido | `PSOAR.hs:68-83` | `isConflict` só verifica a consistência quando os duos se seguem em G → circular inteiro em H contra linear em G → `logicError` |
| 20 | ✅ corrigido | `PFpt.hs:209-216` | `closesCircle` não compara a IR de um circular de 1 gene → `logicError` |
| 21 | ✅ corrigido | `Partition.hs:347-430` | `combine` (SOARComb) funde blocos sem respeitar o emparelhamento → `logicError` |
| 22 | ✅ corrigido | `PFpt.hs:397-405` | FPT não é exato quando um gene casa nos dois sentidos (caps `0`) |
| 23 | ✅ corrigido | `app/Main.hs` (`runOne`) | Uma instância inválida derruba a execução inteira (no paralelo, nenhum arquivo é escrito) |
| 24 | ✅ corrigido | `test/Spec.hs` | `stack test` nunca falha; geradores só fazem um linear com sinal; FPT sem teste |
| 25 | A confirmar | `PSOAR.hs:36-39` | Regras subótimas do SOAR desativadas na árvore de trabalho |
| 26 | Baixa | — | Diversos (ver seção 26) |
| 27 | ✅ corrigido | `PFpt.hs` (`cleanEdges`) | Limpeza de arestas passa do limite de uma testemunha que cobre um circular inteiro → grafo assimétrico → `patternError` |
| 28 | ✅ corrigido | `PFpt.hs` (`sample0`) | Sem sinal, os singletons entravam na amostra inicial só no sentido direto → FPT perdia blocos reversos |
| 29 | ✅ corrigido | `PFpt.hs` (`findLimitsDir`) | Num circular, a testemunha só considerava uma janela máxima (direita primeiro) → FPT perdia blocos que são o circular cortado em outra posição |
| 30 | Em aberto | `PFpt.hs` / `Main.hs` (paralelismo) | Raramente, o `MCFISP-exe -a FPT` para de calcular (CPU perto de zero) sem terminar algumas instâncias, mesmo com tempo limite |

---

## 1. `getBlocksFit` quebra e o emparelhamento está errado [confirmado] — ✅ corrigido

> Corrigido com o algoritmo de Kuhn indexado por bloco de H, propagando `seen` também nos ramos que falham. `checkCommonUnbal` passou a exigir interseção vazia, ignorando os genes dos caps (`geneMapDelete`).

`src/Partition.hs:204-228`

```haskell
getBlocksFit matcher pg ph = foldr selectBlock Map.empty [1 .. length corresp]
  ...
  case foldrM tryFit (Set.empty, blocks_fit_) (corresp !! b_g1) of
```

Há três defeitos:

1. **Índices fora do intervalo.** `corresp` é indexado a partir de 0, mas o laço vai de `1` a `length corresp`. `corresp !! length corresp` lança `Prelude.!!: index too large`, e o bloco 0 nunca é considerado.
2. **Mapas com chaves trocadas.** `blocks_fit` mapeia bloco de G → bloco de H. Mesmo assim, `tryFit` testa `b_g2 `Map.notMember` blocks_fit` e faz `blocks_fit Map.! b_g2`, tratando um bloco de H como se fosse chave de G.
3. **Caminho aumentante errado.** Quando um bloco de H já está ocupado, o bloco de G que o ocupa nunca é realocado. `tryFit` é chamado de novo com um índice de H, e ao final dois blocos de G podem apontar para o mesmo bloco de H.

Além disso, `checkCommonUnbal` (linha 244) usa `zip [1 ..]` para os blocos. Isso bate com as chaves 1-based de G, mas não com os valores 0-based de H.

**Efeito:** `checkCommonUnbal` e `costUnbalanced` falham sempre que são avaliados. Na prática, **todos os algoritmos (SOAR, Greedy, FPT) quebram em qualquer instância desbalanceada**, e o FPT quebra também em instâncias balanceadas quando gera uma partição não comum (ver #3).

```
$ MCFISP-exe -a SOAR -i unbal.txt ...   # G = 1 2 3 4 5, H = 1 2 9 4 5
MCFISP-exe: Prelude.!!: index too large
```

**Correção sugerida:** reescrever com o algoritmo de Kuhn, indexado por bloco de H (como `findUncommon` já faz corretamente) e com índices 0-based:

```haskell
-- | Emparelhamento máximo: bloco de G (índice 0-based) -> bloco de H (índice 0-based)
getBlocksFit :: (Matcher m c1 c2) => m c1 c2 -> Partition c1 -> Partition c2 -> Map.Map Int Int
getBlocksFit matcher pg ph = flipMap $ foldl selectBlock Map.empty [0 .. length corresp - 1]
  where
    corresp = getBlocksCorrespondence_ matcher pg ph
    -- assigned: bloco de H -> bloco de G
    selectBlock assigned b_g = maybe assigned snd (augment Set.empty b_g assigned)

    augment seen b_g assigned = go seen (corresp !! b_g)
      where
        go _ [] = Nothing
        go seen' (b_h : rest)
          | b_h `Set.member` seen' = go seen' rest
          | otherwise =
              let seen'' = Set.insert b_h seen'
               in case Map.lookup b_h assigned of
                    Nothing -> Just (seen'', Map.insert b_h b_g assigned)
                    Just b_g' -> case augment seen'' b_g' assigned of
                      Just (s, assigned') -> Just (s, Map.insert b_h b_g assigned')
                      Nothing -> go seen'' rest
```

Depois disso, troque `zip [1 ..] bs` por `zip [0 ..] bs` em `checkCommonUnbal`. (Para manter a complexidade polinomial, o `seen` de uma tentativa que falhou deveria ser propagado para as irmãs, como em `findUncommon`.)

> **Ponto a verificar:** `checkCommonUnbal` exige `geneMapDifference gmap_g gmap_h == empty`, ou seja, que os genes dos blocos exclusivos de G estejam **contidos** nos genes dos blocos exclusivos de H. Pela definição usual, o que se quer é que os blocos exclusivos de G e de H **não compartilhem** genes (interseção vazia). Se for isso, a condição está invertida.

---

## 2. `mkRGenome`/`mkFGenome` quebram com cromossomo vazio [confirmado] — ✅ corrigido

> `mkRGenome` e `mkFGenome` usam `maximum (0 : filter (/= maxBound) (map abs genes)) + 1`: aceitam cromossomos vazios e consideram os caps, exceto `maxBound`. Excluir `maxBound` é necessário porque ele é usado como cap direito pelo `DB-exe` e pelos testes, e `maxBound + 1` estoura para `minBound` (os testes de `makeSingletons` detectaram isso).

`src/Genomes.hs:565` e `:577`

```haskell
new = maximum (map abs genes_) + 1
```

Um DCJ pode produzir um cromossomo linear que só tem os dois caps. O `DB-exe` grava esse cromossomo como `L  ;`, e ao ler de volta `genes_ == []` faz `maximum []` lançar exceção. Como o valor é preguiçoso, o erro só aparece quando `getNewGene` é usado. Por isso o SOAR quebra (as regras subótimas chamam `makeSingletons`) e o Greedy/Approx não.

```
$ DB-exe RandDB -l 5 -k 40 -n 10 -r 6 -s -d 0 -c -f 0 -o mc
$ MCFISP-exe -a SOAR -i mc ...       # instância 21: "C 3 -4 ... ; L  ;"
MCFISP-exe: Prelude.maximum: empty list
```

**Correção:**

```haskell
new = maximum (0 : map abs genes) + 1   -- usa `genes` (com caps) e aceita lista vazia
```

Usar `genes` em vez de `genes_` também protege o caso em que os caps não são `maxBound`.

---

## 3. FPT: nenhum breakpoint no fechamento de um cromossomo circular [confirmado] — ✅ corrigido

> Corrigido em `breaksOf`: a adjacência que chega ao primeiro gene da testemunha (`startOf`) vira breakpoint, a menos que a testemunha cubra dois circulares inteiros e as adjacências de fechamento casem (`closesCircle`).

`src/PFpt.hs:183-190` (`sampleGraphToPartition.breaksOf`)

Uma adjacência `(i, i')` só vira breakpoint quando os genes vêm de testemunhas diferentes. Em um cromossomo circular inteiramente coberto por uma única testemunha (por exemplo, um circular de 1 gene, ou um circular cujo conteúdo é um segmento linear em H), a adjacência de fechamento `(n, 1)` também tem o mesmo pai. Assim nenhum breakpoint é criado e o bloco vira o cromossomo circular inteiro (`b_whole`). Só que `findLimits` nunca verificou essa adjacência, e o bloco correspondente em H é linear. Como `RFDM` exige `chtypeG == chtypeH`, os blocos não casam, `mkCommonPartition` cai em `checkCommonUnbal` e o programa quebra (#1).

```
G: L 4 ; C 1 2 3        H: L 1 2 3 ; C 4
$ MCFISP-exe -a FPT ...
MCFISP-exe: Prelude.!!: index too large      (Greedy/SOAR/Approx resolvem)
```

**Correção:** uma adjacência só é interna à testemunha se também estiver entre os limites dela. Por exemplo:

```haskell
(i, i') <- adjacencies chr,
let p = parentOf inG (chri, i),
p /= parentOf inG (chri, i') || i' == leftLimit inG p   -- i' começa a testemunha: é corte
```

Aqui `leftLimit inG p` é `fst (w_left_pair p)` em G, ou o extremo esquerdo em H (`w_right_pair` quando `w_inverted`). A outra opção é aceitar o bloco circular inteiro só quando o bloco correspondente em H também for um circular inteiro.

---

## 4. FPT devolve partições subótimas [confirmado] — ✅ corrigido

> `createSampleGraph` agora ignora pares paralelos a (ou em conflito com) uma testemunha que já está na amostra. Com isso, a instância de exemplo passou de 8 para 4 blocos. Nas 120 instâncias de teste (lineares, multicromossomo e com repetições), o FPT dá o mesmo número de blocos que o melhor entre Approx e Greedy.

`src/PFpt.hs:335` (`createSampleGraph`)

`sample0` contém **todos** os singletons, e `createSampleGraph` adiciona uma aresta preta para cada um sem verificar se o par já é paralelo a uma testemunha existente. Quando dois singletons estão no mesmo bloco comum, a segunda `addBlackEdge` chama `cleanEdges`, que apaga as arestas verdes da primeira testemunha. O bloco é então dividido em dois.

```
G: L 1 2 3 4 5 ; C 6 7 8      H: L 1 -3 -2 4 5 ; C 7 8 6
Greedy/SOAR/Approx:  L 0 1 | 2 3 | 4 5 0 ; C 6 7 8            (4 blocos)
FPT ("Exact"):       L 0 1 | 2 | 3 | 4 | 5 0 ; C 6 | 7 | 8    (8 blocos)
```

**Correção:** em `createSampleGraph`, ignore um par se ele for paralelo (ou em conflito) a uma testemunha que já está na amostra, como `edgesFromRareOddPath` já faz:

```haskell
createSampleGraph matcher g h = foldr add sg0
  where
    add idxs sg =
      let bw = findLimits sg idxs
       in if any (\bw' -> checkParallel g h bw bw' || checkConflict bw bw') (sample sg)
            then sg
            else addBlackEdge bw sg
```

---

## 5. SOAR devolve a partição dos genomas modificados [confirmado] — ✅ corrigido

> `part` agora usa `g_`/`h_` (genomas originais). Testado com as regras subótimas ativas: a saída mantém os genes originais e as partições continuam comuns.

`src/PSOAR.hs:36-40`

```haskell
(g, h) = uncurry (suboptimalRulePairs matcher) (suboptimalRuleInterval matcher g_ h_)
part = mkCommonPartition2 matcher g breaksG h breaksH
```

As regras subótimas trocam genes por singletons novos, e a partição é montada sobre esses genomas modificados. Por isso a saída mostra genes que não existem na entrada (inclusive caps trocados):

```
entrada G: L 1 2 3 4 5 ; ...
saída   G: L 11 1 | 9 10 | 4 5 12 ; ...
```

**Correção:** as trocas preservam posições, então basta usar os genomas originais com os breakpoints calculados:

```haskell
part = mkCommonPartition2 matcher g_ breaksG h_ breaksH
```

As trocas são consistentes entre G e H (só ocorrem onde os genes casavam), então a partição continua comum nos genomas originais.

---

## 6. `approxPartition` pode entrar em loop infinito [confirmado] — ✅ corrigido

> `approxPartition` para com `logicError` quando uma iteração não adiciona breakpoint e a partição ainda não é comum (antes repetia a mesma iteração para sempre). Pelo `MCFISP-exe`, o `Approx` passa pelo Balancer e só recebe genomas balanceados.

`src/PApprox.hs:160-169` e `:415-429`

Quando `findBreakpoint` devolve `Nothing`, `compatibilityCompletion` retorna as mesmas partições com `False`. A iteração seguinte de `until done2` encontra os componentes já perfeitos (`done1` imediato) e chama `compatibilityCompletion` sobre a mesma entrada, e isso se repete para sempre. Com uma instância desbalanceada simples (`1 2 3 4 5` × `1 2 9 4 5`), o processo não termina (foi morto pelo `timeout` de 60 s).

**Correção:** detectar quando não há progresso. Por exemplo, fazer `compatibilityCompletion` devolver um terceiro estado ("travado") e, nesse caso, lançar um erro claro ou usar um fallback, como quebrar o próprio bloco sem par em blocos unitários. No mínimo, comparar os breakpoints antes e depois de cada iteração e parar com `error` se não mudaram.

---

## 7. `checkIntersections.updateEntry` mistura esquerda e direita [leitura] — ✅ corrigido

> `updateEntry` agora recebe só o tamanho da nova interseção e o lado (`onRight`), e troca o lado quando a orientação do vértice difere da do representante do componente. O import de `<|>`, que ficou sem uso, foi removido.

`src/PApprox.hs:344-361`

```haskell
updateInterMap ((v, cv), (u, cu), s) interMap_aux =
  let interMap_up1 = updateEntry cu cv (vOri u) lu ru interSize ru interMap_aux
   in updateEntry cv cu (vOri v) lv rv lv interSize interMap_up1
...
else IntMap.insert c' ((select l_old r_new, select r_old l_new), ...)
```

Os parâmetros `l_new`/`r_new` recebem o valor antigo do lado que **não** muda (por exemplo, `ru` para `cu`). Quando a orientação do vértice é diferente da do representante, os lados são trocados e o valor antigo de um lado acaba comparado com o do outro: para `cu` fica `left = min lu ru`. O mínimo de cada lado é contaminado, o que muda a classificação em `separateInterCases` e as posições escolhidas em `addBreaks_ii`.

**Correção:** passe só a nova interseção e o lado em que ela ocorre:

```haskell
-- side = True: interseção à direita do vértice; False: à esquerda
updateEntry c' c'' ori side size m =
  let flipSide = ori /= vOri (cVtx (comps IntMap.! c'))
      right = side /= flipSide
      ((l, r), inter) = m IntMap.! c'
      sel old = Just (maybe size (min size) old)
   in IntMap.insert c' (if right then (l, sel r) else (sel l, r), IntSet.insert c'' inter) m
-- uso: updateEntry cu cv (vOri u) False s . updateEntry cv cu (vOri v) True s
```

---

## 8. `intergenicDeletion` circular com `j == size g` perde a IR de fechamento [confirmado] — ✅ corrigido

> O `if j == size g` foi removido: a IR de fechamento do circular é mantida, e no linear o `slice` final já é vazio. Também foi adicionado um `assert` que proíbe deletar um cromossomo circular inteiro (`i = 1`, `j = n + 1`), caso em que o ramo `i == 1` faria um `slice` de tamanho negativo.

`src/Genomes.hs:754-756`

```haskell
Vec.++ if coerce j == size g then Vec.empty else Vec.slice (coerce j - 1) ...
```

No circular, deletar os genes `i..n-1` (`j = n`) deve manter a IR `n` (entre o gene `n` e o gene 1), mas o `if` a descarta.

```
c = (1)-10-(2)-20-(3)-30-(4)-40-   ;  intergenicDeletion 2 4 5 c
genes [1,4], IRs [5]     -- esperado IRs [5, 40]
```

**Correção:** o `if` é desnecessário mesmo no caso linear, porque `Vec.slice (n-1) 0` é válido. Remova-o:

```haskell
| otherwise =
    Vec.slice 0 (coerce i - 2) vi
      Vec.++ Vec.singleton (R x)
      Vec.++ Vec.slice (coerce j - 1) (Vec.length vi - coerce j + 1) vi
```

---

## 9. Reversão/transposição circular que cobre a origem [confirmado para reversão] — ✅ corrigido com assert

> Os casos que cortam a mesma IR duas vezes são proibidos por `assert`: reversão circular com `i = 1` e `j = n`, e transposição circular com `i = 1` e `k = n + 1`.

`src/Genomes.hs:685-686` (reversão) e `:724-726` (transposição)

- **Reversão** com `i = 1` e `j = n` num circular: `ir_i_idx = n-1` e `j-1 = n-1` são a **mesma** IR. A primeira escrita (`x+y`) é sobrescrita pela segunda (`x_rest+y_rest`), e a soma das IRs muda: no teste foi de 100 para 133.
- **Transposição** com `i = 1` e `k = n+1`: `ir_i_idx` e `k-2` também coincidem.

Os asserts permitem esses casos. Para um circular, essas operações equivalem a reverter o cromossomo inteiro ou a uma rotação, e a IR de fechamento deveria receber `x + y_rest`/`y + x_rest` somados de forma consistente.

**Correção:** proibir os casos nos asserts (`assert (chtype == Linear || not (i == 1 && j == n))`, e o análogo para a transposição), ou tratá-los à parte: reversão total = `intergenicFullReversal` com o ajuste da IR de fechamento.

---

## 10. FPT emparelha genes de orientação incompatível [leitura] — ✅ corrigido

> Os pares candidatos (`sample0` e `getNewEdges`) agora são filtrados por `genesMatch`, que exige que os dois genes casem direto ou em reverso pelo matcher, e `findLimits` tem um `assert` com a mesma condição. No exemplo testado (`RFDM`, gene 2 com os dois sinais) a saída não mudou, porque os pares inválidos eram descartados mais adiante, mas eles não chegam mais a virar testemunhas.

`src/PFpt.hs:318-329` (`getNewEdges`/`getVertexInQ`) e `:347` (`findLimits`)

`getVertexInQ` compara `canonicOri`, isto é, ignora o sinal. Para matchers só diretos (`RRDM`/`RFDM`), um par `(a, -a)` não casa, mas `findLimits` marca o par como `w_inverted = True` mesmo assim (porque `dirH0 = False`) e cria uma testemunha de bloco inexistente.

**Correção:** filtrar os candidatos com o matcher:

```haskell
filter (\idx -> let a = subGenome' gene_idx; b = subGenome' idx
                 in isDirectMatch m a b || isReverseMatch m a b) (allGIdx q)
```

E, em `findLimits`, só aceitar `inverted` se `isReverseMatch` valer.

---

## 11. `DB-exe`: reversões/transposições não atingem os últimos genes [leitura] — ✅ corrigido

> `applyOperations` agora sorteia as posições a partir de `size g` (genes reais em `2 .. size g - 1`). Em 2000 instâncias com uma operação e 8 genes, todas as posições foram alteradas por reversões e por transposições (antes, as posições 7 e 8 nunca mudavam).

`app/DB.hs:302-315`

O genoma tem `db_size + 2` genes (os caps ficam nas posições 1 e `db_size + 2`), então os genes reais estão em `2 .. db_size + 1`. Mas:

- reversão: `i ∈ [2, db_size-2]`, `j ∈ [i+1, db_size-1]`;
- transposição: `k ∈ [j+1, db_size]` (e `k` é exclusivo).

Assim, os genes nas posições `db_size` e `db_size + 1` nunca são reordenados, o que enviesa os bancos gerados. Com `db_size < 4`, os intervalos ficam vazios ou invertidos.

**Correção** (melhor ainda usar `size g`, que já considera indels):

```haskell
i_ <- getRandomR (2, size g - 2)
j_ <- getRandomR (i_ + 1, size g - 1)        -- reversão
-- transposição:
i_ <- getRandomR (2, size g - 2)
j_ <- getRandomR (i_ + 1, size g - 1)
k_ <- getRandomR (j_ + 1, size g)
```

---

## 12. A suíte de testes não compila [confirmado] — ✅ corrigido

> `test/PartitionCheck.hs` foi reescrito para a API com `ChromList`: genomas montados com `toMC` e `toMCF`, breakpoints como `Set (CIdx, Idx)` e blocos extraídos com `subChr`. A propriedade de `combine` falha explicitamente se o resultado não for `CGPbal`. Em `prop_joinIsInverseOfCut`, o corte agora é sorteado em `[1, size - 1]`. A suíte compila e as 26 propriedades passam em três rodadas.

`test/PartitionCheck.hs`

`stack build --test` falha em `PartitionCheck.hs:107/111/115`: os algoritmos agora recebem `ChromList`, e os testes passam `GenesIRsR`/`GenesIRsF`. Outros pontos do mesmo arquivo também estão desatualizados:

- `:27-30` usa `EnumSet` para `Breakpoints` (agora é `Set (CIdx, Idx)`) e passa um `GenesIRsR` para `mkPartitionFromBreakpoints`/`mkPartitionFromBlocks`;
- `:40` e `:51`: `commonPrefix` agora devolve uma tripla e `longestSubstring` recebe `ChromList`;
- `:62` e `:84`: as regras subótimas agora trabalham com `ChromList`;
- `:124`: o padrão `(CGPbal pg ph) <-` falha se `combine` devolver `CGPunbal`.

**Correção:** envolver os genomas com `toMC` (e `toMC`/`ChromList` para os flexíveis), trocar `EnumSet` por `Set.fromList [(1, i) | i <- bps]` e ajustar os padrões.

`test/GenomesCheck.hs:198`: `prop_joinIsInverseOfCut` sorteia `i ∈ [1, size]`. Num cromossomo linear, `cut` com `i = size` faz `subGenome (size+1) size`, que lança `indexError`. O intervalo deve ser `[1, size - 1]`.

---

## 13. `Main`: linha em branco derruba o programa [confirmado] — ✅ corrigido

> Linhas em branco (ou só com espaços) e comentários (`#`, mesmo após espaços) são ignorados por `isContent`. O tempo é gravado como `# Time: 0.000384`, sem aspas.

`app/Main.hs:87`

```haskell
filter ((/= '#') . BS.head) . BS.lines
```

`BS.head` em uma linha vazia lança `Data.ByteString.head: empty ByteString`.

**Correção:**

```haskell
filter (\l -> not (BS.all isSpace l) && BS.head l /= '#') . BS.lines
```

Na mesma função (linha 99), `time` já é um `ByteString`, e `BS.pack . show $ time` grava `# Time: "9.5e-4"` com aspas. Use `"# Time: " <> time`.

---

## 14. `writeMultiC` grava os cromossomos em ordem inversa [confirmado] — ✅ corrigido

> Cada cromossomo é escrito antes dos seguintes (`chr_i <> resto`). Ler e escrever de volta reproduz a ordem e o genoma.

`src/Genomes.hs:646-651`

```haskell
writeMultiC writer mc = foldr go ("", "") [1 .. n]
  where go i (bsg, bsi) = (bsg <> chr_i ..., ...)
```

Com `foldr`, o acumulador `bsg` já contém os cromossomos `i+1..n`, e o cromossomo `i` é concatenado **depois** deles. O resultado é `n, n-1, …, 1`:

```
cut (1,2) em [0 1 2 0]  ->  cromossomos [0 1 0], [0 2 0]
writeMultiC  ->  "L 2 ; L 1 ; "
```

O conteúdo do genoma não muda, mas ler e escrever de volta inverte a ordem, e isso confunde a comparação com a entrada.

**Correção:** `(chr_bsg' <> bsg, chr_bsi' <> bsi)`, ou trocar `foldr` por `foldl`.

---

## 15. Problemas menores e pontos a verificar — ✅ corrigido

> Todos os itens foram tratados:
> - `mkRGenome0`/`mkFGenome0` usam `numIRs` (`n` IRs no circular, `n + 1` no linear com caps).
> - `readRGenome`/`readFGenome` lançam `inputError` quando: o número de IRs não bate com `numIRs`; há IR negativa ou intervalo com `l > u`; um token não é inteiro (`readInputInt` no lugar de `fromJust . BS.readInt`); há gene negativo em entrada sem sinal (ver #18). `readMCGenome` exige o mesmo número de cromossomos nas duas linhas. A validação ficou na leitura, e não em `mkRGenome`, porque `intergenicInsertion` recebe de propósito um trecho linear com `#genes + 1` IRs.
> - `cut` usa `localHead`.
> - `Main` lê H como rígido (convertido com `flexibilizeMC 0`) quando a linha de IRs não tem `:`.
> - `DB-exe`: `-r` menor que `-1` dá `inputError`, `-c -r -1` passa a embaralhar, e com `-z` a inserção usa `ir2 = max 0 (ir - ir1)`, sem remover nucleotídeos.
> - A ajuda do `Main` já citava `SOARComb`.
> - Greedy: `commonPrefixC`/`longestSubstringC` indicam quando o bloco é um circular inteiro que casa com um circular inteiro de H (mesmo teste de `isMatch` usado em `Partition`), e nesse caso nenhum breakpoint é adicionado. `commonPrefix`/`longestSubstring` mantêm a interface antiga.
> - FPT: `swapMVar best_sg $!`, `reverse vertices1` em `getOddPaths`, e `undefinedError` explícito para genomas desbalanceados (verificado antes de montar o grafo; o `!sg0` virou `evaluate`).
> - `combine`: nova propriedade `prop_combineKeepsCommonPartitionWithReplicas` (alfabeto de 4 genes, SOAR e Greedy) passa em 100 testes.


- **`Genomes.hs:568,580` (`mkRGenome0`/`mkFGenome0`):** sempre criam `length genes + 1` IRs. Para `Circular` deveriam ser `length genes` IRs, sem caps.
- **`Genomes.hs:590-600` (`readMCGenome`):** não valida a entrada. `zipWith3` trunca em silêncio se o número de `;` diferir entre a linha de genes e a de IRs, e uma contagem de IRs errada só aparece depois, como `index out of bounds`. Sugiro validar em `mkRGenome`/`mkFGenome` (`#IRs == #genes - 1` para linear com caps, `== #genes` para circular) e lançar `inputError`.
- **`Genomes.hs:904` (`cut`):** usa `head . caps` (parcial, com warning `-Wx-partial`). Use `localHead`.
- **Formato de entrada:** `Main` sempre lê H como flexível (`readMCFGenome`). Os arquivos em `datasets/` e a saída do `DB-exe` sem `-f/-F` são rígidos e dão `ERROR: Input is invalid`. Vale aceitar H rígido (convertendo com `flexibilize 0`) ou documentar a exigência.
- **`app/DB.hs:226` vs `:231`:** `g_final` usa `db_nop > 0` e `h` usa `db_nop /= 0`. Com `-c -r -1`, `nop_target = 0` e o embaralhamento pedido é ignorado. Com `-r` negativo diferente de `-1`, `applyDCJs` recebe `nop < 0` e recursa para sempre.
- **`app/DB.hs:271`:** com `-z`, `ir2 = 0` mesmo quando a IR original é maior que zero, então a inserção remove nucleotídeos.
- **`app/Main.hs:48`:** a ajuda não menciona o algoritmo `SOARComb`.
- **`PGreedy.hs:126-135`:** `commonPrefix'` para em `k < size`, e o bloco sempre recebe um breakpoint (`prevChromIdx g_beg == g_end`). Por isso o Greedy nunca produz um cromossomo circular inteiro sem breakpoint, mesmo quando H tem o mesmo circular (Approx/SOAR produzem). O número de blocos é o mesmo, mas a contagem de breakpoints não.
- **`PFpt.hs:163-170` (`selectSG`):** `swapMVar best_sg (if cost sg < ... )` guarda um thunk que referencia o melhor anterior. Use `swapMVar best_sg $! ...` para evitar o acúmulo.
- **`PFpt.hs:231`:** `vertices1 ++ (v : vertices2)` não inverte `vertices1`, então o caminho fica fora de ordem. Hoje não importa (só o conjunto é usado), mas `branch2` (comentado) usa `evens` e depende da ordem.
- **`PFpt` com genomas desbalanceados:** falha no `assert (checkDegrees 1 2 …)` (`PFpt.hs:96`). É um TODO conhecido, mas o erro poderia ser explícito.
- **`Partition.combine` (a verificar):** juntar dois blocos adjacentes em G e dois em H que casam não garante que o restante continue com emparelhamento perfeito quando há genes repetidos. Nesse caso `mkCommonPartition3` termina em `error logicError`. Não reproduzi, mas vale um teste de propriedade (`combine` preserva `checkCommonBal`).

---

## 16. `randomGenomeWithReplicas` gera genomas inválidos [confirmado] — ✅ corrigido

> `randomGenomeWithReplicas` lança `inputError` quando `low > high` ou `rep * high > n` (pior caso dos genes repetidos). `-r 4 -h 4 -n 14` agora para com `ERROR: Input is invalid.`, e `-r 4 -h 3 -n 14` gera instâncias com 14 genes e 15 IRs.

`src/Genomes.hs:1023-1036` (usado por `DB-exe RepDB`)

Os genes repetidos são `concat (zipWith replicate occs [1 .. rep])`, com `occs ∈ [low, high]`. Depois vêm os singletons `[length l + 1 .. n]`. Se `length l > n` (por exemplo `-r 4 -h 4 -n 14`, que dá até 16 genes repetidos), a lista de singletons fica vazia e o genoma sai com `length l > n` genes, mas só `n + 1` IRs. Todos os algoritmos quebram ao ler essa instância (`invalid slice`).

```
$ DB-exe RepDB -l 2 -h 4 -r 4 -k 40 -n 14 -s -f 10 -o rep
instância 11: 15 genes e 15 IRs  ->  MCFISP-exe: invalid slice (12,4,15)
```

**Correção:** validar os parâmetros (`rep * high <= n`) e lançar `inputError`, ou cortar a lista em `n` genes (`take n`) antes de criar o genoma.

---

## 17. `shuffleGenome` embaralha os caps [confirmado] — ✅ corrigido

`src/Genomes.hs` (`shuffleGenome`, usado por `DB-exe -r -1`)

O vetor embaralhado incluía os caps (`0` e `maxBound`), e depois `mkRGenome` ainda adicionava caps novos. H saía com dois genes a mais e caps no meio do cromossomo (por exemplo `-9223372036854775807` e `0`).

> Agora só os genes internos são embaralhados, os caps originais são mantidos e as IRs são sorteadas com a mesma quantidade e a mesma soma. Em circulares, todos os genes são embaralhados.

---

## 18. Instância com sinal lida sem `-s` [confirmado] — ✅ corrigido (validação de entrada)

As instâncias de `DB-exe RepDB … -s` funcionam em todos os algoritmos quando lidas com `-s`. Lidas sem `-s`, todos (SOAR, Greedy, Approx) quebravam com `logicError` em `mkCommonPartition`, o que parecia um erro do algoritmo.

> `readRGenome`/`readFGenome` lançam `inputError` quando a entrada sem sinal tem genes negativos.

---

# Segunda revisão

Revisão completa feita depois das correções de #1–#18. Quatro frentes em paralelo (Genomes/LocalBase; Partition/Balancer; SOAR/Greedy/Approx; FPT/Main/DB/testes). Cerca de 5 mil instâncias, de bancos do `DB-exe` (todas as combinações de opções) e de um gerador próprio com circulares e alfabetos de 1–4 genes, foram rodadas em todos os algoritmos. As saídas foram conferidas por um validador independente, que remonta os genomas a partir dos blocos e refaz o emparelhamento.

**Bugs antigos:** #1–#18 continuam corrigidos, sem regressão. A suspeita sobre `combine` (seção 15) foi confirmada e virou o #21.

## 19. SOAR: conflito de duos que fecha um circular só em H [confirmado] — ✅ corrigido

> `isConflict` agora também verifica a consistência a partir de H (`followsH` + `consistentG`) e marca como conflito dois duos que cobrem um circular inteiro de 2 genes em só um dos genomas. Os 6 casos conhecidos dão partições válidas.

`src/PSOAR.hs:68-83` (`isConflict`)

A consistência entre dois duos escolhidos só é verificada quando eles se seguem em G (`followsG`). Dois duos podem formar uma cadeia em G e, ao mesmo tempo, fechar um circular em H. Aí H fica com um bloco que é um circular inteiro contra um bloco linear em G, e `mkCommonPartition` falha. Também explica o crash do SOAR que apareceu nas outras frentes (instâncias multicromossomo, com e sem repetição).

```
C -4 -3 4 -3 ;           (-s)
0 0 2 0 ;
C -3 ; C -4 ; C 3 -4 ;
0:0 ; 0:0 ; 2:2 0:0 ;
MCFISP-exe: ERROR: There is something wrong ...   (Greedy/Approx resolvem)
```

Outros casos: `L -6 4 -6 ; 318 97 31 150 ;` contra `L 6 ; C -4 6 ; 270:398 80:114 ; 78:116 24:34 ;` (com `-s`), e `L 1 2 2 8 ; C 2 2` contra `L 1 9 8 2 2 2 2` (sem sinal).

**Correção sugerida** (testada numa cópia: os 5 casos distintos passam, e mil instâncias novas rodaram sem falha):

```haskell
    || followsH idx2 idx2' && not (consistentG due1 due2 idx1 idx1')
    || followsH idx2' idx2 && not (consistentG due1' due2' idx1' idx1)
    || (followsG idx1 idx1' && followsG idx1' idx1) /= (followsH idx2 idx2' && followsH idx2' idx2)
  where
    consistentG d1 d2 i i' =
      followsG i i' && isDirectMatch matcher d1 d2
        || followsG i' i && isReverseMatch matcher d1 d2
```

## 20. FPT: circular de 1 gene com IR diferente [confirmado] — ✅ corrigido

> `closesCircle` exige que a testemunha cubra os dois circulares e compara os cromossomos inteiros com `isMatch`. O exemplo dá `C 6 ; | 27`, como nos outros algoritmos.

`src/PFpt.hs:209-216` (`closesCircle`)

Num circular de 1 gene, `gr == gl`, e `subGenome gr gl` é só o gene: a IR de fechamento nunca é comparada, `closesCircle` dá True e nenhum breakpoint é criado.

```
C 6 ;  27 ;  C 6 ;  39:57 ;
FPT: logicError        Greedy/SOAR/Approx: C 6 ; | 27
```

**Correção:** comparar os cromossomos inteiros, como o Greedy faz agora: `isCircular chrG && isCircular chrH && size chrG == size chrH && isMatch sg_matcher chrG chrH`.

## 21. `combine` (SOARComb) quebra a partição comum [confirmado] — ✅ corrigido

> Uma fusão (em `go1` e em `wrapMerge`) só é aceita se a partição continuar comum e do mesmo tipo da original (`checkCommonBal` ou `checkCommonUnbal`, em `isCommon`). Os dois exemplos passam.

`src/Partition.hs:347-430` (`go1`/`wrapMerge`)

Um par adjacente em G é fundido com um par adjacente em H sempre que os blocos fundidos casam. Só que ninguém verifica se os blocos originais estavam emparelhados entre si. A fusão pode então tomar o único parceiro de outro bloco, e `mkCommonPartition3` termina em `logicError`. Em instâncias geradas com circulares, falhou em 39 (14 delas balanceadas, inclusive fusões pela origem). Nos bancos desbalanceados do `DB-exe`, falhou em 2.

```
L 1 2 1 1 ;            (sem sinal; SOAR dá partição válida)
0 0 0 0 0 ;
L 3 1 1 1 1 ;
0:0 0:0 0:0 0:0 0:0 0:0 ;
SOARComb: logicError
```

A propriedade `prop_combineKeepsCommonPartitionWithReplicas` não pega o caso porque só gera um linear.

**Correção sugerida:** manter o emparelhamento atual (`getBlocksFit`) e só fundir `(ba1, bb1)` com `(ba2, bb2)` se `ba1↔ba2` e `bb1↔bb2` (cruzado no casamento reverso) forem pares dele. Uma alternativa mais simples, porém mais lenta, é aceitar a fusão só se `checkCommonBal || checkCommonUnbal` continuar valendo depois dela.

## 22. FPT não é exato quando um gene casa nos dois sentidos [confirmado] — ✅ corrigido

> `findAllLimits` devolve uma testemunha para cada sentido em que o par casa, e `edgesFromRareOddPath` ramifica nas duas. `findLimits` (usado na amostra inicial) continua preferindo a direta: com sinal, só o cap 0 casa nos dois sentidos, e ele nunca é singleton. O exemplo agora dá `L 0 2 | 2 | 0`. Em 1080 execuções de `DB-exe` (balanceadas, com e sem `-c`, com repetições) e na nova propriedade, o FPT nunca foi pior que outro algoritmo.

`src/PFpt.hs:397-405` (`findLimits`)

Os caps são todos `0` e `0 = -0`, então `isDirectMatch` vale, `dirH0 = True` e a testemunha invertida nunca é construída. O mesmo acontece com qualquer gene que casa nos dois sentidos (o cabeçalho do módulo já diz que reverso sem sinal não é suportado).

```
L 2 2 ;  51 309 193 ;  L 2 -2 ;  157:229 246:360 46:68 ;    (-s)
FPT:        L 0 | 2 | 2 | 0     (3 breakpoints)
SOAR/Approx: L 0 2 | 2 | 0      (G "0 2" com H "-2 0" invertido, IR 51 ∈ [46,68])
```

**Correção:** quando o par casa nos dois sentidos, ramificar nas duas testemunhas (ou fazer `edgesFromRareOddPath`/`createSampleGraph` gerarem as duas orientações).

## 23. `Main`: uma instância inválida derruba todas [confirmado] — ✅ corrigido

> Cada instância roda dentro de um `try`. Um erro é escrito no arquivo da instância (`# Error: <msg>` e `# Time`) e em stderr, as demais instâncias continuam, e o programa termina com código 1 se alguma falhou. Isso vale nos modos paralelo e `--no-par`.

`app/Main.hs` (`runOne`)

Num arquivo de 4 instâncias com a 2ª inválida, o modo paralelo não escreve **nenhum** arquivo, e `--no-par` escreve só `_0001`. Arquivos de execuções anteriores também não são apagados, então resultados velhos podem passar por novos.

**Correção:** envolver cada instância em `try` (`ErrorCall`) e escrever `# Error: <msg>` no arquivo dela.

## 24. Testes: `stack test` nunca falha [confirmado] — ✅ corrigido

> - `Spec.hs` termina com `exitFailure` quando alguma propriedade falha, e inclui `PFptCheck`.
> - Novo gerador `genMultiPair` (`GenomesCheck.hs`): de 1 a 3 cromossomos lineares ou circulares, alfabeto de 1 a 6 genes, com ou sem sinal, H obtido de G por DCJs aleatórios e flexibilizado, opcionalmente desbalanceado (genes exclusivos em G e em H).
> - Novas propriedades `prop_*MultiPartitionIsValid` (Greedy, GreedySin, SOAR, SOARComb, Approx; 300 testes cada): a partição é comum, preserva os genomas de entrada e é balanceada quando a entrada é.
> - `prop_fptIsValidAndNotWorseThanOthers`: FPT em genomas balanceados pequenos termina, é válido e nunca é pior que Greedy, SOAR e Approx.
> - Antes das correções, essas propriedades pegavam #19, #21 e #22. Agora, 33 propriedades passam em 6 execuções seguidas.

`test/Spec.hs`

`main` só faz `print ans`, e o código de saída é sempre 0. **Correção:** `unless ans exitFailure`.

Lacunas de cobertura que esconderam #19–#22:
- as propriedades de partição usam um único linear, com sinal e só reversões (nada de circular, multicromossomo, DCJ, transposição ou desbalanceado);
- o FPT não tem teste: a propriedade está comentada, e `PFptCheck.hs` está vazio e fora do `Spec.hs`.

## 25. Regras subótimas do SOAR desativadas [leitura] — a confirmar (é isso mesmo)

`src/PSOAR.hs:36-39`

Na árvore de trabalho a linha ativa é `(g, h) = (g_, h_)`, enquanto no HEAD as regras `suboptimalRuleInterval`/`suboptimalRulePairs` estão ligadas. Se foi de propósito (por exemplo, porque as regras ainda não suportam multicromossomo), vale deixar isso explícito no comentário. Uma cópia com as regras ligadas rodou 480 instâncias do `DB-exe` sem falha, e a saída manteve os genes originais (#5 continua valendo).

## 26. Problemas menores

- **Gene `0` na entrada** (`Genomes.hs`, `readMCGenome`): é aceito e confundido com os caps. `checkCommonUnbal` ignora esse gene, e então G `L 1 0 2` contra H `L 1 2` sai "válido". **Correção:** em `readRGenome`/`readFGenome`, lançar `inputError` se `0 ∈ genes`.
- **CRLF e tabulação** (`splitElements`/`readInputInt`): a leitura nova rejeita arquivos com `\r\n` e tabs (antes, CRLF sem `;` no fim funcionava). **Correção:** separar por `isSpace` e tratar como vazio o segmento que for só espaço.
- **Genoma vazio** (linhas só com `;`): é aceito e gera partição vazia marcada como `# Exact solution`. **Correção:** `inputError` em `readMCGenome` se não houver cromossomo.
- **`randomGenomeWithReplicas` com `-l 0`** (`Genomes.hs:1048`): um gene repetido pode receber 0 ocorrências, e aí os singletons `[length l + 1 ..]` reaproveitam rótulos de repetidos (`-l 0 -h 1` gera um gene duas vezes). **Correção:** numerar os singletons a partir de `rep + 1` (`take (n - length l)`) e/ou rejeitar `low < 1`.
- **`DB-exe`** (`app/DB.hs`): `-p` fora de [0,100] e `-d` negativo entram em loop infinito, e `-n 1 -r 1` falha no `assert` (Genomes.hs:700). **Correção:** validar com `inputError`, junto do teste de `-r`.
- **SOAR em circular de 1 gene** (`PSOAR.hs:101-104`): `getDuo` usa `subGenomeBySize i 2`, então a adjacência (1,1) nunca é escolhida, e o SOAR sempre corta `C 1` mesmo quando H tem o mesmo circular. (✅ no SOARComb: o `combine` agora fecha o circular que ficou com um único bloco, ver o item 11 de `CompatibilityWithText.md`; o `-a SOAR` sem `combine` continua cortando.) Também corta circulares e regiões exclusivas sem necessidade (`8 | 9 | 7` contra `8 9 7` do Greedy). O resultado é válido, só subótimo.
- ✅ **Greedy em empate** (`PGreedy.hs`, `longestSubstringC`): `maxWith` comparava só `d`, e com repetições podia escolher o casamento que não fecha o circular. Corrigido com `maxWith (\(_, _, d, c) -> (d, c))` nos dois lugares; ver o item 8 de `CompatibilityWithText.md` (o exemplo `C 1 ; C 1` não reproduzia; o de lá reproduz).
- **Caps em excesso no Balancer** (`Balancer.hs:134-135`): os caps são tratados como genes comuns, e quando G e H têm números diferentes de lineares, a escolha "mais à esquerda" pode tirar os dois caps do mesmo cromossomo. Exemplo: Approx com 5 breakpoints contra 3 do SOAR/Greedy em G `L 1 2 ; L 3` × H `L 1 2 3`. É válido, só subótimo.
- **`ApproxLB` em instâncias desbalanceadas** (`Balancer.hs:169-184`): não é um limite inferior de verdade (ficou acima de uma solução válida em dezenas de instâncias), como o próprio comentário do código admite. A ajuda do `Main` diz "lower bound" sem essa ressalva. Em instâncias balanceadas não houve nenhuma violação.
- ✅ **Custo** (`Partition.hs:190-198`): `costBalanced`/`costUnbalanced` contavam blocos, mas os comentários diziam breakpoints. Com isso o FPT não diferenciava um circular inteiro de um circular cortado uma vez. Agora os dois contam os breakpoints de G e de H (ver o item 1 de `CompatibilityWithText.md`).
- **Tempo do FPT** (`app/Main.hs`): fixo em 600 s e não configurável, então o caminho `# Partial solution` não é fácil de testar.

---

## 27. FPT: `cleanEdges` passa do limite da testemunha num circular [confirmado] — ✅ corrigido

`src/PFpt.hs` (`cleanEdges`)

Ao remover as arestas de uma testemunha a partir de um vértice, `cleanEdges` anda em cada genoma afastando-se do par preto. A caminhada só parava antes do par preto ou num vértice sem aresta da testemunha, e não no limite da testemunha. Quando a testemunha cobre um cromossomo circular inteiro, passar do limite leva ao outro extremo da mesma testemunha, cujas arestas também eram apagadas. Isso acontecia só em um dos genomas, então o grafo de amostra ficava assimétrico e `finalizeCorrespondence` terminava em `patternError`.

```
C 1 1 1 1 ;  0 0 1 0 ;  C 1 ; C 1 1 1 ;  0:0 ; 1:1 0:0 0:0 ;     (-s)
FPT: ERROR: Pattern shouldn't be possible.   (Greedy/SOARComb/Approx resolvem)
```

O bug apareceu quando a poda do FPT passou a usar breakpoints (item 1 de `CompatibilityWithText.md`). Com circulares essa poda é mais fraca, e o FPT passou a avaliar folhas que antes eram descartadas.

> **Correção:** `move` para no limite da testemunha na direção da caminhada (o par direito ou o esquerdo, trocados no lado de H quando a testemunha é invertida). O exemplo dá custo 4, e a propriedade do FPT passa em 5 execuções seguidas.


## 28. FPT sem sinal: singletons fixados no sentido direto [confirmado] — ✅ corrigido

`src/PFpt.hs` (`fptPartition`, `sample0`)

A amostra inicial do FPT contém os genes que aparecem uma só vez em cada genoma (singletons), e a testemunha de cada um era construída com `findLimits`, que prefere o sentido direto. Com sinal, só os telômeros casam nos dois sentidos, e eles nunca são singletons; por isso a correção do #22 (ramificar nos dois sentidos durante a busca) bastava. Sem sinal, **todo** gene casa nos dois sentidos, e os singletons entravam na amostra com testemunhas diretas, o que impedia blocos reversos.

```
L 1 2 3 ;  1 2 3 4 ;  L 3 2 1 ;  4:4 3:3 2:2 1:1 ;    (sem -s)
FPT antes:  L 0 | 1 | 2 | 3 | 0   (8 breakpoints)
FPT agora:  L 0 1 2 3 0           (0, casamento reverso inteiro)
```

**Correção:** os singletons que casam nos dois sentidos ficam fora da amostra inicial (`genesMatchOneWay`). O vértice isolado de cada um é um caminho ímpar raro, e a regra de ramificação 1 o acrescenta com as testemunhas dos dois sentidos (`findAllLimits`). Com sinal, nenhum singleton casa nos dois sentidos, então nada muda: em 200 instâncias com sinal, o custo do FPT ficou igual ao anterior. Em 400 instâncias aleatórias sem sinal (lineares e circulares), o FPT nunca foi pior que o SOAR, o Approx ou o Guloso. `stack test` passou. O cabeçalho do módulo, que dizia que o caso sem sinal não era suportado, foi atualizado.

## 29. FPT: só uma janela máxima por testemunha num circular [confirmado] — ✅ corrigido

Uma testemunha estende o par de genes para a direita e para a esquerda enquanto os duos casam, mas não pode cobrir mais que um circular inteiro (`max_steps = min |G_c| |H_c| − 1`). `findLimitsDir` estendia primeiro para a direita até esse limite e só usava o que sobrava para a esquerda. Quando o circular inteiro casa, existem várias janelas máximas (rotações) que contêm o par, e só uma era gerada. O ótimo pode precisar de outra, ou seja, do circular cortado em outra posição.

```
-s
G: L 1 2 1 -1 -2 ; C -2 2 -1     (IRs 0)
H: L 1 2 1 -2 2 1 1 2            (IRs 0)
FPT antes:  custo 10, marcado como exato
FPT agora:  custo 8   (G L 0 1 2 1 | -1 | -2 | 0 ; C -1 -2 2 cortado no fechamento,
                       com -1 -2 2 casando com -2 2 1 de H invertido)
```

O bloco ótimo `-1 -2 2` (posições 3, 1, 2 do circular) só sai do par (G3, H6). A partir do par (G1, H5), por exemplo, a extensão para a direita gastava os dois passos e dava `-2 2 -1` ↔ `1 -2 2`.

**Correção:** `findLimitsDir` calcula a extensão máxima de cada lado (cada uma limitada por `max_steps`). Se a soma passa de `max_steps`, devolve todas as janelas com `max_steps` passos; senão, a única janela. `findAllLimits` junta as janelas dos dois sentidos, e a regra de ramificação 1 tenta todas. Na amostra inicial, um singleton com mais de uma testemunha fica de fora (`uniqueWitness`), como os que casam nos dois sentidos (#28), e é tratado pela ramificação.

**Verificação:** uma força bruta independente (todos os conjuntos de breakpoints, com circulares inteiros e circulares cortados uma vez) foi comparada com o FPT em 1500 instâncias pequenas com circulares (1000 com sinal, 500 sem sinal). O FPT novo deu o ótimo em todas; o antigo errava em 7. Em outras 800 instâncias com circulares, o FPT novo nunca foi pior que Guloso, SOARComb e Approx, e o tempo total não mudou (34 s nas 400 de `s3`). `stack test` passou.


## 30. FPT para de calcular sem terminar algumas instâncias [observado uma vez] — em aberto

Na calibração do experimento do FPT (2026-11-24), `MCFISP-exe -a FPT -s -t 60 -i N50_SREV25_I0_L50_F0.txt` (8 instâncias geradas com `gen_db_fpt.sh`, rodando ao mesmo tempo que outro processo igual) gravou 3 saídas em cerca de 3 s e ficou mais de 80 minutos sem gravar as outras 5, dormindo (55 s de CPU no total). O mesmo arquivo terminou em 60 s em todas as repetições: sozinho, com `--no-par`, com `-t 10` e com quatro processos simultâneos. A causa não foi encontrada. Suspeitas: a interação do `race` (tempo limite) com o pool do `parallel-io` ou o cancelamento de uma thread que não aloca.

**Contorno:** o `run_fpt.sh` limita cada chamada do `MCFISP-exe` com `timeout`; um arquivo incompleto fica na pasta temporária e é refeito na próxima execução do script.
