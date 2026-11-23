# Códigos da tese

Códigos usados nos experimentos da tese, separados do repositório de trabalho do doutorado.

## Estrutura

| pasta | conteúdo |
|---|---|
| `cyc/` | heurísticas de empacotamento do grafo de adjacências por ciclos (`dec`), distâncias de rearranjo (`dist`) e conversão para permutações (`toperm`); veja `cyc/README.md` |
| `scripts/` | scripts que rodam os experimentos da tese |

## Compilação

```sh
cd cyc
make            # gera cyc/dec, cyc/dist, cyc/part e cyc/toperm
```

Requer `cmake`, `g++` com OpenMP e as dependências listadas em `cyc/shell.nix` (com Nix: `nix-shell cyc/shell.nix`).

## Experimentos de empacotamento e distância

| script | função |
|---|---|
| `scripts/run_dec_tese.sh` | roda `dec` (heurísticas `ga` e `rand`, decomposições `bfs`, `map` ou `1c2c3c`) na base da tese |
| `scripts/run_dist_tese.sh` | roda `dist` na base da tese, com mapeamentos aleatórios (`map`) ou com os empacotamentos de `run_dec_tese.sh` (`dec`) |

Os scripts entram em `cyc/` e usam caminhos relativos a essa pasta; as pastas de saída padrão (`ANS_TESE/`, `TESE_OUT/`) ficam na raiz deste repositório e não são versionadas. As variáveis de ambiente aceitas (pastas de entrada e saída, número de iterações, valores dos parâmetros das instâncias) estão descritas no cabeçalho de cada script.

## Partição de genomas

Os códigos da parte de partição de genomas ainda serão adicionados a este repositório. As opções `PART` e `FIX` de `run_dec_tese.sh` e `run_dist_tese.sh`, as opções `-f`/`-p` de `dec` e o executável `cyc/part` dependem dessa parte.
