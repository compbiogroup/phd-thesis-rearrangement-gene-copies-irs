# Base de dados do experimento do FPT

Instâncias menores, geradas com `scripts/gen_db_fpt.sh` para avaliar o algoritmo FPT de partição de genomas em instâncias que ele resolve de forma exata (`scripts/run_fpt_small.sh`). O gerador é o mesmo da base `db/` (`DB-exe` de `partition/`), com os parâmetros proporcionais ao número de genes n.

A organização dos arquivos, o formato e o significado dos parâmetros são os mesmos de `db/` (veja `db/README.md`):

```
db_fpt/N<n>/<mod>/N<n>_<mod><op>_I<indels>_L<alfabeto>_F0.txt
```

| parâmetro | valores |
|---|---|
| `n` | 25, 50 |
| `mod` | `SDCJ`, `SREV`, `SREV_TRANS`, `TRANS` |
| `op` | n/4 e n/2 (6 e 12 para n = 25, 12 e 25 para n = 50) |
| `indels` | 0 e n/4 (0 e 6 para n = 25, 0 e 12 para n = 50) |
| `alfabeto` | n/2 e n (12 e 25 para n = 25, 25 e 50 para n = 50) |
| `flex` | 0 (relação de igualdade, para a qual o FPT é exato) |

São 64 arquivos com 100 instâncias cada. Como a geração é aleatória, os arquivos estão no repositório para que os resultados possam ser reproduzidos.
