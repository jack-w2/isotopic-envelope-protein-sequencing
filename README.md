# Isotopic-envelope protein sequencing

Companion source code for the MSc thesis

> Jacek Wolf, *Protein sequencing in mass spectrometry via detection of isotopic envelopes*, University of Warsaw, Faculty of Mathematics, Informatics and Mechanics, 2026. Supervisor: dr Michał Startek.

This is a research prototype used to produce the results in the thesis, not a supported sequencing tool. The search starts from the empty sequence and scores candidate one-residue extensions by Masserstein deconvolution of IsoSpec isotopic envelopes.

Repository: https://github.com/jack-w2/isotopic-envelope-protein-sequencing

## Requirements

Python 3.10 or newer. From the repository root:

```bash
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

The pipeline depends on [IsoSpec](https://github.com/lcbm-group/IsoSpec) (`IsoSpecPy`) and [Masserstein](https://github.com/mciach/masserstein).

## Reproducing a single run

```bash
python create_spectrum.py -f MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGLVEQCCTSLCSLYQLENYCN -l insulin.csv
```

If `-f` is omitted, the default sequence in `config.toml` is used (human insulin). `-T` sets the iteration budget (default `110`, as in the autonomous evaluation). `--teacher-force` switches to the calibration regime, which forces the true prefix at every step.

## Batch runs

Autonomous evaluation (`T = 110`, free search):

```bash
python multiple_runs_launcher.py --outdir tests/logs_run -T 110
```

Calibration dataset (`T = 26`, teacher-forced):

```bash
python multiple_runs_launcher.py --outdir tests/logs_run_25 -T 26 --teacher-force
```

`tests/proteins.fasta` holds the 100-protein UniProt set used in the thesis. Isoleucine is replaced by leucine on load. MOTS-c (`A0A0C5B5G6`) is in the FASTA but produced no log in the original run and is excluded from the reported statistics.

## Datasets committed with the thesis

| Path | Role in the thesis |
| --- | --- |
| `tests/logs/` | Autonomous evaluation corpus (`T = 110`), 99 proteins |
| `tests/logs_25/` | Teacher-forced calibration corpus (26 steps), 99 proteins |
| `tests/ablation/` | Seed-matched ablation, first noise draw (`T = 30`) |
| `tests/ablation_draws/draw2/` | Second noise draw of the same ablation |
| `tests/proteins.fasta` | Target sequences |

Do not overwrite these directories if you want the logged numbers to stay identical to the thesis. New runs should go to a separate `--outdir`.

## Other scripts

- `linear_model.py` — ordinary-least-squares fit of the scoring constants \(c_1,\ldots,c_4\) from `tests/logs_25/`
- `cross_validation.py` — algorithm-vs-random performance summary on `tests/logs/`
- `ablation_experiment.py` / `ablation_launcher.py` / `ablation_launcher_extra_draws.py` — supplementary ablation of \(\alpha\) and of the isobaric-resolution mechanism

## Layout

| File | Role |
| --- | --- |
| `create_spectrum.py` | Driver: spectrum simulation, scoring, search, CSV logging |
| `sequence.py` | Sequence object and conversion to a molecular formula |
| `molecular_formula.py` | C/H/N/O/S arithmetic |
| `analyse_spectrum.py` | Wrapper around `masserstein.estimate_proportions` |
| `priority_queue.py` | Best-first queue with a visited set keyed on formula |
| `config.toml` | Default target, intensity profile, noise, calibration constants |
| `amino_acids.csv` | One-letter code to molecular formula |
| `amino_acids_replacements.csv` | Isobaric pairs N↔GG, Q↔AG |

## Licence

MIT. See `LICENSE`.
