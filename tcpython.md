# tcpython Notes

This document captures how Thermo-Calc (`tc_python`) was run for this repo.

## Environment
- Conda env: `tcpython`
- Thermo-Calc install: `TC25B_HOME=/home/vela/thermocalc/2025b`
- Repo root on `PYTHONPATH` so `strength_model.py` is importable.

## Command Used
From the repo root:

```bash
conda run -n tcpython --cwd /home/vela/projects/ult/redesign_2_2_26 \
  python analysis/run_tc_properties.py
```

If `strength_model` is not found, add:

```bash
export PYTHONPATH=/home/vela/projects/ult/redesign_2_2_26
```

## Notes
- In the Codex sandbox, Thermo-Calc startup hung. Running outside the sandbox was required.
- The script writes `trajectory_search/properties_list.csv`.
- The script writes `trajectory_search_grid2/properties_list.csv`.
