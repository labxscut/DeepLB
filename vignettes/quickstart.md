# Quickstart Vignette

This vignette demonstrates compatibility-safe ways to launch existing DeepLB
workflows.

## 1) Existing shell command (unchanged)
```bash
bash Scripts/DeepLB_pipeline.sh -r /path/to/DeepLB -t lihc -g TH -u part3 -q 0.1 -k hyper -v 1
```

## 2) Python module wrapper
```bash
python -m deeplb -- -r /path/to/DeepLB -t lihc -g TH -u part3 -q 0.1 -k hyper -v 1
```

## 3) Dry-run with wrapper
```bash
python -m deeplb --dry-run -- -r /path/to/DeepLB -t lihc -g TH -u part3 -q 0.1 -k hyper -v 1
```

All commands above still execute through `Scripts/DeepLB_pipeline.sh`.
