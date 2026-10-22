# RNA paired-end workflow: DeepImmuno and EasyFuse options

## Check DeepImmuno from a run config

Given a TOML configuration such as:

```toml
[tools.executables]
deepimmuno_dir = "/opt/neoag/tools/DeepImmuno"
deepimmuno_python = "/opt/neoag/envs/deepimmuno/bin/python"
```

run:

```bash
neoag check-tools --config /path/to/run.toml
```

Inspect the `deepimmuno` row. The availability check resolves `deepimmuno_dir` to
`deepimmuno-cnn.py`. The same directory and `deepimmuno_python` are used by
the execution runner.

## Run expression and splice analysis without EasyFuse

Set the following in the RNA sample manifest:

```yaml
inputs:
  easyfuse_enabled: false
```

The generated RNA profile omits EasyFuse discovery, fusion union, fusion
cross-validation, and fusion peptide generation. RNA QC, STAR alignment,
expression quantification, junction extraction, SNAF, SpliceMutr, splice
normalization, and downstream splice-peptide analysis remain available.
Fusion evidence is unassessed in this mode; it is not a negative fusion result.

## Reuse an existing EasyFuse result

To bypass discovery and the EasyFuse reference while retaining fusion
downstream analysis:

```yaml
inputs:
  easyfuse_pass_csv: /path/to/SAMPLE/fusions.pass.csv
  easyfuse_all_csv: /path/to/SAMPLE/fusions.csv # optional
```

`fusions.pass.csv` is required for this reuse mode. `fusions.csv` is optional;
when supplied it is used to rescue high-confidence internal-caller rows that
were not classified as EasyFuse PASS. The workflow records that no discovery
stage was run and starts at candidate union.

## Stage the official EasyFuse reference

The official EasyFuse documentation currently publishes the reference as an
FTP archive. The server reported a compressed size of 80,190,986,066 bytes on
2026-10-22. A resumable site path is:

```text
/mnt/neoag_100T/liup/neoag_refs/easyfuse/easyfuse_ref_v4.tar.gz
```

Download directly on the compute server:

```bash
aria2c -c -x 8 -s 8 --file-allocation=none \
  -d /mnt/neoag_100T/liup/neoag_refs/easyfuse \
  ftp://easyfuse.tron-mainz.de/easyfuse_ref_v4.tar.gz
```

Or transfer an archive downloaded on another machine:

```bash
rsync -ah --partial --append-verify --info=progress2 \
  easyfuse_ref_v4.tar.gz \
  na@10.200.50.134:/mnt/neoag_100T/liup/neoag_refs/easyfuse/
```

Extract and configure it:

```bash
cd /mnt/neoag_100T/liup/neoag_refs/easyfuse
pigz -dc easyfuse_ref_v4.tar.gz | tar -xf -
export NEOAG_EASYFUSE_REF=/mnt/neoag_100T/liup/neoag_refs/easyfuse/easyfuse_ref_v4
test -f "$NEOAG_EASYFUSE_REF/BEFORE_EXECUTING_EASYFUSE"
```

Keep the archive until extraction and the EasyFuse preflight checks succeed.
Do not interpret a failed or skipped reference download as biological evidence.
