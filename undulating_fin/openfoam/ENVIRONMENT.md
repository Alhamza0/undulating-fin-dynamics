# OpenFOAM environment

Reproducible setup for the CFD half of the project, as installed and
verified on 2026-08-31.

| Item | Value |
|---|---|
| Host | Windows 11, WSL 2 |
| Distro | `Ubuntu-24.04` (noble) |
| OpenFOAM | v2412 (ESI / openfoam.com), `openfoam2412-default` |
| Run user | `foam` (uid 1000) — **not root**, see below |
| Resources | 14 cores, 7 GB RAM, 955 GB disk |

## Setup

```bash
# 1. Distro (--no-launch avoids the interactive first-run prompt)
wsl --install Ubuntu-24.04 --no-launch

# 2. OpenFOAM from the official ESI apt repository
wsl -d Ubuntu-24.04 --user root -- bash -c '
  apt-get update && apt-get install -y curl ca-certificates gnupg
  curl -sS https://dl.openfoam.com/add-debian-repo.sh | bash
  apt-get install -y openfoam2412-default'

# 3. A non-root user — REQUIRED, see below
wsl -d Ubuntu-24.04 --user root -- bash -c '
  useradd -m -s /bin/bash foam && usermod -aG sudo foam
  echo "foam ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/foam
  chmod 440 /etc/sudoers.d/foam
  printf "[user]\ndefault=foam\n\n[interop]\nappendWindowsPath=false\n" > /etc/wsl.conf
  echo "source /usr/lib/openfoam/openfoam2412/etc/bashrc" >> /home/foam/.bashrc
  chown foam:foam /home/foam/.bashrc'

wsl --terminate Ubuntu-24.04     # pick up /etc/wsl.conf
```

Run anything with a **login** shell so `.bashrc` sources the OpenFOAM
environment:

```bash
wsl -d Ubuntu-24.04 -- bash -lc 'blockMesh ...'
```

## ⚠️ Coded boundary conditions will not run as root

The fin motion is prescribed through a `codedFixedValue` boundary
condition on `pointDisplacement`, which OpenFOAM compiles into a shared
library at run time and loads with `dlopen`. It refuses to do this for a
privileged user:

```
--> FOAM FATAL IO ERROR:
This code should not be executed by someone with administrator rights
for security reasons. It generates a shared library which is loaded
using dlopen
    From static void Foam::dynamicCode::checkSecurity(...)
```

A distro installed with `wsl --install --no-launch` has **no user
account and defaults to root**, so this hits immediately. Hence the
`foam` user in step 3. Symptom if it regresses: solver exits 1 before
the first time step and `dynamicCode/` is empty.

## Verification performed

| Check | Result |
|---|---|
| `blockMesh`, `snappyHexMesh`, `pimpleFoam`, `checkMesh`, `surfaceFeatureExtract`, `decomposePar` present | pass |
| Moving-mesh solve (`pimpleFoam` on the `movingCone` tutorial, `dynamicMotionSolverFvMesh`) | pass |
| `codedFixedValue` run-time compilation as `foam` | pass — `.so` produced under `dynamicCode/` |
| same, as `root` | fails with the security error above (expected) |
| `surfaceCheck` on a generated fin STL | pass — *"Surface is closed. All edges connected to two faces"*, no illegal triangles |
